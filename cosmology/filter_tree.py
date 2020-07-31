import matplotlib
matplotlib.use('agg')
from galaxy_analysis.plot.plot_styles import *
import matplotlib.pyplot as plt

import numpy as np
import yt
import ytree
import sys

from galaxy_analysis.cosmology.nearest import nn_search



def select_root_halos(a, criteria, fields):
    """
    Select desired halos given a set of criteria

    Parameters
    -----------
    a        : ytree object
    criteria : sequence of conditional strings that can be
               applied to the tree nodes
    fields   : fields used in the criteria

    Returns
    --------
    roots    : List of root halos in tree node that satisfy the
               criteria
    """
    selection = a.select_halos(criteria, fields = fields)
    roots = [x.find_root() for x in selection]

    return roots


def get_neighbors(a, all_halos, n = 1, selection = None):
    """
    Get the properties of the nearest neighbor to all of the halos.
    """

    all_neighbors = []

    if selection is None:
        selection = a['id']==a['id']

    all_ids = a['id'][selection]
    x = a['x'][selection].to('kpc')
    y = a['y'][selection].to('kpc')
    z = a['z'][selection].to('kpc')

    nn_sorted, nn_distances = nn_search(x=x,y=y,z=z,
                                        return_distances=True, n = n)
    print(np.min(nn_distances), np.max(nn_distances))

    for halos in all_halos:

        nhalos = len(halos)
        nn_dict = {}
        fields = ['mass','virial_radius','id','x','y','z']
        for field in ['distance'] + fields:
            nn_dict[field] = np.zeros(nhalos)

        for i,h in enumerate(halos):
            host_index = np.argwhere(a['id'] == h['id'])[0]
            nn_index   = nn_sorted[host_index]

            for field in ['mass','virial_radius','id']:
                nn_dict[field][i] = a[field][nn_index]

            for field in ['x','y','z']:
                nn_dict[field][i] = a[field][nn_index].to('kpc').value / a.box_size.to('kpc').value

            nn_dict['distance'][i] = nn_distances[host_index] # in kpc

        all_neighbors.append(nn_dict)

    return all_neighbors

def get_root_halos(a, criteria, fields,
                      unique = True,
                      mutually_exclusive = ''):
    """
    Loop over a list of criteria and return root halos
    that satisfy each criteria.

    Parameters
    -----------
    a        : ytree object
    criteria : sequence of conditional strings that can be
               applied to the tree nodes
    fields   : fields used in the criteria

    unique   : bool, optional
               ensure that the root nodes within each criteria selection
               are unique. Default : True

    mutually_exlcusive : bool, optional
        Method to ensure that halos across criteria are mutually exclusive.
        By default this is '', and this is not used. If set to 'sequential',
        this is meant to be used with criteria that are (in effect) nested
        such that the first is the most restrictive and the last the most
        general, such that it is likely that all roots in each subsequent
        selection are also included in the previous. Default : 'sequential'

    Returns
    ---------
    final_roots : list
        list of the list of root nodes that satisfy the
        criteria
    """

    num_criteria = len(criteria)

    all_roots = []
    for c in criteria:
        all_roots.append(select_root_halos(a, c, fields))

    # now make sure these are unique

    final_roots = [None]*num_criteria
    if mutually_exclusive == 'sequential':

        final_roots[0] = all_roots[0] # always in this mode

        for i in np.arange(1, num_criteria):
            final_roots[i] = [all_roots[i][ii] for ii in np.arange(np.size(all_roots[i])) if (not(all_roots[i][ii] in all_roots[i-1]))]

    elif mutually_exclusive == '':
        final_roots = all_roots
    else:
        print("Mutually exclusive selection type " + mutually_exclusive + " not yet supported")
        raise(ValueError)

    if unique:
        #
        # make sure roots halos are unique within each list
        # np.unique does not work on TreeNode objects (unfortunately)
        # so doing this in a really dumb way
        #
        final_roots_ids = [None]*num_criteria
        for i in np.arange(num_criteria):
            final_roots_ids[i] = [x['id'] for x in final_roots[i]]

        final_roots = [None]*num_criteria
        for i in np.arange(num_criteria):
            final_roots[i] = [ a.query(ii) for ii in np.arange(a.size) if a.query(ii)['id'] in final_roots_ids[i]]


    return final_roots


def print_halo_info(a, treenodes, neighbors = None, header = None, outfile = None):
    """
    Print some information about the list of treenodes

    Parameters
    -----------
    a         :   ytree object from which the nodes come from
    treenodes : list of tree nodes
    neighbors : output from get_neighbors
    """

    if outfile is None:
        writefunc = lambda x : print(x, end='')
    else:
        f = open(outfile,'w')
        writefunc = f.write

    writefunc(header + '\n')

    if neighbors is None:
        for tn in treenodes:
            writefunc("%4i %5.5E %5.5E %.6f %.6f %.6f\n"%(tn['id'],tn['mass'],tn['virial_radius'],tn['x']/a.box_size,tn['y']/a.box_size,tn['z']/a.box_size))
    else:
        for i,tn in enumerate(treenodes):

            flag = "iso"
            if neighbors['distance'][i] < 1.5*tn['virial_radius']:
                flag = "near"

            writefunc("%4i %5.5E %5.5E %.6f %.6f %.6f %5.5E %4s %5.5E %5.5E %4i %.6f %.6f %.6f\n"%(tn['id'],tn['mass'],tn['virial_radius'],
                                                          tn['x']/a.box_size,tn['y']/a.box_size,tn['z']/a.box_size,
                                                          neighbors['distance'][i], flag, neighbors['mass'][i], neighbors['virial_radius'][i],
                                                          neighbors['id'][i],neighbors['x'][i],neighbors['y'][i],neighbors['z'][i]))


    if not (outfile is None):
        f.close()

    return

def plot_halos(a,treenodes, outfile = None):
    """
    For a list of tree nodes, plot their growth history
    """

    if outfile is None:
        outfile = "halo_growth_history.png"


    fig, ax = plt.subplots(1,2)
    fig.set_size_inches(12,6)

    for i,tn in enumerate(treenodes):
        x = tn['prog','redshift']
        y = tn['prog','mass']
        color = "C%01i"%(i)

        ax[0].scatter(x,y, color = color, s = 40)
        ax[0].plot(x,y,color=color,lw=3,label=tn['id'])

        y = tn['prog','virial_radius']
        ax[1].scatter(x,y,color=color,s=40)
        ax[1].plot(x,y,color=color,lw=3,label=tn['id'])

    if len(treenodes) > 9:
        ax[1].legend(loc='best',ncol=2)
    else:
        ax[0].legend(loc='best')
    ax[0].set_ylim(1.0E8,8.0E9)
    ax[0].set_ylabel(r"Rockstar Halo Mass (M$_{\odot}$)")
    ax[0].semilogy()

    ax[1].set_ylim(1, 100)
    ax[1].set_ylabel(r"Rockstar Virial Radius (kpccm)")
    ax[1].semilogy()



    for axes in ax:
        axes.set_xlabel("z")
        axes.set_xlim(15.0,4.5)

    ax[0].plot(ax[0].get_xlim(), [1.0E9,1.0E9], lw = 2, ls = '--' ,color = 'black')

    plt.minorticks_on()
    plt.tight_layout()
    fig.savefig(outfile) #, bbox_inches='tight', pad_inches=0.0)

    return


def my_halo_filter(tree_file):
    mass_string = "1E9"

    all_criteria = [ '(tree["tree", "redshift"] > 12.0) *'
                     '(tree["tree", "redshift"] < 20.0) *'
                     '(tree["tree", "mass"] > ' + mass_string + ')',

                     '(tree["tree", "redshift"] > 11.0) *'
                     '(tree["tree", "redshift"] < 12.0) *'
                     '(tree["tree", "mass"] > ' + mass_string + ')',

                     '(tree["tree", "redshift"] > 10.0) *'
                     '(tree["tree", "redshift"] < 11.0) *'
                     '(tree["tree", "mass"] > ' + mass_string + ')',

                     '(tree["tree", "redshift"] > 9.0) *'
                     '(tree["tree", "redshift"] < 10.0) *'
                     '(tree["tree", "mass"] > ' + mass_string + ')',

                     '(tree["tree", "redshift"] > 8.0) *'
                     '(tree["tree", "redshift"] < 9.0) *'
                     '(tree["tree", "mass"] > ' + mass_string + ')',

                     '(tree["tree", "redshift"] > 7.0) *'
                     '(tree["tree", "redshift"] < 8.0) *'
                     '(tree["tree", "mass"] > ' + mass_string + ')',

                     '(tree["tree", "redshift"] > 6.0) *'
                     '(tree["tree", "redshift"] < 7.0) *'
                     '(tree["tree", "mass"] > ' + mass_string + ')', 

                     '(tree["tree", "redshift"] > 5.0) *'
                     '(tree["tree", "redshift"] < 6.0) *'
                     '(tree["tree", "mass"] > ' + mass_string + ')']

    fields = ['redshift','mass']


    a = ytree.load(tree_file)

    all_halos      = get_root_halos(a, all_criteria, fields, unique=True, mutually_exclusive='sequential')
    selection      = a['mass'] > 1.0E8
    all_neighbors  = get_neighbors(a, all_halos, n = 1, selection=selection)

    outfiles, outplots = [], []
    for z in np.arange(12,4.5,-1):
        outfiles.append(mass_string+"_above_z%i.dat"%(z))
        outplots.append(mass_string+"_above_z%i.png"%(z))

    for i in np.arange(len(all_criteria)):
        print_halo_info(a, all_halos[i], header = all_criteria[i], neighbors=all_neighbors[i],outfile = outfiles[i])
        plot_halos(a,all_halos[i], outfile = outplots[i])

    return 

if __name__ == "__main__":


    my_halo_filter(str(sys.argv[1]))
