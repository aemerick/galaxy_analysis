import matplotlib
matplotlib.use('agg')
from galaxy_analysis.plot.plot_styles import *
import matplotlib.pyplot as plt

import numpy as np
import yt
import ytree
import sys




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


def print_halo_info(a, treenodes, header = None, outfile = None):
    """
    Print some information about the list of treenodes

    Parameters
    -----------
    a         :   ytree object from which the nodes come from
    treenodes : list of tree nodes
    """

    if outfile is None:
        writefunc = lambda x : print(x, end='')
    else:
        f = open(outfile,'w')
        writefunc = f.write

    writefunc(header + '\n')

    for tn in treenodes:
        writefunc("%4i %5.5E %5.5E %.6f %.6f %.6f\n"%(tn['id'],tn['mass'],tn['virial_radius'],tn['x']/a.box_size,tn['y']/a.box_size,tn['z']/a.box_size))


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

if __name__ == "__main__":

    #
    #
    #
    all_criteria = [ '(tree["tree", "redshift"] > 8.0) *'
                     '(tree["tree", "redshift"] < 12.5) *'
                     '(tree["tree", "mass"] > 1.0E9)',

                     '(tree["tree", "redshift"] > 7.0) *'
                     '(tree["tree", "redshift"] < 8.1) *'
                     '(tree["tree", "mass"] > 1.0E9)',

                     '(tree["tree", "redshift"] > 6.0) *'
                     '(tree["tree", "redshift"] < 7.1) *'
                     '(tree["tree", "mass"] > 1.0E9)']

    fields = ['redshift','mass']


    a = ytree.load(str(sys.argv[1]))

    all_halos = get_root_halos(a, all_criteria, fields, unique=True, mutually_exclusive='sequential')

    outfiles = ['1E9_above_z8.dat','1E9_above_z7.dat', '1E9_above_z6.dat']
    outplots = ['1E9_above_z8.png','1E9_above_z7.png', '1E9_above_z6.png']
    for i in np.arange(len(all_criteria)):
        print_halo_info(a,all_halos[i], header = all_criteria[i], outfile = outfiles[i])
        plot_halos(a,all_halos[i], outfile = outplots[i])
