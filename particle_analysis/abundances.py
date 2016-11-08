import yt.mods as yt
import numpy as np
import matplotlib.pyplot as plt
from collections import Iterable, OrderedDict

import glob
import os
import h5py

# -- internal --
from galaxy_analysis.utilities import convert_abundances

def compute_aratio(ds, data, ratios, particle_type = 11):
    """
    For a given data set (ds) and data structure, generates abundances ratios
    desired (e.g. 'Fe/H') as normalized to solar value. This returns a
    dictionary of the abundance ratios, where the ratios are given in normalized
    log units, e.g. [Fe / H] = log (Fe/H star) - log (Fe/H solar). This is
    the standard observational definition
    """

    birth_mass = data['birth_mass'].value * yt.units.Msun
    ptype      = data['particle_type']

    if not isinstance(ratios, Iterable):
        ratios = [ratios]
    #
    # split string
    #
    aratios = OrderedDict()

    for ratio in ratios:
        if '/' in ratio:
            ele1, ele2 = ratio.rsplit('/')
        else:
            print "Must provide abundance ratio string of style: Fe/H"
            return

        enzo_name1 = ('io','particle_' + ele1 + '_fraction')
        enzo_name2 = ('io','particle_' + ele2 + '_fraction')

        mass1 = data[enzo_name1] * birth_mass
        mass2 = data[enzo_name2] * birth_mass

        aratios[ratio] = convert_abundances.abundance_ratio( (ele1, mass1), 
                                                             (ele2, mass2), 
                                                             'mass' )
    return aratios

def plot_abundances(h5file = 'abundances.h5', dir = './abundances/', plot_type = 'standard', color_by_age=False):
    """
    Given an hdf5 file of stored abundances generated from Enzo particle data,
    plot all abundance ratios of certain plot type. Currently only supports the standard
    plot type, but this can be expanded upon easily. Saves images using data dump name in
    same directory as hdf5 file

    standard: log [ x /H ] vs. log [Fe / H]
    """

    hf = h5py.File(dir + h5file, 'r')

    xlabels = {'standard': r'log [Fe / H]'}
    ylabels = {'standard': r'log [X / Fe]'}

    if plot_type == 'standard':
        denom1 = 'Fe'
        denom2 = 'H'

    for dsname in hf.keys():

        # make one plot for each file
        t  = hf[dsname]['Time'].value
        ns = hf[dsname]['Nstars'].value
        ms = hf[dsname]['Mstars'].value

        # always going to be N - 1

        abund = hf[dsname]['abundances']
        elements = [x for x in abund.keys() if (x!= denom1) and (x!=denom2)]
        nabundances = len(elements)

        nrow, ncol = 2 , 4

        fig, ax = plt.subplots(nrow,ncol)
        fig.set_size_inches(4*ncol,4*nrow)
        fig.suptitle("Time = %.1f Myr - Nstars = %.2E - Mstars = %.2E Msun"%(t, ns, ms))

        i,j = 0,0

        for ele in elements:

            if color_by_age:
                c = ax[(i,j)].scatter( abund[ele][denom2], abund[ele][denom1], s = 15, alpha = 0.5,
                                   c = t - hf[dsname]['creation_time'].value, label=ele, cmap='viridis')
                plt.colorbar()
            else:
                ax[(i,j)].scatter( abund[ele][denom2], abund[ele][denom1], s =15, alpha =0.75,
                                                       color = 'black', label = ele)

            ax[(i,j)].set_xlabel(xlabels[plot_type])
            ax[(i,j)].set_ylabel(ylabels[plot_type])
            ax[(i,j)].set_xlim(-10, 1)
            ax[(i,j)].set_ylim(-4,4)
            ax[(i,j)].minorticks_on()

            ax[(i,j)].legend(loc='upper right')


            j = j + 1
            if j >= ncol:
                j = 0
                i = i + 1

        plt.tight_layout()
        plt.savefig(dir + dsname + '_abundances.png')
        plt.close()

    return

def generate_abundances(outfile = 'abundances.h5', dir = './abundances/', overwrite = False):
    """
    Function to generate hdf5 output file containing abundance ratios for all metal
    species in all main sequence stars in all data sets in the given directory.

    Abundances can be plotted with:
        > abundances.plot_abundances()
    """
    #
    # do this for all
    #
    if not os.path.exists(dir):
        os.makedirs(dir)

    if not os.path.isfile(dir + outfile):
        hf = h5py.File(dir + outfile, 'w')
    elif overwrite:
        print "WARNING FILE EXISTS: Overwrite set to false - finishing"
        return
    else:
        hf = h5py.File(dir + outfile, 'w')

    ds_list = np.sort( glob.glob('./DD????/DD????') )
    times   = np.zeros(np.size(ds_list))

    # get elements present:
    ds              = yt.load(ds_list[-1])
    fields          = ds.field_list
    element_fields  = [x[1] for x in fields if ('particle_' in x[1] and\
                                                '_fraction' in x[1]) and\
                                               ('all' in x[0])]
    elements = np.sort([x.rsplit('_')[1] for x in element_fields])
    metals   = [x for x in elements if (x != 'H' and x != 'He')]
    ratios   = [ x +'/H' for x in metals]

    if 'Mg' in metals:
        ratios = ratios + [ x + '/Mg' for x in metals]

    if 'Fe' in metals:
        ratios = ratios + [ x + '/Fe' for x in metals]

    for i, dsname in enumerate(ds_list):
        ds   = yt.load(dsname)
        data = ds.all_data()

        groupname = dsname.rsplit('/')[1]
        g = hf.create_group(groupname)
        g.create_dataset('Time'  , data = ds.current_time.convert_to_units('Myr').value)

        if ('io', 'particle_type') in ds.field_list:

            aratios = compute_aratio(ds, data, ratios)

            g.create_dataset('Nstars', data = np.size(data['particle_mass'][ data['particle_type'] == 11]))
            g.create_dataset('Mstars', data = np.sum( data['particle_mass'][ data['particle_type'] == 11].convert_to_units('Msun').value))
            g.create_dataset('creation_time', data = data['creation_time'][data['particle_type'] == 11].convert_to_units('Myr').value)
            g.create_dataset('birth_mass', data = data['birth_mass'][data['particle_type'] == 11].value)

            sg = hf.create_group(groupname + '/abundances')
            for abundance in aratios.keys():
                sg.create_dataset( abundance, data = aratios[abundance])
        else:
            g.create_dataset('Nstars', data = 0.0)
            g.create_dataset('Mstars', data = 0.0)
            sg = hf.create_group(groupname + '/abundances')
            for abundance in aratios.keys():
                sg.create_dataset( abundance, data = 0.0)


    hf.close()

    return

if __name__=='__main__':

    generate_abundances()
    plot_abundances(plot_type = 'standard')
