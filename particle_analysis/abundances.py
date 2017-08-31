#import yt.mods as yt
import yt
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from collections import Iterable, OrderedDict

import glob
import os
import h5py

# -- internal --
from galaxy_analysis.utilities import convert_abundances
from galaxy_analysis.utilities import utilities

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

    if particle_type == 'all':
        select = (ptype == ptype)
    else:
        select = (ptype == particle_type)

    birth_mass = birth_mass[select]

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

        mass1 = data[enzo_name1][select] * birth_mass
        mass2 = data[enzo_name2][select] * birth_mass

        aratios[ratio] = convert_abundances.abundance_ratio( (ele1, mass1),
                                                             (ele2, mass2),
                                                             'mass' )
    return aratios

def plot_time_evolution(h5file = 'abundances.h5', dir = './abundances/',
                        plot_type = 'standard', ds_list = None,
                        show_std = True, show_quartile = True):

    hf = h5py.File(dir + h5file, 'r')

    if ds_list is None:
        ds_list = hf.keys()

    if plot_type == 'standard':
        denom = 'Fe'

    for dsname in ds_list:

        # make one plot for each file
        t  = hf[dsname]['Time'].value
        ns = hf[dsname]['Nstars'].value
        ms = hf[dsname]['Mstars'].value

        abund    = hf[dsname]['abundances']
        elements = [x for x in abund.keys() if ( (x != 'H') and (x != 'He'))]
        nabundances = len(elements)

        outname = dir + dsname + '_abundances_time_evolution.png'

        nrow, ncol = utilities.rowcoldict[nabundances]

        fig, ax = plt.subplots(nrow,ncol)
        fig.set_size_inches(4*ncol, 4*nrow)

        i,j = 0,0
        for ele in elements:
            index = (i,j)
            if ele == denom:
                ele2 = 'H'
            else:
                ele2 = denom

            g = hf[dsname]['statistics']['10Myr'][ele][ele2]
            bins = hf[dsname]['statistics']['10Myr']['bins']
            x = (bins[1:] + bins[:-1])*0.5

            # plot average
            ax[index].plot(x, g['mean'], lw = 3, color = 'black', label = ele)

            if show_std:
                mean = np.array(g['mean'])
                std  = np.array(g['std'])
                std[std == 9999] = 0.0
                y1, y2 = mean - std, mean + std
                ax[index].fill_between(x, y1, y2, color = 'grey', alpha = 0.75, lw = 3)
            if show_quartile:
                ax[index].fill_between(x, g['Q1'], g['Q3'], color = 'grey', alpha = 0.33, lw = 3)

            ax[index].set_ylabel(r'log( [' + ele + '/' + ele2 + '])')
            ax[index].set_xlabel(r'Time (Myr)')
            ax[index].minorticks_on()

            ax[index].set_ylim(-3, 3)
            if ele == denom:
                ax[index].set_ylim(-10, 0)

            j = j + 1
            if j >= ncol:
                j = 0
                i = i + 1

        plt.tight_layout()

        plt.savefig(outname)
        plt.close()

    return

def plot_abundances(h5file = 'abundances.h5', dir = './abundances/', plot_type = 'standard', color_by_age=False,
                    ds_list = None):
    """
    Given an hdf5 file of stored abundances generated from Enzo particle data,
    plot all abundance ratios of certain plot type. Currently only supports the standard
    plot type, but this can be expanded upon easily. Saves images using data dump name in
    same directory as hdf5 file

    standard: log [ x /H ] vs. log [Fe / H]
    """

    hf = h5py.File(dir + h5file, 'r')

    # add new labels for plot types here, (e.g. X / Mg vs. Mg / H)
    xlabels = {'standard': r'log [Fe / H]'}
    ylabels = {'standard': r'log [X / Fe]'}

    if plot_type == 'standard':
        denom1 = 'Fe'
        denom2 = 'H'

    if ds_list is None: # do all
        ds_list = hf.keys()


    for dsname in ds_list:

        # make one plot for each file
        t  = hf[dsname]['Time'].value
        ns = hf[dsname]['Nstars'].value
        ms = hf[dsname]['Mstars'].value

        # always going to be N - 1

        abund = hf[dsname]['abundances']
        elements = [x for x in abund.keys() if (x!= denom1) and (x!=denom2)]
        nabundances = len(elements)

        outname = dir + dsname + '_abundances.png'

        nrow, ncol = utilities.rowcoldict[nabundances]

        fig, ax = plt.subplots(nrow,ncol)
        fig.set_size_inches(4*ncol,4*nrow)
        fig.suptitle("Time = %.1f Myr - Nstars = %.2E - Mstars = %.2E Msun"%(t, ns, ms))

        i,j = 0,0

        for ele in elements:

            if color_by_age:
                age = np.array(t - hf[dsname]['creation_time'].value)
                c = ax[(i,j)].scatter( np.array(abund[denom1][denom2].value), np.array(abund[ele][denom1].value), s = 7.5, alpha = 0.25,
                                   c = age, label=ele, cmap='algae')
            else:
                ax[(i,j)].scatter( abund[denom1][denom2].value, abund[ele][denom1].value, s =15, alpha =0.75,
                                                       color = 'black', label = ele)

            ax[(i,j)].set_xlabel(xlabels[plot_type])
            ax[(i,j)].set_ylabel(ylabels[plot_type])

            if ele in ['Eu','Ba','Y']:
                ax[(i,j)].set_ylim(-8,0)
            else:
                ax[(i,j)].set_ylim(-4,4)

            ax[(i,j)].set_xlim(-15, 0)
            ax[(i,j)].minorticks_on()

            ax[(i,j)].legend(loc='upper right')


            j = j + 1
            if j >= ncol:
                j = 0
                i = i + 1

        plt.tight_layout()
        cbar = fig.colorbar(c)
        cbar.ax.set_label('Age (Myr)')
        plt.savefig(outname)
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

    if not os.path.isfile(dir + outfile) or overwrite:
        hf = h5py.File(dir + outfile, 'w')
    else:
        hf = h5py.File(dir + outfile, 'a')

    ds_list = np.sort( glob.glob('./DD????/DD????') )
    times   = np.zeros(np.size(ds_list))

    # get elements present:
    ds              = yt.load(ds_list[-1])
    fields          = ds.field_list
    elements = utilities.species_from_fields(fields, include_primordial=True)
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

        if groupname in hf and not overwrite:
            continue # skip this one, it already exists

        g = hf.create_group(groupname)
        g.create_dataset('Time'  , data = ds.current_time.convert_to_units('Myr').value)

        if ('io', 'particle_type') in ds.field_list:

            aratios = compute_aratio(ds, data, ratios)
            MS = data['particle_type'] == 11

            g.create_dataset('Nstars', data = np.size(data['particle_mass'][ MS]))
            g.create_dataset('Mstars', data = np.sum( data['particle_mass'][ MS].convert_to_units('Msun').value))
            g.create_dataset('creation_time', data = data['creation_time'][MS].convert_to_units('Myr').value)
            g.create_dataset('birth_mass', data = data['birth_mass'][MS].value)

            sg = hf.create_group(groupname + '/abundances')
            for abundance in aratios.keys():
                sg.create_dataset( abundance, data = aratios[abundance])

            # now compute statistics on all of the data, and store them
            statgroup = hf.create_group(groupname + '/statistics')
            all = statgroup.create_group('all')
            for abundance in aratios.keys():
                stats = utilities.compute_stats(aratios[abundance], return_dict = True)
                g = all.create_group(abundance)
                for k in stats.keys():
                    g.create_dataset(k, data = stats[k])

            #
            # now do it in 10 Myr bins
            #
            aratios = compute_aratio(ds, data, ratios, particle_type = 'all')

            g = statgroup.create_group('10Myr')
            dt = 10
            t  = ds.current_time.convert_to_units('Myr').value
            tmax = np.around(t, decimals = -len(str(dt)) + 1)
            if tmax < t:
                tmax = tmax + dt
            tbins = np.arange(0.0, tmax + 0.5*dt, dt)

            index = np.digitize(data['creation_time'].convert_to_units('Myr').value, tbins)
            hist, bins  = np.histogram(data['creation_time'].convert_to_units('Myr').value, bins = tbins)
            g.create_dataset('bins', data = tbins)
            g.create_dataset('hist', data = hist)

            stats_array_dict = {}
            for abundance in aratios.keys():
                stats_array_dict[abundance] = {}
                for k in stats.keys():
                    stats_array_dict[abundance][k] = np.zeros(np.size(tbins) - 1)

            for i in np.arange(np.size(tbins)-1):
                for abundance in aratios.keys():
                    if i == 0:
                        sub_g = g.create_group(abundance)
                    if hist[i] > 0:
                        stats = utilities.compute_stats(aratios[abundance][index == i+1], return_dict = True)
                        for k in stats.keys():
                            stats_array_dict[abundance][k][i] = stats[k]
                    else:
                        for k in stats.keys():
                            stats_array_dict[abundance][k][i] = 9999

            for abundance in aratios.keys():
                g = hf[groupname + '/statistics/10Myr/' + abundance]
                for k in stats.keys():
                    g.create_dataset(k, data = stats_array_dict[abundance][k])

            # ------------ can do a correlation across time bins here too --------- 
            # Pick some time t_o, for the ith bin past t_o, do correlation between
            # those two populations of stars
            # x  = np.array([stars in t_o bin] + [stars in t_i bin])
            # corr[i] = np.correlate(x,x, mode = 'full')
            # allow to plot correlation as a function of time.


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

    plot_abundances(plot_type = 'standard', color_by_age = True)
    plot_time_evolution()
