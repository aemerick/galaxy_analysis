import yt
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc, cm
from collections import Iterable, OrderedDict

COMPUTE_ACF = False # global disable for now - not working and time consuming

fsize = 17
rc('text', usetex=False)
rc('font', size=fsize)#, ftype=42)
line_width = 3
point_size = 30
plasma     = cm.get_cmap('plasma')

import glob
import os
import h5py
from scipy.interpolate import interp1d

# -- parallel --
from multiprocessing import Pool
from contextlib import closing
import itertools

# -- internal --
from galaxy_analysis.plot.plot_styles import *
from galaxy_analysis.utilities import convert_abundances
from galaxy_analysis.utilities import utilities
from galaxy_analysis.analysis import Galaxy
from astroML.density_estimation import bayesian_blocks, knuth_bin_width

def compute_mass_fractions(ds, data, elements, particle_type = 11):

    ptype = data['particle_type']
    if particle_type == 'all':
        select = (ptype ==ptype)
    else:
        select = (ptype == particle_type)


    mass_fractions = OrderedDict()

    for e in elements:
        mass_fractions[e] = data['particle_' + e + '_fraction'][select]

    return mass_fractions


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

        if ele1 == ele2:
            continue

        aratios[ratio] = data[('io','particle_' + ele1 + '_over_' + ele2)][select]
#
#        enzo_name1 = ('io','particle_' + ele1 + '_fraction')
#        enzo_name2 = ('io','particle_' + ele2 + '_fraction')
#
#        mass1 = data[enzo_name1][select] * birth_mass
#        mass2 = data[enzo_name2][select] * birth_mass
#
#        aratios[ratio] = convert_abundances.abundance_ratio( (ele1, mass1),
#                                                             (ele2, mass2),
#                                                             'mass' )
    return aratios

def plot_acf(h5file = 'star_abundances.h5', dir = './abundances/', ds_list = None):
    """
    Plots the time evolution of all abundance ratios [X/Fe]. Optionally
    can include bands giving the standard deviation and 1st and 3rd 
    quartile of the distributions.
    """

    hf = h5py.File(dir + h5file, 'r')

    if ds_list is None:
        ds_list = hf.keys()


    for dsname in ds_list:

        # make one plot for each file
        t  = hf[dsname]['Time'].value
        ns = hf[dsname]['Nstars'].value
        ms = hf[dsname]['Mstars'].value

        abund    = hf[dsname]['abundances']
        elements = utilities.sort_by_anum([x for x in abund.keys() if ( (x != 'H') and (x != 'He'))])
        nabundances = len(elements)
        outname = dir + dsname + '_abundances_acf.png'
        nrow, ncol = utilities.rowcoldict[nabundances]

        fig, ax = plt.subplots(nrow,ncol)
        fig.set_size_inches(4*ncol, 4*nrow)

        i,j = 0,0
        denom = 'Fe'
        for ele in elements:
            index = (i,j)
            if ele == 'Fe':
                ele2 = 'H'
            else:
                ele2 = denom

            g   = hf[dsname]['statistics']
            acf = g['all_particles'][ele][ele2]['acf']
            x   = g['all_particles'][ele][ele2]['acf_bins']
            print x.value, acf.value
    #        acf = g['0Myr'][ele][ele2]['acf']
            #print g['0Myr'][ele][ele2]['interp_mean'], np.array(acf)
    #        x   = g['0Myr']['bins'][-len(acf)-1:]
            x   = 0.5 * (x[1:] + x[:-1])
            x   = x - x[0]

            ax[index].plot(x, acf, lw = 3, color = 'black', label = ele)

            ax[index].set_ylabel(r'ACF log([' + ele + '/' + ele2 + '])')
            ax[index].set_xlabel(r'Time (Myr)')
            ax[index].minorticks_on()

            ax[index].set_ylim(-0.1, 1.2)

            j = j + 1
            if j >= ncol:
                j = 0
                i = i + 1

        plt.tight_layout()

        plt.savefig(outname)
        plt.close()

    return


def plot_time_evolution(h5file = 'star_abundances.h5', dir = './abundances/',
                        plot_type = 'standard', ds_list = None,
                        show_std = True, show_quartile = True):
    """
    Plots the time evolution of all abundance ratios [X/Fe]. Optionally
    can include bands giving the standard deviation and 1st and 3rd 
    quartile of the distributions.
    """

    hf = h5py.File(dir + h5file, 'r')

    if ds_list is None:
        ds_list = hf.keys()

    if plot_type == 'standard':
        denom = 'Fe'

    for dsname in [ ds_list[-1] ]:

        # make one plot for each file
        t  = hf[dsname]['Time'].value
        t  = t # - t[0]
        ns = hf[dsname]['Nstars'].value
        ms = hf[dsname]['Mstars'].value

        abund    = hf[dsname]['abundances']
        elements = utilities.sort_by_anum([x for x in abund.keys() if ( (x != 'H') and (x != 'He') and (not 'alpha' in x))])
        elements = elements + ['alpha']
        nabundances = len(elements)

        for dt in [1, 10]:
            outname = dir + dsname + '_abundances_time_%i_evolution.png'%(dt)

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

                g = hf[dsname]['statistics']['%iMyr'%(dt)][ele][ele2]
                bins = hf[dsname]['statistics']['%iMyr'%(dt)]['bins']
                x = (bins[1:] + bins[:-1])*0.5

                # plot average
                mean = np.array(g['mean'])
                ax[index].plot(x, mean, lw = 3, color = 'black', label = ele)

                if show_std:
                    std  = np.array(g['std'])
                    y1, y2 = mean - std, mean + std
                    ax[index].fill_between(x, y1, y2, color = 'grey', alpha = 0.5, lw = 1.5)
                if show_quartile:
                    ax[index].fill_between(x, g['Q1'], g['Q3'], color = 'black', alpha = 0.5, lw = 1.5)

                ax[index].set_ylabel(r'log( [' + ele + '/' + ele2 + '])')
                ax[index].set_xlabel(r'Time (Myr)')
                ax[index].minorticks_on()

                ax[index].set_ylim(-3, 3)
                if ele == denom:
                    ax[index].set_ylim(-10, 0)

                xmin = np.min( np.array(x))
                ax[index].set_xlim( xmin , bins[-1])

                j = j + 1
                if j >= ncol:
                    j = 0
                    i = i + 1

            plt.tight_layout()

            plt.savefig(outname)
            plt.close()

    return

def plot_abundances(h5file = 'star_abundances.h5', dir = './abundances/', plot_type = 'standard', color_by_age=False,
                    ds_list = None, show_average = False):
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
        elements = utilities.sort_by_anum([x for x in abund.keys() if (x!= denom1) and (x!=denom2) and (not 'alpha' in x)])
        elements = elements + ['alpha']
        nabundances = len(elements)

        outname = dir + dsname + '_abundances.png'

        nrow, ncol = utilities.rowcoldict[nabundances]

        fig, ax = plt.subplots(nrow,ncol)
        fig.set_size_inches(4*ncol,4*nrow)
        fig.suptitle("Time = %.1f Myr - Nstars = %.2E - Mstars = %.2E Msun"%(t, ns, ms))

        i,j = 0,0

        for ele in elements:
            index = (i,j)

            xpoints = np.array( abund[denom1][denom2].value)
            ypoints = np.array( abund[ele][denom1].value)
            if color_by_age:
                age = np.array(t - hf[dsname]['creation_time'].value)
                c = ax[index].scatter( xpoints, ypoints, s = 7.5, alpha = 0.25,
                                   c = age, label=ele, cmap='plasma')
            else:
                ax[index].scatter( xpoints, ypoints, s =15, alpha =0.75,
                                                       color = 'black', label = ele)

            if show_average:
                # need to bin the data
                xbins = np.arange(-15, 4.05, 0.25) # larger than needed bins
                yavg  = np.ones(np.size(xbins)-1) * 9999
                for ii in np.arange(0, len(xbins) - 1):
                    yavg[ii] = np.average( ypoints[ (xpoints > xbins[ii]) * (xpoints <= xbins[ii+1])])
                xcent = 0.5 * (xbins[1:] + xbins[:-1])
                yavg[yavg==9999] = None
#                if color_by_age:
#                    ax[index].plot( xcent, yhist, lw = 3, alpha = 0.75, c = age, cmap = 'plasma')
#                else:
                ax[index].step(xcent, yavg, lw = 3, color = 'black', where = 'post')

            ax[index].set_xlabel(xlabels[plot_type])
            ax[index].set_ylabel(r'log([' + ele + '/' + denom1 +'])')

            if ele in ['Eu','Ba','Y']:
                ax[index].set_ylim(-8,0)
            else:
                ax[index].set_ylim(-4,4)

            ax[index].set_xlim(-15, 0)
            ax[index].minorticks_on()

#            ax[index].legend(loc='upper right')


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

def single_MDF(x, bins = None, norm = 'peak', ax = None, label = True, **kwargs):

    if bins is None:
        bins = np.linspace(np.min(x), np.max(x), 50)

    fig = None
    if ax is None:
        fig, ax = plt.subplots()

    hist, bins = np.histogram(x, bins = bins) #knuth_bin_width(points,True,disp=False)[1]) #bayesian_blocks(points))

    if norm == 'fraction':
        A = 1.0 / (1.0 * np.sum(hist))
        if label:
            ax.set_ylabel(r'N / N$_{\rm total}$')
    elif norm == 'peak':
        max = np.max( [np.max(hist), 1.0] )
        A = 1.0 / (1.0 * max)
        if label:
            ax.set_ylabel(r'N / N$_{\rm max}$')
        ax.set_ylim(0,1)
    elif norm == 'PDF':
        A = 1.0 / (bins[1:] - bins[:-1])
        if label:
            ax.set_ylabel(r'PDF')

    y = hist * A
    plot_histogram(ax, bins, y, **kwargs)

    if not (fig is None):
        fig.set_size_inches(8,8)

        return fig, ax

    else:
        return


def plot_MDF(h5file = 'star_abundances.h5', dir = './abundances/', plot_type = 'standard',
             ds_list = None, show_average = False, norm = 'peak'):
    """
    Plot the MDF of each element (for now, hard coded as all over Fe and [Fe/H]). Choose
    normalization as: 
        1) 'fraction' : (N / N_total)
        2) 'peak'     : (N / N_max)  - normalizes peak to 1.0
        3) 'PDF'      : (dN / dx)    - independent of bin sizing
    where N is number in a given bin, N_total is total number of stars, and N_max is
    maximum number in any given bin.
    """

    hf = h5py.File(dir + h5file, 'r')

    # add new labels for plot types here, (e.g. X / Mg vs. Mg / H)
    xlabels = {'standard': r'log [Fe / H]'}

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
        elements = utilities.sort_by_anum([x for x in abund.keys() if (x!= denom1) and (x!=denom2) and (not 'alpha' in x)])
        elements = elements + ['alpha']
        nabundances = len(elements)

        outname = dir + dsname + '_MDF.png'

        nrow, ncol = utilities.rowcoldict[nabundances]

        fig, ax = plt.subplots(nrow,ncol)
        fig.set_size_inches(4*ncol,4*nrow)
        fig.suptitle("Time = %.1f Myr - Nstars = %.2E - Mstars = %.2E Msun"%(t, ns, ms))

        i,j = 0,0
        for ele in elements:
            index = (i,j)
            points = np.array(abund[ele][denom1].value)

            bins = np.arange(-3,3.1,0.1)
            single_MDF(points, bins = bins, norm = 'peak', ax = ax[index])

            j = j + 1
            if j >= ncol:
                j = 0
                i = i + 1

            ax[index].set_xlabel(r'['+ ele + '/' + denom1 + ']')

            ax[index].minorticks_on()

        # plot Fe over H
        bins = np.arange(-8,-1,0.1)
        single_MDF(abund['Fe']['H'].value, bins = bins, norm = 'peak', ax = ax[(i,j)])
        ax[(i,j)].set_xlabel('[Fe/H]')
        ax[(i,j)].minorticks_on()


        plt.tight_layout()
        plt.savefig(outname)
        plt.close()


    return


def generate_abundances(ds_list = None, outfile = 'star_abundances.h5', dir = './abundances/', overwrite = False):
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

    if ds_list is None:
        ds_list = np.sort( glob.glob('./DD????/DD????') )
        times   = np.zeros(np.size(ds_list))
    elif (not (type(ds_list) is list)):
        # assume a single string passed
        ds_list = [ds_list]

    # get elements present:
    ds              = yt.load(ds_list[-1])
    fields          = ds.field_list
    elements = utilities.species_from_fields(fields, include_primordial=True)
    metals   = [x for x in elements if (x != 'H' and x != 'He')]
    metals   = metals + ['alpha', 'alpha_5'] # add these two by hand for aggregate metal abundances
    ratios   = [ x +'/H' for x in metals]

    if 'Mg' in metals:
        ratios = ratios + [ x + '/Mg' for x in metals]

    if 'Fe' in metals:
        ratios = ratios + [ x + '/Fe' for x in metals]

    if 'O' in metals:
        ratios = ratios + [ x + '/O' for x in metals]

    if 'C' in metals:
        ratios = ratios + [ x + '/C' for x in metals]

    if 'Ba' in metals:
        ratios = ratios + [ x + '/Ba' for x in metals]

#    if 'alpha' in metals:
#        ratios = ratios + [ x + '/alpha' for x in metals]

    for i, dsname in enumerate(ds_list):
        ds   = yt.load(dsname)
        data = ds.all_data()

        groupname = dsname.rsplit('/')[1]

        if groupname in hf and not overwrite:
            continue # skip this one, it already exists

        if ('io','particle_type') in ds.field_list:
            g = hf.create_group(groupname)
            g.create_dataset('Time'  , data = ds.current_time.convert_to_units('Myr').value)

#        if ('io', 'particle_type') in ds.field_list:

            #
            # Compute and store abundance ratios and relevant properties for all MS stars
            #
            aratios         = compute_aratio(ds, data, ratios) # by default, only does MS stars
            mass_fractions  = compute_mass_fractions(ds, data, elements)

            MS = data['particle_type'] == 11

            Nstars = np.size(data['particle_mass'][MS])
            g.create_dataset('Nstars', data = Nstars)
            g.create_dataset('Mstars', data = np.sum( data['particle_mass'][ MS].convert_to_units('Msun').value))
            g.create_dataset('creation_time', data = data['creation_time'][MS].convert_to_units('Myr').value)
            g.create_dataset('birth_mass', data = data['birth_mass'][MS].value)
            g.create_dataset('metallicity', data = data['metallicity_fraction'][MS].value)
            spatial = g.create_group('kinematics')

            r  = np.zeros(Nstars)
            vr = np.zeros(Nstars)
            for i, xname in enumerate(['x','y','z']):
                x  = (data['particle_position_' + xname][MS] - ds.domain_center[i]).convert_to_units('pc').value
                vx = (data['particle_velocity_' + xname][MS]).convert_to_units('km/s').value
                r  += x**2
                vr += vx**2
                spatial.create_dataset( xname, data = x)
            spatial.create_dataset('r', data = np.sqrt(r))
            spatial.create_dataset('vr', data = np.sqrt(vr))

#
            mf = hf.create_group(groupname + '/mass_fractions')
            for e in elements:
                mf.create_dataset(  e, data = mass_fractions[e])
            mf_statgroup = hf.create_group(groupname + '/mass_fraction_statistics')
            all = mf_statgroup.create_group('all_MS')
            for e in elements:
                stats = utilities.compute_stats( mass_fractions[e], return_dict = True)
                g     = all.create_group(e)
                for k in stats.keys():
                    g.create_dataset(k, data = stats[k])

#
            sg = hf.create_group(groupname + '/abundances')
            for abundance in aratios.keys():
                sg.create_dataset( abundance, data = aratios[abundance])

            # now compute statistics on the MS stars, and store them
            #
            statgroup = hf.create_group(groupname + '/statistics')
            all = statgroup.create_group('all_MS')
            for abundance in aratios.keys():
                stats = utilities.compute_stats(aratios[abundance], return_dict = True)
                g = all.create_group(abundance)
                for k in stats.keys():
                    g.create_dataset(k, data = stats[k])

            #
            # Now, do this for all particles, regardless of type.
            # Aka... ignore observational / physical reality and treat them all as tracers
            #
            aratios = compute_aratio(ds, data, ratios, particle_type = 'all')
            tracers = statgroup.create_group('all_particles')
            for abundance in aratios.keys():
                stats = utilities.compute_stats(aratios[abundance], return_dict = True)
                g     = tracers.create_group(abundance)

                if COMPUTE_ACF: # hide this for now - not working
                    t    = data['creation_time'].convert_to_units('Myr').value
                    t_n  = t - np.min(t)
                    dt   = 1.0

                    bins = np.arange(0.0, np.ceil(np.max(t_n)) + dt, dt)
                    y    = aratios[abundance]
                    y    = y + np.min(y)*2.0
                    dy   = np.abs(0.001 * y) # error should be irrelevant, but must be non-zero
                    dy[dy == 0.0] = 0.00001
                    acf, acf_error, acf_bins = utilities.acf(t_n, y, dy = dy, bins = bins)

                    stats['acf'] = acf
                    stats['acf_error'] = acf_error
                    stats['acf_bins']  = acf_bins

                for k in stats.keys():
                    g.create_dataset(k, data = stats[k])

            mass_fractions      = compute_mass_fractions(ds, data, elements, particle_type = 'all')
            tracers = mf_statgroup.create_group('all_particles')
            for e in elements:
                stats = utilities.compute_stats(mass_fractions[e], return_dict = True)
#
# left off here
#

            g = mf_statgroup.create_group("cumulative")
            t = ds.current_time.convert_to_units('Myr').value
            tmax = np.ceil(t)
            tbins = np.arange(0.0, tmax + 0.1, 0.5)
            hist,bins = np.histogram(data['creation_time'].convert_to_units('Myr').value, bins = tbins)
            g.create_dataset('bins', data = tbins)
            g.create_dataset('hist', data = np.array(hist))
            t_form = data['creation_time'].convert_to_units('Myr').value
            lifetime = data[('io','particle_model_lifetime')].convert_to_units('Myr').value
            age = t - t_form

            mf_stats_array_dict = {}
            for e in elements:
                mf_stats_array_dict[e] = {}
                for k in stats.keys():
                    mf_stats_array_dict[e][k] = np.zeros(np.size(tbins)-1)

            for i in np.arange(np.size(tbins)-1):

                age = tbins[i] - t_form
                selection = (age >= 0.0)*(age <= lifetime)
                for e in elements:
                    if i == 0:
                        sub_g = g.create_group(e)

                    if np.size(age[selection]) > 1:
                        stats = utilities.compute_stats(mass_fractions[e][selection], return_dict = True) # +1 b/c index starts at 1
                        for k in stats.keys():
                            mf_stats_array_dict[e][k][i] = stats[k]
                    else:
                        for k in stats.keys():
                            mf_stats_array_dict[e][k][i] = None

            for e in elements:
                g = hf[groupname + '/mass_fraction_statistics/cumulative/' + e]
                for k in mf_stats_array_dict[e].keys():
                    g.create_dataset(k, data = mf_stats_array_dict[e][k])

            for dt in [0.1, 1, 10]:
                g = mf_statgroup.create_group('%iMyr'%(dt))
                t  = ds.current_time.convert_to_units('Myr').value
                tmax = np.around(t, decimals = -len(str(dt)) + 1)
                if tmax < t:
                    tmax = tmax + dt
                tbins = np.arange(0.0, tmax + 0.5*dt, dt)

                index = np.digitize(data['creation_time'].convert_to_units('Myr').value, tbins)
                hist, bins  = np.histogram(data['creation_time'].convert_to_units('Myr').value, bins = tbins)
                g.create_dataset('bins', data = tbins)
                g.create_dataset('hist', data = np.array(hist))

                mf_stats_array_dict = {}
                for e in elements:
                    mf_stats_array_dict[e] = {}
                    for k in stats.keys():
                        mf_stats_array_dict[e][k] = np.zeros(np.size(tbins) - 1)

                for i in np.arange(np.size(tbins)-1):
                    for e in elements:
                        if i == 0:
                            sub_g = g.create_group(e)
                        if hist[i] > 0:
                            stats = utilities.compute_stats(mass_fractions[e][index == i+1], return_dict = True) # +1 b/c index starts at$
                            for k in stats.keys():
                                mf_stats_array_dict[e][k][i] = stats[k]
                        else:
                            for k in stats.keys():
                                mf_stats_array_dict[e][k][i] = None

                for e in elements:
                    # - - - - - Produce a gap-less, interpolated mean to compute the ACF
                    if False: # don't do this anymore
                        first        = np.where( np.logical_not(np.isnan( mf_stats_array_dict[e]['mean'] )))[0][0]
                        mean         = mf_stats_array_dict[e]['mean'][first:]
                        select       = np.logical_not(np.isnan(mean))
                        clean_mean   = mean[select]
                        tcent        = 0.5 * (tbins[1:] + tbins[:-1])
                        tcent        = tcent[first:]
                        clean_t      = tcent[select]
                        f_interp     = interp1d(clean_t, clean_mean)
                        interp_mean  = mean
                        interp_mean[np.logical_not(select)] = f_interp( tcent[np.logical_not(select)] )
                        mf_stats_array_dict[e]['interp_mean'] = interp_mean
                        mf_stats_array_dict[e]['acf'] = utilities.acf(interp_mean, nlags = len(tcent))

                    g = hf[groupname + '/mass_fraction_statistics/%iMyr/'%(dt) + e]
                    for k in mf_stats_array_dict[e].keys():
                        g.create_dataset(k, data = mf_stats_array_dict[e][k])


            #
            # now do it in time bins to get time evolution
            #

            # First, lets do the observational version, where we compute the total
            #  MDF at each point in time (using all stars) and compute median and spread, etc.
            #  next we will do the instantaneous (binned) version of this
            g = statgroup.create_group("cumulative")
            t = ds.current_time.convert_to_units('Myr').value
            tmax = np.ceil(t)
            tbins = np.arange(0.0, tmax + 0.1, 0.5) # can go arbitrarily small here
            hist, bins = np.histogram(data['creation_time'].convert_to_units('Myr').value, bins = tbins)
            g.create_dataset('bins', data = tbins)
            g.create_dataset('hist', data = np.array(hist))

            t_form   = data['creation_time'].convert_to_units('Myr').value
            # unfortunately we can't use dynamical_time because we are doing this for a single data output
            # and want to get WD and SN remnant stars binned appropriately, but their dynamical_time values change
            # when they form...
            lifetime = data[('io','particle_model_lifetime')].convert_to_units('Myr').value
            age      = t - t_form

            stats_array_dict = {}
            for abundance in aratios.keys():
                stats_array_dict[abundance] = {}
                for k in stats.keys():
                    stats_array_dict[abundance][k] = np.zeros(np.size(tbins) - 1)
            for i in np.arange(np.size(tbins)-1):

                age = tbins[i] - t_form
                selection = (age >= 0.0)*(age <= lifetime)
                for abundance in aratios.keys():
                    if i == 0:
                        sub_g = g.create_group(abundance)

                    if np.size(age[selection]) > 1:
                        stats = utilities.compute_stats(aratios[abundance][selection], return_dict = True) # +1 b/c index starts at 1
                        for k in stats.keys():
                            stats_array_dict[abundance][k][i] = stats[k]
                    else:
                        for k in stats.keys():
                            stats_array_dict[abundance][k][i] = None

            for abundance in aratios.keys():
                g = hf[groupname + '/statistics/cumulative/' + abundance]
                for k in stats_array_dict[abundance].keys():
                    g.create_dataset(k, data = stats_array_dict[abundance][k])

            # now bin by times (using various dt) to get instantaneous median and spread in SF
            # at any given point in time. This is NOT an observational quantity, but rather a theoretical
            # bit of information to understand how much formed stars vary in abundance ratio at any
            # given point in time (i.e. this is the stellar analog to the gas version of these plots)
            for dt in [0.1, 1, 10]:
                g = statgroup.create_group('%iMyr'%(dt))
                t  = ds.current_time.convert_to_units('Myr').value
                tmax = np.around(t, decimals = -len(str(dt)) + 1)
                if tmax < t:
                    tmax = tmax + dt
                tbins = np.arange(0.0, tmax + 0.5*dt, dt)

                index = np.digitize(data['creation_time'].convert_to_units('Myr').value, tbins)
                hist, bins  = np.histogram(data['creation_time'].convert_to_units('Myr').value, bins = tbins)
                g.create_dataset('bins', data = tbins)
                g.create_dataset('hist', data = np.array(hist))

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
                            stats = utilities.compute_stats(aratios[abundance][index == i+1], return_dict = True) # +1 b/c index starts at 1
                            for k in stats.keys():
                                stats_array_dict[abundance][k][i] = stats[k]
                        else:
                            for k in stats.keys():
                                stats_array_dict[abundance][k][i] = None

                for abundance in aratios.keys():
                    # - - - - - Produce a gap-less, interpolated mean to compute the ACF
                    if False: # don't do this anymore
                        first        = np.where( np.logical_not(np.isnan( stats_array_dict[abundance]['mean'] )))[0][0]
                        mean         = stats_array_dict[abundance]['mean'][first:]
                        select       = np.logical_not(np.isnan(mean))
                        clean_mean   = mean[select]
                        tcent        = 0.5 * (tbins[1:] + tbins[:-1])
                        tcent        = tcent[first:]
                        clean_t      = tcent[select]
                        f_interp     = interp1d(clean_t, clean_mean)
                        interp_mean  = mean
                        interp_mean[np.logical_not(select)] = f_interp( tcent[np.logical_not(select)] )
                        stats_array_dict[abundance]['interp_mean'] = interp_mean
                        stats_array_dict[abundance]['acf'] = utilities.acf(interp_mean, nlags = len(tcent))

                    g = hf[groupname + '/statistics/%iMyr/'%(dt) + abundance]
                    for k in stats_array_dict[abundance].keys():
                        g.create_dataset(k, data = stats_array_dict[abundance][k])

            # ------------ can do a correlation across time bins here too --------- 
            # Pick some time t_o, for the ith bin past t_o, do correlation between
            # those two populations of stars
            # x  = np.array([stars in t_o bin] + [stars in t_i bin])
            # corr[i] = np.correlate(x,x, mode = 'full')
            # allow to plot correlation as a function of time.


        else:
            continue
#
#            g.create_dataset('Nstars', data = 0.0)
#            g.create_dataset('Mstars', data = 0.0)
#            sg = hf.create_group(groupname + '/abundances')
#            for abundance in aratios.keys():
#                sg.create_dataset( abundance, data = 0.0)


    hf.close()

    return

if __name__=='__main__':

    ds_list = np.sort(glob.glob("DD????/DD????"))
    ds_list = [ds_list[-1]] # just do most recent file
    dsnames = [x.split('/DD')[0] for x in ds_list]


    gal = Galaxy(dsnames[-1])
    del(gal)

    generate_abundances(ds_list = ds_list, overwrite = True)

    plot_MDF(ds_list = dsnames)
    plot_abundances(ds_list = dsnames, plot_type = 'standard', color_by_age = True, show_average = True)
    plot_time_evolution(ds_list = dsnames, show_quartile=True)


#    plot_acf() - computeation of ACF currently broken (see global shutoff parameter)

