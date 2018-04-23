SKIP_ABUNDANCES = False

import yt
import numpy as np
from galaxy_analysis.plot.plot_styles import *
fsize = 14

import matplotlib.pyplot as plt
from collections import Iterable, OrderedDict

import glob
import os
import h5py
import deepdish as dd

# parallel
from multiprocessing import Pool
from contextlib import closing
import itertools


# --- internal ---
from galaxy_analysis import Galaxy
from galaxy_analysis.utilities import utilities as utilities
from galaxy_analysis.utilities import functions
from galaxy_analysis.static_data import ISM


#
# Brainstorming:
#    - Functions to go through all star forming regions in the
#      galaxy, defined as regions where:
#        1) n > n_thresh
#        2) T < T_thresh
#        3) div(v) < 0
#
#      Make a derived field that is a mask for the above.
#      Need to define velocity field gradient. 
#
#
#    structure dictionaries as:
#
#
#          Time         --- Fractions  ---> (all fractions)
#         /            /
#    DD -- mask_type -- --- Abundances ---> (all abundances)
#                      \
#                       --- General    ---> (n, T, P, G_o, t_dyn dists)
#
#
#   mask_types: -
#                 1) Star forming regions
#                 2) Disk CNM
#                 3) Disk WNM
#                 4) Disk HIM
#                 5) Halo



def _star_forming_region(galaxy):
    """
    """
    mask        = galaxy.disk[('gas','is_star_forming')]
    data_source = galaxy.disk
    return data_source, mask

def _disk(galaxy):
    mask        = np.ones(np.shape(galaxy.disk['x'])) # all cells
    data_source = galaxy.disk
    return data_source, mask


def _CNM(galaxy):
    return _ISM_region(galaxy,'CNM')
def _WIM(galaxy):
    return _ISM_region(galaxy,'WIM')
def _molecular(galaxy):
    return _ISM_region(galaxy,'Molecular')
def _WNM(galaxy):
    return _ISM_region(galaxy,'WNM')
def _HIM(galaxy):
    return _ISM_region(galaxy,"HIM")

def _halo(galaxy):

    mask = galaxy.halo_sphere['spherical_radius'] > galaxy.disk.radius
    data_source = galaxy.halo_sphere

    return data_source, mask

def _ISM_region(galaxy, name):
    mask        = None
    data_source = galaxy.disk.cut_region(ISM[name])
    return data_source, mask


_mask_color = {'star_forming' : 'gold',
               'CNM'          : 'blue',
               'WNM'          : 'green',
               'WIM'          : 'orange',
               'Molecular'    : 'purple',
               'HIM'          : 'red',
               'halo'         : 'black',
               'Disk'         : 'cyan'}

_mask_ls    = {'star_forming' : '-',
               'CNM'          : '-',
               'WNM'          : '-',
               'HIM'          : '-',
               'halo'         : '-',
               'WIM'          : '-',
               'Molecular'    : '-',
               'Disk'         : '-'}


def compute_stats_all_masks(galaxy, fraction_fields = None, 
                                    abundance_fields = None):

    # define the standard masks, then compute things for all of them
    all_masks   = {'star_forming': _star_forming_region,
                   'CNM': _CNM,
                   'WNM': _WNM,
                   'HIM': _HIM,
                   'halo': _halo,
                   'WIM' : _WIM,
                   'Molecular' : _molecular,
                   'Disk'   : _disk}

    data = {}
    for m in all_masks.keys():
        data_source, mask = all_masks[m](galaxy)
        data[m] = compute_abundance_stats(galaxy.ds,
                                          data_source, mask, fraction_fields,
                                          abundance_fields)

    return data

def compute_abundance_stats(ds, data_source, mask = None,
                                fraction_fields = None,
                                abundance_fields = None,
                                mask_abundances = False,
                                unpoluted_threshold = -18):
    """
    Computes the mass and volume weighted distributions of metal
    species mass fractions AND metal species abundance ratios (in
    the form [X/Y]) for all available fields (default behavior).
    Alternatively a limited list of fraction or abundance fields can
    be provided.

    A better way to compute the distributions of the abundance ratios
    is likely to mask out cells that have not differed from the initial
    values, assuming all cells begin at very low  [X/H], as would occur
    when the initial mass fractions are set to a tiny number. This would
    make abundance ratios between two metal species easier to interpret.
    If 'mask_abundances' is True, this masking is done, removing cells
    from the computation if [X/H] is below 'unpoluted_threshold'. This
    is the default behavior, should lead to neary empty distributions
    for early times in the galaxy's evolution.
    """

    # if None, do this for everything
    if mask is None:
        mask = np.ones(np.shape(data_source['x'])) # all cells
        # mask = (mask == 1)

    # make the masked data set

    # want to compute distributions of mass fractions
    # for each individual species in each cell. In addition to
    # common abundance ratios (do X / Fe, X / H, and X / Mg)

    fbins = np.logspace(-20,  0, 401)
    abins = np.linspace(-10,  4, 501)

    rbins = np.arange(0.0, 620.0, 20.0) * yt.units.pc # 20 parsec bins?

    data_dict = {}
    data_dict['volume_fraction']  = {}
    data_dict['mass_fraction']    = {}
    data_dict['radial_profile']   = {}

    data_dict['volume_fraction']['bins'] = fbins
    data_dict['mass_fraction']['bins'] = fbins
    data_dict['volume_fraction']['abins'] = abins
    data_dict['mass_fraction']['abins'] = abins
    data_dict['radial_profile']['rbins'] = rbins

    mask = np.array(mask).astype(bool)

    mask_empty = not any(mask)

    cv = data_source['cell_volume'][mask]
    cm = data_source['cell_mass'][mask]
    total_volume = np.sum(cv.value) * cv.unit_quantity * 1.0     # total volume of masked cells
    total_mass   = np.sum(cm.value) * cm.unit_quantity * 1.0     # total mass   of masked cells

    all_fields = None
    if not (fraction_fields is None):
        all_fields = fraction_fields

    if not (abundance_fields is None):
        all_fields = all_fields + abundance_fields

    for field in all_fields:

        if 'over' in field:
            bins = abins

            if SKIP_ABUNDANCES:
                continue # WARNING FOR DEBUGGING PURPOSES

        else:
            bins = fbins
        centers = 0.5 * (bins[1:] + bins[:-1])

        #  for the abundance ratio fields, we want to get rid
        #  of the things that have 'primordial' abundances
        #  to make interpretation cleaner -
        if mask_abundances and 'over' in field:
            ele            = field.split('_over_')[0]
#            over_H         = ele + '_over_H'
            test_field     = ('gas','O_Fraction')
            unpoluted_mask = data_source[test_field] > unpoluted_threshold

            mask = mask * unpoluted_threshold
            mask = np.array(mask).astype(bool)
            mask_empty = not any(mask)

        # if there is no data, fill the arrays with zeros
        #print field, mask_empty
        #print len(data_source[field][mask])

        if mask_empty or (len(data_source[field][mask]) < 3):
            data_dict['volume_fraction'][field] = {'hist':np.zeros(np.size(bins)-1)}
            data_dict['mass_fraction'][field]   = {'hist':np.zeros(np.size(bins)-1)}
            stats = utilities.compute_weighted_stats(np.zeros(10), np.ones(10), return_dict=True) # arbitrary - just trying to get keys
            data_dict['radial_profile'][field] = {}

            for k in stats:
                data_dict['volume_fraction'][field][k] = None
                data_dict['mass_fraction'][field][k]   = None
                data_dict['radial_profile'][field][k]  = [None] * (np.size(rbins) - 1)

            data_dict['mass_fraction'][field]['mode'] = None
            data_dict['volume_fraction'][field]['mode'] = None
            data_dict['radial_profile'][field]['mode'] = [None] * (np.size(rbins) - 1)

        else:

            fdata = data_source[field][mask]
            r_cyl = data_source['cylindrical_radius'][mask]

            # compute the histograms of the data
            mass_hist, temp = np.histogram(fdata.value, weights = cm.value, bins = bins) / total_mass.value
            vol_hist, temp  = np.histogram(fdata.value, weights = cv.value, bins = bins) / total_volume.value

            # now compute descriptive statistics for each weighting
            stats = utilities.compute_weighted_stats(fdata, cv, return_dict = True)
            data_dict['volume_fraction'][field] = {'hist': vol_hist}
            data_dict['volume_fraction'][field]['mode'] = centers[np.argmax(vol_hist)]

            stats2 = utilities.compute_weighted_stats(fdata, cm, return_dict = True)
            data_dict['mass_fraction'][field]   = {'hist': mass_hist}
            data_dict['mass_fraction'][field]['mode'] = centers[np.argmax(mass_hist)]

            # save these into the dictionary
            data_dict['radial_profile'][field] = {} # and set up rprof
            for k in stats:
                data_dict['volume_fraction'][field][k] = stats[k]
                data_dict['mass_fraction'][field][k]   = stats2[k]

                # set up radial profile bins
                data_dict['radial_profile'][field][k]  = np.zeros(np.size(rbins)-1)
            data_dict['radial_profile'][field]['mode'] = np.zeros(np.size(rbins)-1)

            # now compute the radial profile - mass weighted ONLY
            for i in np.arange(np.size(rbins)-1):
                selection = (r_cyl >= rbins[i]) * (r_cyl < rbins[i+1])
                x         = fdata[selection] # sub-select field
                w         = cm[selection]    #

                if np.size(x) >= 3:
                    hist, temp = np.histogram(x, weights = w, bins = bins)

                    stats     = utilities.compute_weighted_stats(x, w, return_dict = True)
                    for k in stats:
                        data_dict['radial_profile'][field][k][i] = stats[k]
                    data_dict['radial_profile'][field]['mode'][i] = centers[np.argmax(hist)]

                else:
                    for k in stats:
                        data_dict['radial_profile'][field][k][i] = None
                    data_dict['radial_profile'][field]['mode'][i] = None

    #
    # general properties
    #
    data_dict['general'] = {'total_volume' : total_volume.convert_to_units('pc**3'),
                            'total_mass'   : total_mass.convert_to_units('Msun')}

    # data_dict['general']['number_density'] = masked_data


    return data_dict

def _parallel_loop(dsname, fraction_fields):

    groupname = dsname.rsplit('/')[1]
    print "starting computation on ", groupname
    gal = Galaxy(groupname)

#
#
#
#
    # limit ourselves to Fe and H demoninators for now to speed up computation
    abundance_fields = utilities.abundance_ratios_from_fields(gal.ds.derived_field_list,
                                                              select_denom = ['Fe','H','Ba','Mg','O','N'])
    species = utilities.species_from_fields(gal.ds.field_list,include_primordial=True)
    metal_species = utilities.species_from_fields(gal.ds.field_list)

#
#    don't recompute things unless overwrite is set
#
#    if (not overwrite) and (existing_keys is None):
#        print "Cannot check for overwrites without knowing existing data"
#        raise ValueError
#
#
#        if (groupname in existing_keys) and (not overwrite):
#           return {groupname + '_skip_this_one' : {}} # g['skip_this_one_' + dsname]

    dictionary = {groupname : {}}
    #hf[groupname] = {} # make an empty key for this group
    g = dictionary[groupname]

    g['general'] = {}
    g['general']['Time']    = gal.ds.current_time.convert_to_units('Myr').value
    # generalized function to loop through all mask types and compute stats
    gas_data  = compute_stats_all_masks(gal, fraction_fields  = fraction_fields,
                                      abundance_fields = abundance_fields)

    # make the stats group and the histogram group
    for k in gas_data.keys():
        g[k] = gas_data[k]

    del(gal)

    print "ending computation on ", groupname

    return dictionary

def _parallel_loop_star(combined):
    return _parallel_loop(*combined)

def generate_all_stats(outfile = 'gas_abundances.h5',
                        dir = './abundances/', overwrite=False, nproc = 1):
    """
    For all data files, generate gas abundance statistics for all
    element fractions and abundance ratios (as defined below). This is
    an expensive operation.
    """

    if not os.path.exists(dir):
        os.makedirs(dir)

    hdf5_filename = dir + outfile

    if not os.path.isfile(hdf5_filename) or overwrite:
        hf = h5py.File(hdf5_filename, 'w')
        hf.close()

    hf = dd.io.load(hdf5_filename)

    ds_list = np.sort( glob.glob('./DD???0/DD???0') +\
                       glob.glob('./DD???5/DD???5'))

    print "WARNING: Only doing limited number of outputs for ease of use"

    for i, dsname in enumerate(ds_list):
        ds = yt.load(dsname)
        if ds.parameters['NumberOfParticles'] > 0:
            start_index = i
            del(ds)
            break
        del(ds)

    times = np.zeros(np.size(ds_list))
    ds_list = ds_list[start_index:]
    times   = times[start_index:]

    # get the fields
    ds = yt.load(ds_list[0])
    metals = utilities.species_from_fields(ds.field_list)

    # additional fraction fields for the non-equillibrium chemistry species
    fraction_fields = ['H_p0_fraction','H_p1_fraction','He_p0_fraction',
                       'He_p1_fraction','He_p2_fraction','H2_fraction']
    for m in metals:
        fraction_fields += [m + '_Fraction']

    # make the abundance ratio fields
    #   - should do everything over H
    #   - should do everything over Fe
    #   - should do everything over Mg

    # loop through all data files
    if nproc == 1:
        for i, dsname in enumerate(ds_list):
            print i, dsname
            groupname = dsname.rsplit('/')[1]
            gal = Galaxy(groupname)

            # if first loop, define the fields
            if i == 0:
                # limit ourselves to Fe and H demoninators for now to speed up computation
                abundance_fields = utilities.abundance_ratios_from_fields(gal.ds.derived_field_list,
                                                                          select_denom = ['Fe','H','O','Mg','N'])
                species = utilities.species_from_fields(gal.ds.field_list,include_primordial=True)
                metal_species = utilities.species_from_fields(gal.ds.field_list)

            # don't recompute things unless overwrite is set
            if groupname in hf.keys() and (not overwrite):
                continue

            hf[groupname] = {} # make an empty key for this group
            g = hf[groupname]

            g['general'] = {}
            g['general']['Time']    = gal.ds.current_time.convert_to_units('Myr').value

            # generalized function to loop through all mask types and compute stats
            gas_data  = compute_stats_all_masks(gal, fraction_fields  = fraction_fields,
                                              abundance_fields = abundance_fields)

            # make the stats group and the histogram group
            for k in gas_data.keys():
                g[k] = gas_data[k]

            del(gal)
    else: # parallel

        # select out data sets that already exist in output
        if not overwrite:
            ds_list = [x for x in ds_list if ( not any( [x.rsplit('/')[1] in y for y in hf.keys() ]))]

        # construct the pool, and map the results to a holder
        #   pool splits computation among processors

        #
        # do this in a loop, so we can save progressively once a full set of processors is complete
        #   saves you if there is a crash (kinda). This way we do better memory management if
        #   operating on many large datasets.
        #
        for sub_list in itertools.izip_longest(*(iter(ds_list),) * nproc):

            sub_list = list(sub_list)
            sub_list = [s for s in sub_list if s is not None] # remove None values
            reduced_nproc = np.min( [len(sub_list), nproc] )  # only run on needed processors

            pool = Pool(reduced_nproc)
            results = pool.map_async(_parallel_loop_star,
                                      itertools.izip(sub_list, itertools.repeat(fraction_fields)))
            pool.close() # no more processes
            pool.join()  # wait and join running processes

            # gather results and add to output
            for r in results.get():
                hf[r.keys()[0]] = r[r.keys()[0]]
            del(results)


        # define these fields (which are defined in the Pool funtion but don't want to
        # deal with passing this back) ------
        if len(ds_list) > 0:
            gal = Galaxy(ds_list[0].split('/')[1])
            abundance_fields = utilities.abundance_ratios_from_fields(gal.ds.derived_field_list,
                                                                      select_denom = ['Fe','H','O','Mg','N'])
            species = utilities.species_from_fields(gal.ds.field_list,include_primordial=True)
            metal_species = utilities.species_from_fields(gal.ds.field_list)
            del(gal)

    # save field names
    if not ('species' in hf.keys()):
        hf['species']          = species
    if not ('metal_species' in hf.keys()):
        hf['metal_species']    = metal_species
    if not ('abundance_fields' in hf.keys()):
        hf['abundance_fields'] = abundance_fields

    dd.io.save(hdf5_filename, hf)

    return

def plot_gas_fractions(dir = './abundances/', fname = 'gas_abundances.h5', overwrite = True,
                       fraction_type = 'volume'):
    """
    Panel plot showing distributions of gas mass fractions for regions
    in galaxy
    """

    # for each dataset, plot the distributions of gas fractions in each phase

    all_data = dd.io.load(dir + fname)

    plot_fields = all_data['metal_species']
    nplots      = len(plot_fields)
    nrow, ncol  = utilities.rowcoldict[nplots]

    ds_loop = np.sort([x for x in all_data.keys() if 'DD' in x])
    ds_loop = [ds_loop[-1]]

    for dsname in ds_loop:
        if (len(glob.glob(dsname + '*_fraction*.png')) > 0) and not overwrite:
            continue # do not remake plots

        data  = all_data[dsname]

        fig, ax = plt.subplots(nrow, ncol)
        fig.set_size_inches(nrow*4,ncol*4)
        all_phases  = ['CNM','WNM','HIM','star_forming','halo']

        for j,phase in enumerate(all_phases):
            print j, phase
            phase_data = data[phase][fraction_type + '_fraction']
            logbins = np.log10(phase_data['bins'])

            axi, axj = 0,0
            for field in plot_fields:
                axind = (axi,axj)
                ax[axind].plot(logbins[:-1], np.log10(phase_data[field + '_Fraction']['hist']), drawstyle='steps-post', lw = 3,
                          color = _mask_color[phase], ls = _mask_ls[phase], label=phase)

                ax[axind].set_ylim(-5, 0)
                ax[axind].set_xlim(-15, -2)
                ax[axind].set_ylabel('log [' + fraction_type + ' Fraction]')
                ax[axind].set_xlabel('log [' + field + ' Fraction By Number)')

                axj = axj + 1
                if axj >= ncol:
                    axj = 0
                    axi = axi + 1

        ax[(0,0)].legend(loc='best')

        plt.tight_layout()
        plt.minorticks_on()
        fig.savefig(dsname + '_metal_' + fraction_type +'_fractions.png')

        plt.close()

    return


def plot_abundances(plot_type = 'standard', dir = './abundances/', fname = 'gas_abundances.h5', overwrite = True,
                    fraction_type = 'volume'):
    """
    Panel plot showing distributions of gas mass fractions for regions
    in galaxy
    """

    # for each dataset, plot the distributions of gas fractions in each phase

    all_data   = dd.io.load(dir + fname)
    all_fields = all_data['abundance_fields']

    if plot_type == 'standard' or plot_type == 'Fe':
        all_fields  = all_data['abundance_fields']

        # plot all over Fe and Fe over H
        plot_fields = [x for x in all_fields if '_over_Fe' in x]
        plot_fields = plot_fields + ['Fe_over_H']

    elif len(plot_type) <= 2:
        # assume it is a species name
        plot_fields = [x for x in all_fields if ('_over_' + plot_type) in x]
        if (plot_type + '_over_H') in all_fields:
            plot_fields = plot_fields + [plot_type + '_over_H']

    nplots     = len(plot_fields)
    nrow, ncol = utilities.rowcoldict[nplots]

    # only do most recent
    ds_loop = np.sort([x for x in all_data.keys() if 'DD' in x])
    ds_loop = [ds_loop[-1]]

    for dsname in ds_loop:
        if (len(glob.glob(dsname + '*_abundances*.png')) > 0) and not overwrite:
            continue # do not remake plots

        data  = all_data[dsname]

        # hard code for now
        fig, ax = plt.subplots(nrow, ncol)
        fig.set_size_inches(4*nrow, 4*ncol)
        all_phases  = ['CNM','WNM','HIM','star_forming','halo']

        for j,phase in enumerate(all_phases):
#            print j, phase
            phase_data = data[phase][fraction_type + '_fraction'] # mass or volume fraction
            logbins =  phase_data['abins']

            axi, axj = 0, 0
            for field in plot_fields:
                axind = (axi,axj)

                ax[axind].plot(logbins[:-1], np.log10(phase_data[field]['hist']), drawstyle='steps-post', lw = 3,
                          color = _mask_color[phase], ls = _mask_ls[phase], label=phase)

                ax[axind].set_ylim(-4, 0)

                if '_over_H' in field:
                    ax[axind].set_xlim(-8,1)
                else:
                    ax[axind].set_xlim(-3, 3)

                ax[axind].set_ylabel('log [' + fraction_type + ' Fraction]')
                ax[axind].set_xlabel('[' + field + '] Abundance')

                axj = axj + 1
                if axj >= ncol:
                    axj = 0
                    axi = axi + 1


        ax[(0,0)].legend(loc='best')

        plt.tight_layout()
        plt.minorticks_on()
        fig.savefig(dsname + '_' + fraction_type +'_abundances.png')

        plt.close()

    return

def plot_time_evolution(dir = './abundances/', fname = 'gas_abundances.h5',
                        abundance = False, show_std = False, show_quartile = False,
                        fraction_type = 'mass', plot_type = 'Fe'):

    # for each dataset, plot the distributions of gas fractions in each phase

    all_data   = dd.io.load(dir + fname)
    all_fields = all_data['abundance_fields']

    share_axis = False
    if abundance:

        if plot_type == 'standard':
            all_fields  = all_data['abundance_fields']

            # plot all over Fe and Fe over H
            plot_fields = utilities.sort_by_anum([x for x in all_fields if '_over_Fe' in x])
            plot_fields = plot_fields + ['Fe_over_H']

        elif len(plot_type) <= 2:
            # assume it is a species name
            plot_fields = utilities.sort_by_anum([x for x in all_fields if ('_over_' + plot_type) in x])


#           No longer doing this - plot separately with everything else to allow for share axis
#            if (plot_type + '_over_H') in all_fields:
#                plot_fields = plot_fields + [plot_type + '_over_H']

            share_axis = True
    else:
        plot_fields = utilities.sort_by_anum(all_data['metal_species'])
        plot_fields = [x + '_Fraction' for x in plot_fields]


    nplots     = len(plot_fields)
    nrow, ncol = utilities.rowcoldict[nplots]

    # all_ds represents list of dataset names and kwargs to pull from
    all_ds = np.sort([x for x in all_data.keys() if 'DD' in x])

    # gather times for all data sets
    times = utilities.extract_nested_dict_asarray(all_data, ['general','Time'], loop_keys = all_ds)
    # normalize time evolution
    times = times - times[0]

    # set up plot - define the phase names we want to loop over
    if share_axis:
        fig, ax = plt.subplots(nrow,ncol, sharex=True, sharey=True)
    else:
        fig, ax = plt.subplots(nrow, ncol)

    fig.set_size_inches(4*nrow, 4*ncol)
    fig.subplots_adjust(hspace=0.0, wspace = 0.0)

    all_phases  = ['CNM','WNM','HIM','star_forming','halo'] # these are all that are available

    axi, axj = 0, 0
    xmin, xmax = 1.0E99, -1.0E99
    for field in plot_fields:
        axind = (axi,axj)

        # plot all phases for a fixed plot - add std or quartile if desired
        ymax = -1.0E99
        for j,phase in enumerate(all_phases):
            key_list = [phase, fraction_type + '_fraction', field]

            avg = utilities.extract_nested_dict_asarray(all_data, key_list + ['mean'], loop_keys = all_ds)
            ax[axind].plot(times, avg, color = _mask_color[phase], ls = _mask_ls[phase], label=phase, lw = line_width)
            ymax = np.max( [ymax, np.max(avg)] )

            if show_std:
                std = utilities.extract_nested_dict_asarray(all_data, key_list + ['std'], loop_keys = all_ds)
                ax[axind].fill_between(times, avg - std, avg + std, color = _mask_color[phase], alpha = 0.5, lw = 2.0)

            if show_quartile:
                Q1  = utilities.extract_nested_dict_asarray(all_data, key_list + ['Q1'], loop_keys = all_ds)
                Q3  = utilities.extract_nested_dict_asarray(all_data, key_list + ['Q3'], loop_keys = all_ds)
                # filter out "None" that can appear
                select = (Q1>-np.inf)*(Q3>-np.inf)*(times>-np.inf)
#                print times[select]
#                print Q1[select]
#                print Q3[select]
#                print _mask_color[phase], phase
                xmin = np.min( [xmin, np.min(times[select])])
                xmax = np.max( [xmax, np.max(times[select])])
                ax[axind].fill_between(times[select], np.array(Q1[select]).astype(float),
                                           np.array(Q3[select]).astype(float), color = _mask_color[phase], alpha = 0.5, lw = 2.0)

        # set axis limits and labels based on whether or not we are plotting
        # abundance ratios or number fractions
        if abundance:
            label = field.split('_over_')
            if share_axis:
                xy = (0.125, 0.8)
                ax[axind].annotate(label[0], xy = xy, xytext=xy,
                             xycoords = 'axes fraction', textcoords = 'axes fraction')

                label = r'log([X/' + label[1] + '])'
            else:
                label = r'log([' + label[0] +'/' + label[1] + '])'

            # set axis limits and labels according to shared restrictions
            if (share_axis and axj == 0) or (not share_axis):
                ax[axind].set_xlim(xmin, xmax)
                ax[axind].set_ylabel(label)

            if (share_axis and axi == nrow) or (not share_axis):
                ax[axind].set_xlabel(r'Time (Myr)')

            if (share_axis and axi == 0) or (not share_axis):
                dx = 0.0 # some slop to get axis labels to not clash
                if share_axis:
                    dx = 0.15
                if '_over_H' in field:
                    ax[axind].set_ylim(-8 - dx, 0 + dx)
                else:
                    ax[axind].set_ylim(-3 - dx, 3 + dx)

        else:

            if share_axis:
                xy = (0.125, 0.8)
                ax[axind].annotate(field.split('_Frac')[0], xy = xy,
                                      xytext=xy, xycoords = 'axes fraction',
                                      textcoords = 'axes fraction')

                label = r'log([X Fraction])'
            else:
                label = r'log([' + field + '])'

            if (share_axis and axj ==0) or (not share_axis):
                ax[axind].set_xlim(xmin, xmax)
                ax[axind].set_ylabel(label)

            if (share_axis and axi ==nrow) or (not share_axis):
                ax[axind].set_xlabel(r'Time (Myr)')

            if (share_axis and axi == 0) or (not share_axis):
                ax[axind].set_ylim(1.0E-12, 3.0E-2)
                ax[axind].semilog_y()


        ax[axind].minorticks_on()
        axj = axj + 1
        if axj >= ncol:
            axj = 0
            axi = axi + 1

    ax[(0,0)].legend(loc='best', ncol=2)

    if not share_axis:
        plt.tight_layout()

    plt.minorticks_on()
    if abundance:
        outname = 'time_evolution_' + fraction_type + '_abundances_over_' + plot_type + '.png'
    else:
        outname = 'time_evolution_' + fraction_type + '_fractions.png'

    if show_quartile:
        outname = 'quartile_' + outname

    fig.savefig(outname)

    plt.close()

    return

def collate_to_time_array(filepath = None):
    if filepath is None:
        filepath = './gas_abundances.h5'
    data = dd.io.load(filepath)
    # make a new array to hold times
    if (not ('time_evolution') in data.keys()):
        data['time_evolution'] = {}
    return

if __name__ == '__main__':

    generate_all_stats(overwrite=False, nproc = 28)

#    plot_time_evolution(abundance=True, plot_type = 'Fe')
#    plot_time_evolution(abundance=True, plot_type = 'H')
#    plot_time_evolution(abundance=True, plot_type = 'Fe', show_quartile = True)
#    plot_time_evolution(abundance=True, plot_type = 'H',  show_quartile = True)


#    plot_gas_fractions(overwrite=True, fraction_type = 'volume')
#    plot_gas_fractions(overwrite=True, fraction_type = 'mass')

#    plot_abundances(plot_type = 'standard', fraction_type = 'mass')
