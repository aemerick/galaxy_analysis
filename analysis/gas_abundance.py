import yt
import numpy as np
from matplotlib import rc

fsize = 17
rc('text', usetex=False)
rc('font', size=fsize)#, ftype=42)
line_width = 3
point_size = 30

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
from collections import Iterable, OrderedDict

import glob
import os
import h5py
import deepdish as dd

# --- internal ---
from galaxy_analysis import Galaxy
from galaxy_analysis.utilities import utilities as utilities
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
#   mask_types: - 1) Star forming regions
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


def _CNM(galaxy):
    return _ISM_region(galaxy,'CNM')
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
               'HIM'          : 'red',
               'halo'         : 'black'}

_mask_ls    = {'star_forming' : '-',
               'CNM'          : '-',
               'WNM'          : '-',
               'HIM'          : '-',
               'halo'         : '-'}


def compute_stats_all_masks(galaxy, fraction_fields = None, 
                                    abundance_fields = None):

    # define the standard masks, then compute things for all of them
    all_masks   = {'star_forming': _star_forming_region,
                   'CNM': _CNM,
                   'WNM': _WNM,
                   'HIM': _HIM,
                   'halo': _halo}

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
                                mask_abundances = True,
                                unpoluted_threshold = -10):
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
        mask = np.ones(np.shape(data_source['Density']))
        mask = (mask == 1)

    # make the masked data set

    # want to compute distributions of mass fractions
    # for each individual species in each cell. In addition to
    # common abundance ratios (do X / Fe, X / H, and X / Mg)

    fbins = np.logspace(-20, 0, 401)
    abins = np.linspace(-10,10, 401)

    data_dict = {}
    data_dict['volume_fraction']  = {}
    data_dict['mass_fraction'] = {}
    data_dict['abundance'] = {}

    data_dict['volume_fraction']['bins'] = fbins
    data_dict['mass_fraction']['bins'] = fbins
    data_dict['volume_fraction']['abins'] = abins
    data_dict['mass_fraction']['abins'] = abins

    mask = np.array(mask).astype(bool)

    mask_empty = not any(mask)

    cv = data_source['cell_volume'][mask]
    cm = data_source['cell_mass'][mask]
    total_volume = np.sum(cv) * 1.0     # total volume of masked cells
    total_mass   = np.sum(cm) * 1.0     # total mass   of masked cells

    all_fields = None
    if not (fraction_fields is None):
        all_fields = fraction_fields

    if not (abundance_fields is None):
        all_fields = all_fields + abundance_fields

    for field in all_fields:

        if 'over' in field:
            bins = abins
        else:
            bins = fbins

        #  for the abundance ratio fields, we want to get rid
        #  of the things that have 'primordial' abundances
        #  to make interpretation cleaner
        if mask_abundances and 'over' in field:
            ele            = field.split('_over_')[0]
            over_H         = ele + '_over_H'
            unpoluted_mask = data_source[over_H] > unpoluted_threshold

            mask = mask * unpoluted_threshold
            mask = np.array(mask).astype(bool)
            mask_empty = not any(mask)

        # if there is no data, fill the arrays with zeros
        if mask_empty or len(data_source['Density']) < 1:
            data_dict['volume_fraction'][field] = {'hist':np.zeros(np.size(bins)-1)}
            data_dict['volume_fraction'][field].update(utilities.compute_stats(np.zeros(np.size(bins)-1)), return_dict=True)
            data_dict['mass_fraction'][field] = {'hist':np.zeros(np.size(bins)-1)}
            data_dict['mass_fraction'][field].update(utilities.compute_stats(np.zeros(np.size(bins)-1)), return_dict=True)
        else:

            fdata = data_source[field][mask]

            # compute volume fraction
            hist  = np.zeros(np.size(bins) -1)
            hist2 = np.zeros(np.size(bins) -1)


            for i in np.arange(np.size(bins) - 1):
                hist[i]  = np.sum( cv[ (fdata < bins[i+1]) * (fdata >= bins[i]) ] ) / total_volume
                hist2[i] = np.sum( cm[ (fdata < bins[i+1]) * (fdata >= bins[i]) ] ) / total_mass

            stats = utilities.compute_stats(hist, return_dict = True)
            data_dict['volume_fraction'][field] = {'hist':hist}
            data_dict['volume_fraction'][field].update(stats)

            stats2 = utilities.compute_stats(hist2, return_dict = True)
            data_dict['mass_fraction'][field]   = {'hist':hist}
            data_dict['mass_fraction'][field].update(stats)

    #
    # general properties
    #

    # data_dict['general']['number_density'] = masked_data


    return data_dict

def generate_all_stats(outfile = 'gas_abundances.h5',
                        dir = './abundances/', overwrite=False):
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

    ds_list = np.sort( glob.glob('./DD????/DD????'))
    times = np.zeros(np.size(ds_list))

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
    for i, dsname in enumerate(ds_list):
        print i, dsname
        groupname = dsname.rsplit('/')[1]
        gal = Galaxy(groupname)

        # if first loop, define the fields
        if i == 0:
            abundance_fields = utilities.abundance_ratios_from_fields(gal.ds.derived_field_list)
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

    # save field names
    hf['species']          = species
    hf['metal_species']    = metal_species
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

    for dsname in [x for x in all_data.keys() if 'DD' in x]:
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

    for dsname in [x for x in all_data.keys() if 'DD' in x]:
        if (len(glob.glob(dsname + '*_abundances*.png')) > 0) and not overwrite:
            continue # do not remake plots

        data  = all_data[dsname]

        # hard code for now
        fig, ax = plt.subplots(nrow, ncol)
        fig.set_size_inches(4*nrow, 4*ncol)
        all_phases  = ['CNM','WNM','HIM','star_forming','halo']

        for j,phase in enumerate(all_phases):
            print j, phase
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

def plot_time_evolution(filepath = None, abundance = False, show_std = False, show_quartile = False):

    # for each dataset, plot the distributions of gas fractions in each phase

    all_data   = dd.io.load(dir + fname)
    all_fields = all_data['abundance_fields']

    if abundance:
        plot_type = 'standard'

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
    else:
        plot_fields = all_data['metal_species']


    nplots     = len(plot_fields)
    nrow, ncol = utilities.rowcoldict[nplots]

    all_ds = [x for x in all_data.keys() if 'DD' in x]
    times = utilities.extract_nested_dict_asarray(all_data, ['Time'], loop_keys = all_ds)

    fig, ax = plt.subplots(nrow, ncol)
    fig.set_size_inches(4*nrow, 4*ncol)
    all_phases  = ['CNM','WNM','HIM','star_forming','halo']

    axi, axj = 0, 0
    for field in plot_fields:
        axind = (axi,axj)

        for j,phase in enumerate(all_phases):
            key_list = [phase, fraction_type + '_fraction', field]
            avg = utilities.extract_nested_dict_asarray(all_data, key_list + ['avg'], loop_keys = all_ds)
            ax[axind].plot(times, avg, color = _mask_color[phase], ls = _mask_ls[phase], label=phase)


            if show_std:
                std = utilities.extract_nested_dict_asarray(all_data, key_list + ['std'], loop_keys = all_ds)
                ax[axind].fill_between(times, avg - std, avg + std, color = _mask_color[phase], alpha = 0.5, lw = 2.0)

            if show_quartile:
                Q1  = utilities.extract_nested_dict_asarray(all_data, key_list + ['Q1'], loop_keys = all_ds)
                Q3  = utilities.extract_nested_dict_asarray(all_data, key_list + ['Q3'], loop_keys = all_ds)

                ax[axind].fill_between(times, Q1, Q3, color = _mask_color[phase], alpha = 0.5, lw = 2.0)

 
        if abundance:
            ax[axind].set_ylim(-4, 0)

            if '_over_H' in field:
                ax[axind].set_xlim(-8,1)
            else:
                ax[axind].set_xlim(-3, 3)

            ax[axind].set_ylabel('log [' + fraction_type + ' Fraction]')
            ax[axind].set_xlabel('[' + field + '] Abundance')
        else:
            ax[axind].set_ylim(-5, 0)
            ax[axind].set_xlim(-15, -2)
            ax[axind].set_ylabel('log [' + fraction_type + ' Fraction]')
            ax[axind].set_xlabel('log [' + field + ' Fraction By Number)')

        axi = axi + 1
        if axj >= ncol:
            axj = 0
            axi = axi + 1

    ax[(0,0)].legend(loc='best')

    plt.tight_layout()
    plt.minorticks_on()
    if abundances:
        outname = 'time_evolution_' + fraction_type + '_abundances.png'
    else:
        outname = 'time_evolution_' + fraction_type + '_fractions.png'

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

    generate_all_stats(overwrite=False)
    plot_time_evolution(abundance=False)


    plot_gas_fractions(overwrite=True, fraction_type = 'volume')
    plot_gas_fractions(overwrite=True, fraction_type = 'mass')

    plot_abundances(plot_type = 'standard', fraction_type = 'mass')
