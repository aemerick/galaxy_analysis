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

    all_masks   = {'star_forming': _star_forming_region,
                   'CNM': _CNM,
                   'WNM': _WNM,
                   'HIM': _HIM,
                   'halo': _halo}

    data = {}
    for m in all_masks.keys():
        data_source, mask = all_masks[m](galaxy)
        print m, mask
        data[m] = compute_abundance_stats(galaxy.ds,
                                          data_source, mask, fraction_fields,
                                          abundance_fields)

    return data

def compute_abundance_stats(ds, data_source, mask = None,
                                fraction_fields = None, abundance_fields = None):

    # if None, do this for everything
    if mask is None:
        mask = np.ones(np.shape(data_source['Density']))
        mask = (mask == 1)

    # make the masked data set

    # want to compute distributions of mass fractions
    # for each individual species in each cell. In addition to
    # common abundance ratios (do X / Fe, X / H, and X / Mg)

    fbins = np.logspace(-20, 0, 401)

    data_dict = {}
    data_dict['fraction']  = {}
    data_dict['abundance'] = {}

    data_dict['fraction']['bins'] = fbins

    mask = np.array(mask).astype(bool)

    mask_empty = not any(mask)

    cv = data_source['cell_volume'][mask]
    total_volume = np.sum(cv) * 1.0
    for field in fraction_fields:

        if mask_empty or len(data_source['Density']) < 1:
            data_dict['fraction'][field] = [np.zeros(np.size(fbins)-1),
                                            np.zeros(5)]
        else:

            fdata = data_source[field][mask]

            # compute volume fraction
            hist = np.zeros(np.size(fbins) -1)

            for i in np.arange(np.size(fbins) - 1):
                hist[i] = np.sum( cv[ (fdata < fbins[i+1]) * (fdata >= fbins[i]) ] ) / total_volume
#            hist, bins = np.histogram(fdata, bins = fbins)

            stats = utilities.compute_stats(fdata)

            data_dict['fraction'][field] = [hist, stats]

    #
    # general properties
    #

    # data_dict['general']['number_density'] = masked_data


    return data_dict

def generate_all_stats(outfile = 'gas_abundances.h5',
                        dir = './abundances/', overwrite=False):

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

    fraction_fields = ['H_p0_fraction','H_p1_fraction','He_p0_fraction',
                       'He_p1_fraction','He_p2_fraction','H2_fraction']
    for m in metals:
        fraction_fields += [m + '_Fraction']

    # make the abundance ratio fields

    abundance_fields = None

    for i, dsname in enumerate(ds_list):
        print i, dsname
        groupname = dsname.rsplit('/')[1]
        gal = Galaxy(groupname)

        if groupname in hf.keys() and (not overwrite):
            continue

        hf[groupname] = {} # make an empty key for this group
        g = hf[groupname]

#        g = hf.create_group(groupname)
        g['general'] = {}
        g['general']['Time'] = gal.ds.current_time.convert_to_units('Myr').value

        gas_data  = compute_stats_all_masks(gal, fraction_fields  = fraction_fields,
                                          abundance_fields = abundance_fields)


        # make the stats group and the histogram group
        for k in gas_data.keys():
            g[k] = gas_data[k]

        del(gal)

    dd.io.save(hdf5_filename, hf)
#    hf.close()

    return

def plot_gas_fractions(dir = './abundances/', fname = 'gas_abundances.h5', overwrite = False):
    """
    Panel plot showing distributions of gas mass fractions for regions
    in galaxy
    """

    # for each dataset, plot the distributions of gas fractions in each phase

    all_data = dd.io.load(dir + fname)

    for dsname in all_data.keys():
        if (len(glob.glob(dsname + '*_fraction*.png')) > 0) and not overwrite:
            continue # do not remake plots

        data  = all_data[dsname]

        #metal_fractions = data[data.keys()[0]]['fraction'].keys()
        #metal_fractions = [x for x in metal_fractions if ( (not 'H' in x) and (not x == 'bins'))]
        #n_metal         = len(metal_fractions)

        # hard code for now
        fig, ax = plt.subplots(2,5)
        all_phases  = ['CNM','WNM','HIM','star_forming','halo']
        plot_fields = ['C','N','O','Mg','Si','S','Ca','Mn','Fe','Ni']
        index       = [(0,0),(0,1),(0,2),(0,3),(0,4),(1,0),(1,1),(1,2),(1,3),(1,4)]

        for j,phase in enumerate(all_phases):
            print j, phase
            phase_data = data[phase]['fraction']
            logbins = np.log10(phase_data['bins'])
            for i, field in enumerate(plot_fields):
                axi = index[i]
                ax[axi].plot(logbins[:-1], np.log10(phase_data[field + '_Fraction'][0]), drawstyle='steps-post', lw = 3,
                          color = _mask_color[phase], ls = _mask_ls[phase], label=phase)

                ax[axi].set_ylim(-5, 0)
                ax[axi].set_xlim(-15, -2)
                ax[axi].set_ylabel('log [Volume Fraction]')
                ax[axi].set_xlabel('log [' + field + ' Fraction By Number)')

        ax[(0,0)].legend(loc='best')


        fig.set_size_inches(40,16)
        plt.tight_layout()
        plt.minorticks_on()
        fig.savefig(dsname + '_metal_fractions.png')

        plt.close()

    return

if __name__ == '__main__':

    generate_all_stats(overwrite=False)

    plot_gas_fractions(overwrite=True)

    plot_abundances(plot_type = 'standard')
