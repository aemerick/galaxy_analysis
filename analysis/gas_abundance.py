import yt
import numpy as np
import matplotlib.pyplot as plt
from collections import Iterable, OrderedDict

import glob
import os
import h5py

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
    data_source = galaxy.halo_spherical_region

    return data_source, mask

def _ISM_region(galaxy, name):
    mask        = None
    data_source = galaxy.disk.cut_region(ISM[name])
    return data_source, mask

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

    # make the masked data set

    # want to compute distributions of mass fractions
    # for each individual species in each cell. In addition to
    # common abundance ratios (do X / Fe, X / H, and X / Mg)

    fbins = np.logspace(-10, 1, 100)

    data_dict = {}
    data_dict['fraction']  = {}
    data_dict['abundance'] = {}
    data_dict['general']    = {}


    data_dict['fraction']['bins'] = fbins

    mask_empty = np.sum(mask)

    for field in fraction_fields:

        if mask_empty or len(data_source['Density']) < 1:
            data_dict['fraction'][field] = [np.zeros(np.size(fbins)),
                                            np.zeros(5)]
        else:
            
            fdata = data_source[field][mask]

            bins, hist = np.histogram(fdata, bins = fbins)

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

    if not os.path.isfile(dir + outfile) or overwrite:
        hf = h5py.File(dir + outfile, 'w')
    else:
        hf = h5py.File(dir + outfile, 'a')

    ds_list = np.sort( glob.glob('./DD????/DD????'))
    times = np.zeros(np.size(ds_list))

    # get the fields 
    ds = yt.load(ds_list[0])
    metals = utilities.species_from_fields(ds.field_list)

    fraction_fields = ['H_p0_fraction','H_p1_fraction','He_p0_fraction',
                       'He_p1_fraction','He_p2_fraction','H2_fraction']
    for m in metals:
        fraction_fields += [m + '_fraction']

    # make the abundance ratio fields

    abundance_fields = None

    for i, dsname in enumerate(ds_list):

        groupname = dsname.rsplit('/')[1]
        gal = Galaxy(groupname)

        if groupname in hf and not overwrite:
            continue

        g = hf.create_group(groupname)
        g.create_dataset('Time'  , data = gal.ds.current_time.convert_to_units('Myr').value)

        gas_data = compute_stats_all_masks(gal, fraction_fields  = fraction_fields,
                                          abundance_fields = abundance_fields)

        print gas_data
        for k in gas_data.keys():
            g.create_dataset(k, data = gas_data[k])

        del(ds)

    hf.close()

    return

def plot_gas_fractions():
    """
    Panel plot showing distributions of gas mass fractions for regions
    in galaxy
    """

    return

if __name__ == '__main__':

    generate_all_stats(overwrite=True)

    plot_gas_fractions()

    plot_abundances(plot_type = 'standard')
