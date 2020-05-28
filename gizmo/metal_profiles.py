import numpy as np
import deepdish as dd
import sys,os, glob

# PFH python:
#from pfh_python import colors_sps.colors_table as sps
import gizmo_analysis as gizmo
import utilities as gizmo_ut
from utilities.basic.binning import BinClass

FIREDIR = "/scratch/projects/xsede/GalaxiesOnFIRE/metal_diffusion"

def compute_profile_dict(part, part0, fields = None):
    if fields is None:
        fields = ['mass','mass.metals',
                  'mass.si','mass.o','mass.c','mass.fe']

    # profiling utilitiy
    SpeciesProfile = gizmo_ut.particle.SpeciesProfileClass(
                     scaling = 'linear', limits=[r_min,r_max], width=dr,
                     dimension_number=3)

    profile = {}
    norm_dict = {}
    initial_dict = {}
    for field in fields:
        profile[field] = SpeciesProfile.get_profiles(part, species=['gas','star'], property_name=field,
                                                        property_statistic='sum',
                                                        # property_select = select_dict,
                                                        weight_by_mass = False)
        # compute the initial amount of these fields
        initial_dict[field] = np.sum(part0['gas'].prop(field))
        #
        # compute the normalization (total amount produced by stars in the sim)
        # so we need to also subtract out the amount initially present
        profile['norm'][field] = np.sum(part['star'].prop(field)) + np.sum(part['gas'].prop(field)) -\
                              initial_dict[field]


    return profile

def compute_galaxy_stats(part):

    #
    # load galaxy and compute half light radius
    #

    return

def grab_halo_info(wdir= FIREDIR):

    #
    # look for halo info of 2 types:
    #

    if (os.isdir(wdir + '/halo/ahf_all-SGK/')):
        filename = wdir + '/halo/ahf_all-SGK/catalog/snapshot_600.z0.000.AHF_halos')

    elif (os.isdir(wdir + '/halo/ahf/')):
        filename = wdir + '/halo/ahf/halo_00000.dat'

    if not (os.isfile(filename)):
        print("Cannot find halo file: ", filename)
        raise RuntimeError

    full_halo_data = np.genfromtxt(filename,names=True)

    halo_data = {'Mvir' :   full_halo_data['Mvir4'][0],
             'Rvir' :   full_halo_data['Rvir12'][0],
             'M_gas':   full_halo_data['M_gas45'][0],
             'M_star':  full_halo_data['M_star65'][0]}

    return halo_data


def compute_datasets(data_dirs, wdir = FIREDIR, outdir = '.'):

    for galdir in data_dirs:

        sim_index = 600

        part = gizmo.io.Read.read_snapshots(['star','gas'], # types of particles to load. Gas and/or stars
                                        'index',        # what the next value describes (check docstring for more)
                                        sim_index,      # simulation output index (nth output)
                                        assign_host_principal_axes=True,    # yes! compute the disk of the galaxy
                                        simulation_directory = wdir + '/' + galdir)

        #
        # (and for the sake of this analysis, I'm also grabbing the 1st snapshot, aka the initial conditions)
        #    no stars in this one, so just need to get gas
        part0 = gizmo.io.Read.read_snapshots(['gas'], 'index', 0,
                                        # assign_host_principal_axes=True, don't need to do this for this one
                                        simulation_directory = wdir + '/' + galdir)


        dataset['profiles']  = compute_profile_dict(part, fields = None):
        dataset['halo_data'] = grab_halo_info(widr + '/' + galdir)

        dd.io.save(outdir + galdir + '_profile.h5', dataset)

    return


if __name__ == "__main__":

    data_dirs = ['m11d_res710','m11e_res7100','m11h_res7100','m11i_res7100',
                 'm11b_res260','m11h_res880','m11q_res880']

    data_dirs = ['m11q_res880']

    compute_datasets(data_dirs, outdir = '.')
