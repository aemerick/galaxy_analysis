import matplotlib
matplotlib.use("agg")
import numpy as np
import deepdish as dd
import sys,os, glob

# PFH python:
###from pfh_python import colors_sps.colors_table as sps
import gizmo_analysis as gizmo
import utilities as gizmo_ut
from utilities.basic.binning import BinClass

FIREDIR = "/scratch/projects/xsede/GalaxiesOnFIRE/metal_diffusion"

def compute_profile_dict(part, part0, r_min, r_max, dr, fields = None):
    if fields is None:
        fields = ["mass","mass.metals",
                  "mass.si","mass.o","mass.c","mass.fe"]

    # profiling utilitiy
    SpeciesProfile = gizmo_ut.particle.SpeciesProfileClass(
                     scaling = "linear", limits=[r_min,r_max], width=dr,
                     dimension_number=3)

    profile = {}
    profile["norm"] = {}
    profile["initial"] = {}
    for field in fields:
        savefield = field.replace(".","_") # make hdf5 happier at the end

        profile[savefield] = SpeciesProfile.get_profiles(part, species=["gas","star"], property_name=field,
                                                        property_statistic="sum",
                                                        # property_select = select_dict,
                                                        weight_by_mass = False)
        # compute the initial amount of these fields
        if (not (part0 is None)):
            profile["initial"][savefield] = np.sum(part0["gas"].prop(field))
        else:
            profile["initial"][savefield] = 0.0

        #
        # compute the normalization (total amount produced by stars in the sim)
        # so we need to also subtract out the amount initially present
        profile["norm"][savefield] = np.sum(part["star"].prop(field)) + np.sum(part["gas"].prop(field)) -\
                                 profile["initial"][savefield]


    return profile

def compute_galaxy_stats(part, Rvir):


    rselect = 0.1 * Rvir

    #
    # load galaxy and compute half light radius
    #
    rgas  = part["gas"].prop("host.distance.spherical.radius")[:,0]
    rstar = part["star"].prop("host.distance.spherical.radius")[:,0]

    galprop = {  "m_star" : np.sum(part["star"].prop("mass")[ rstar < rselect]),
                 "m_gas"  : np.sum(part["gas"].prop("mass")[ rgas < rselect]),
                 "m_hi"   : np.sum(part["gas"].prop("mass")[rgas<rselect] * part["gas"].prop("hydrogen.neutral.fraction")[rgas<rselect])
              }

    return galprop

def grab_halo_info(wdir= FIREDIR):

    #
    # look for halo info of 2 types:
    #

    if (os.path.exists(wdir + "/halo/ahf_all-SGK/")):
        filename = wdir + "/halo/ahf_all-SGK/catalog/snapshot_600.z0.000.AHF_halos"
    elif (os.path.exists(wdir + "/halo/ahf/")):
        filename = wdir + "/halo/ahf/halo_00000.dat"
    else:
        print("Cannot find halo dir for ", wdir)
        raise RuntimeError

    if not (os.path.isfile(filename)):
        print("Cannot find halo file: ", filename)
        raise RuntimeError

    full_halo_data = np.genfromtxt(filename,names=True)

    print(full_halo_data.dtype.names)

    Mvirkey = [x for x in full_halo_data.dtype.names if 'Mvir' in x][0]
    Rvirkey = [x for x in full_halo_data.dtype.names if 'Rvir' in x][0]
    Mstarkey = [x for x in full_halo_data.dtype.names if 'M_star' in x][0]
    Mgaskey = [x for x in full_halo_data.dtype.names if "M_gas" in x][0]

    halo_data = {"Mvir" :   full_halo_data[Mvirkey][0],
             "Rvir" :   full_halo_data[Rvirkey][0],
             "M_gas":   full_halo_data[Mgaskey][0],
             "M_star":  full_halo_data[Mstarkey][0]}

    return halo_data


def compute_datasets(data_dirs, wdir = FIREDIR, outdir = "."):

    dataset = {}

    for galdir in data_dirs:
        dataset["halo_data"] = grab_halo_info(wdir + "/" + galdir)

        sim_index = 600

        part = gizmo.io.Read.read_snapshots(["star","gas"], # types of particles to load. Gas and/or stars
                                        "index",        # what the next value describes (check docstring for more)
                                        sim_index,      # simulation output index (nth output)
                                        assign_host_principal_axes=True,    # yes! compute the disk of the galaxy
                                        simulation_directory = wdir + "/" + galdir)

        #
        # (and for the sake of this analysis, I"m also grabbing the 1st snapshot, aka the initial conditions)
        #    no stars in this one, so just need to get gas
        try:
            part0 = gizmo.io.Read.read_snapshots(["gas"], "index", 0,
                                        # assign_host_principal_axes=True, don"t need to do this for this one
                                        simulation_directory = wdir + "/" + galdir)
        except:
            part0 = None

        r_min, r_max = 0.0, 1.0 # in units of R_vir
        r_min = r_min * dataset["halo_data"]["Rvir"] # now in kpc
        r_max = r_max * dataset["halo_data"]["Rvir"]
        nbins = 100
        dr    = (r_max-r_min)/(float(nbins))
        dataset["profiles"]  = compute_profile_dict(part, part0,r_min,r_max,dr)
        dataset["stats"] = compute_galaxy_stats(part, r_max)


        outpath = outdir + "/" + galdir + "_profiles.h5"
        print("saving to ", outpath)

        dd.io.save(outpath, dataset)

    return


if __name__ == "__main__":

    possible_data_dirs = ["m11d_res7100","m11e_res7100","m11h_res7100","m11i_res7100",
                 "m11b_res260","m11h_res880","m11q_res880"]

    
    data_dirs = [x for x in possible_data_dirs if (not (os.path.isfile('./'+x+"_profiles.h5")))]

    print("Computing for ", data_dirs)
    compute_datasets(data_dirs, outdir = ".")
