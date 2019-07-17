import yt
import numpy as np
from galaxy_analysis.plot.plot_styles import *
fsize = 22

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
#from galaxy_analysis.utilities import functions
#from galaxy_analysis.static_data import ISM


MAX_NUM = 700


GLOBAL_DR = 20.0 * yt.units.pc  # fix this


# function to do this for a single data set
def stellar_environment(ds, data, dead_only = True, write_to_file = True,
                        output_file = 'stellar_environment.dat', dR = GLOBAL_DR):
    """
    Goes through dataset and computes properties of local ISM
    conditions for ALL stars. Either writes this to file,
    using ``return_type == "file"'' or as a dictionary
    with kwargs for each particle ID number.
    """

    #
    t_now  = ds.current_time.convert_to_units('Myr').value
    min_dx = ds.length_unit.to('pc') / ( ds.domain_dimensions[0] * 2.0**ds.max_level)

    if dR is None:
        dR     = (ds.parameters['IndividualStarFeedbackStencilSize'] + 0.5) * min_dx

    # gather properties for all stars
    pid      = data['particle_index'].value
    ptype    = data['particle_type'].value
    M_o      = data['birth_mass'].value        # initial mass of particle
    M_p      = data['particle_mass'].convert_to_units('Msun').value
    t_o      = data['creation_time'].convert_to_units('Myr').value
    lifetime = data[('io','particle_model_lifetime')].convert_to_units('Myr').value # MS star lifetime
    r_cyl    = data['particle_position_cylindrical_radius'].convert_to_units('pc').value
    z_cyl    = data['particle_position_z'].convert_to_units('pc').value - ds.domain_center[2].to('pc').value
    age      = t_now - t_o
    dynamical_time = data['dynamical_time'].convert_to_units('Myr').value


    # coordinates
    x  = data['x'].convert_to_units('pc')
    y  = data['y'].convert_to_units('pc')
    z  = data['z'].convert_to_units('pc')
    px = data['particle_position_x'].convert_to_units('pc')
    py = data['particle_position_y'].convert_to_units('pc')
    pz = data['particle_position_z'].convert_to_units('pc')

    # now compute the environment properties for all stars
    prop = {}

    if dead_only:
        loop_indexes = np.where(  ( (lifetime-age) > 0) * ( (lifetime-age) <= 1.0) )[0]

    else:
        loop_indexes = np.arange(np.size(pid)) # all stars

    print(np.size(pid), np.size(loop_indexes))

    for i in loop_indexes:

        #if (dead_only and (lifetime - age > 1.0)): # skip for stars still alive
        #    continue

        ID = int(pid[i])

        r = np.sqrt(  (x-px[i])**2 + (y-py[i])**2 + (z-pz[i])**2 )

        select = r <= dR
        #print i, pid[i], np.size(r[select])
        if np.size(r[select]) == 0:
            select = r <= np.min(r)
            #print '---', i, pid[i], np.size(r[select]), np.min(r), dR

        M      = (data['cell_mass'].to('Msun'))[select]
        V      = (data['cell_volume'].to('cm**(3)'))[select]
        n      = (data['number_density'])[select]
        T      = (data['temperature'])[select]

        prop[ID] = {}
        prop[ ID ]['env'] = {'n_min' : np.min(n), 'n_max' : np.max(n),
                                 'n_v_avg' : np.sum(n*V)/np.sum(V), 'n_med' : np.median(n),
                                 'T_m_avg' : np.sum(T*M)/np.sum(M), 'T_v_avg' : np.sum(T*V)/np.sum(V),
                                 'M_tot' : np.sum(M)}

        prop[ID]['c_prop'] = {'M_o' : M_o[i], 'lifetime' : lifetime[i],
                                               't_o' : t_o[i]}

        prop[ID]['s_prop'] = {'age' : age[i], 'r_cyl' : r_cyl[i], 'z_cyl' : z_cyl[i],
                                           'M' : M_p[i], 'ptype' : ptype[i], 'dyn_time' : dynamical_time[i]}


#        if write_to_file:
#            file.write("%i %i"%(ID,ptype[i]))
#
#            results = [M_o[i], r_cyl[i], z_cyl[i], t_now, t_o[i], lifetime[i], age[i],
#                             prop[ID]['env']['n_min'], prop[ID]['env']['n_max'],
#                             prop[ID]['env']['n_v_avg'], prop[ID]['env']['n_med'], prop[ID]['env']['T_m_avg'], 
#                             prop[ID]['env']['T_v_avg'], prop[ID]['env']['M_tot']]
#
#            for val in results:
#                file.write(" %4.4E"%(val))
#
#            file.write("\n")

    return prop


def _parallel_loop(dsname):

    groupname = dsname.rsplit('/')[1]
    print("starting computation on ", groupname)
    gal = Galaxy(groupname)

    dictionary = {groupname : {}}

    g = dictionary[groupname]

    g['Time']    = gal.ds.current_time.convert_to_units('Myr').value
    # generalized function to loop through all mask types and compute stats
    data  = stellar_environment(gal.ds, gal.df)

    for k in list(data.keys()):
        g[k] = data[k]

    del(gal)

    print("ending computation on ", groupname)

    return dictionary

def compute_stats_all_datasets(overwrite = False, 
                               dir = './', outfile = 'stellar_environment.h5',
                               write_to_text = True, text_file = 'stellar_environment.dat',
                               nproc = 24, dR = None):


    if not (dR is None):
        outfile   = "%4.4f"%(dR.value) + outfile
        text_file = "%4.4f"%(dR.value) + text_file

    hdf5_filename = dir + outfile

    if not os.path.isfile(hdf5_filename) or overwrite:
        hf = h5py.File(hdf5_filename, 'w')
        hf.close()

    hf = dd.io.load(hdf5_filename)

    ds_list = np.sort( glob.glob('./DD????/DD????'))
    for i, dsname in enumerate(ds_list):
        ds = yt.load(dsname)
        if ds.parameters['NumberOfParticles'] > 0:
            start_index = i
            del(ds)
            break
        del(ds)

    times = np.zeros(np.size(ds_list))
    ds_list = np.array(ds_list[start_index:])
    times   = np.array(times[start_index:])

    ds_list = ds_list[:np.min([np.size(ds_list),MAX_NUM])]
    times   = times[:np.min([np.size(ds_list),MAX_NUM])]

    if write_to_text:
        if not os.path.exists(text_file):
            file = open(text_file,'w')
            file.write("#PID ptype M_o M r z t t_o lifetime age dyn_time n_min n_max n_v_avg n_med T_m_avg T_v_avg M_tot\n")
        else:
            file = open(text_file,'a') 


    def _write(_file, prop, ID, t_now):
        _file.write("%i %i"%(ID, prop['s_prop']['ptype']))

        results = [prop['c_prop']['M_o'], prop['s_prop']['M'],
                   prop['s_prop']['r_cyl'],
                   prop['s_prop']['z_cyl'], t_now, prop['c_prop']['t_o'], prop['c_prop']['lifetime'], prop['s_prop']['age'],
                   prop['s_prop']['dyn_time'], prop['env']['n_min'], prop['env']['n_max'],
                   prop['env']['n_v_avg'], prop['env']['n_med'], prop['env']['T_m_avg'], 
                   prop['env']['T_v_avg'], prop['env']['M_tot']]

        for val in results:
            file.write(" %4.4E"%(val))

        file.write("\n")
        return

####
    if nproc == 1:
        for i, dsname in enumerate(ds_list):
            #print i, dsname
            groupname = dsname.rsplit('/')[1]
            gal = Galaxy(groupname)

            hf[groupname] = {}
            g = hf[groupname]
            g['Time']    = gal.ds.current_time.convert_to_units('Myr').value
            # generalized function to loop through all mask types and compute stats
            data  = stellar_environment(gal.ds, gal.df, dR = dR)

            for k in list(data.keys()):
                g[k] = data[k]

            if write_to_text:
                for PID in list(g.keys()):
                    if PID == 'Time':
                        continue
                    _write(file, g[PID], PID, g['Time'])

            del(gal)

            print("ending computation on ", groupname)

    else: # parallel

        # select out data sets that already exist in output
        if not overwrite:
            ds_list = [x for x in ds_list if ( not any( [x.rsplit('/')[1] in y for y in list(hf.keys()) ]))]

        # construct the pool, and map the results to a holder
        #   pool splits computation among processors

        #
        # do this in a loop, so we can save progressively once a full set of processors is complete
        #   saves you if there is a crash (kinda). This way we do better memory management if
        #   operating on many large datasets.
        #

        for sub_list in itertools.zip_longest(*(iter(ds_list),) * nproc):

            sub_list = list(sub_list)
            sub_list = [s for s in sub_list if s is not None] # remove None values
            reduced_nproc = np.min( [len(sub_list), nproc] )  # only run on needed processors

            pool = Pool(reduced_nproc)
            results = pool.map_async(_parallel_loop, sub_list)

            pool.close() # no more processes
            pool.join()  # wait and join running processes

            # gather results and add to output
            for r in results.get():
                hf[list(r.keys())[0]] = r[list(r.keys())[0]]

            if write_to_text:
                for dsname in sub_list:
                    k = dsname.split('/')[1]

                    for PID in list(hf[k].keys()):
                        if PID == 'Time':
                            continue

                        _write(file, hf[k][PID], PID, hf[k]['Time'])

            del(results)


    if write_to_text:
        file.close()

    dd.io.save(hdf5_filename, hf)

    return


if __name__ == "__main__":
    # do things here
    compute_stats_all_datasets(nproc = 28, dR = GLOBAL_DR)


