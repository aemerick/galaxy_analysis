import numpy as np
import deepdish as dd
import h5py
import glob
import os
import yt

def generate_dataset(wdir = '.', overwrite = False, filename = 'orbit.h5'):

    times = np.array([])
    pid   = np.array([])
    orbit_data = {}
    if os.path.isfile(wdir + '/' + filename):
        with h5py.File(wdir + '/' + filename, 'r') as hf:
            times = hf['times']
            pid   = hf['particles'].keys() # get all particle ID's

        orbit_data = dd.io.load(wdir + '/' + filename)
    else:
        orbit_data['times'] = times
        orbit_data['pid']   = {}

#    hf   = h5py.File( wdir + '/' + filename, 'w')
#
#    if not ('times' in hf.keys()):
#        hf.create_dataset
    orbit_data   = dd.io.load(wdir + '/' + filename)

    # identify the data outputs in wdir
    data_files   = np.sort(glob.glob(wdir +'/' + 'DD????/DD????'))
    dnames       = [x[-6:] for x in data_files]

    # check last file - see if there are any new particles
    ds   = yt.load( data_files[-1] )
    data = ds.all_data()

    latest_pid    = data['particle_index']
    new_particles = np.setdiff1d(latest_pid, pid)
    for p in new_particles:
        orbit_data['particles'][p] = {}
        for k in ['x','y','z','vx','vy','vz','pt','M']:
            orbit_data['particles'][p][k] = np.ones(np.size(dsnames))* (-999)

        orbit_data['particles'][p]['t_o'] = data['creation_time'][data['particle_index'] == p].convert_to_units('Myr').value
        orbit_data['particles'][p]['M_o'] = data['birth_mass'][data['particle_index'] == p]

    # if file already existed, go through old PID and append slots if needed
    if len(pid) > 0:
        if len(data_files) > len(times):
            num_new = len(data_files) - len(times)
            for p in pid:
                for k in orbit_data['particles'][p].keys():
                    orbit_data['particles'][p][k] =\
                          orbit_data['particles'][p][k].append( np.ones(num_new) * (-999))

    # now loop and do everything
    all_times = np.zeros(np.size(data_files))

    if np.size(times) > 1:
        all_times[:np.size(times)] = times

    for d,i in enumerate(dnames):
        ds  = yt.load( data_files[i] )
        t   = ds.current_time.convert_to_units("Myr").value

        if all_times[i] == t:
            continue

        all_times[i] = t

        for p,di in enumerate(data['particle_index']):
            j = 0
            for coord in ['x','y','z']:
                orbit_data['particles'][p][coord][i] = (data['particle_position_' + coord][di] - ds.domain_center[j]).convert_to_units('pc').value
                orbit_data['particles'][p]['v' + coord][i] = (data['particle_' + coord + '_velocity'][di].convert_to_units('km/s').value)
                j = j + 1

            orbit_data['M']   = data['particle_mass'][di].convert_to_units('Msun').value
            orbit_data['pt']  = data['particle_type'][di]

    orbit_data['times'] = all_times

    dd.io.save(wdir + '/' + filename, orbit_data)

    return
