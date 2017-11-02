from galaxy_analysis.plot.plot_styles import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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
        orbit_data['particles']   = {}

#    hf   = h5py.File( wdir + '/' + filename, 'w')
#
#    if not ('times' in hf.keys()):
#        hf.create_dataset

    # identify the data outputs in wdir
    data_files   = np.sort(glob.glob(wdir +'/' + 'DD????/DD????'))
    dnames       = [x[-6:] for x in data_files]

    # check last file - see if there are any new particles
    ds   = yt.load( data_files[-1] )
    data = ds.all_data()

    latest_pid    = data['particle_index'].value
    new_particles = np.setdiff1d(latest_pid, pid)
    for p in new_particles:
        orbit_data['particles'][p] = {}
        for k in ['x','y','z','vx','vy','vz','pt','M']:
            orbit_data['particles'][p][k] = np.ones(np.size(dnames))* (-999)

        orbit_data['particles'][p]['t_o'] = data['creation_time'][data['particle_index'] == p].convert_to_units('Myr').value
        orbit_data['particles'][p]['M_o'] = data['birth_mass'][data['particle_index'] == p].value

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

    for i,d in enumerate(dnames):
        ds  = yt.load( data_files[i] )
        data = ds.all_data()
        t   = ds.current_time.convert_to_units("Myr").value

        if all_times[i] == t:
            continue

        all_times[i] = t

        for di,p in enumerate(data['particle_index'].value):
            j = 0
            for coord in ['x','y','z']:
                orbit_data['particles'][p][coord][i] = 1.0*(data['particle_position_' + coord][di] - ds.domain_center[j]).convert_to_units('pc').value
                orbit_data['particles'][p]['v' + coord][i] = 1.0*(data['particle_velocity_'+coord][di].convert_to_units('km/s').value)
                j = j + 1

            orbit_data['particles'][p]['M'][i]   = 1.0*data['particle_mass'][di].convert_to_units('Msun').value
            orbit_data['particles'][p]['pt'][i]  = 1.0*data['particle_type'][di].value

    orbit_data['times'] = all_times

    dd.io.save(wdir + '/' + filename, orbit_data)

    return


def plot_orbit_evolution(wdir = '.', filename = 'orbit.h5',
                         dt = 1.0):

    data  = dd.io.load(wdir + '/' + filename)

    times = data['times']
    pid   = data['particles'].keys()

    plot_times = np.arange(np.min(times), np.max(times) + dt*0.5, dt)
    n_particles = np.size(pid)
    for plot_i, t in enumerate(plot_times):

        all_x = np.array([-1E99] * n_particles)
        all_y = np.array([-1E99] * n_particles)
        all_z = np.array([-1E99] * n_particles)
        # need to generate x,y,z coordinates to plot
        # if each particle is alive
        for i,p in enumerate(pid):
            if data['particles'][p]['t_o'] <= t:
                x,y,z = data['particles'][p]['x'], data['particles'][p]['y'], data['particles'][p]['z']

                select = x != -999
                all_x[i] =  np.interp( t, times[select], x[select])
                all_y[i] =  np.interp( t, times[select], y[select])
                all_z[i] =  np.interp( t, times[select], z[select])

        fig = plt.figure()
        ax  = fig.add_subplot(111,projection='3d')
        fig.set_size_inches(6,6)
        ax.scatter(all_x, all_y, all_z, s = 10, alpha = 0.75, color = 'black')
        ax.set_xlim(-500,500); ax.set_ylim(-500,500); ax.set_zlim(-500,500)
        plt.minorticks_on()
        ax.set_xlabel(r'x (pc)'); ax.set_ylabel(r'y (pc)'); ax.set_zlabel(r'z (pc)')
        fig.savefig("orbit/%00004i_orbit.png"%(plot_i))
        plt.close()

    return

if __name__ == "__main__":
#    generate_dataset()
    plot_orbit_evolution()
