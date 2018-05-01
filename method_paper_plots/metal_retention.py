from galaxy_analysis.plot.plot_styles import *
from galaxy_analysis.utilities import utilities
#----------------------------------------------


import matplotlib.pyplot as plt
import numpy as np
import glob as glob
import deepdish as dd

TMAX = 500.0

line_width = 3.0

# would be nice to start making gather functions
# for all of these plot functions to not have to
# do any more looping over ALL data sets to gather
# wdir = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11_30km/'
def plot_metal_retention(workdir = './', outdir = './'):

    labels     = {'Halo' : 'CGM' , 'Disk' : 'Disk', 'Outside Halo' : 'Outside Halo'} 
    lstyle     = {'Halo' : '-', 'Disk' : '--', 'Outside Halo' : ':'}

    gather_keys = {'Disk' : ['gas_meta_data', 'masses', 'Disk', 'Total Tracked Metals'],
                   'Halo' : ['gas_meta_data', 'masses', 'Halo', 'Total Tracked Metals'],
                   'FB'   : ['gas_meta_data', 'masses', 'FullBox', 'Total Tracked Metals'],
                   'Outside Box' : ['gas_meta_data', 'masses', 'OutsideBox', 'Total Tracked Metals']}

    all_data = {}

    data_list, times = utilities.select_data_by_time(dir = workdir,
                                                     tmin=0.0,tmax= 650.0)
    all_data['times'] = times
    for k in gather_keys.keys():
        all_data[k] = utilities.extract_nested_dict_asarray(None, gather_keys[k], data_list, False)

    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    total = all_data['FB'] + all_data['Outside Box']
    disk_frac = all_data['Disk'] / total
    halo_frac = all_data['Halo'] / total
    outside_halo_frac = (all_data['FB'] - all_data['Halo'] - all_data['Disk'] + all_data['Outside Box']) / total

    t = all_data['times'] - all_data['times'][0]

    ax.plot(t, halo_frac, lw = line_width, ls = lstyle['Halo'], color = 'black', label = labels['Halo'])
    ax.plot(t, disk_frac, lw = line_width,  ls = lstyle['Disk'], color = 'black', label = labels['Disk'])
    ax.plot(t, outside_halo_frac, lw = line_width, ls = lstyle['Outside Halo'], color = 'black', label = labels['Outside Halo'])

    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(r'Fraction of Metals')
    ax.set_xlim(0.0, TMAX)
    ax.set_ylim(0.0, 1.0)
    ax.legend(loc = 'best')
    plt.minorticks_on()
    plt.tight_layout()
    fig.savefig(outdir + 'metal_retention.png')
    plt.close()

    return

if __name__ == "__main__":

    plot_metal_retention()
