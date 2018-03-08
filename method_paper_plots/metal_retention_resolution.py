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
wdir = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11_30km/'
def plot_metal_retention_resolution():

    labels     = {'3pc_hsn' : '3.6 pc - SNx2', '3pc' : '3.6 pc', 'final_sndriving' : 'Fiducial', '6pc_hsn' : '7.2 pc'}
    lstyle     = {'3pc_hsn' : '--', '3pc' : ':', 'final_sndriving' : '-', '6pc_hsn' : '-.'}
#    for l in lstyle:
#        lstyle[l] = '-'

    dirs   = {}

    for k in labels.keys():
        dirs[k] = wdir + k + '/'

    gather_keys = {'Disk TM' : ['gas_meta_data', 'masses', 'Disk', 'Total Tracked Metals'],
                   'Halo TM' : ['gas_meta_data', 'masses', 'Halo', 'Total Tracked Metals'],
                   'FB TM'   : ['gas_meta_data', 'masses', 'FullBox', 'Total Tracked Metals'],
                   'Outside Box TM' : ['gas_meta_data', 'masses', 'OutsideBox', 'Total Tracked Metals']}

    all_data = {}

    for sim in labels.keys():
        data_list, times = utilities.select_data_by_time(dir = dirs[sim],
                                                         tmin=0.0,tmax= 650.0)
        all_data[sim] = {}
        all_data[sim]['times'] = times
        for k in gather_keys.keys():

            all_data[sim][k] = utilities.extract_nested_dict_asarray(None, gather_keys[k], data_list, False)

    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    for sim in ['final_sndriving','3pc','3pc_hsn','6pc_hsn']:
        total = all_data[sim]['FB TM'] + all_data[sim]['Outside Box TM']


        disk_frac = all_data[sim]['Disk TM'] / total
        halo_frac = all_data[sim]['Halo TM'] / total
        outside_halo_frac = (all_data[sim]['FB TM'] - all_data[sim]['Halo TM'] - all_data[sim]['Disk TM'] + all_data[sim]['Outside Box TM']) / total

        t = all_data[sim]['times'] - all_data[sim]['times'][0]

        ax.plot(t, disk_frac, lw = line_width,  ls = lstyle[sim], color = 'navy')
        ax.plot(t, outside_halo_frac, lw = line_width, ls = lstyle[sim], color = 'C1')

    # plot this so it shows up on diagram
    ax.plot([-10,-1],[0,0],lw=line_width,ls='-',color='navy', label = 'Disk')
    ax.plot([-10,-1],[0,0],lw=line_width,ls='-',color='C1', label = 'Outside Halo')

    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(r'Fraction of Metals')
    ax.set_xlim(0.0, TMAX)
    ax.set_ylim(0.0, 1.0)
    ax.legend(loc = 'best')
    plt.minorticks_on()
    plt.tight_layout()
    fig.savefig('metal_retention_resolution.png')
    plt.close()

    return

if __name__ == "__main__":

    plot_metal_retention_resolution()
