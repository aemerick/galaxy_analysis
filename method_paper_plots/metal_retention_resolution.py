from galaxy_analysis.plot.plot_styles import *
from galaxy_analysis.utilities import utilities
#----------------------------------------------
rc('font',size=22)

import matplotlib.pyplot as plt
import numpy as np
import glob as glob
import deepdish as dd
import sys

TMAX = 500.0

line_width = 3.0


# would be nice to start making gather functions
# for all of these plot functions to not have to
# do any more looping over ALL data sets to gather
def plot_metal_retention_resolution(work_dir = './', output_dir = None, comparison = None):

    if output_dir is None:
        output_dir = work_dir

    if comparison is None:
        labels = {'3pcH2' : '3.6 pc' , '6pcH2' : '7.2 pc', 'Fiducial' : 'Fiducial'}
        lstyle = {'3pcH2' : '--', '6pcH2' : '-.', 'Fiducial' : '-'}
        dirs   = {'3pcH2' : '../3pc_H2/' , '6pcH2' : '../6pc_H2/', 'Fiducial' : work_dir}

    else:
	for k in comparison.keys():
            dirs[k]   = work_dir + comparison[0]
            labels[k] = comparison[1]
            lstyle[k] = comparison[2]


    gather_keys = {'Disk TM' : ['gas_meta_data', 'masses', 'Disk', 'Total Tracked Metals'],
                   'Halo TM' : ['gas_meta_data', 'masses', 'Halo', 'Total Tracked Metals'],
                   'FB TM'   : ['gas_meta_data', 'masses', 'FullBox', 'Total Tracked Metals'],
                   'Outside Box TM' : ['gas_meta_data', 'masses', 'OutsideBox', 'Total Tracked Metals']}

    all_data = {}

    for sim in labels.keys():
        data_list, times = utilities.select_data_by_time(dir = dirs[sim],
                                                         tmin=0.0,tmax= 1000.0)
        all_data[sim] = {}
        all_data[sim]['times'] = times
        for k in gather_keys.keys():

            all_data[sim][k] = utilities.extract_nested_dict_asarray(None, gather_keys[k], data_list, False)

    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    for sim in all_data.keys():
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
    fig.savefig(work_dir + output_dir + 'metal_retention_resolution.png')
    plt.close()

    return

if __name__ == "__main__":

    work_dir = './'
    if len(sys.argv) > 1:
        work_dir = sys.argv[1]
    output_dir = None
    if len(sys.argv) > 2:
        output_dir = sys.argv[2]

    plot_metal_retention_resolution(work_dir = work_dir, output_dir = output_dir)
