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
def plot_mass_loading_comparison(work_dir = './', output_dir = None, comparison = None):

    if output_dir is None:
        output_dir = work_dir

    if comparison is None:
        labels = {'3pcH2' : '3.6 pc' , '6pcH2' : '7.2 pc', 'Fiducial' : 'Fiducial'}
        lstyle = {'3pcH2' : '--', '6pcH2' : '-.', 'Fiducial' : '-'}
        dirs   = {'3pcH2' : '../3pc_H2/' , '6pcH2' : '../6pc_H2/', 'Fiducial' : work_dir}

    else:
        dirs = {}
        labels = {}
        lstyle = {}
        for k in comparison.keys():
            dirs[k]   = work_dir + comparison[k][0]
            labels[k] = comparison[k][1]
            lstyle[k] = comparison[k][2]


    gather_keys = {'mass_outflow' : ['gas_profiles','outflow','sphere']}

    all_data = {}

    for sim in labels.keys():
        data_list, times = utilities.select_data_by_time(dir = dirs[sim],
                                                         tmin=0.0,tmax= 1000.0)
        all_data[sim] = {}
        all_data[sim]['times'] = times
        for k in gather_keys.keys():

            all_data[sim][k] = utilities.extract_nested_dict_asarray(None, gather_keys[k], data_list, False)

        all_data[sim]['SFR'] = utilities.extract_nested_dict_asarray(None, ['time_data','SFR'], data_list[-1], False)
        all_data[sim]['SFR_time'] = utilities.extract_nested_dict_asarray(None, ['time_data','time'], data_list[-1], False)
        all_data[sim]['SFR_fit'] = lambda x: np.interp(x, all_data['sim']['SFR_time'] / 1.0E6, all_data['sim']['SFR'])

    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    for sim in all_data.keys():
        t = all_data[sim]['times'] - all_data[sim]['times'][0]

        Mdot  = all_data[sim][('gas','cell_mass')][:,1]
        SFR   = all_data[sim]['SFR_fit'](t)

        ax.plot(t, Mdot / SFR, lw = line_width,  ls = lstyle[sim], color = color[sim])

    # plot this so it shows up on diagram
#    ax.plot([-10,-1],[0,0],lw=line_width,ls='-',color='navy', label = 'Disk')
#    ax.plot([-10,-1],[0,0],lw=line_width,ls='-',color='C1', label = 'Outside Halo')

    ax.semilogy()
    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(r'Mass Loading Factor at 0.25 R$_{\rm vir}$')
    ax.set_xlim(0.0, TMAX)
    ax.legend(loc = 'best')
    plt.minorticks_on()
    plt.tight_layout()
    fig.savefig(work_dir + output_dir + 'mass_loading_comparison.png')
    plt.close()

    return

if __name__ == "__main__":

    work_dir = './'
    if len(sys.argv) > 1:
        work_dir = sys.argv[1]
    output_dir = None
    if len(sys.argv) > 2:
        output_dir = sys.argv[2]

    plot_mass_loading_comparison(work_dir = work_dir, output_dir = output_dir)
