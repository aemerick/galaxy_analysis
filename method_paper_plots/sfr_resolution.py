from galaxy_analysis.plot.plot_styles import *
import deepdish as dd
import matplotlib.pyplot as plt
import numpy as np
import glob
import sys

rc('font', size = 22)

def sfr_resolution(work_dir = './', output_dir = None, comparison = None, ylim = None):

    if output_dir is None:
        output_dir = work_dir

    if comparison is None:
        labels = {'3pcH2' : '3.6 pc' , '6pcH2' : '7.2 pc', 'Fiducial' : 'Fiducial'}
        lstyle = {'3pcH2' : '--', '6pcH2' : '-.', 'Fiducial' : '-'}
        dirs   = {'3pcH2' : '../3pc_H2/' , '6pcH2' : '../6pc_H2/', 'Fiducial' : work_dir}
        comparison = {}
        for k in labels.keys():
            comparison[k] = (dirs[k],labels[k],lstyle[k])
    else:
        dirs = {}
        labels = {}
        lstyle = {}
	for k in comparison.keys():
            dirs[k]   = work_dir + comparison[k][0]
            labels[k] = comparison[k][1]
            lstyle[k] = comparison[k][2]

    all_data = {}
    for k in labels.keys():
        data_list = np.sort(glob.glob(dirs[k] + 'DD*_galaxy*.h5'))
        all_data[k] = dd.io.load( data_list[-1] )

    #data_list, times = utilities.select_data_by_time(dir = work_dir,
    #                                                 tmin=0.0,tmax=1000)


    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    for k in comparison.keys(): # labels.keys():
        x = (all_data[k]['time_data']['time'] - all_data[k]['time_data']['time'][0])/1.0E6

        plot_histogram(ax, x,
                          all_data[k]['time_data']['SFR'], lw = line_width, ls = lstyle[k], label = labels[k])

    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(r'SFR (M$_{\odot}$ yr$^{-1}$)')
    ax.semilogy()
    if ylim is None:
        ylim = (1.0E-5, 2.0E-3)
        
    ax.set_ylim(ylim)
    ax.set_xlim(0.0, 500.0)
    ax.legend(loc='best')

    plt.tight_layout()
    fig.savefig(work_dir + output_dir + 'sfr_resolution_study.png')

    return


if __name__ == '__main__':
    work_dir = './'
    if len(sys.argv) > 1:
        work_dir = sys.argv[1]

    sfr_resolution(work_dir)
