from galaxy_analysis.plot.plot_styles import *
import deepdish as dd
import matplotlib.pyplot as plt
import numpy as np
import glob

labels = {'3pc_hsn' : '3.6 pc - SNx2', '3pc' : '3.6 pc', 'final_sndriving' : 'Fiducial', '6pc_hsn' : '7.2 pc'}
lstyle     = {'3pc_hsn' : '--', '3pc' : ':', 'final_sndriving' : '-', '6pc_hsn' : '-.'}

#for l in labels:
#    lstyle[l] = '-'
#
#
filepath = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11_30km/'

wdirs = {}
for k in labels:
    wdirs[k] = filepath + k + '/'

all_data = {}
for k in labels.keys():
    data_list = np.sort(glob.glob(wdirs[k] + 'DD*_galaxy*.h5'))
    all_data[k] = dd.io.load( data_list[-1] )

#data_list, times = utilities.select_data_by_time(dir = work_dir,
#                                                 tmin=0.0,tmax=1000)


fig, ax = plt.subplots()
fig.set_size_inches(8,8)

for k in ['final_sndriving','3pc','3pc_hsn','6pc_hsn']: # labels.keys():
    x = (all_data[k]['time_data']['time'] - all_data[k]['time_data']['time'][0])/1.0E6

    plot_histogram(ax, x,
                          all_data[k]['time_data']['SFR'], lw = line_width, ls = lstyle[k], label = labels[k])

ax.set_xlabel(r'Time (Myr)')
ax.set_ylabel(r'SFR (M$_{\odot}$ yr$^{-1}$)')
ax.semilogy()
ax.set_ylim(1.0E-5, 2.0E-3)
ax.set_xlim(0.0, 500.0)
ax.legend(loc='best')

plt.tight_layout()
fig.savefig('sfr_resolution_study.png')
