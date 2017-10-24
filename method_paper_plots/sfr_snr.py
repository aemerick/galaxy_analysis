from galaxy_analysis.plot.plot_styles import *
from galaxy_analysis.analysis.compute_time_average import compute_time_average
from galaxy_analysis.utilities import utilities
import deepdish as dd
import numpy as np
import matplotlib.pyplot as plt
import glob

filepath = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11_30km/final_sndriving'



work_dir      = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11_30km/final_sndriving/'
#work_dir      = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11/corrected_sndriving/'

data_list = np.sort(glob.glob('DD*_galaxy*.h5'))
data = dd.io.load( data_list[-1] )
#data_list, times = utilities.select_data_by_time(dir = work_dir,
#                                                 tmin=0.0,tmax=1000)

all_data = {}
times = data['time_data']['time'] / 1.0E6
sfr   = (data['time_data']['SFR'])
snr   = (data['time_data']['SNII_snr'])
#sfr = np.ones(np.size(data_list))
#snr = np.ones(np.size(data_list))
#for i,k in enumerate(data_list):
#    sfr[i] = dd.io.load(k, '/meta_data/SFR')
#    snr[i] = dd.io.load(k, '/time_data/SNII_snr')[-1]

times = times - times[0]

fig, ax = plt.subplots()
fig.set_size_inches(8,8)

plot_histogram(ax, times, sfr, lw = line_width, color = 'black', label = 'SFR')


ax2 = ax.twinx()
plot_histogram(ax2, times, snr, lw = line_width, color = 'black', label = 'SNR', ls = '--')
ax2.semilogy()
ax2.set_ylabel(r'SNR (yr$^{-1}$)')

ax.set_xlabel(r'Time (Myr)')
ax.set_ylabel(r'SFR (M$_{\odot}$ yr$^{-1}$)')
ax.semilogy()

ax.set_xlim(np.min(times),np.max(times))
#ax.legend(loc='best')
plt.tight_layout()
plt.minorticks_on()

ax.set_ylim(0.5E-6,5.0E-4)

ax2.set_ylim(ax.get_ylim())

fig.savefig('sfr_snr.png')
