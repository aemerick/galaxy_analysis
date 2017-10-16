from galaxy_analysis.plot.plot_styles import *
import matplotlib.pyplot as plt

import deepdish as dd
from galaxy_analysis.utilities import utilities
import numpy as np

from scipy.interpolate import interp1d

outflow_field = 'cell_mass'
work_dir      = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11_30km/final_sndriving/'
#work_dir      = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11/corrected_sndriving/'



data_list, times = utilities.select_data_by_time(dir = work_dir,
                                                 tmin=0.0,tmax=1000)

x = dd.io.load(data_list[0], '/gas_profiles/outflow/sphere')
x = x['centers_rvir']

_temp_dict = {}
all_data = {}
sfr = np.ones(np.size(data_list))
for i,k in enumerate(data_list):
    _temp_dict[k] = dd.io.load(k,
                         '/gas_profiles/outflow/sphere')
    sfr[i] = dd.io.load(k, '/meta_data/SFR')
    all_data[k] =  _temp_dict[k][ ('gas',outflow_field) ] # Becuase i'm dumb
for i in np.arange(np.size(sfr)):
    if sfr[i] > 0.0:
        break
sfr[:i] = sfr[i]
f = interp1d(times[sfr>0], sfr[sfr>0])
sfr[sfr == 0.0] = f(times[sfr==0.0])


#
# Mass Outflow Plot
#
fig, ax = plt.subplots()
fig.set_size_inches(8,8)

binned_y = np.array( [all_data[k] for k in all_data.keys()] )
for i,loc in enumerate([0.1, 0.25, 0.5, 1.0]):
    loc_bin = np.argmin(np.abs(x-loc)) # get the right position bin
    y = binned_y[:,loc_bin]

    bin_edges=utilities.bin_edges(np.round(times))
    # rebin with 10 Myr bins, rather than previous 1 Myr bins
    newx,rebiny=utilities.simple_rebin(1.0*bin_edges,1.0*y,
                   np.arange(np.min(bin_edges),np.max(bin_edges)+2,10), 'average')


    plot_histogram(ax, newx-newx[0], rebiny, lw = line_width, color = plasma(i/4.0),
                   label = r"%0.1f R$_{\rm vir}$"%(loc))

ax.set_xlabel(r'Time (Myr)')
ax.set_ylabel(r'Outflow Rate (M$_{\odot}$ yr$^{-1}$)')
ax.semilogy()
ax.set_xlim(0.0, np.max(newx-newx[0]))
ax.set_ylim(7E-6, 0.01)

plt.tight_layout()
ax.legend(loc='best')
plt.minorticks_on()
fig.savefig('total_mass_outflow.png')
plt.close()

#
# Mass Loading plot
#
fig, ax = plt.subplots()
fig.set_size_inches(8,8)

binned_y = np.array( [all_data[k] for k in all_data.keys()] )
for i,loc in enumerate([0.1, 0.25, 0.5, 1.0]):
    loc_bin = np.argmin(np.abs(x-loc)) # get the right position bin
    y = binned_y[:,loc_bin] / sfr

    bin_edges=utilities.bin_edges(np.round(times))
    # rebin with 10 Myr bins, rather than previous 1 Myr bins
    newx,rebiny=utilities.simple_rebin(1.0*bin_edges,1.0*y,
                   np.arange(np.min(bin_edges),np.max(bin_edges)+2,10), 'average')

    plot_histogram(ax, newx-newx[0], rebiny, lw = line_width, color = plasma(i/4.0),
                   label = r"%0.1f R$_{\rm vir}$"%(loc))

ax.set_xlabel(r'Time (Myr)')
ax.set_ylabel(r'Mass Loading Factor')
ax.semilogy()
ax.set_xlim(0.0, np.max(newx-newx[0]))
ax.set_ylim(0.1,300)

plt.tight_layout()
#ax.legend(loc='best')
plt.minorticks_on()
fig.savefig('total_mass_loading.png')
plt.close()

