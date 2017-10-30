from galaxy_analysis.plot.plot_styles import *
import matplotlib.pyplot as plt

import deepdish as dd
from galaxy_analysis.utilities import utilities
import numpy as np

from scipy.interpolate import interp1d

work_dir      = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11_30km/final_sndriving/'
#work_dir      = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11/corrected_sndriving/'


data_list, times = utilities.select_data_by_time(dir = work_dir,
                                                 tmin=0.0,tmax=1000)

x = dd.io.load(data_list[0], '/observables')

# pre-load all the data
_all_data = {}
for i,k in enumerate(data_list):
    _all_data[k] = dd.io.load(k, '/observables')

# gather all data so it can be readily plotted
all_data = {}
for k in _all_data[data_list[0]].keys():
    all_data[k] = utilities.extract_nested_dict_asarray(_all_data, [k], self_contained = True)

#
# Plot SD_SFR vs. atomic
#
fig, ax = plt.subplots()

ax.scatter( np.log10(all_data['SD_HI_sf']), np.log10(all_data['SD_SFR_sf'], label = 'HI Observed')
ax.set_xlim(-0.5,5.2)
ax.set_ylim(-5.0, 3.5)
ax.set_ylabel(r'log($\Sigma_{\rm SFR}$ / (M$_{\odot}$ yr$^{-1}$ kpc$^{-2}$))')
ax.set_xlabel(r'log($\Sigma_{\rm HI}$ / (M$_{\odot}$ pc$^{-2}$))')
plt.minorticks_on()
fig.set_size_inches(8,8)
plt.tight_layout()
fig.savefig('atomic_schmidt_law_evolution.png')
plt.close()

#
#
#

fig, ax = plt.subplots()

ax.scatter( np.log10(all_data['SD_total_obs_sf']), np.log10(all_data['SD_SFR_sf'], label = 'HI Observed')
ax.set_xlim(-0.5,5.2)
ax.set_ylim(-5.0, 3.5)
ax.set_ylabel(r'log($\Sigma_{\rm SFR}$ / (M$_{\odot}$ yr$^{-1}$ kpc$^{-2}$))')
ax.set_xlabel(r'log($\Sigma_{\rm gas, obs}$ / (M$_{\odot}$ pc$^{-2}$))')

# from Shi et. al. 2011 - Table 4
x = np.linspace(ax.get_xlim(), 100)
y = -3.90 + 1.38 * x
std = 0.112
ax.plot(x, y, lw = 3, ls = '-', color = 'black')
ax.plot(x, y - std, lw = 3, ls = '--', color = 'black')
ax.plot(x, y + std, lw = 3, ls = '--', color = 'black')

plt.minorticks_on()
fig.set_size_inches(8,8)
plt.tight_layout()
fig.savefig('all_gas_schmidt_law_evolution.png')
plt.close()

#
#
#

fig, ax = plt.subplots()

x = np.log10(all_data['SD_gas_obs_sf'])
ax.scatter( x , np.log10(all_data['SD_SFR_sf'] - x - 6)
ax.set_xlim(-0.5,5.2)
ax.set_ylim(-11.2, -6)
ax.set_ylabel(r'log(($\Sigma_{\rm SFR}$ / $\Sigma_{\rm gas}$) / yr$^{-1}$)')
ax.set_xlabel(r'log($\Sigma_{\rm HI}$ / (M$_{\odot}$ pc$^{-2}$))')

# from Shi et. al. 2011 - Table 4
x = np.linspace(-0.5, 5.2, 100)
y = -9.85 + 0.35 * x
std = 0.092
ax.plot(x, y, lw = 3, ls = '-', color = 'black')
ax.plot(x, y + std, lw = 3, ls = '--', color = 'black')
ax.plot(x, y - std, lw = 3, ls = '--', color = 'black')

plt.minorticks_on()
fig.set_size_inches(8,8)
plt.tight_layout()
fig.savefig('all_gass_efficiency_evolution.png')
plt.close()

