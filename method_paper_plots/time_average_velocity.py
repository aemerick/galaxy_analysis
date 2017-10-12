from galaxy_analysis.plot.plot_styles import *
from galaxy_analysis.analysis.compute_time_average import compute_time_average
from galaxy_analysis.utilities import utilities

import numpy as np
import matplotlib.pyplot as plt

filepath = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11_30km/final_sndriving'

phase_colors = {'cold' : 'C0', 'warm' : 'C1', 'hot' : 'C3'}
labels = {'cold' : 'Cold' , 'warm' : 'Warm', 'hot' : 'Hot'}

fig, ax  = plt.subplots()
fig.set_size_inches(8,8)

sum = None
for phase in ['cold','warm','hot']:

    x,avg,min,max,std = compute_time_average(['gas_profiles','velocity','halo',phase], tmin = 100, tmax = 200,
                                             dir = filepath, x_field = 'vbins')
    x, avg = utilities.simple_rebin(x, avg, 10) # re-bin in 10 km/s

    plot_histogram(ax, x, avg, color = phase_colors[phase], lw = line_width,
                               ls = '-', label = labels[phase])
    if sum is None:
        sum = 1.0 * avg
    else:
        sum += avg

plot_histogram(ax, x, sum, color = 'black', lw = line_width, ls = '-', label = 'Total')

ax.set_xlabel(r'Outflow Velocity (km s$^{-1})$')
ax.set_ylabel(r'Mass (M$_{\odot}$)')
ax.semilogy()
ax.set_xlim(0.0 ,np.max(x))
ax.set_ylim(1.0, 1.0E5)
plt.minorticks_on()
ax.legend(loc='best')

fig.savefig('velocity_distribution_time_average.png')
plt.close()

cum_sum = np.cumsum(sum)
percent = cum_sum / (cum_sum[-1]) * 100

for q in np.arange(0,100,5):
    bin = len(percent[percent <= q])
    print q, " percentile ", bin, x[bin]
