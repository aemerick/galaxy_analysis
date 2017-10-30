from galaxy_analysis.plot.plot_styles import *
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import deepdish as dd
from galaxy_analysis.utilities import utilities
import numpy as np

work_dir      = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11_30km/final_sndriving/'
#work_dir      = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11/corrected_sndriving/'

data_list, times = utilities.select_data_by_time(dir = work_dir,
                                                 tmin=0.0,tmax=1000)
# pre-load all the data
_all_data = {}
for i,k in enumerate(data_list):
    try:
        _all_data[k] = dd.io.load(k, '/observables')
    except:
        print 'skipping data entry ', i, k

# gather all data so it can be readily plotted
all_data = {}
for k in _all_data[ _all_data.keys()[0] ].keys():
    all_data[k] = utilities.extract_nested_dict_asarray(_all_data, [k], self_contained = True)

#print all_data['r_sf']

times = times - times[0]
#
# Plot SD_SFR vs. atomic
#
fig, ax = plt.subplots()
ax.scatter( np.log10(all_data['SD_HI_sf']), np.log10(all_data['SD_SFR_sf']), label = 'HI Observed', 
            c = times, cmap = 'viridis')
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
# Do this again, but overplot on the actual image
#
# defined pixel positions previously for Roychowdhury Fig 4
#
x0 = (0,109.5)
x1 = (1,229.0)
y0 = (-4, 529.2)
y1 = (-3, 459.0)
fig = plt.figure()
ax  = fig.add_axes([0.,0.,1.,1.,])
ax.set_xticks([])
ax.set_yticks([])
p_x, p_y = utilities.map_to_pixels(x0,x1,y0,y1)

img=mpimg.imread('Roychowdhury_2014_f8.png')
img_size = np.shape(img)
dpi = 80.0
fig.set_size_inches(img_size[1]/dpi, img_size[0]/dpi)
ax.imshow(img)
ax.scatter( p_x(np.log10(all_data['SD_gas_sf'])), p_y(np.log10(all_data['SD_SFR_sf'])),
               c = times, cmap = 'viridis', s = 10, alpha = 0.75)
ax.scatter( p_x(np.log10(all_data['SD_gas_sf_obs'])), p_y(np.log10(all_data['SD_SFR_sf'])),
               c = times, cmap = 'plasma', s = 10, alpha = 0.75)
ax.set_xticks([])
ax.set_yticks([])
#plt.tight_layout()
fig.savefig('overplot_roychowdhury_2014_f8.png')
plt.close()

########################

fig, ax = plt.subplots()

ax.scatter( np.log10(all_data['SD_gas_sf_obs']), np.log10(all_data['SD_SFR_sf']), label = 'HI Observed', c = times, cmap = 'viridis')
ax.set_xlim(-0.5,5.2)
ax.set_ylim(-5.0, 3.5)
ax.set_ylabel(r'log($\Sigma_{\rm SFR}$ / (M$_{\odot}$ yr$^{-1}$ kpc$^{-2}$))')
ax.set_xlabel(r'log($\Sigma_{\rm gas, obs}$ / (M$_{\odot}$ pc$^{-2}$))')

# from Shi et. al. 2011 - Table 4
x = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 100)
y = -3.90 + 1.38 * x
std = 0.49
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

x = np.log10(all_data['SD_gas_sf_obs'])
ax.scatter( x , np.log10(all_data['SD_SFR_sf']) - x - 6, c = times, cmap = 'viridis')
ax.set_xlim(-0.5,5.2)
ax.set_ylim(-11.2, -6)
ax.set_ylabel(r'log(($\Sigma_{\rm SFR}$ / $\Sigma_{\rm gas}$) / yr$^{-1}$)')
ax.set_xlabel(r'log($\Sigma_{\rm HI}$ / (M$_{\odot}$ pc$^{-2}$))')

# from Shi et. al. 2011 - Table 4
x = np.linspace(-0.5, 5.2, 100)
y = -9.85 + 0.35 * x
std = 0.49
ax.plot(x, y, lw = 3, ls = '-', color = 'black')
ax.plot(x, y + std, lw = 3, ls = '--', color = 'black')
ax.plot(x, y - std, lw = 3, ls = '--', color = 'black')

plt.minorticks_on()
fig.set_size_inches(8,8)
plt.tight_layout()
fig.savefig('all_gass_efficiency_evolution.png')
plt.close()

#
#
#
#
#
x0 = (0,208.9)
x1 = (1,340.0)
y0 = (-11, 646.3)
y1 = (-7, 155.5)
fig = plt.figure()
ax  = fig.add_axes([0.,0.,1.,1.,])
ax.set_xticks([])
ax.set_yticks([])
p_x, p_y = utilities.map_to_pixels(x0,x1,y0,y1)

img=mpimg.imread('roychowdhury_2017_f2.png')
img_size = np.shape(img)
dpi = 80.0
fig.set_size_inches(img_size[1]/dpi, img_size[0]/dpi)
ax.imshow(img)

xdata = np.log10(all_data['SD_gas_sf'])
ydata = np.log10(all_data['SD_SFR_sf']) - xdata - 6
ax.scatter( p_x(xdata), p_y(ydata),
               c = times, cmap = 'viridis', s = 10, alpha = 0.75)
ax.scatter( p_x(np.log10(all_data['SD_gas_sf_obs'])), p_y(ydata),
               c = times, cmap = 'plasma', s = 10, alpha = 0.75)
ax.set_xticks([])
ax.set_yticks([])
#plt.tight_layout()
fig.savefig('overplot_roychowdhury_2017_f2.png')
plt.close()

