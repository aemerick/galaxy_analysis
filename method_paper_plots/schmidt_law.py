from galaxy_analysis.plot.plot_styles import *
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import deepdish as dd
from galaxy_analysis.utilities import utilities
import numpy as np
psize = 40.0

plot_obs = True     # plot just HI * 1.36 and no H2
plot_theory = False # plot total gas mass (HI + H2 + He + HII)

obs_cmap = 'viridis' # 'plasma'

sim_cmap = 'viridis'

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

#print all_data['A_sf']
#print all_data['A']

#print np.median(all_data['A_sf'])
#print np.max(all_data['A_sf'])
#print np.min(all_data['A_sf'])

times = times - times[0]
#
# Plot SD_SFR vs. atomic
#
fig, ax = plt.subplots()
ax.scatter( np.log10(all_data['SD_HI_sf']), np.log10(all_data['SD_SFR_sf']), label = 'HI Observed',
            c = times, cmap = sim_cmap)
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

if plot_theory:
    ax.scatter( p_x(np.log10(all_data['SD_gas_sf'])), p_y(np.log10(all_data['SD_SFR_sf'])),
               c = times, cmap = sim_cmap, s = psize, alpha = 0.75)
if plot_obs:
    ax.scatter( p_x(np.log10(all_data['SD_gas_sf_obs'])), p_y(np.log10(all_data['SD_SFR_sf'])),
               c = times, cmap = obs_cmap, s = psize, alpha = 0.75)
ax.set_xticks([])
ax.set_yticks([])
#plt.tight_layout()
fig.savefig('overplot_roychowdhury_2014_f8.png')
plt.close()

########################

fig, ax = plt.subplots()

if plot_obs:
    ax.scatter( np.log10(all_data['SD_gas_sf_obs']), np.log10(all_data['SD_SFR_sf']), label = 'This Work', c = times, cmap = obs_cmap, s = psize)
if plot_theory:
    ax.scatter( np.log10(all_data['SD_gas_sf']), np.log10(all_data['SD_SFR_sf']),c = times, cmap = sim_cmap, label = 'Total Gas', s = psize)
#ax.set_xlim(-0.5,5.2)
ax.set_xlim(-0.75, 1.5)
ax.set_ylim(-5.8, -2.25)
ax.set_ylabel(r'log($\Sigma_{\rm SFR}$ / (M$_{\odot}$ yr$^{-1}$ kpc$^{-2}$))')
ax.set_xlabel(r'log($\Sigma_{\rm gas}$ / (M$_{\odot}$ pc$^{-2}$))')
#ax.legend(loc = 'best')

diff = np.log10(all_data['SD_gas_sf']) - np.log10(all_data['SD_gas_sf_obs'])
#print diff
#print np.min(diff), np.max(diff), np.median(diff), np.average(diff)
ratio = (np.log10(all_data['SD_gas_sf_obs']) / np.log10(all_data['SD_gas_sf']))
#print np.min(ratio), np.max(ratio), np.average(ratio)

# from Shi et. al. 2011 - Table 4
#x = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 100)
#y = -3.90 + 1.38 * x
#std = 0.49
#ax.plot(x, y, lw = 3, ls = '-', color = 'black')
#ax.plot(x, y - std, lw = 3, ls = '--', color = 'black')
#ax.plot(x, y + std, lw = 3, ls = '--', color = 'black')

# From Towchudhury 2014
x = np.linspace(ax.get_xlim()[0], 2.0, 100)
y = 0.91 * ( x -np.log10( 1.36)) - 3.844
#yp = (0.91 + 0.00) * ( x - np.log10(1.36)) - 3.99
#ym = (0.91 -0.00) * ( x -np.log10(1.36)) - (3.84 - 0.19)
ax.plot(x, y, lw = 4, ls = '--', color = 'C0', label = 'Roychowdhury et. al. 2014 best fit')
#ax.plot(x, yp, lw = 3, ls = '--', color ='blue')
#ax.plot(x, ym, lw = 3, ls = '--', color = 'blue')

# from Kennicutt & Evans 2012
#y = 2.0 * x + np.log10(1.7E-5)
#ax.plot(x, y, lw = 3, ls = '-', color = 'red')

# Teich et. al. 2016
teich_HI = np.array([0.65,0.21,0.50,0.43,0.67,0.19,0.32,0.37,0.50,-0.10,0.74,0.13]) + np.log10(1.34)
teich_HI_error = np.array([0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.13,0.07,0.07,0.07,0.07])
teich_SFR = np.array([-2.75,-3.13,-2.99,-2.61,-3.48,-2.97,-3.27,-2.91,-2.54,-3.28,-2.76,-3.31])
teich_SFR_error = np.array([0.16,0.16,0.16,0.16,0.17,0.16,0.20,0.16,0.1,0.16,0.16,0.16])

ax.errorbar(teich_HI, teich_SFR, xerr=teich_HI_error, yerr=teich_SFR_error, color = 'black', lw = 3, fmt = 's', label = 'Teich et. al. 2016')
#ax.scatter(teich_HI, teich_SFR, color = 'black', marker = 's', s = 40.0, label = 'Teich et. al. 2016')


# Wang et. al. 2017
wang_HI  = np.array([-1,-1,-1,6.53,5.44,2.94,4.74,1.83,7.73,-1,5.09,-1,-1,8.45,3.8,-1,5.08,-1,-1,-1,-1,3.91,-1,-1,-1,4.41,4.64,-1,-1,-1,4.15,-1,-1,-1,-1,-1,-1,4.08,-1,-1,-1,-1,2.56,1.60,-1,1.77,-1,-1,-1,6.97,4.02,-1,4.25,-1,1.99,-1,3.74,-1,2.65,6.77,-1,-1,2.47,-1,-1,3.72,-1,-1,-1,-1,-1,3.67,-1,8.19,-1,-1,3.92,8.03,-1,7.02,2.12,6.22])
wang_SFR = np.array([-2.92,-3.16,-3.67,-3.14,-2.86,-2.38,-3.37,-2.97,-2.74,-3.45,-3.48,0.0,-3.02,-2.40,-3.15,-4.35,-3.00,-2.06,0.0,0.0,0.0,-3.71,-2.66,0.0,0.0,0.0,-2.32,0.0,-3.35,-1.19,-2.60,0.0,0.0,0.0,0.0,0.0,-3.12,0.0,0.0,-3.06,0.0,-3.57,-2.74,-3.81,-3.97,0.0,-3.99,-3.05,-2.97,-2.85,-4.33,-3.39,-2.11,-3.55,-3.07,-2.89,-1.79,-2.26,-3.20,0.0,0.0,0.0,-4.07,0.0,-2.79,-1.46,0.0,0.0,0.0,0.0,0.0,-3.78,0.0,0.0,-1.59,0.0,-3.30,-2.59,-2.83,-3.33,-3.14,-2.64])
if np.size(wang_HI) != np.size(wang_SFR):
    apples = oranges
select = (wang_HI > 0.0) * (wang_SFR < 0.0)
wang_HI = wang_HI[select]
wang_SFR = wang_SFR[select]
#ax.scatter(wang_HI, wang_SFR, color = 'C0', marker = 'd', label = 'Wang et. al. 2017')

plt.minorticks_on()
fig.set_size_inches(8,8)
plt.tight_layout()
ax.legend(loc='lower right')
leg = ax.get_legend()
leg.legendHandles[1].set_color(viridis(0.75))
fig.savefig('all_gas_schmidt_law_evolution.png')
plt.close()

#
#
#

fig, ax = plt.subplots()

x = np.log10(all_data['SD_gas_sf_obs'])
if plot_obs:
    ax.scatter( x , np.log10(all_data['SD_SFR_sf']) - x - 6, c = times, cmap=obs_cmap, s = psize, label = 'This Work')
if plot_theory:
    ax.scatter( np.log10(all_data['SD_gas_sf']) , np.log10(all_data['SD_SFR_sf']) - np.log10(all_data['SD_gas_sf']) - 6, c = times, cmap = sim_cmap, s = psize)
ax.set_xlim(-0.75, 1.5)
ax.set_ylim(-11, -8.5)
ax.set_ylabel(r'log(($\Sigma_{\rm SFR}$ / $\Sigma_{\rm gas}$) / yr$^{-1}$)')
ax.set_xlabel(r'log($\Sigma_{\rm gas}$ / (M$_{\odot}$ pc$^{-2}$))')

# from Shi et. al. 2011 - Table 4
x = np.linspace(-2, 2, 100)
y = -9.85 + 0.35 * x
std = 0.49
ax.plot(x, y, lw = 3, ls = '-', color = 'black')
ax.plot(x, y + std, lw = 3, ls = '--', color = 'black')
ax.plot(x, y - std, lw = 3, ls = '--', color = 'black')

# from Teich et. al. 2016
teich_SFE = np.array([-9.40, -9.33, -9.49, -9.03, -10.14, -9.16, -9.59, -9.28, -9.04, -9.18, -9.50, -9.44])
teich_SFE_error = np.array([0.33,0.33,0.31,0.33,0.30,0.32,0.22,0.70,0.26,0.31,0.33,0.32])

ax.errorbar(teich_HI, teich_SFE, xerr=teich_HI_error, yerr=teich_SFE_error, color = 'black', lw = 3, fmt = 's')
ax.scatter(teich_HI, teich_SFE, marker = 's', color = 'black', s = 50.0, label = 'Teich et. al. 2016')

# 1.0E6 converts from kpc^-2 to pc^-2 before making ratio since HI in pc^-2
# ax.scatter(wang_HI, (wang_SFR/(1.0E6)) /wang_HI , color = 'C0', marker = 'd', label = 'Wang et. al. 2017')


plt.minorticks_on()
fig.set_size_inches(8,8)
plt.tight_layout()
ax.legend(loc='upper left')
leg = ax.get_legend()
leg.legendHandles[0].set_color(viridis(0.75))
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
if plot_theory:
    ax.scatter( p_x(xdata), p_y(ydata),
               c = times, cmap = sim_cmap, s = psize)
if plot_obs:
    ax.scatter( p_x(np.log10(all_data['SD_gas_sf_obs'])), p_y(ydata),
               c = times, cmap = obs_cmap, s = psize)
ax.set_xticks([])
ax.set_yticks([])
#plt.tight_layout()
fig.savefig('overplot_roychowdhury_2017_f2.png')
plt.close()

