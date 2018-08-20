from galaxy_analysis.plot.plot_styles import *
import numpy
import matplotlib.pyplot as plt

rc('font', size = 22)

names = ['Fiducial','Shortrad','No RT']

dirs = {'Fiducial' : './sn_H2atten_H2sh', 'Shortrad' : './sn_H2atten_H2_shortrad',
        'No RT' : './sn_H2atten_H2_noion'}

data = {}
for n in names:
    data[n] = np.genfromtxt(dirs[n] + '/number_density_cut_mass.dat',names=True)


fig, ax = plt.subplots()
fig.set_size_inches(8,8)

colors = {}
colors['Fiducial'] = 'C0'
colors['Shortrad'] = 'C1'
colors['No RT'] = 'C3'

for i,n in enumerate(names):
    ax.plot(data[n]['Time'] - 119.0, data[n]['Mass']/1.0E4, color = colors[n],
             lw = 3, ls = '-', label = n)

ax.set_xlabel(r'Time (Myr)')
ax.set_ylabel(r'Mass Above SF Threshold (10$^4$ M$_{\odot}$)')
#ax.semilogy()
#ax.set_ylim(1.0E2,1.0E5)
ax.set_ylim(0.0, 6E4 / 1.0E4)
ax.set_xlim(0.0, 100.0)
ax.legend(loc='best')
plt.minorticks_on()
plt.tight_layout()
fig.savefig('mass_density_cut.png')
plt.close()

