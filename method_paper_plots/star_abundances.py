from galaxy_analysis.plot.plot_styles import *
import matplotlib.pyplot as plt

#import deepdish as dd
import yt
from galaxy_analysis.utilities import utilities
import numpy as np

#
from galaxy_analysis.analysis import Galaxy


workdir = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11_30km/final_sndriving/'

gal = Galaxy('DD0179', wdir = workdir)

#
#
#
fig,ax = plt.subplots()
fig.set_size_inches(8,7)

ptype     = gal.df['particle_type']
fe_over_h = gal.df[('io','particle_Fe_over_H')]
alpha     = gal.df[('io','particle_alpha_over_Fe')]
age       = (gal.ds.current_time - gal.df[('io','creation_time')]).convert_to_units('Myr')

age = age - np.min(age)

p = ax.scatter(fe_over_h[ptype==11], alpha[ptype==11],
              s = point_size, lw = 2, c = age[ptype==11], cmap = 'plasma_r', alpha = 0.75)
p.set_clim([0.0, np.max(age)])
cb = fig.colorbar(p)
cb.set_label(r'Stellar Age (Myr)')

ax.set_xlim(-9,-1)
ax.set_ylim(-1.75,1.75)

ax.set_xlabel(r'[Fe/H]')
ax.set_ylabel(r'[$\rm \alpha$/Fe]')

plt.minorticks_on()
plt.tight_layout()
fig.savefig('alpha_over_fe.png')


