import yt
from galaxy_analysis.analysis import Galaxy
import numpy as np
from galaxy_analysis.static_data import ISM
from galaxy_analysis.plot.plot_styles import *
from scipy.stats import binned_statistic
import matplotlib.pyplot as plt
import sys

if len(sys.argv) > 1:
    dsname = str(sys.argv[1])
else:
    dsname = 'DD0401'

gal = Galaxy(dsname)
#
#
#
#scipy.stats.binned_statistic(x, values, statistic='mean', bins=10, range=None)[source]
#

n = gal.disk['number_density']
M = gal.disk['cell_mass'].convert_to_units('Msun')
M_tot = np.sum(M)

#
n_bins = np.arange(-6, 3.01, 0.05) # 0.05 dex bins in n

def compute_pdf(n, M, norm = None):
    stat, bin_edges, dummy = binned_statistic(np.log10(n), M, statistic = 'sum', bins = n_bins)
    bin_sizes = 0.5 * (10.0**n_bins[1:] + 10.0**(n_bins[:-1]))

#    bin_sizes = 0.5 * (bins

    if norm is None:
        norm = 1.0 * np.sum(M)

    frac = stat / (norm)
    PDF  = frac / bin_sizes
    return PDF

fig, ax = plt.subplots()
fig.set_size_inches(8,8)


PDF = compute_pdf(n,M)
plot_histogram(ax, n_bins, PDF, lw = 3, color = 'black', label = 'Total')


for field in ['CNM','WNM','WIM','HIM']:
    n = gal.disk.cut_region(ISM[field])['number_density']
    M = gal.disk.cut_region(ISM[field])['cell_mass'].convert_to_units('Msun')
    
    PDF = compute_pdf(n,M, norm = M_tot)
    plot_histogram(ax, n_bins, PDF, lw = 3, color = color_dict[field], label = field)


ax.set_xlabel(r'n (cm$^{-3}$)')
ax.set_ylabel(r'PDF')
ax.set_xlim(-6,3)
ax.set_ylim(4.0E-6, 11.0)
ax.legend(loc='best')
ax.semilogy()
plt.tight_layout()
fig.savefig(dsname + '_density_PDF.png')


fig, ax = plt.subplots()
fig.set_size_inches(8,8)


n = gal.disk['number_density']
V = gal.disk['cell_volume']
V_tot = np.sum(V) * 1.0

PDF = compute_pdf(n,V, norm = V_tot)
plot_histogram(ax, n_bins, PDF, lw = 3, color = 'black', label = 'Total')


for field in ['CNM','WNM','WIM','HIM']:
    n = gal.disk.cut_region(ISM[field])['number_density']
    V = gal.disk.cut_region(ISM[field])['cell_volume']
    
    PDF = compute_pdf(n,V, norm = V_tot)
    plot_histogram(ax, n_bins, PDF, lw = 3, color = color_dict[field], label = field)


ax.set_xlabel(r'n (cm$^{-3}$)')
ax.set_ylabel(r'PDF')
ax.set_xlim(-6,3)
ax.set_ylim(1.0E-8, 1.0E4)
ax.legend(loc='best')
ax.semilogy()
plt.tight_layout()
fig.savefig(dsname + '_density_PDF_volume_weighted.png')

