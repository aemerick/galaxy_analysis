from galaxy_analysis.plot.plot_styles import *
import matplotlib.pyplot as plt
import glob
import deepdish as dd
import yt
from galaxy_analysis.utilities import utilities
import numpy as np
from matplotlib.ticker import NullFormatter
#
from galaxy_analysis.analysis import Galaxy

from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py

# grab the most recent file
workdir = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11_30km/final_sndriving/'
#workdir = '/home/emerick/work/enzo_runs/pleiades/starIC/run11_30km/final_sndriving/'
data_files = np.sort(glob.glob(workdir + 'DD????'))
name = data_files[-1].split('final_sndriving/')[1]
gal = Galaxy(name, wdir = workdir)

#
#
#
def plot_alpha_vs_fe():
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

    return

def plot_alpha_vs_fe_with_histograms():

    sep = 0.02
    left, width = 0.125, 0.65
    bottom, height = 0.1, 0.65
    left_h = left + width + sep
    bottom_h = bottom + height + sep

    rect_scatter = [left,bottom,width,height]
#    rect_colorbar =
#    rect_histx   = [left, bottom_h, width, 0.95 - bottom_h - (left-bottom)]
#    rect_histy   = [left_h, bottom, 0.95 - left_h, height]

#    fig,ax = plt.subplots()
    fig = plt.figure(1, figsize=(8,8))
#    fig.set_size_inches(8,8)

    ax_scatter = plt.axes(rect_scatter)
#    ax_hist_x  = plt.axes(rect_histx)
#    ax_hist_y  = plt.axes(rect_histy)
#    ax_color   = plt.axes(rect_colorbar)

    ptype     = gal.df['particle_type']
    fe_over_h = gal.df[('io','particle_Fe_over_H')]
    alpha     = gal.df[('io','particle_alpha_over_Fe')]
    age       = (gal.ds.current_time - gal.df[('io','creation_time')]).convert_to_units('Myr')

    age = age - np.min(age)

    # scatter plot
    p = ax_scatter.scatter(fe_over_h[ptype==11], alpha[ptype==11],
                  s = point_size, lw = 2, c = age[ptype==11], cmap = 'plasma_r', alpha = 0.75)
    p.set_clim([0.0, np.max(age)])

    cb = fig.colorbar(p, ax = ax_scatter, orientation = 'horizontal', pad = 0.125, fraction = 0.046,
                         aspect = 40)
    cb.set_label(r'Stellar Age (Myr)')
#
#
#
    ax_scatter.set_xlim(-9,-1)
    ax_scatter.set_ylim(-1.75,1.75)
    ax_scatter.tick_params(axis='x',which='minor',bottom='on')
    ax_scatter.tick_params(axis='y',which='minor',bottom='on')

    ax_scatter.set_xlabel(r'[Fe/H]')
    ax_scatter.set_ylabel(r'[$\rm \alpha$/Fe]')
    plt.minorticks_on()

    #
    # find main plot and construct histograms
    #
    divider = make_axes_locatable(ax_scatter)
    left, bottom, width, height  = divider.get_position()
#    width, height = divider.get_horizontal(), divider.get_vertical()
    sep = 0.01
    thickness = np.min( np.array([0.95 - left - width - sep, 0.95 - bottom - height - sep]))
    rect_histx = [left, bottom + height + sep, width, thickness]
    rect_histy = [left + width + sep, bottom, thickness, height]
    ax_hist_x  = plt.axes(rect_histx)
    ax_hist_y  = plt.axes(rect_histy)


    nbins = 100
    hist,bins = np.histogram(fe_over_h, bins = nbins)
    weights   = np.ones(np.size(fe_over_h)) * (1.0 / (1.0*np.max(hist)))
    ax_hist_x.hist(fe_over_h, color = 'C0', bins = nbins, weights = weights)
#    plot_histogram(ax_hist_x, bins, hist / (1.0*np.max(hist)), color = 'black')
    plt.minorticks_on()
#    hist,bins = np.histogram(alpha, bins = 24)
#    plot_histogram(ax_hist_y, bins, hist / (1.0*np.max(hist)), color = 'black', orientation = 'horizontal')
    nbins = 50
    hist,bins = np.histogram(alpha, bins = nbins)
    weights = np.ones(np.size(fe_over_h)) * (1.0 / (1.0*np.max(hist)))
    ax_hist_y.hist(alpha, orientation='horizontal', color = 'C0', bins = nbins, weights = weights)

    ax_hist_x.xaxis.set_major_formatter(NullFormatter())
    ax_hist_y.yaxis.set_major_formatter(NullFormatter())
    ax_hist_x.set_xlim(ax_scatter.get_xlim())
    ax_hist_y.set_ylim(ax_scatter.get_ylim())
    ticks = [0.0,0.25,0.5,0.75,1.0]
    ax_hist_x.set_yticks(ticks)
    ax_hist_y.set_xticks(ticks)
    ax_hist_y.set_xticklabels(ticks, rotation = 270)

    plt.minorticks_on()
#    plt.tight_layout()
    fig.savefig('alpha_over_fe_hist.png')

    return

def plot_time_evolution():
    """
    Make a panel plot of the time evolution of all elemental abundance ratios
    with respect to both H and Fe (as separate plots)
    """
    filename = workdir + '/abundances/abundances/abundances.h5'

    hdf5_data   = h5py.File(filename, 'r')
    dfiles = hdf5_data.keys()
    dfile  = dfiles[-1]  # do this with most recent data file

    data = dd.io.load(filename, '/' + str(dfile))
    elements = utilities.sort_by_anum([x for x in data['abundances'].keys() if (not 'alpha' in x)])
    elements = elements + ['alpha']


    for time_type in ['cumulative','10Myr']:
        for base in ['H','Fe']:
            fig, ax = plt.subplots(4,4, sharex = True, sharey = True)
            fig.set_size_inches(6*4+1,6*4+1)
            fig.subplots_adjust(hspace=0.0, wspace = 0.0)


            i,j = 0,0
            for e in elements:
                if (base == e):
                    continue
                print "plotting " + e + "/" + base + " time evolution"
                index = (i,j)

                t  = data['statistics'][time_type]['bins']
                y  = data['statistics'][time_type][e][base]['median']
                Q1 = data['statistics'][time_type][e][base]['Q1']
                Q3 = data['statistics'][time_type][e][base]['Q3']
                select = (y*0 == 0) # remove nan values

                t  = t[select]
                t  = t - t[0]

                ax[index].plot( t, y[select], lw = line_width, ls = '-', color = 'black', label = r' ' + e +' ')
                ax[index].fill_between(t, Q1[select], Q3[select], color = 'black', alpha = 0.5, lw = 0.5 * line_width)
                ax[index].set_xlim(0.0, np.max(t))
                ax[index].plot( [0.0,1000.0], [0.0,0.0], ls = ':', color = 'black', lw = line_width)
                ax[index].legend(loc = 'upper right')

                j = j + 1
                if j >= 4:
                    j = 0
                    i = i + 1

            for i in np.arange(4):
                ax[(3,i)].set_xlabel(r'Time (Myr)')
                ax[(i,0)].set_ylabel(r'[X/' + base +']')

                if base == 'H':
                    ax[(i,0)].set_ylim(-12.25, 0.125)
                elif base == 'Fe':
                    ax[(i,0)].set_ylim(-3.25, 3.25)

#            for j in np.arange(3):
#                ax[(j,i)].set_xticklabels([])
#                ax[(i,j+1)].set_yticklabels([])
#            ax[(3,i)].set_xticklabels(np.arange(0,np.max(t)+20,20))
#            if base == 'Fe':
#                ax[(i,0)].set_yticklabels([-3,-2,-1,0,1,2,3,])
#            else:
#                ax[(i,0)].set_yticklabels([-12, -10, -8, -6, -4, -2, 0])

            plt.minorticks_on()
            fig.savefig('stellar_x_over_' + base + '_' + time_type +'_evolution.png')
            plt.close()

    return


if __name__ == '__main__':
    plot_time_evolution()

#    plot_alpha_vs_fe_with_histograms()

#    plot_alpha_vs_fe()


