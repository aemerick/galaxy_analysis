
from galaxy_analysis.plot.plot_styles import *
import matplotlib.pyplot as plt
"-----------------------------------------"

import deepdish as dd
from galaxy_analysis.utilities import utilities
from galaxy_analysis.plot.plot_styles import color_dict
import numpy as np
import sys

TMAX = 500.0 # maximum (normalized) plot time


def get_fractions(ftype, data = None, tmin = None, tmax = None,
                  average = False,
                  phases = ['Molecular','CNM','WNM','WIM','HIM']):

    if (not (ftype == 'mass')) and (not (ftype == 'volume')):
        print  "ftype must be either 'mass' or 'volume'"
        raise ValueError

    if data is None:
        if tmin is None or tmax is None:
            print "need to specify data set or time range"
            raise ValueError

        if average:
            # compute averages for each field
            fractions = {}
            min, max, std = {}, {}, {}
            for k in phases:
                x, fractions[k], min[k], max[k], std[k] =\
                     compute_time_average(['gas_meta_data', ftype +'_fractions'.k],
                                          tmin = tmin, tmax = tmax, xbins = None)

            return fractions, min, max, std

        else:
            data_list, times = utilities.select_data_by_time(tmin = tmin, tmax = tmax)

            all_fractions = [None]*len(data_list)
            for i, dname in enumerate(data_list):
                all_fractions[i] = dd.io.load(dname, '/gas_meta_data/' + ftype + '_fractions/')

            # transfer to plotable format
            combined_fractions = {}
            for k in phases:
                combined_fractions[k] = [all_fractions[i][k] for i in np.arange(len(data_list))]

            return times[(times>tmin)*(times<tmax)], combined_fractions
    else:
        fractions = data['gas_meta_data'][ftype +'_fractions']
        return fractions


def plot_mass_fraction(t, y, std = None, outdir = './',
                       phases = ['Molecular','CNM','WNM','WIM','HIM'] ):

    fig, ax = plt.subplots()

    for k in phases:
        ax.plot(t-t[0], y[k], lw = line_width, label = k, color = color_dict[k])

    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(r'ISM Mass Fraction')
    ax.legend(loc='upper right', ncol=2)
    plt.minorticks_on()
    fig.set_size_inches(8,8)
    plt.tight_layout()
    ax.set_ylim(0.0,0.9)
    ax.set_xlim(0.0, np.min([TMAX,np.max(t-t[0])])  )
    fig.savefig('phase_mass_fraction_evolution.png')
    ax.set_ylim(1.0E-5, 1.0)
    ax.semilogy()
    ax.legend(loc='lower right', ncol=2)
    fig.savefig(outdir + 'phase_mass_fraction_evolution_log.png')
    plt.close()

    return

def plot_volume_fraction(t, y, std = None, outdir = './',
                         phases = ['Molecular','CNM','WNM','WIM','HIM']):

    fig, ax = plt.subplots()

    for k in phases:
        ax.plot(t-t[0], y[k], lw = line_width, label = k, color = color_dict[k])

    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(r'ISM Volume Fraction')
    ax.legend(loc='upper right', ncol=2)
    ax.set_xlim(0.0, np.min([TMAX,np.max(t-t[0])]) )
    plt.minorticks_on()
    fig.set_size_inches(8,8)
    ax.set_ylim(0.0, 0.9)
    plt.tight_layout()
    fig.savefig('phase_volume_fraction_evolution.png')
    ax.set_ylim(1.0E-5, 1.0)
    ax.semilogy()
    ax.legend(loc='lower right', ncol = 2)
    fig.savefig(outdir + 'phase_volume_fraction_evolution_log.png')
    plt.close()

    return

def plot_fractions(outdir = './',
                   phases = ['Molecular','CNM','WNM','WIM','HIM']):
    times, mass = get_fractions(tmin = 50, tmax = 1260,
                                ftype = 'mass', phases = phases)
    plot_mass_fraction(times, mass, outdir = outdir, phases = phases)
    times, volume = get_fractions(tmin = 50, tmax = 1260,
                                  ftype = 'volume', phases = phases)
    plot_volume_fraction(times, volume, outdir = outdir, phases = phases)
    return

if __name__=='__main__':

    if len(sys.argv) > 1:
        outdir = sys.argv[1]

    phases = ['CNM','WNM','WIM','HIM']

    plot_fractions(outdir = outdir, phases = phases)
