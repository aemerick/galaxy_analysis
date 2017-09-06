from galaxy_analysis.plot.plot_styles import *
import matplotlib.pyplot as plt
"-----------------------------------------"

import deepdish as dd
from galaxy_analysis.utilities import utilities
import numpy as np


_phases = ['Molecular','CNM','WNM','WIM','HIM']

def get_fractions(ftype, data = None, tmin = None, tmax = None, average = False):

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
            for k in _phases:
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
            for k in _phases:
                combined_fractions[k] = [all_fractions[i][k] for i in np.arange(len(data_list))]

            return times, combined_fractions
    else:
        fractions = data['gas_meta_data'][ftype +'_fractions']
        return fractions


def plot_mass_fraction(t, y, std = None):

    fig, ax = plt.subplots()

    for k in _phases:
        ax.plot(t, y[k], lw = 3, label = k)

    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(r'ISM Mass Fraction')

    fig.set_size_inches(8,8)
    plt.tight_layout()
    fig.savefig('phase_mass_fraction_evolution.png')

    return

#def plot_volume_fraction(t, y, std = None):
#    return


if __name__ == '__main__':

    times, mass = get_fractions(tmin = 100, tmax = 126, ftype = 'mass')
    plot_mass_fraction(times, mass)
