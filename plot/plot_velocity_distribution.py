from galaxy_analysis.plot.plot_styles import *
import matplotlib.pyplot as plt
"-----------------------------------------"
import sys

from joblib import Parallel, delayed
import multiprocessing

import deepdish as dd
from galaxy_analysis.utilities import utilities
import numpy as np
import glob as glob

from galaxy_analysis.static_data import ISM, CUT_REGION

_mask_colors = {'Molecular' : plasma(0.0), 'CNM' : plasma(1.0/5.0),
                 'WNM' : plasma (2.0/5.0), 'WIM' : plasma(3.0/5.0),
                 'HIM' : plasma(4.0/5.0), 'stars' : plasma(1.0),
                 'cold' : 'blue','warm':'orange','hot':'red'}

_mask_ls    = {'star_forming' : '-',
               'CNM'          : '-',
               'WNM'          : '-',
               'HIM'          : '-',
               'halo'         : '-'}

def plot_velocity_distribution(data_name, data = None, rebin = True):

    if data is None:
        data = dd.io.load(data_name, '/gas_profiles/velocity')

    fig, ax = plt.subplots(1,2)
    fig.set_size_inches(16,8)

    sum = None
    vbins = data['vbins']
    if hasattr(vbins,'value'):
        vbins = vbins.value

    for k in [x for x in ISM.keys() if (not (x=='star_forming'))]:
        x, y = utilities.simple_rebin(vbins, (data['halo'][k]), 10) # rebin to 10 km/s

        plot_histogram(ax[0], x, y,
                       lw = line_width, color = _mask_colors[k], label = k)

        if sum is None:
            sum = np.zeros(np.size(y))
        sum = sum + y

    for k in ['cold','warm','hot']:
         x,y = utilities.simple_rebin(vbins, (data['halo'][k]), 10) # rebin to 10 km/s
         plot_histogram(ax[1], x,y,
                        lw = line_width, color = _mask_colors[k], label = k)


    plt.minorticks_on()
    for a in ax:
        a.semilogy()
        a.set_xlabel('Velocity (km/s)')
        a.set_ylabel('Mass (Msun)')
        a.legend(loc='best')
        a.set_ylim(1.0, 3.0E5)
        plot_histogram(a, x, sum, lw = line_width, color = 'black', label = 'Total')

    name = data_name.split('_')

    plt.savefig(name[0] + '_velocity_distribution.png')

    return


if __name__ == "__main__":

    all_data = np.sort(glob.glob('DD*galaxy*.h5'))

    if len(sys.argv) == 1:
        # assume not parallel
        plot_velocity_distribution(all_data[-1])
    elif len(sys.argv) == 3:

        i,j = int(sys.argv[1]), int(sys.argv[2])
        n_jobs = multiprocessing.cpu_count()

        Parallel(n_jobs = n_jobs)(\
                delayed(plot_velocity_distribution)(x) for x in all_data[i:j])

