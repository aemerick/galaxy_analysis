from matplotlib import rc

fsize = 17
rc('text', usetex=False)
rc('font', size=fsize)#, ftype=42)
line_width = 3
point_size = 30

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt

import deepdish as dd
import numpy as np

import glob

from galaxy_analysis.utilities import utilities as util

def compute_time_average(field_path,
                         dir = '.', tmin = None, tmax = None,
                         times = None, data_list = None,
                         x_field = 'xbins'):

    if data_list is None:
        data_list = np.sort(glob.glob(dir + '/DD????_galaxy_data.h5'))

    # get list of times from all data sets
    if times is None:
        times = np.zeros(len(data_list))
        for i, d in enumerate(data_list):
            times[i] = dd.io.load(d, '/meta_data/Time')

    # now select which data sets to average
    if tmin is None and tmax is None:
        avg_data_list = data_list
    elif tmin is None:
        avg_data_list = data_list[ times < tmax ]
    elif tmax is None:
        avg_data_list = data_list[  times >= tmin ]
    else:
        avg_data_list = data_list[(times<tmax)*(times>=tmin)]

    # now load all fields, avg, compute min and max
    temp = dd.io.load(avg_data_list[0])
    sum  = util.extract_nested_dict(temp, field_path)
    min  =  np.inf * np.ones(np.size(sum))
    max  = -1 * min
    std  = 0.0 * sum

    del(temp)

    s0 = 0
    s2 = 0.0 * std
    avg = 0.0 * sum
    M2  = 0.0 * sum
    for d in avg_data_list:
        data = dd.io.load(d)
        y    = util.extract_nested_dict(data, field_path)
        sum += y

        min  = np.min( [y,min], axis=0)
        max  = np.max( [y,max], axis=0)

        s0    += 1
        delta  = y - avg
        avg   += delta /(1.0 * s0)
        M2    += delta*(y - avg)

    std = np.sqrt( M2 / (1.0*(s0 - 1)))

    avg = sum / (1.0 * s0)

    field_path[-1] = x_field

    x = util.extract_nested_dict(data, field_path)

    return x, avg, min, max, std

def plot_time_average(x, y, std = None, min = None, max = None):

    fig, ax = plt.subplots()

    fill = False
    if (not std is None):
        fillmin = y - std
        fillmax = y + std
        fill = True
    elif ( not min is None) and (not max is None):
        fillmin = min
        fillmax = max
        fill = True

    if fill:
        ax.fill_between(x, fillmin, fillmax, facecolor = 'black', interpolate=True,
                           alpha = 0.5)

        ax.plot(x, fillmin, color = 'black', lw = 1.5, ls = '-')
        ax.plot(x, fillmax, color = 'black', lw = 1.5, ls = '-')

    ax.plot(x, avg, color = 'black', lw = 3, ls = '--')

    return fig, ax


if __name__=="__main__":

    x, avg, min, max, std = compute_time_average(['gas_profiles','surface_density','disk',
                                           ('enzo','HI_Density')],
                                            tmin = 50.0, tmax = 1000.0)

    # x is bins, want bin centers for plotting
    x = 0.5*(x[1:] + x[:-1])
    print avg
    print std
    print min
    print max
    fig, ax = plot_time_average(x/1000.0, avg, std=std) #min = min, max = max)

    ax.set_xlabel(r'R (kpc)')
    ax.set_ylabel(r'$\Sigma_{\rm HI}}$ (M$_{\odot}$ pc$^{-2}$)')

    ax.set_ylim(1.0E-3, 2.0)
    ax.set_xlim(0.0, 2.25)

    fig.set_size_inches(8,8)
    plt.minorticks_on()
    plt.tight_layout()

    ax.semilogy()

    fig.savefig('test_time_average.png')
    plt.close()
