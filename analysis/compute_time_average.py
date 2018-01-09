from galaxy_analysis.plot.plot_styles import *
import matplotlib.pyplot as plt

import deepdish as dd
import numpy as np

import glob

from galaxy_analysis.utilities import utilities as util

def compute_time_average(field_path,
                         dir = '.', tmin = None, tmax = None,
                         times = None, data_list = None,
                         self_contained = False, sc_data = None, index = None,
                         x_field = 'xbins', return_quartiles = False):
    """
    Computes the time average of some quantity pre-computed for a given
    set of simulations dumps using the galaxy analysis framework. The quantity
    can be a single value or an array (e.g. a profile); the assumption is that
    binning in consistent if it is the latter. If an 'x_field' is also provided,
    searches for this value in the 'field_path' and returns this as well. Otherwise
    does nothing if x_field is None (as should be the case for a single value).

    Additionally computes the min, max, and standard deviation over the averaging
    time.

    If self_contained is true, then the data is all contained in a single
    file, rather than one file per output. This is currently a hacky way
    to handle this, but it works. If self_contained is true, then
    data_list must be the filename of the global hdf5 file containing all
    data.
    """

    if data_list is None and self_contained:
        print "Must provide name of file as `data_list' if self_contained is True"
        raise ValueError

    if data_list is None:
        data_list = np.sort(glob.glob(dir + '/DD????_galaxy_data.h5'))

    if self_contained:
        if sc_data is None:
            sc_data   = dd.io.load(data_list)
            data_list = np.sort(sc_data.keys())

    # get list of times from all data sets
    if times is None:
        times = np.zeros(len(data_list))

        if self_contained:
            for i, d in enumerate(data_list):
                times[i] = sc_data[d]['general']['Time']
        else:
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
    if self_contained:
        temp = sc_data[avg_data_list[0]]
    else:
        temp = dd.io.load(avg_data_list[0])

    sum  = util.extract_nested_dict(temp, field_path)
    if not (index is None):
        sum = sum[index]
    min  =  np.inf * np.ones(np.size(sum))
    max  = -1 * min
    std  = 0.0 * sum

    del(temp)

    s0 = 0
    s2 = 0.0 * std
    avg = 0.0 * sum
    M2  = 0.0 * sum
    store_profiles = [None] * np.size(avg_data_list)

    for i,d in enumerate(avg_data_list):
        if self_contained:
            data = sc_data[ d ]
        else:
            data = dd.io.load(d)

        y    = util.extract_nested_dict(data, field_path)
        if not (index is None):
            y = y[index]

        sum += y
        store_profiles[i] = y

        min  = np.min( [y,min], axis=0)
        max  = np.max( [y,max], axis=0)

        s0    += 1
        delta  = y - avg
        avg   += delta /(1.0 * s0)
        M2    += delta*(y - avg)

    std = np.sqrt( M2 / (1.0*(s0 - 1)))

    avg = sum / (1.0 * s0)

    q1  = np.percentile(np.array(store_profiles), 25, axis = 0)
    q2  = np.percentile(np.array(store_profiles), 50, axis = 0)
    q3  = np.percentile(np.array(store_profiles), 75, axis = 0)

    if np.size(avg) > 0 and (not (x_field is None)):
        try:
            temp_field_path     = list(field_path)
            temp_field_path[-1] = x_field
            print temp_field_path
            x = util.extract_nested_dict(data, temp_field_path)
        except:
            try:
                temp_field_path = list(field_path)
                del(temp_field_path[-1])
                temp_field_path[-1] = x_field
                x = util.extract_nested_dict(data, temp_field_path)
            except:
                print field_path
                print "x_field not found in current layer or above layer", x_field
                raise ValueError

    else:
        x = None

    if return_quartiles:
        return x, avg, min, max, std, q1, q2, q3
    else:
        return x, avg, min, max, std

def plot_time_average(x, y, std = None, min = None, max = None,
                            facecolor = 'none', color = 'black', hatch = None, hatchcolor=None,
                            label = None, fig = None, ax = None):

    if fig is None and ax is None:
        fig, ax = plt.subplots()

    fill = False
    if (not std is None):
        fillmin = y - std
        fillmax = y + std
        fillmin[fillmin<0] = 0.0
        fill = True
    elif ( not min is None) and (not max is None):
        fillmin = min
        fillmax = max
        fill = True

    if fill:
        ax.fill_between(x, fillmin, fillmax, facecolor = facecolor, interpolate=True,
                           hatch = hatch, edgecolor = hatchcolor, lw = 0.0)

        ax.plot(x, fillmin, color = color, lw = line_width*0.5, ls = '-')
        ax.plot(x, fillmax, color = color, lw = line_width*0.5, ls = '-')

    ax.plot(x, y, color = color, lw = line_width, ls = '--', label = label)

    return fig, ax


if __name__=="__main__":

    # Example useage, where will we compute the time averaged HI
    # surface density profile, plotting the average along with the
    # standard deviation

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
