import numpy as np
import yt
from yt import derived_field
from galaxy_analysis.plot.plot_styles import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.optimize import brentq
import deepdish as dd
from galaxy_analysis.analysis import compute_time_average as cta

from multiprocessing import Pool
from contextlib import closing
import itertools
import os
import glob

import sys

@derived_field(name = "mag_z", units = "pc")
def _mag_z(field, data):
    return np.abs(data['cylindrical_z']).convert_to_units("pc")

def compute_scale_height(region, rbins):
#    rho  = region['Density']
#    magz = region['mag_z']
    scale_height = np.zeros(np.size(rbins)-1)
    for i in np.arange(1,np.size(rbins)):
        cr_string = "(obj['cylindrical_radius'].in_units('pc') > %.8f) * (obj['cylindrical_radius'].in_units('pc') < %.8f)"%(rbins[i-1],rbins[i])
        cut_region = region.cut_region(cr_string)
        prof = yt.create_profile(cut_region, 'mag_z', "Density", weight_field = 'cell_volume',
                                 logs = {'mag_z':False}, n_bins = 100,
                                 extrema={'mag_z':(0.0*yt.units.pc,900*yt.units.pc)})

        x = prof.x.convert_to_units('pc').value
        y = prof['Density'].convert_to_units('g/cm**3').value
        interp_y = interp1d(x,y)
        root_find = lambda value : interp_y(value) - (y[0]/np.e)
        try:
            scale_height[i-1] = brentq(root_find, x[0], x[-1])
        except:
            scale_height[i-1] = 1000.0

#        scale_height[i-1] = x[np.argmin(np.abs(y-y[0]/(np.e)))]

    return scale_height

_RBINS = np.arange(0,601,20)*yt.units.pc

def _parallel_loop(dsname):

#    for dsname in ds_list:
    print dsname,
    d  = {dsname : {}}
    ds = yt.load(dsname + '/' + dsname)

    rbins  = _RBINS
    region = ds.disk([0.5,0.5,0.5],[0,0,1], np.max(rbins), 2.0*yt.units.kpc)

    d[dsname]['scale_height'] = compute_scale_height(region, rbins.value)
    d[dsname]['times'] = ds.current_time.convert_to_units('Myr')
    d[dsname]['phases'] = {}
    for phase in ['Molecular','CNM','WNM','WIM','HIM']:
        phase_region = region.cut_region(ISM[phase])
        d[dsname]['phases'][phase] = compute_scale_height(phase_region, rbins.value)
    del(region)
    del(ds)

    return d


def compute_all_data(nproc = 28):

    if os.path.isfile('scale_height_data.h5'):
        all_data = dd.io.load('scale_height_data.h5')
    else:
        all_data = {}
        all_data['times'] = None

    ds_list = np.sort(glob.glob('DD????/DD????'))

    if all_data['times'] is None:
        all_data['times'] = np.zeros(np.size(ds_list))
        it = 0
    else:
        old_times = 1.0 * all_data['times']
        all_data['times'] = np.zeros(np.size(ds_list))
        all_data['times'][:np.size(old_times)] = old_times
        it = np.size(old_times)

    already_computed = np.sort([x for x in all_data.keys() if 'DD' in x])
    ds_list = [x.split('/')[0] for x in ds_list if (not (x.split('/')[0] in already_computed))]

    for sub_list in itertools.izip_longest(*(iter(ds_list),) * nproc):
        sub_list = list(sub_list)
        sub_list = [s for s in sub_list if s is not None] # remove None values
        reduced_nproc = np.min( [len(sub_list), nproc] )  # only run on needed processors

        pool    = Pool(reduced_nproc)
        results = pool.map_async(_parallel_loop, sub_list)
        pool.close() # no more processes
        pool.join()  # wait and join running processes

        # gather results and add to output
        for r in results.get():
            print r.keys()[0], r[r.keys()[0]].keys(), it

            all_data[r.keys()[0]] = {}
            for k in r[r.keys()[0]].keys():
                all_data[r.keys()[0]][k] = r[r.keys()[0]][k]
#            all_data[r.keys()[0]]['scale_height'] = r[r.keys()[0]]['scale_height']
#            all_data[r.keys()[0]]['phases'] = r[r.keys()[0]]['phases']
#            all_data['times'][it] = r[r.keys()[0]]['times']
            it = it +1
        del(results)

    all_data['xbins'] = _RBINS

    dd.io.save('scale_height_data.h5', all_data)

    return all_data

def plot_all_data(t_o = 46, dt = 20, t = [150,300,500]):
    """
    Plot data that must already be generated and saved as scale_height_data.h5.

    t_o is the initial time (first star formation)
    t is the list of additional times to plot using a time average of total
    width dt, centered on each t. t is relative to t_o
    """

    data = dd.io.load('scale_height_data.h5')
    data_list = np.array([x for x in data.keys() if 'DD' in x])

    fig,ax = plt.subplots()
    fig.set_size_inches(8,8)

    colors = ['C0','C1','C2','C3','C4','C5']

    i = 0


    plot_histogram(ax, data['xbins'], data['DD0001']['scale_height'],
                      lw = line_width, color = 'black', label = 'Initial Conditions', ls = '--')

    ds_o = data.keys()[np.argmin(np.abs(data['times'] - t_o))] # name of initial dataset

    plot_histogram(ax, data['xbins'], data[ds_o]['scale_height'],
                      lw = line_width, color = plasma(0.0), label = '0 Myr', ls = '-')

    i = 0
    ls = ['--','-.',':']*3
    for t_center in t:
        tmin = t_center - 0.5*dt + t_o
        tmax = t_center + 0.5*dt + t_o
        x, avg, min, max, std = cta.compute_time_average(['scale_height'], dir = '.', data_list = data_list, sc_data=data,  #data_list = 'scale_height_data.h5',
                                                         self_contained = True, x_field = None,
                                                         times = data['times'], tmin = tmin, tmax = tmax)

        plot_histogram(ax, data['xbins'], avg,
                                   lw = line_width, color = plasma((i+1) / (1.0*len(t)+1)), ls = ls[i], label = '%3i Myr'%(t_center+t_o))

#    x, avg, min, max, std = cta.compute_time_average(['scale_height'], dir = '.', data_list = data_list, sc_data = data, #data_list = 'scale_height_data.h5',
#                                                     self_contained = True, x_field = None,
#                                                     times = data['times'], tmin = 326.0, tmax = 366.0)

#    plot_histogram(ax, data['xbins'], avg,
#                               lw = line_width, color = plasma(0.5), label = '300 Myr', ls = '-.')

#    x, avg, min, max, std = cta.compute_time_average(['scale_height'], dir = '.', data_list = data_list, sc_data = data, # data_list = 'scale_height_data.h5',
#                                                     self_contained = True, x_field = None,
#                                                     times = data['times'], tmin = 526.0, tmax = 566.0)
#    plot_histogram(ax, data['xbins'], avg,
#                               lw = line_width, color = plasma(0.75), label = '500 Myr', ls = ':')


    ax.set_xlabel(r'R (pc)')
    ax.set_ylabel(r'Scale Height (pc)')
#    ax.semilogy()
    plt.minorticks_on()
    plt.tight_layout()
    ax.set_xlim(0, 600)
    ax.set_ylim(0, 300)
    ax.legend(loc = 'upper right')
    fig.savefig('scale_height.png')

    return

if __name__ == "__main__":

    if len(sys.argv) == 1:
#    compute_all_data()

        plot_all_data()
    else:
        if len(sys.argv < 4):
            print "if passing arguments need to pass at least 3, t_o, dt, and then at least 1 time to plot (prefereably more)"
            raise ValueError
        # assume all args are provided:
        plot_all_data( t_o = sys.argv[1], dt = sys.argv[2], t = sys.argv[3:])
