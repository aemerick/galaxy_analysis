import deepdish as dd
import glob as glob
import numpy as np
mpl.use('Agg')
from matplotlib import rc

fsize = 17
rc('text', usetex=False)
rc('font', size=fsize)#, ftype=42)
line_width = 3
point_size = 30

import matplotlib.pyplot as plt


def process_yield_data(dname, norm_data = None):
    """
    Takes a data set and returns a dictionary containing the yields
    for fullbox, halo, stars, and ISM of the galaxy. If a norm_data
    is passed, these values are reduced by those in norm_data (i.e.
    representing those produced between norm_data and data).
    """

    # this is bad, should do this straight from data set
    ele = ['C','N','O','Mg','Si','Ca','Mn','Fe','Ni','Y','Ba','Eu']

    yield_dict = {}
    yield_dict['elements'] = ele

    data      = dd.io.load(dname, '/gas_meta_data/masses')

    if not (norm_data is None):
        norm_data = dd.io.load(norm_data, '/gas_meta_data/masses')

    for k in data.keys():
        yield_dict[k] = np.zeros(len(ele))

        if not (norm_data is None):
            for i,e in enumerate(ele):
                yield_dict[k][i] = data[k][e] - norm_data[k][e]
        else:
            for i,e in enumerate(ele):
                yield_dict[k][i] = data[k][e]

    yield_dict['halo'] = yield_dict['FullBox'] - yield_dict['Disk'] - yield_dict['stars']

    return yield_dict

def plot_yield_data(data, xnames, fig = None, ax = None,
                    color = 'black', label = None, marker = None, s = 25,
                    norm_data = None):

    if fig is None and ax is None:
        fig,ax = plt.subplots()

    if norm_data is None:
        norm_data = np.ones(len(data))

    x = np.arange(1, len(data),+1)
    ax.scatter(x, data / norm_data, c = color, s = s, label = label,
                                    marker = marker)

    ax.set_xlim(ax.get_xlim())
    plt.minorticks_on()

    ax.set_xticks(x)
    ax.set_xticklabels(xnames)
    
    return fig, ax


if __name__ == "__main__":

    data_dict = process_yield_data('DD0483_galaxy_data.h5',norm_data='DD0039_galaxy_data.h5')
    fig, ax = plot_yield_data(data_dict['Disk'], norm_data = data_dict['FullBox'], c = 'black', label = "Disk")
    fig, ax = plot_yield_data(data_dict['Stars'], norm_data = data_dict['FullBox'], c = 'organge', marker = '*', fig = fig, ax = ax, label = "Stars")
    fig, ax = plot_yield_data(data_dict['Halo'], norm_data = data_dict['FullBox'], c = 'purple', marker = 's', label = "Halo", fig = fig, ax = ax)
    fig.save('element_yield_analysis.png')
