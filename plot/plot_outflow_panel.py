from galaxy_analysis.plot.plot_styles import *
import matplotlib.pyplot as plt
"-----------------------------------------"

import deepdish as dd
from galaxy_analysis.utilities import utilities
import numpy as np


def plot_outflow_panel(dir = '.', tmin = None, tmax = None, mass_loading = False):
    """
    Given a directory and desired time range (in Myr), plots the outflow rates or
    mass loading factors for all species over that time at all radii sampled.
    This plot will be VERY large and hard to understand, but will be useful to
    quickly look at many things simultaneously for analysis purposes.
    """

    # plot in this order:
    fields        = ["cell_mass", "H_total_mass", "H2_mass", "metal_mass"]
    fields_labels = {'cell_mass' : 'Total', 'H_total_mass' : "H",
                     'H2_mass' : r'H$_{2}$', "metal_mass" : 'Total Metal'}

    # obtain ds list based off time range selection
    data_list, times = utilities.select_data_by_time(tmin=tmin, tmax = tmax)
    n_non_metal = len(fields)
    # this is a hacky way of getting the metal fields - should just save this to file
    metal_fields = dd.io.load(data_list[0], '/gas_meta_data/masses/FullBox')
    exclude      = ['H','HI','HII','H2','Metals','Total','He','HeI','HeII','HeIII']
    metal_fields = utilities.sort_by_anum([x for x in metal_fields if (not any([y in x for y in exclude]))])

    for k in metal_fields:
        new_name = k + '_Mass'
        fields_labels[new_name] = k
        fields = fields + [new_name]

    nplots = len(fields)
    nrow, ncol = utilities.rowcoldict[nplots]
    fig,ax = plt.subplots(nrow,ncol)

#    key_select = ['gas_profiles', 'outflow', 'sphere']

    x = dd.io.load(data_list[0], '/gas_profiles/outflow/sphere')
    x = x['centers_rvir']

    # load everything to memory at start - limits number of reads from disk
    all_data_dict = {}
    if mass_loading:
        norm = np.ones(np.size(data_list))
        for i,k in enumerate(data_list):
            all_data_dict[k] = dd.io.load(k, '/gas_profiles/outflow/sphere')
            norm[i]          = dd.io.load(k, '/meta_data/SFR_100')
    else:
        for k in data_list:
            all_data_dict[k] = dd.io.load(k, '/gas_profiles/outflow/sphere')
        norm = np.ones(len(all_data_dict.keys()))
    times[norm == 0.0] = None

    axi, axj = 0,0
    for field in fields:
        axind = (axi, axj)

        load_field = ('gas',field) # because I'm dumb and keyed them like this
        binned_y   = utilities.extract_nested_dict_asarray(all_data_dict, load_field)

        # plot at 0.1, 0.25, 0.5, and 1 Rvir for now:
        for loc in [0.1, 0.25, 0.5, 1.0]:
            loc_bin = np.argmin(np.abs(x-loc))

#            #
            y = binned_y[:,loc_bin]/norm # already checking for divide by zero in times above
#            y = [norm == 0.0] = None
            ax[axind].plot(times, y, lw = line_width, color = plasma(loc),
                                label = r"%0.1f R$_{\rm vir}$"%(loc))

        ax[axind].set_xlabel(r'Time (Myr)')
        if  mass_loading:
            ax[axind].set_ylabel(r' ' + fields_labels[field] + ' Mass Loading Factor')
        else:
            ax[axind].set_ylabel(r' ' + fields_labels[field] + ' Outflow Rate (M$_{\odot}$ yr$^{-1}$')

        ax[axind].semilogy()

        if axi * nrow + axj < n_non_metal:
            if mass_loading:
                ax[axind].set_ylim(1.0E-6, 1000.0)
            else:
                ax[axind].set_ylim(1.0E-10, 1.0)
        else:
            if mass_loading:
                ax[axind].set_ylim(1.0E-8, 100.0)
            else:
                ax[axind].set_ylim(1.0E-12, 1.0E-2)

        axj = axj + 1
        if axj >= ncol:
            axj = 0
            axi = axi + 1



    ax[(0,0)].legend(loc='best')

    fig.set_size_inches(6*nrow,6*ncol)
    plt.tight_layout()
    plt.minorticks_on()
    if mass_loading:
        outname = 'mass_loading_factor_'
    else:
        outname = 'mass_outflow_rate_'
    fig.savefig(outname + 'panel_plot_evolution.png')
    plt.close()

    return



if __name__ == "__main__":

    plot_outflow_panel(tmin = 0.0, tmax = np.inf)
    plot_outflow_panel(tmin = 0.0, tmax = np.inf, mass_loading = True)
