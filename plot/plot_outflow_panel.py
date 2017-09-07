from galaxy_analysis.plot.plot_styles import *
import matplotlib.pyplot as plt
"-----------------------------------------"

import deepdish as dd
from galaxy_analysis.utilities import utilities
import numpy as np


def plot_outflow_panel(dir = '.', tmin = None, tmax = None, mass_loading = False):

    # plot in this order:
    fields        = ["cell_mass", "H_total_mass", "H2_mass", "metal_mass"]
    fields_labels = {'cell_mass' : 'Total', 'H_total_mass' : "H",
                     'H2_mass' : r'H$_{2}$', "metal_mass" : 'Total Metal'}

    # obtain ds list based off time range selection
    data_list, times = utilities.select_data_by_time(tmin=tmin, tmax = tmax)

    # this is a hacky way of getting the metal fields - should just save this to file
    metal_fields = dd.io.load(data_list[0], 'gas_meta_data/masses/FullBox')
    exclude      = ['H','HI','HII','H2','Metals','Total','He','HeI','HeII','HeIII']
    metal_fields = [x for x in metal_fields if (not any([y in x for y in exclude]))]

    for k in metal_fields:
        fields_labels[k] = k
    fields = fields + metal_fields

    nplots = len(plot_fields)
    nrow, ncol = utilities.rowcoldict[nplots]

    key_select = ['gas_profiles', 'outflow', 'sphere']

    x = dd.io.load(data_list[0], 'gas_profiles/outflow/sphere/centers_rvir')

    # load everything to memory at start - limits number of reads from disk
    all_data_dict = {}
    for k in data_list:
        all_data_dict[k] = dd.io.load(k, 'gas_profiles/outflow/sphere')

    axi, axj = 0,0
    for field in fields:
        axind = (axi, axj)

        load_field = ('gas',field) # because I'm dumb and keyed them like this
        binned_y   = utilities.extract_nested_dict_asarray(all_data_dict, load_field)
        norm       = 1.0
        if mass_loading:
            norm = 1.0 # need to set SFR here

        # plot at 0.1, 0.25, 0.5, and 1 Rvir for now:
        for loc in [0.1, 0.25, 0.5, 1.0]:
            loc_bin = np.argmin(np.abs(x-loc))
            ax[axind].plot(x, binned_y[:,loc_bin]/norm, lw = line_width, color = plasma(loc),
                                label = r"%0.1f R$_{\rm vir}$"%(loc))

        ax[axind].set_xlabel(r'Time (Myr)')
        if  mass_loading:
            ax[axind].set_ylabel(r' ' + fields_labels[field] + ' Mass Loading Factor')
        else:
            ax[axind].set_ylabel(r' ' + fields_labels[field] + ' Outflow Rate (M$_{\odot}$ yr$^{-1}$')

        ax[axind].semilogy()

    ax[(0,0)].legend(loc='best')

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

    plot_outfow_panel(tmin = 0.0, tmax = np.inf)
