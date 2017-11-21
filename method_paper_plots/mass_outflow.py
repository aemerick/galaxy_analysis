from galaxy_analysis.plot.plot_styles import *
import matplotlib.pyplot as plt

import deepdish as dd
from galaxy_analysis.utilities import utilities
import numpy as np

from scipy.interpolate import interp1d

outflow_field = 'cell_mass'
work_dir      = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11_30km/final_sndriving/'
#work_dir      = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11/corrected_sndriving/'


data_list, times = utilities.select_data_by_time(dir = work_dir,
                                                 tmin=0.0,tmax=1000)

xpos = dd.io.load(data_list[0], '/gas_profiles/outflow/sphere')
xpos = xpos['centers_rvir']
temp = dd.io.load(data_list[0], '/gas_meta_data/masses/CNM')
elements = utilities.sort_by_anum([x for x in temp.keys() if len(x) <= 2 and (not (x in ['H','He','H2','HI','HII','HeI']))])


def obtain_mprod(fieldname):
    all_data = {}

    for i,k in enumerate(data_list):
        x = dd.io.load(k, '/gas_meta_data/masses')
        all_data[k] = x['FullBox'][fieldname] + x['OutsideBox'][fieldname]

    all_data = np.array([all_data[k] for k in all_data.keys()])

    return all_data

def obtain_outflow_rates(outflow_field):
    _temp_dict = {}
    all_data = {}
    for i,k in enumerate(data_list):
        _temp_dict[k] = dd.io.load(k,
                             '/gas_profiles/outflow/sphere')
        all_data[k] =  _temp_dict[k][ ('gas',outflow_field) ] # Becuase i'm dumb

    all_data =np.array( [all_data[k] for k in all_data.keys()])
    return all_data

def obtain_sfr(smooth_sfr = True):
    sfr = np.ones(np.size(data_list))

    for i,k in enumerate(data_list):
        sfr[i] = dd.io.load(k, '/meta_data/SFR')

    if smooth_sfr:
        for i in np.arange(np.size(sfr)):
            if sfr[i] > 0.0:
                break

        sfr[:i] = sfr[i]
        f = interp1d(times[sfr>0], sfr[sfr>0])
        sfr[sfr == 0.0] = f(times[sfr==0.0])

    return sfr

def plot_species_outflow_panel(method = 'fraction'):
    """
    Default to plotting the following:
        For each species, x, plots   (dM_x/dt / M_x_prod) / SFR
        or the fractional mass outflow rate per unit SFR.
    """
    fig, ax = plt.subplots(4, 4, sharex=True, sharey=True)
    fig.set_size_inches(16,16)
    fig.subplots_adjust(hspace = 0.0, wspace = 0.0)

    SFR = obtain_sfr()

    axi,axj = 0,0
    
    for e in elements:
        index = (axi,axj)

        field_name = e + '_Mass'
        all_data = obtain_outflow_rates(field_name)
        M_prod   = obtain_mprod(e)

        binned_y = 1.0 * all_data # np.array( [all_data[k] for k in all_data.keys()] )
        for i, loc in enumerate([0.1,0.25,0.5,1.0]):
            loc_bin = np.argmin(np.abs(xpos-loc)) # get the right position bin
            y = binned_y[:,loc_bin]

            bin_edges=utilities.bin_edges(np.round(times))
            # rebin with 10 Myr bins, rather than previous 1 Myr bins

            # convert y from outflow rate to outflow rate / mass produced / sfr
            if method == 'fraction':
                y = (y/M_prod) * 1.0E6 # convert to 1/Myr
            else:
                y = (y/M_prod)/ SFR

            newx,rebiny=utilities.simple_rebin(1.0*bin_edges,1.0*y,
                           np.arange(np.min(bin_edges),np.max(bin_edges)+2, 5.0), 'average')


            plot_histogram(ax[index], newx-newx[0], rebiny, lw = line_width, color = plasma(i/4.0),
                           label = r"%0.1f R$_{\rm vir}$"%(loc))

            if i == 0:
                xy = (np.max(newx-newx[0])*0.8, np.max(rebiny)*10.0)

        ax[index].annotate(e, xy=xy,xytext=xy)

        axj = axj + 1
        if axj >= 4:
            axj = 0
            axi = axi + 1

    for i in np.arange(4):
        ax[(3,i)].set_xlabel(r'Time (Myr)')
#        ax[(i,0)].set_ylabel(r'log(X Fractional Outflow per SFR [M$_\odot$ yr$^{-1}$]$^{-1}$)')
        ax[(i,0)].set_xlim(newx[0]-newx[0], newx[-1]-newx[0])

        if method == 'fraction':
            ax[(0,i)].set_ylim(1.0E-5, 1.0E-1)
        else:
            ax[(0,i)].set_ylim(1.0E-7, 1.0E-3)

        ax[(0,i)].semilogy()

    if method == 'fraction':
        ax[(2,0)].set_ylabel(r'log(X Fractional Outflow Rate [Myr$^{-1}$])')
    else:
        ax[(2,0)].set_ylabel(r'log(X Fractional Outflow per SFR [M$_\odot$ yr$^{-1}$]$^{-1}$)')

    plt.minorticks_on()
    if method == 'fraction':
        fig.savefig('X_Fractional_Outflow_panel.png')
    else:
        fig.savefig('X_Fractional_Outflow_loading_panel.png')
    plt.close()

    return

def basic_outflow_loading_plot():
    #
    # Mass Outflow Plot
    #
    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    all_data = obtain_outflow_rates('cell_mass')
    sfr      = obtain_sfr()

    binned_y = 1.0 * all_data # np.array( [all_data[k] for k in all_data.keys()] )
    for i,loc in enumerate([0.1, 0.25, 0.5, 1.0]):
        loc_bin = np.argmin(np.abs(xpos-loc)) # get the right position bin
        y = binned_y[:,loc_bin]

        bin_edges=utilities.bin_edges(np.round(times))
        # rebin with 10 Myr bins, rather than previous 1 Myr bins
        newx,rebiny=utilities.simple_rebin(1.0*bin_edges,1.0*y,
                       np.arange(np.min(bin_edges),np.max(bin_edges)+2,10), 'average')


        plot_histogram(ax, newx-newx[0], rebiny, lw = line_width, color = plasma(i/4.0),
                       label = r"%0.1f R$_{\rm vir}$"%(loc))

    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(r'Outflow Rate (M$_{\odot}$ yr$^{-1}$)')
    ax.semilogy()
    ax.set_xlim(0.0, np.max(newx-newx[0]))
    ax.set_ylim(7E-6, 0.01)

    plt.tight_layout()
    ax.legend(loc='best')
    plt.minorticks_on()
    fig.savefig('total_mass_outflow.png')
    plt.close()

    #
    # Mass Loading plot
    #
    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

#    binned_y = np.array( [all_data[k] for k in all_data.keys()] )
    for i,loc in enumerate([0.1, 0.25, 0.5, 1.0]):
        loc_bin = np.argmin(np.abs(xpos-loc)) # get the right position bin
        y = binned_y[:,loc_bin] / sfr

        bin_edges=utilities.bin_edges(np.round(times))
        # rebin with 10 Myr bins, rather than previous 1 Myr bins
        newx,rebiny=utilities.simple_rebin(1.0*bin_edges,1.0*y,
                       np.arange(np.min(bin_edges),np.max(bin_edges)+2,10), 'average')

        plot_histogram(ax, newx-newx[0], rebiny, lw = line_width, color = plasma(i/4.0),
                       label = r"%0.1f R$_{\rm vir}$"%(loc))

    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(r'Mass Loading Factor')
    ax.semilogy()
    ax.set_xlim(0.0, np.max(newx-newx[0]))
    ax.set_ylim(0.1,300)

    plt.tight_layout()
    #ax.legend(loc='best')
    plt.minorticks_on()
    fig.savefig('total_mass_loading.png')
    plt.close()

    return


if __name__ == "__main__":

    plot_species_outflow_panel()
#    plot_basic_outflow_and_loading()

