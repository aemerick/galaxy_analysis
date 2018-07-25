from galaxy_analysis.plot.plot_styles import *
import matplotlib.pyplot as plt

import deepdish as dd
from galaxy_analysis.utilities import utilities
import numpy as np
import sys
from scipy.interpolate import interp1d
rc('font',size=22)

TMAX = 500.0




def obtain_mprod(data_list, fieldname):
    all_data = {}

    for i,k in enumerate(data_list):
        x = dd.io.load(k, '/gas_meta_data/masses')
        all_data[k] = x['FullBox'][fieldname] + x['OutsideBox'][fieldname]

    all_data = np.array([all_data[k] for k in np.sort(all_data.keys())])

    return all_data

def obtain_outflow_rates(data_list, outflow_field, elements):
    _temp_dict = {}
    all_data = {}

    for i,k in enumerate(data_list):
        _temp_dict[k] = dd.io.load(k,
                             '/gas_profiles/outflow/sphere')
        if outflow_field == 'Total Tracked Metals':
            all_data[k] = np.zeros(np.shape(_temp_dict[k][ ('gas', elements[0] + '_Mass')]))
            for x in elements:
                all_data[k] = all_data[k] + np.array( _temp_dict[k][ ('gas', x + '_Mass')])
        else:
            all_data[k] =  _temp_dict[k][ ('gas',outflow_field) ] # Becuase i'm dumb

    all_data =np.array( [all_data[k] for k in np.sort(all_data.keys())])
    print np.shape(all_data)
    return all_data

def obtain_stellar_mass(data_list):
    m = np.ones(np.size(data_list))
    for i,k in enumerate(data_list):
        m[i] = dd.io.load(k, '/particle_meta_data/total_birth_mass')
    return m

def obtain_metal_mass(data_list):
    m = np.ones(np.size(data_list))
    for i,k in enumerate(data_list):
        m[i] = dd.io.load(k, '/gas_meta_data/masses/FullBox')['Total Tracked Metals'] +\
               dd.io.load(k, '/gas_meta_data/masses/OutsideBox')['Total Tracked Metals']
    return m

def obtain_sfr(data_list, times, smooth_sfr = True):
    sfr = np.ones(np.size(data_list))

    for i,k in enumerate(data_list):
        sfr[i] = dd.io.load(k, '/meta_data/SFR')

    if smooth_sfr:
        print 'WARNING: Smoothing SFR - Assuming 1 Myr even samples in time'
        # smooth into 100 Myr bins, centered on each
        # sample point
        di = 250
        old_sfr = 1.0 * sfr
        for i in np.arange(np.size(sfr)):
            imin = np.max( [0, i - di])
            imax = np.min( [np.size(sfr), i + di - 1])

            if imin == 0:
                imax = 100
            if imax == np.size(sfr):
                imin = imax - di*2

            sfr[i] = np.average( sfr[imin:imax] )

    if True:
        for i in np.arange(np.size(sfr)):
            if sfr[i] > 0.0:
                break

        sfr[:i] = sfr[i]
        f = interp1d(times[sfr>0], sfr[sfr>0])
        sfr[sfr == 0.0] = f(times[sfr==0.0])

    return sfr

def plot_species_outflow_panel(work_dir = './', t_min = 0.0, t_max = 1000.0,
                               method = 'fraction', outdir = './'):
    """
    Default to plotting the following:
        For each species, x, plots   (dM_x/dt / M_x_prod) / SFR
        or the fractional mass outflow rate per unit SFR.
    """

    data_list, times = utilities.select_data_by_time(dir = work_dir,
                                                     tmin=0.0,tmax=1000)

    xpos = dd.io.load(data_list[0], '/gas_profiles/outflow/sphere')
    xpos = xpos['centers_rvir']
    temp = dd.io.load(data_list[0], '/gas_meta_data/masses/CNM')
    elements = utilities.sort_by_anum([x for x in temp.keys() if len(x) <= 2 and (not (x in ['H','He','H2','HI','HII','HeI']))])

    fig, ax = plt.subplots(4, 4, sharex=True, sharey=True)
    fig.set_size_inches(16,16)
    fig.subplots_adjust(hspace = 0.0, wspace = 0.0)

    SFR = obtain_sfr(data_list, times)

    axi,axj = 0,0

    for e in elements:
        index = (axi,axj)

        field_name = e + '_Mass'
        all_data = obtain_outflow_rates(data_list, field_name, elements)
        M_prod   = obtain_mprod(data_list, e)

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
        ax[(i,0)].set_xlim(newx[0]-newx[0], np.min([TMAX,newx[-1]-newx[0]]))

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
        fig.savefig(outdir + 'X_Fractional_Outflow_panel.png')
    else:
        fig.savefig(outdir + 'X_Fractional_Outflow_loading_panel.png')
    plt.close()

    return

def plot_basic_outflow_and_loading(work_dir = './', t_min = 0.0, t_max = 1000.0,
                                   outdir = './'):

    data_list, times = utilities.select_data_by_time(dir = work_dir,
                                                     tmin=t_min,tmax=t_max)

    xpos = dd.io.load(data_list[0], '/gas_profiles/outflow/sphere')
    xpos = xpos['centers_rvir']
    temp = dd.io.load(data_list[0], '/gas_meta_data/masses/CNM')
    elements = utilities.sort_by_anum([x for x in temp.keys() if len(x) <= 2 and (not (x in ['H','He','H2','HI','HII','HeI']))])

    #
    # Mass Outflow Plot
    #
    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    all_data   = obtain_outflow_rates(data_list, 'cell_mass', elements)
    metal_data = obtain_outflow_rates(data_list, 'Total Tracked Metals', elements)
    metal_mass = obtain_metal_mass(data_list)
    stellar_mass = obtain_stellar_mass(data_list) # total mass in stars produced at a time, NOT M_* of galaxy
    sfr        = obtain_sfr(data_list, times)

    binned_y = 1.0 * all_data # np.array( [all_data[k] for k in all_data.keys()] )
    binned_metal_mass = 1.0 * metal_data
    for i,loc in enumerate([0.1, 0.25, 0.5, 1.0]):
        loc_bin = np.argmin(np.abs(xpos-loc)) # get the right position bin
        y = binned_y[:,loc_bin]
        print i, loc, loc_bin
        bin_edges=utilities.bin_edges(np.round(times))
        # rebin with 10 Myr bins, rather than previous 1 Myr bins
        newx,rebiny=utilities.simple_rebin(1.0*bin_edges,1.0*y,
                       np.arange(np.min(bin_edges),np.max(bin_edges)+2,10), 'average')


        plot_histogram(ax, newx-newx[0], rebiny, lw = line_width, color = plasma(i/4.0),
                       label = r"%0.1f R$_{\rm vir}$"%(loc))

    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(r'Outflow Rate (M$_{\odot}$ yr$^{-1}$)')
    ax.semilogy()
    ax.set_xlim(0.0, np.min([TMAX,np.max(newx-newx[0])]))
    ax.set_ylim(1.0E-5, 0.03)

    plt.tight_layout()
    ax.legend(loc='best')
    plt.minorticks_on()
    fig.savefig(outdir + 'total_mass_outflow.png')
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
    ax.set_xlim(0.0, np.min([TMAX,np.max(newx-newx[0])]))
    ax.set_ylim(0.1,300)

    plt.tight_layout()
    #ax.legend(loc='best')
    plt.minorticks_on()
    fig.savefig(outdir + 'total_mass_loading.png')
    plt.close()

    #
    # Metal Mass Loading
    #
    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    for i, loc in enumerate([0.1,0.25,0.5,1.0]):
        loc_bin = np.argmin(np.abs(xpos-loc)) # get the right position bin
        print np.shape(binned_y), np.shape(binned_metal_mass), np.size(sfr), np.size(metal_mass), np.size(stellar_mass)
        y = binned_metal_mass[:,loc_bin] / (sfr  * (metal_mass / stellar_mass)) # metal mass loading factor
        print loc, np.min(y), np.max(y), np.average(y)
        print np.min(metal_mass), np.max(metal_mass), np.min(stellar_mass),np.max(stellar_mass)

        bin_edges=utilities.bin_edges(np.round(times))
        # rebin with 10 Myr bins, rather than previous 1 Myr bins
        newx,rebiny=utilities.simple_rebin(1.0*bin_edges,1.0*y,
                       np.arange(np.min(bin_edges),np.max(bin_edges)+2,10), 'average')

        plot_histogram(ax, newx-newx[0], rebiny, lw = line_width, color = plasma(i/4.0),
                       label = r"%0.1f R$_{\rm vir}$"%(loc))

    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(r'Metal Mass Loading Factor')
    ax.semilogy()
    ax.set_xlim(0.0, np.min([TMAX,np.max(newx-newx[0])]))
    ax.set_ylim(0.04, 20)

    plt.tight_layout()
    ax.legend(loc='upper right')
    plt.minorticks_on()
    fig.savefig(outdir + 'metal_mass_loading.png')
    plt.close()

    #
    # Metal Mass Loading
    #
    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    for i, loc in enumerate([0.1,0.25,0.5,1.0]):
        loc_bin = np.argmin(np.abs(xpos-loc)) # get the right position bin
        print np.shape(binned_y), np.shape(binned_metal_mass), np.size(sfr), np.size(metal_mass), np.size(stellar_mass)
        y = binned_metal_mass[:,loc_bin] / (sfr) # metal mass loading factor
        print loc, np.min(y), np.max(y), np.average(y)
        print np.min(metal_mass), np.max(metal_mass), np.min(stellar_mass),np.max(stellar_mass)

        bin_edges=utilities.bin_edges(np.round(times))
        # rebin with 10 Myr bins, rather than previous 1 Myr bins
        newx,rebiny=utilities.simple_rebin(1.0*bin_edges,1.0*y,
                       np.arange(np.min(bin_edges),np.max(bin_edges)+2,10), 'average')


        print np.sum(rebiny) / (1.0*np.size(rebiny)), '-------------'
        plot_histogram(ax, newx-newx[0], rebiny, lw = line_width, color = plasma(i/4.0),
                       label = r"%0.1f R$_{\rm vir}$"%(loc))

    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(r'Metal Mass Loading Factor')
    ax.semilogy()
    ax.set_xlim(0.0, np.min([TMAX,np.max(newx-newx[0])] ))
    ax.set_ylim(1.0E-6, 1.0)

    plt.tight_layout()
    ax.legend(loc='best')
    plt.minorticks_on()
    fig.savefig(outdir + 'metal_mass_loading_sfr.png')
    plt.close()

    return


if __name__ == "__main__":
    work_dir      = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11_30km/final_sndriving/'

#    plot_species_outflow_panel()
    if len(sys.argv) > 1:
        work_dir = sys.argv[1]

    plot_basic_outflow_and_loading(work_dir = work_dir)
