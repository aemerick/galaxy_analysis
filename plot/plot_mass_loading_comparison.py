from galaxy_analysis.plot.plot_styles import *
from galaxy_analysis.utilities import utilities
#----------------------------------------------
rc('font',size=22)

import matplotlib.pyplot as plt
import numpy as np
import glob as glob
import deepdish as dd
import sys

TMAX = 500.0

line_width = 3.0


# would be nice to start making gather functions
# for all of these plot functions to not have to
# do any more looping over ALL data sets to gather
def plot_mass_loading_comparison(work_dir = './', output_dir = None,
                                 comparison = None, color = None, rbin = 1, ylim = None):

    if output_dir is None:
        output_dir = work_dir

    if comparison is None:
        labels = {'3pcH2' : '3.6 pc' , '6pcH2' : '7.2 pc', 'Fiducial' : 'Fiducial'}
        lstyle = {'3pcH2' : '--', '6pcH2' : '-.', 'Fiducial' : '-'}
        dirs   = {'3pcH2' : '../3pc_H2/' , '6pcH2' : '../6pc_H2/', 'Fiducial' : work_dir}

    else:
        dirs = {}
        labels = {}
        lstyle = {}
        for k in list(comparison.keys()):
            dirs[k]   = work_dir + comparison[k][0]
            labels[k] = comparison[k][1]
            lstyle[k] = comparison[k][2]

    if color is None:
        color = {}
        for i,k in enumerate(comparison.keys()):
            color[k] = 'C%1i'%(i)

    gather_keys = {'mass_outflow' : ['gas_profiles','outflow','sphere']}

    all_data = {}

    for sim in list(comparison.keys()):
        data_list, times = utilities.select_data_by_time(dir = dirs[sim],
                                                         tmin=0.0,tmax= 1000.0)
        all_data[sim] = {}
        all_data[sim]['times'] = times
        for k in list(gather_keys.keys()):

            all_data[sim][k] = utilities.extract_nested_dict_aslist(None, gather_keys[k], data_list, False)

        all_data[sim]['SFR'] = utilities.extract_nested_dict_asarray(None, ['time_data','SFR_1'], [data_list[-1]], False)[0]
        all_data[sim]['SFR_time'] = utilities.extract_nested_dict_asarray(None, ['time_data','time_1'], [data_list[-1]], False)[0]
        # take instantaneous SFR data and average over 100 Myr window - generate fit to make SFR continuous
        sfr_average = np.zeros(np.size(all_data[sim]['SFR_time']))
        for i,sfrval in enumerate(all_data[sim]['SFR_time']):
            t_o = sfrval - 50.0E6
            t_max = sfrval + 50.0E6
            sfr_average[i] = np.average(  all_data[sim]['SFR'][  (all_data[sim]['SFR_time'] > t_o) * (all_data[sim]['SFR_time'] < t_max)])
            #all_data[sim]['SFR_time']-all_data[sim]['SFR_time'][0]
        all_data[sim]['sfr_average'] = 1.0 * sfr_average
        print("---", np.size(all_data[sim]['SFR']), np.size(all_data[sim]['SFR_time']))
        print("---", np.shape(all_data[sim]['SFR_time']), np.shape(sfr_average))
        all_data[sim]['SFR_fit'] = lambda x: np.interp(x,
                                                      (all_data[sim]['SFR_time'] - all_data[sim]['SFR_time'][0])/1.0E6 ,
                                                       all_data[sim]['sfr_average']
                                                       )


    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    for sim in list(comparison.keys()):
        #t = all_data[sim]['times'] - all_data[sim]['times'][0]
        t = all_data[sim]['times'] - all_data[sim]['times'][0]

        outflow  = all_data[sim]['mass_outflow'] #[('gas','cell_mass')]
        Mdot     = np.array([x[('gas','cell_mass')] for x in outflow])
        SFR_func      = all_data[sim]['SFR_fit']
        yplot = Mdot[:,rbin]# / SFR_func(t)

        bin_edges=utilities.bin_edges(np.round(t*1.0))
# rebin with 10 Myr bins, rather than previous 1 Myr bins
        newx,rebiny=utilities.simple_rebin(1.0*bin_edges,1.0 * Mdot[:,rbin],
                    np.arange(np.min(bin_edges),np.max(bin_edges)+2,50), 'average')

#newx,rebiny=utilities.simple_rebin(1.0*bin_edges,1.0*y,
#               np.arange(np.min(bin_edges),np.max(bin_edges)+2,10), 'average')

        ax.plot(t-t[0],yplot, lw = line_width, ls = '-', color = color[sim]) #, label = labels[sim])
        ax.plot(t-t[0],Mdot[:,5], lw = line_width, ls = '--', color = color[sim])#, label = sim)
        #xc  = 0.5  *(newx[1:] + newx[:-1])
        #ax.plot(xc-xc[0], rebiny, lw = line_width,  ls = ':', color = color[sim], label = sim)
        #plot_histogram(ax, newx-newx[0], rebiny, lw = line_width,  ls = ':', color = color[sim])#, label = sim)

        print(np.size(Mdot[:,2]), np.size(rebiny))

    # plot this so it shows up on diagram
    ax.plot([-10,-1],[0,0],lw=line_width,ls='-',color='black', label = r'0.25 R$_{\rm vir}$')
    ax.plot([-10,-1],[0,0],lw=line_width,ls='--',color='black', label = r'1.0 R$_{\rm vir}$')

    ax.semilogy()
    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(r"Mass Outflow Rate [M$_{\odot}$ yr$^{-1}$]")
    # Loading Factor at 0.25 R$_{\rm vir}$')
    ax.set_xlim(0.0, TMAX)
    if not (ylim is None):
        ax.set_ylim(ylim)
    ax.legend(loc = 'best')
    plt.minorticks_on()
    plt.tight_layout()
    fig.savefig(work_dir + output_dir + 'mass_outflow_comparison.png')
    plt.close()
#############################
#############################
# do same but with mass loading factor

    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    for sim in list(comparison.keys()):
        #t = all_data[sim]['times'] - all_data[sim]['times'][0]
        t = all_data[sim]['times'] - all_data[sim]['times'][0]

        outflow  = all_data[sim]['mass_outflow'] #[('gas','cell_mass')]
        Mdot     = np.array([x[('gas','cell_mass')] for x in outflow])
        SFR_func      = all_data[sim]['SFR_fit']
        print(np.size(Mdot[:,rbin]), np.size(t))
        yplot = Mdot[:,rbin] / SFR_func(t)

        bin_edges=utilities.bin_edges(np.round(t*1.0))
# rebin with 10 Myr bins, rather than previous 1 Myr bins
        newx,rebiny=utilities.simple_rebin(1.0*bin_edges,1.0 * Mdot[:,rbin],
                    np.arange(np.min(bin_edges),np.max(bin_edges)+2,50), 'average')

#newx,rebiny=utilities.simple_rebin(1.0*bin_edges,1.0*y,
#               np.arange(np.min(bin_edges),np.max(bin_edges)+2,10), 'average')

        ax.plot(t-t[0],yplot, lw = line_width, ls = '-', color = color[sim], label = labels[sim])
        ax.plot(t-t[0],Mdot[:,5]/SFR_func(t), lw = line_width, ls = '--', color = color[sim])#, label = sim)
        #xc  = 0.5  *(newx[1:] + newx[:-1])
        #ax.plot(xc-xc[0], rebiny, lw = line_width,  ls = ':', color = color[sim], label = sim)
        #plot_histogram(ax, newx-newx[0], rebiny, lw = line_width,  ls = ':', color = color[sim])#, label = sim)

    ax.semilogy()
    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(r"Mass Loading Factor")
    # Loading Factor at 0.25 R$_{\rm vir}$')
    ax.set_xlim(0.0, TMAX)
    #if not (ylim is None):
    ax.set_ylim(0.01,1000.0)
    #ax.legend(loc = 'best')
    plt.minorticks_on()
    plt.tight_layout()
    fig.savefig(work_dir + output_dir + 'mass_loading_comparison.png')
    plt.close()


    return

if __name__ == "__main__":

    work_dir = './'
    if len(sys.argv) > 1:
        work_dir = sys.argv[1]
    output_dir = None
    if len(sys.argv) > 2:
        output_dir = sys.argv[2]

    plot_mass_loading_comparison(work_dir = work_dir, output_dir = output_dir)
