import deepdish as dd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import glob
import numpy as np


workdir = '/mnt/ceph/users/emerick/enzo_runs/stampede/leo_p/fullphysics/fullres'

DATA_PATHS = { 'lm'        : workdir + '/run11/200cc',
               'lm_noRT'   : workdir + '/run11/200cc_noRT',
               'lm_nowind' : workdir + '/run11/200cc_nowind',
               'lm_xx'     : workdir + '/run11/othercloudy',
               'mm'        : workdir + '/run13',
               'mm_3pc'    : workdir + '/run13/3parsec',
               'hm'        : workdir + '/run15'}

ls_dict = { 'lm' : '-', 'lm_noRT' : '--', 'lm_nowind' : ':', 'lm_xx' : '--',
             'mm' : '-', 'mm_3pc' : '-.', 'hm' : '-'}


purple  = '#7d03a8'
magenta = '#cb4679'
blue    = '#0c0887'
orange  = '#fdc328'

color_dict = {'lm' : purple, 'lm_noRT' : purple, 'lm_nowind' : purple,
              'mm' : magenta  , 'mm_3pc'  : magenta  , 'hm' : orange, 'lm_xx' : blue}


linewidth = 4.5

ALL_DATA = {}
for s in DATA_PATHS.keys():
    ALL_DATA[s] = np.sort(glob.glob(DATA_PATHS[s] + '/DD*.h5'))

def time_first_star(data = None, t = None, sfr = None):
    if data is None and (t is None or sfr is None):
        print "If supplying no data set, must supply both time array and sfr array"
        raise ValueError

    if not data is None:
        t   = data['time_data']['time']
        sfr = data['time_data']['SFR']

    t_first = np.min( t[sfr > 0])

    return t_first

def plot_mass_loading(sim_names = None, species = 'total', z = 100.0, mass_loading = False):
    """
    Given a dictionary that contains list of simulation names and filepaths,
    go through all simulation data analysis outputs in that path for each
    simulation and plot the time evolution of the outflow rate. 

    Inputs
    ------ 
        sim_names    : dictionary, optional
        species      : String. Elemental species (or 'total' or 'metal' for all gas or
                       all metals) of which to plot the outflow rate. optional. default: total
        z            : Fixed height above / below disk to plot the rate (in pc). Default 100
        mass_loading : Normalize by SFR to get the mass loading factor instead. Default False


    """

    if sim_names is None:
        sim_names = DATA_PATHS.keys()

    # do some checking to match desired species with dict names
    if species == 'total':
        species = 'total'
    elif species == 'H':
        species = 'H_total_mass'
    elif species == 'metal' or species == 'Metal' or species == 'Metal_Mass':
        species = 'metal_mass'
    elif species == 'He':
        species = 'He_total_mass'
    elif not '_Mass' in species:
        species = species + '_Mass'


    fig, ax = plt.subplots()

    # loop over simulations
    for s in sim_names:
        # always load most recent file to check
        data = dd.io.load(ALL_DATA[s][-1])
        t_first = time_first_star(data) / 1.0E6 # now in Myr

        # make arrays for plotting
        t  = np.zeros(np.size(ALL_DATA[s]))
        ml = np.zeros(np.size(ALL_DATA[s]))

        # now go through every data analysis dump
        for i,x in enumerate(ALL_DATA[s]):
            xdata = dd.io.load(ALL_DATA[i], '/gas_profiles/outflow_rate')
            t[i]  = dd.io.load(ALL_DATA[i], '/meta_data/Time')

            # find the bin whose center is closest to the desired z
            zbin  = np.argmin( np.abs( (xdata['zbins'][1:] + xdata['zbins'][:-1])*0.5 - z) )

            if species == 'total':
                ml[i] = xdata['H_total_mass'][zbin] + xdata['He_total_mass'][zbin] + xdata['metal_mass'][zbin]
            else:
                ml[i] = xdata[species][zbin]

        # normalize t and plot only after SF occurs
        t         = t - t_first
        t_plot    = t[t >= 0.0]
        mass_plot = mass[t >= 0.0]

        norm = 1.0
        if mass_loading: # maybe compute as the ~ 20 Myr average on either side of the data point?
            norm = 1.0 / SFR  

        ax.plot(t_plot, mass_plot*norm, lw = linewidth, ls = ls_dict[s], color = color_dict[s], label = s, drawstyle = 'steps-post')

    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(species + r' Mass Outflow Rate (M$_{\odot}$)')
    ax.semilogy()

    ax.set_xlim(0.0, 300.0)
    ax.legend(loc='best')

    fig.set_size_inches(8,8)
    plt.tight_layout()
    plt.minorticks_on()

    plottype = 'outflow'
    if mass_loading:
        plottype = 'loading'

    outname = species + '_mass_' + plottype + '_z%2f.png'%(z)

    fig.savefig(outname)
    plt.close()

    return            

def plot_mass(sim_names = None, species = 'HI'):


    if sim_names is None:
        sim_names = DATA_PATHS.keys()

    fig,ax = plt.subplots()

    for s in sim_names:
        # always load most recent file in every case
        data = dd.io.load(ALL_DATA[s][-1])

        t_first = time_first_star(data)

        if species == 'HI' or species == 'M_HI':
            fname = 'M_HI'
        elif species == 'total' or species == 'Total' or species == 'gas' or species == 'M_gas':
            fname = 'M_gas'
        elif species == 'star' or species == 'stars' or species == 'stellar' or species == 'M_star':
            fname = 'M_star'

        t = np.zeros(np.size(ALL_DATA[s]))
        mass = np.zeros(np.size(ALL_DATA[s]))
        for i,x in enumerate(ALL_DATA[s]):
            xdata = dd.io.load(x)
            t[i] = xdata['meta_data']['Time']
            mass[i] = xdata['meta_data'][fname]
        
        t = t - t_first
        t_plot   = t[t >= 0.0]
        mass_plot = mass[t >= 0.0]

        ax.plot(t_plot, mass_plot, lw = linewidth, ls = ls_dict[s], color = color_dict[s], label = s, drawstyle = 'steps-post')

    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(species + r' Mass (M$_{\odot}$)')
    ax.semilogy()

    ax.set_xlim(0.0,200.0)
    ax.legend(loc='best')

    fig.set_size_inches(8,8)
    plt.tight_layout()
    plt.minorticks_on()

    fig.savefig(species + '_mass_comparison.png')
    plt.close()

    return

def plot_sfr(sim_names = None):

    if sim_names is None:
        sim_names = DATA_PATHS.keys()

    fig,ax = plt.subplots()

    for s in sim_names:
        # always load most recent file in every case
        data = dd.io.load(ALL_DATA[s][-1])

        t    = data['time_data']['time']
        sfr  = data['time_data']['SFR']

        t_first = np.min( t[sfr > 0])

        t = t - t_first

        t_plot   = t[t >= 0.0] / 1.0E6
        sfr_plot = sfr[t >= 0.0]

        ax.plot(t_plot, sfr_plot, lw = linewidth, ls = ls_dict[s], color = color_dict[s], label = s, drawstyle = 'steps-post')

    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(r'SFR (M$_{\odot}$ yr$^{-1}$)')
    ax.set_xlim(0.0, 75.0)
    ax.semilogy()

    ax.legend(loc='best')

    fig.set_size_inches(8,8)
    plt.tight_layout()
    plt.minorticks_on()

    fig.savefig('sfr_comparison.png')
    plt.close()

    return


def plot_snr(sim_names = None):

    if sim_names is None:
        sim_names = DATA_PATHS.keys()

    fig,ax = plt.subplots()

    for s in sim_names:
        # always load most recent file in every case
        data = dd.io.load(ALL_DATA[s][-1])

        t    = data['time_data']['time']
        sfr  = data['time_data']['SFR']
        snr  = data['time_data']['SNII_snr']

        t_first = np.min( t[sfr > 0])

        t = t - t_first

        t_plot   = t[t >= 0.0] / 1.0E6
        snr_plot = snr[t >= 0.0]

        ax.plot(t_plot, snr_plot, lw = linewidth, ls = ls_dict[s], color = color_dict[s], label = s, drawstyle = 'steps-post')

    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(r'SNR (yr$^{-1}$)')
    ax.set_xlim(0.0,100.0)
    ax.semilogy()
    plt.minorticks_on()

    ax.legend(loc='best')

    fig.set_size_inches(8,8)
    plt.tight_layout()

    fig.savefig('snr_comparison.png')
    plt.close()

    return


if __name__ == '__main__':

    all_s = ['lm','lm_noRT','lm_nowind','mm','mm_3pc','hm', 'lm_xx']

    for species in ['total', 'metals', 'C', 'Fe', 'H']:

        plot_mass_loading(sim_names = all_s, species = species, z = 100)
        plot_mass_loading(sim_names = all_s, species = species, z = 500)

    plot_snr(sim_names = all_s) 
    plot_sfr(sim_names = all_s)
    plot_mass(sim_names = all_s, species = 'HI')
    plot_mass(sim_names = all_s, species = 'stars')
    plot_mass(sim_names = all_s, species = 'total')
    
