import deepdish as dd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import glob
import numpy as np

# plot style defines
from galaxy_analysis.plot import plot_styles as ps


global_tmin = 0.0
global_tmax = 110.0

#
#
#

workdir  = '/mnt/ceph/users/emerick/enzo_runs'
stampede  = workdir + '/stampede/leo_p/fullphysics/fullres'
pleiades  = workdir + '/pleiades/'
local     = workdir + '/leo_p/fullres'

workdir = '/mnt/ceph/users/emerick/enzo_runs/stampede/leo_p/fullphysics/fullres'


DATA_PATHS = { 'lm'        : stampede + '/run11/200cc',
               'lm_noRT'   : stampede + '/run11/200cc_noRT',
               'lm_nowind' : stampede + '/run11/200cc_nowind',
               'lm_xx'     : stampede + '/run11/othercloudy',
               'mm'        : stampede + '/run13',
               'mm_3pc'    : stampede + '/run13/3parsec',
               'hm'        : stampede + '/run15'}

PERT_DATA_PATHS = {
                   'lm_p'      : local    + '/run11/200cc/perturb',
                   'lm_p_noOT' : local    + '/run11/200cc/perturb_iononly',
                   'lm_p_sn'   : stampede + '/run11/200cc/perturb_snonly',
                   'mm_p'      : pleiades + '/run13/no_wind',
                   'hm_p'      : pleiades + '/run15/no_wind',
                   'hm_p_sn'   : pleiades + '/run15/sn_only',
                   'hm_p_noOT' : pleiades + '/run15/no_otrad',
                   'hm_p_noRT' : pleiades + '/run15/no_ion',
                   'hm_p_3pc'  : pleiades + '/run15/3parsec',
                   'vhm_pert'  : stampede + '/run21/perturb',
                   'vhm_nopert': stampede + '/run21'}

STAR_IC = {
           'lm_ICs'     : pleiades + '/starIC/run11/lowsf',
#           'lm_ICs_3pc' : pleiades + '/starIC/run11/lowsf/3parsec',
           'mm_ICs'     : pleiades + '/starIC/run13/lowsf',
#           'mm_ICs_3pc' : pleiades + '/starIC/run13/lowsf/3parsec',
           'hm_ICs'     : pleiades + '/starIC/run15/lowsf',
           'hm_IC_compact' : pleiades + '/starIC/run15_compact/no_wind',
           'hm_p_noRT' : pleiades + '/starIC/run15/no_RT',
           'hm_p_noOT' : pleiades + '/starIC/run15/no_otrad',
           'hm_p_sn'   : pleiades + '/starIC/run15/sn_only',
           'hm_IC_compact_msf' : pleiades + '/starIC/run15_compact/msf',
           'hm_IC_compact_hsf' : pleiades + '/starIC/run15_compact/highsf'}
#           'hm_ICs_3pc' : pleiades + '/starIC/run15/lowsf/3parsec',
#           'lm_nostar'  : local + '/run11/200cc/perturb',
#           'lm_nopert'  : stampede + '/run11/200cc',
#           'mm_nostar'  : pleiades + '/run13/no_wind',
#           'mm_nopert'  : stampede + '/run13',
#           'hm_nostar'  : pleiades + '/run15/no_wind'}
#           'hm_nopert'  : stampede + '/run15'}



ls_dict = { 'lm' : '-', 'lm_noRT' : '--', 'lm_nowind' : ':', 'lm_xx' : '--',
             'mm' : '-', 'mm_3pc' : '-.', 'hm' : '-',
             'lm_p' : '-', 'lm_p_noOT' : '-.', 'lm_p_sn' : ':',
             'hm_p' : '-', 'hm_p_sn': ':', 'hm_p_noOT' : '-.', 'hm_p_noRT' : '--',
             'vhm_pert' : '-', 'vhm_nopert' : ':', 'hm_p_3pc' : '-', 'mm_p' : '-'}

add_to_ls = {'lm_ICs' : '-', 'mm_ICs' : '-', 'hm_ICs' : '-',
           'lm_nostar' : '--', 'mm_nostar' : '--', 'hm_nostar' : '--',
           'lm_nopert' : '-.', 'mm_nopert' : '-.', 'hm_nopert' : '-.',
           'lm_ICs_3pc' : ':', 'mm_ICs_3pc' : ':', 'hm_ICs_3pc' : ':', 'hm_IC_compact' : '-',
           'hm_IC_compact_msf' : '--', 'hm_IC_compact_hsf' : '-.'}

for k in add_to_ls:
    ls_dict[k] = add_to_ls[k]

color_dict = {'lm' : ps.purple, 'lm_noRT' : ps.purple, 'lm_nowind' : ps.purple,
              'mm' : ps.magenta  , 'mm_3pc' : ps.magenta  , 'hm' : ps.orange, 'lm_xx' : ps.blue}

for k in PERT_DATA_PATHS.keys() + STAR_IC.keys():

    if '3pc' in k:
        color_dict[k] = 'black'
    elif '_compact' in k:
        color_dict[k] = ps.blue
    elif 'lm' in k:
        color_dict[k] = ps.purple
    elif 'mm' in k:
        color_dict[k] = ps.magenta
    elif 'vhm' in k:
        color_dict[k] = ps.blue
    elif 'hm' in k:
        color_dict[k] = ps.orange


ALL_DATA = {}
for s in DATA_PATHS.keys():
    ALL_DATA[s] = np.sort(glob.glob(DATA_PATHS[s] + '/DD*.h5'))

for s in PERT_DATA_PATHS.keys():
    ALL_DATA[s] = np.sort(glob.glob(PERT_DATA_PATHS[s] + '/DD*.h5'))

for s in STAR_IC.keys():
    ALL_DATA[s] = np.sort(glob.glob(STAR_IC[s] + '/DD*.h5'))

def time_first_star(data = None, t = None, sfr = None):

    if data is None and (t is None or sfr is None):
        print "If supplying no data set, must supply both time array and sfr array"
        raise ValueError

    if not data is None:

        if ('t_first_star' in data['particle_meta_data']):
            return data['particle_meta_data']['t_first_star']

        t   = data['time_data']['time']
        sfr = data['time_data']['SFR']

    t_first = np.min( t[sfr > 0])

    return t_first

def plot_stellar_abundance(sim_names = None, species = 'metallicity'):
    """
    Plots average stellar abundance as a function of time for all simulations.
    This is done by passing either
    """
    if sim_names is None:
        sim_names = DATA_PATHS.keys()
    fig,ax = plt.subplots()

    for s in sim_names:
        # always load most recent file in every case
        data = dd.io.load(ALL_DATA[s][-1])

        t_first = time_first_star(data)

        t = np.zeros(np.size(ALL_DATA[s]))
        mass = np.zeros(np.size(ALL_DATA[s]))
        for i,x in enumerate(ALL_DATA[s]):
            xdata = dd.io.load(x)
            t[i] = xdata['meta_data']['Time']
            z[i] = xdata['particle_meta_data']['metallicity_stars'][2]

        t = t - t_first
        t_plot   = t[t >= 0.0]
        z_plot = z[t >= 0.0]

        ax.plot(t_plot, z_plot, lw = ps.lw, ls = ls_dict[s], color = color_dict[s], label = s, drawstyle = 'steps-post')


    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel('Average Stellar Abundance: ' + species)
    ax.semilogy()

    ax.set_xlim(global_tmin, global_tmax)
    ax.legend(loc='best')

    fig.set_size_inches(8,8)
    plt.tight_layout()
    plt.minorticks_on()

    fig.savefig('average_stellar_ ' + species + '_comparison.png')
    plt.close()

    return



    return

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
        species = ('gas','H_total_mass')
    elif species == 'metal' or species == 'Metal' or species == 'Metal_Mass' or species == 'metals' or species == 'Metals':
        species = ('gas','metal_mass')
    elif species == 'He':
        species = ('gas','He_total_mass')
    elif not '_Mass' in species:
        species = ('gas', species + '_Mass')


    fig, ax = plt.subplots()

    # loop over simulations
    for s in sim_names:
        # always load most recent file to check
        data = dd.io.load(ALL_DATA[s][-1])
        t_first = time_first_star(data)

        # make arrays for plotting
        t  = np.zeros(np.size(ALL_DATA[s]))
        ml = np.zeros(np.size(ALL_DATA[s]))

        # now go through every data analysis dump
        for i,x in enumerate(ALL_DATA[s]):
            try:
                xdata = dd.io.load(ALL_DATA[s][i], '/gas_profiles/outflow/disk')
            except:
                print "outflow rates load failed for " + ALL_DATA[s][i]
                continue

            t[i]  = dd.io.load(ALL_DATA[s][i], '/meta_data/Time')

            # find the bin whose center is closest to the desired z
            zbin  = np.argmin( np.abs( (xdata['xbins'][1:].value + xdata['xbins'][:-1].value)*0.5 - z ) )

            if species == 'total':
                ml[i] = xdata[('gas','H_total_mass')][zbin] +\
                        xdata[('gas','He_total_mass')][zbin] +\
                        xdata[('gas','metal_mass')][zbin]
            else:
                ml[i] = xdata[species][zbin]

        # normalize t and plot only after SF occurs
        t         = t - t_first
        t_plot    = t[t >= 0.0]
        mass_plot = ml[t >= 0.0]

        norm = 1.0
        if mass_loading: # maybe compute as the ~ 20 Myr average on either side of the data point?
            norm = 1.0 / SFR  

        ax.plot(t_plot, mass_plot*norm, lw = ps.lw, ls = ls_dict[s], color = color_dict[s], label = s, drawstyle = 'steps-post')

    sname = species
    if len(species) > 1:
        sname = species[1]

    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(sname + r' Mass Outflow Rate (M$_{\odot}$ yr$^{-1}$)')
    ax.semilogy()

    ax.set_xlim(global_tmin, global_tmax)
    ax.legend(loc='best')

    fig.set_size_inches(8,8)
    plt.tight_layout()
    plt.minorticks_on()

    plottype = 'outflow'
    if mass_loading:
        plottype = 'loading'
    else:
        ax.set_ylim(1.0E-11, 1.0E-1)

    outname = sname + '_mass_' + plottype + '_z%2f.png'%(z)

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

        ax.plot(t_plot, mass_plot, lw = ps.lw, ls = ls_dict[s], color = color_dict[s], label = s, drawstyle = 'steps-post')

    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(species + r' Mass (M$_{\odot}$)')
    ax.semilogy()

    ax.set_xlim(global_tmin, global_tmax)
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

        t    = data['time_data']['time'] / 1.0E6
        sfr  = data['time_data']['SFR']

        t_first = time_first_star(data)

        t = t - t_first

        t_plot   = t[t >= 0.0]
        sfr_plot = sfr[t >= 0.0]

        ax.plot(t_plot, sfr_plot, lw = ps.lw, ls = ls_dict[s], color = color_dict[s], label = s, drawstyle = 'steps-post')

    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(r'SFR (M$_{\odot}$ yr$^{-1}$)')
    ax.set_xlim(global_tmin, global_tmax)
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

        t    = data['time_data']['time'] / 1.0E6
        sfr  = data['time_data']['SFR']
        snr  = data['time_data']['SNII_snr']

        t_first = np.min( t[sfr > 0])

        t = t - t_first

        t_plot   = t[t >= 0.0]
        snr_plot = snr[t >= 0.0]

        ax.plot(t_plot, snr_plot, lw = ps.lw, ls = ls_dict[s], color = color_dict[s], label = s, drawstyle = 'steps-post')

    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(r'SNR (yr$^{-1}$)')
    ax.set_xlim(global_tmin, global_tmax)
    ax.semilogy()
    plt.minorticks_on()

    ax.legend(loc='best')

    fig.set_size_inches(8,8)
    plt.tight_layout()

    fig.savefig('snr_comparison.png')
    plt.close()

    return


if __name__ == '__main__':

#    all_s = ['lm','lm_noRT','lm_nowind','mm','mm_3pc','hm']

    all_s = STAR_IC.keys()

#    all_s = PERT_DATA_PATHS.keys()

    plot_mass(sim_names = all_s, species = 'HI')
    plot_mass(sim_names = all_s, species = 'stars')
    plot_mass(sim_names = all_s, species = 'total')


    plot_sfr(sim_names = all_s)
    plot_snr(sim_names = all_s)

    for species in ['total', 'metals', 'C', 'Fe', 'H']:

        plot_mass_loading(sim_names = all_s, species = species, z = 100)
#        plot_mass_loading(sim_names = all_s, species = species, z = 500)

