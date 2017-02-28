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

def plot_mass(sim_names = None, species = 'HI'):


    if sim_names is None:
        sim_names = DATA_PATHS.keys()

    fig,ax = plt.subplots()

    for s in sim_names:
        # always load most recent file in every case
        data = dd.io.load(ALL_DATA[s][-1])

        t    = data['time_data']['time']
        sfr  = data['time_data']['SFR']

        t_first = np.min( t[sfr > 0]) / 1.0E6

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

    plot_snr(sim_names = all_s) 
    plot_sfr(sim_names = all_s)
    plot_mass(sim_names = all_s, species = 'HI')
    plot_mass(sim_names = all_s, species = 'stars')
    plot_mass(sim_names = all_s, species = 'total')
    
