import deepdish as dd
from matplotlib import rc

fsize = 17
rc('text', usetex=False)
rc('font', size=fsize)#, ftype=42)
line_width = 3
point_size = 30

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
import glob
import numpy as np

# plot style defines
from galaxy_analysis.plot import plot_styles as ps


global_tmin = 0.0
global_tmax = 1000.0

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


run11_IC = { 'r11_smooth' : stampede + '/run11/200cc',
             'r11_pert'   : local + '/run11/200cc/perturb',
             'r11_lowsf'  : pleiades + '/starIC/run11/lowsf',
             'r11_lbox'   : pleiades + '/starIC/run11_largebox/no_wind',
             'r11 fiducial' : pleiades + '/starIC/run11_largebox/sndriving',
             'r11 2xr_s'    : pleiades + '/starIC/run11_2rs/sndriving',
             'r11 30km'     : pleiades + '/starIC/run11_30km/sndriving' ,
             'r11 40km'     : pleiades + '/starIC/run11_40km/sndriving'}

color_dict['r11 fiducial'] = ps.orange
color_dict['r11 2xr_s'] = ps.magenta
color_dict['r11 30km']  = ps.purple
color_dict['r11 40km']  = ps.black
ls_dict['r11 fiducial'] = '-'
ls_dict['r11 2xr_s'] = '-'
ls_dict['r11 30km'] = '-'
ls_dict['r11 40km'] = '-'

color_dict['r11_smooth'] = ps.blue
color_dict['r11_pert'  ] = ps.blue
ls_dict['r11_smooth'] = '-'
ls_dict['r11_pert'  ] = '--'

color_dict['r11_lbox'] = ps.black
ls_dict['r11_lbox']    = '-'

color_dict['r11_lowsf'] = ps.purple
color_dict['r11_lbox']  = ps.purple
ls_dict['r11_lowsf']    = '-'
ls_dict['r11_lowsf']    = ':'


#
# run11_feedback models
#

run11_feedback = { 'r11 - No Wind' : pleiades + '/starIC/run11_largebox/no_wind',
                   'r11 - No Wind 2': pleiades + '/starIC/run11_largebox/stampede/no_wind',
                   'r11 - No Ionizing'  : pleiades + '/starIC/run11_largebox/no_ion',
                   'r11 - No FUV/LW'    : pleiades + '/starIC/run11_largebox/no_otrad',
                   'r11 - SN only'      : pleiades + '/starIC/run11_largebox/sn_only',
                   'r11 - Full Physics' : pleiades + '/starIC/run11_largebox/fullphys'}

color_dict['r11 - No Wind'] = ps.black
color_dict['r11 - No Wind 2'] = ps.black
color_dict['r11 - No Ionizing'] = ps.magenta
color_dict['r11 - No FUV/LW']   = ps.blue
color_dict['r11 - SN only']     = ps.orange
color_dict['r11 - Full Physics'] = 'green'

for k in run11_feedback:
    ls_dict[k] = '-'
ls_dict['r11 - No Wind 2'] = '--'

run11_stampede_feedback = { 'r11 - st - No Wind' : pleiades + '/starIC/run11_largebox/stampede/no_wind',
                   'r11 - st - No Ionizing'  : pleiades + '/starIC/run11_largebox/stampede/no_ion',
                   'r11 - st - No FUV/LW'    : pleiades + '/starIC/run11_largebox/stampede/no_otrad',
                   'r11 - st - SN only'      : pleiades + '/starIC/run11_largebox/stampede/sn_only'}

color_dict['r11 - st - No Wind'] = ps.black
color_dict['r11 - st - No Ionizing'] = ps.magenta
color_dict['r11 - st - No FUV/LW']   = ps.blue
color_dict['r11 - st - SN only']     = ps.orange
for k in run11_stampede_feedback:
    ls_dict[k] = '-'


run15_feedback = {'r15 - No Wind' : pleiades + '/starIC/run15/lowsf',
                  'r15 - No Ionizing' : pleiades + '/starIC/run15/no_RT',
                  'r15 - No FUV/LW'   : pleiades + '/starIC/run15/no_otrad',
                  'r15 - SN only'     : pleiades + '/starIC/run15/sn_only',
                  'r15 - No Wind Compact' : pleiades + '/starIC/run15_compact/no_wind'}

color_dict['r15 - No Wind'] = ps.black
color_dict['r15 - No Ionizing'] = ps.magenta
color_dict['r15 - No FUV/LW'] = ps.blue
color_dict['r15 - SN only'] = ps.orange
color_dict['r15 - No Wind Compact'] = 'green'

for k in run15_feedback:
    ls_dict[k] = '-'



run15_IC = {'r15_smooth'   : stampede + '/run15',
            'r15_pert'     : pleiades + '/run15/no_wind',
            'r15_lowsf'    : pleiades + '/starIC/run15/lowsf',
            'r15_lowsf_c'  : pleiades + '/starIC/run15_compact/no_wind',
            'r15_medsf_c'  : pleiades + '/starIC/run15_compact/msf',
            'r15_highsf_c' : pleiades + '/starIC/run15_compact/highsf',
            'r15_lowsf_large' : pleiades + '/starIC/run15_largebox/no_wind/compact'}


color_dict['r15_smooth'] = ps.blue
color_dict['r15_pert']   = ps.blue
ls_dict['r15_smooth'] = '-'
ls_dict['r15_pert'] = '--'

color_dict['r15_lowsf'] = ps.black
color_dict['r15_lowsf_c'] = ps.orange
color_dict['r15_medsf_c'] = ps.orange
color_dict['r15_highsf_c'] = ps.orange
ls_dict['r15_lowsf'] = '-'
ls_dict['r15_lowsf_c'] = '-'
ls_dict['r15_medsf_c'] = '--'
ls_dict['r15_highsf_c'] = ':'

color_dict['r15_lowsf_large'] = ps.magenta
ls_dict['r15_lowsf_large'] = '-'


comparison_sim = {'Hu et. al. 2017' : pleiades + '/hu',
                  'Forbes et. al. 1 kpc' : pleiades + '/forbes/starIC',
                  'run11 2 x r_s'  : pleiades + '/starIC/run11_2rs/no_wind',
                  'run11 sndriving' : pleiades + '/starIC/run11_largebox/sndriving',
                  'Hu sndriving' : pleiades +'/hu/sndriving',
                  'Hu shortrad'  : pleiades +'/hu/sndriving_shortrad',
#                  'run11 no radpressure' : pleiades + '/starIC/run11_largebox/no_radpressure',
                  'run11 original' : pleiades + '/starIC/run11_largebox/no_wind',
                  'run15 lbox' : pleiades + '/starIC/run15_largebox/better_IC'}

ls_dict['Hu et. al. 2017'] = '-'
color_dict['Hu et. al. 2017'] = ps.black
ls_dict['Hu sndriving'] = '-'
color_dict['Hu sndriving'] = ps.black
ls_dict['Hu shortrad'] = '-'
color_dict['Hu shortrad'] = ps.blue

ls_dict['Forbes et. al. 1 kpc'] = '-'
color_dict['Forbes et. al. 1 kpc'] = ps.blue
color_dict['run15 lbox'] = 'red'

ls_dict['run11 2 x r_s'] = '-'
ls_dict['run11 no radpressure'] = ':'
ls_dict['run11 original'] = '-'
ls_dict['run11 sndriving'] = ':'
ls_dict['run15 lbox'] = '-'

color_dict['run11 sndriving'] = ps.orange
color_dict['run11 2 x r_s'] = 'green'
color_dict['run11 no radpressure'] = ps.orange
color_dict['run11 original'] = ps.orange


ALL_DATA = {}
for s in DATA_PATHS.keys():
    ALL_DATA[s] = np.sort(glob.glob(DATA_PATHS[s] + '/DD*.h5'))

for s in PERT_DATA_PATHS.keys():
    ALL_DATA[s] = np.sort(glob.glob(PERT_DATA_PATHS[s] + '/DD*.h5'))

for s in STAR_IC.keys():
    ALL_DATA[s] = np.sort(glob.glob(STAR_IC[s] + '/DD*.h5'))

for s in run11_IC.keys():
    ALL_DATA[s] = np.sort(glob.glob(run11_IC[s] + '/DD*.h5'))

for s in run15_IC.keys():
    ALL_DATA[s] = np.sort(glob.glob(run15_IC[s] + '/DD*.h5'))

for s in run11_feedback.keys():
    ALL_DATA[s] = np.sort(glob.glob(run11_feedback[s] + '/DD*.h5'))

for s in run11_stampede_feedback.keys():
    ALL_DATA[s] = np.sort(glob.glob(run11_stampede_feedback[s] + '/DD*.h5'))

for s in run15_feedback.keys():
    ALL_DATA[s] = np.sort(glob.glob(run15_feedback[s] + '/DD*.h5'))

for s in comparison_sim.keys():
    ALL_DATA[s] = np.sort(glob.glob(comparison_sim[s] + '/DD*.h5'))

def time_first_star(data = None, t = None, sfr = None):

    if data is None and (t is None or sfr is None):
        print "If supplying no data set, must supply both time array and sfr array"
        raise ValueError

    if not data is None:

        if ('t_first_star' in data['particle_meta_data']):
            return data['particle_meta_data']['t_first_star']

        if ('SFR_1' in data['time_data']):
            t   = data['time_data']['time_1']
            sfr = data['time_data']['SFR_1']
        else:
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
        print s
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




def plot_mass_loading(sim_names = None, species = 'total', z = 0.25, mass_loading = False):
    """
    Given a dictionary that contains list of simulation names and filepaths,
    go through all simulation data analysis outputs in that path for each
    simulation and plot the time evolution of the outflow rate. 

    Inputs
    ------
        sim_names    : dictionary, optional
        species      : String. Elemental species (or 'total' or 'metal' for all gas or
                       all metals) of which to plot the outflow rate. optional. default: total
        z            : Fixed fraction of virial radius to plot the outflow from. Default 0.25
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
                xdata = dd.io.load(ALL_DATA[s][i], '/gas_profiles/outflow/sphere')
            except:
                print "outflow rates load failed for " + ALL_DATA[s][i]
                continue

            t[i]  = dd.io.load(ALL_DATA[s][i], '/meta_data/Time')

            # find the bin whose center is closest to the desired z
            zbin  = np.argmin( np.abs( (xdata['centers_rvir'][1:] + xdata['centers_rvir'][:-1])*0.5 - z ) )

            if species == 'total':
                try:
                    ml[i] = xdata[('gas','cell_mass')][zbin]
                except:
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

    ymin = 1.0E99
    ymax = -1.0E99

    for s in sim_names:
        # always load most recent file in every case
        print s
        data = dd.io.load(ALL_DATA[s][-1])

        t_first = time_first_star(data)

        if species == 'HI' or species == 'M_HI':
            fname = 'M_HI'
        elif species == 'total' or species == 'Total' or species == 'gas' or species == 'M_gas':
            fname = 'M_gas'
        elif species == 'star' or species == 'stars' or species == 'stellar' or species == 'M_star':
            fname = 'M_star'
        else:
            if species == 'Metals' or species == 'Metal' or species == 'metals':
                species = 'metal'

            fname = 'M_' + species

        t = np.zeros(np.size(ALL_DATA[s]))
        mass = np.zeros(np.size(ALL_DATA[s]))
        for i,x in enumerate(ALL_DATA[s]):
            xdata = dd.io.load(x, '/meta_data')
            t[i] = xdata['Time']

            if fname in xdata: # new style
                mass[i] = xdata[fname]
            else: # for backwards compatability
                tempdata = dd.io.load(x, '/gas_profiles/accumulation/disk')
                if not (species == 'metal'):
                    fname = ('gas', species + '_Mass')
                else:
                    fname = ('gas', species + '_mass')

                if not (fname in tempdata):
                    print "Breaking for field ", fname, " in plot mass evolution"
                    break

                mass[i] = np.sum(tempdata[fname])

        t = t - t_first
        t_plot   = t[t >= 0.0]
        mass_plot = mass[t >= 0.0]

        ax.plot(t_plot, mass_plot, lw = ps.lw, ls = ls_dict[s], color = color_dict[s], label = s, drawstyle = 'steps-post')

        ymin = np.min(  [ymin, np.min(mass_plot)])
        ymax = np.max(  [ymax, np.max(mass_plot)])

    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(species + r' Mass (M$_{\odot}$)')
    ax.semilogy()

    ymin = np.min( [ymin, 1.0E-2*ymax])

    ymin = np.max( [ymin, 1.0E-10*ymin])

    ax.set_ylim(ymin, ymax)

    ax.set_xlim(global_tmin, global_tmax)
    ax.legend(loc='best')

    fig.set_size_inches(8,8)
    plt.tight_layout()
    plt.minorticks_on()

    fig.savefig(species + '_mass_comparison.png')
    plt.close()

    return

def plot_sfr(sim_names = None, sampling = None):

    if sim_names is None:
        sim_names = DATA_PATHS.keys()

    fig,ax = plt.subplots()

    for s in sim_names:
        # always load most recent file in every case
        data = dd.io.load(ALL_DATA[s][-1])

        if sampling is None or (sampling == 10):  # load default (10 Myr)
            t    = data['time_data']['time'] / 1.0E6
            sfr  = data['time_data']['SFR']
        elif sampling == 1:
            t = data['time_data']['time_1'] / 1.0E6
            sfr = data['time_data']['SFR_1']
        elif sampling == 100:
            t = data['time_data']['time_100'] / 1.0E6
            sfr = data['time_data']['SFR_100']

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

    outstr = 'sfr_comparison'
    if not (sampling is None):
        if sampling == 1:
            outstr += '_1'
        elif sampling == 100:
            outstr += '_100'

    fig.savefig(outstr  + '.png')
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

#    all_s = STAR_IC.keys()
#    all_s = run15_IC.keys()
#    all_s = run11_stampede_feedback.keys()
#    all_s = comparison_sim.keys()
    all_s = ['Hu sndriving','Hu shortrad']
    all_s = ['r11 fiducial', 'r11 2xr_s', 'r11 30km', 'r11 40km']

#    all_s = PERT_DATA_PATHS.keys()

    plot_mass(sim_names = all_s, species = 'HI')
    plot_mass(sim_names = all_s, species = 'stars')
    plot_mass(sim_names = all_s, species = 'total')
    plot_mass(sim_names = all_s, species = 'metals')
    plot_mass(sim_names = all_s, species = 'Fe')
    plot_mass(sim_names = all_s, species = 'O')
    plot_mass(sim_names = all_s, species = 'C')
    plot_mass(sim_names = all_s, species = 'Mg')


    plot_sfr(sim_names = all_s)
    plot_sfr(sim_names = all_s, sampling = 1)
    plot_sfr(sim_names = all_s, sampling = 100)
    plot_snr(sim_names = all_s)

    for species in ['total', 'metals', 'O', 'C', 'Mg', 'Fe']:
        plot_mass_loading(sim_names = all_s, species = species, z = 0.25)
        plot_mass_loading(sim_names = all_s, species = species, z = 1.0)
#        plot_mass_loading(sim_names = all_s, species = species, z = 500)

