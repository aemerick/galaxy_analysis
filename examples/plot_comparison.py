#
# This is a great example of bad code
#
import deepdish as dd
from matplotlib import rc,cm
plasma = cm.get_cmap('plasma')

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
global_tmax = 500.0
#global_tmax = 270.0
#
#
#

workdir  = '/mnt/ceph/users/emerick/enzo_runs'
stampede  = workdir + '/stampede/leo_p/fullphysics/fullres'
stampede2  = workdir + '/stampede/leo_p/starIC'
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

for k in list(PERT_DATA_PATHS.keys()) + list(STAR_IC.keys()):

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
             '17' : pleiades + '/starIC/run11/corrected_sndriving',
             '25'    : pleiades + '/starIC/run11_2rs/corrected_sndriving',
#             'r11 ibug'     : pleiades + '/starIC/run11_2rs/ion_bugfix',
             '30'     : pleiades + '/starIC/run11_30km/final_sndriving' ,
             '30_3pc' : pleiades + '/starIC/run11_30km/3pc',
             '30_3pch' : pleiades + '/starIC/run11_30km/3pc_hsn',
             '30_3pcmh' : pleiades + '/starIC/run11_30km/3pc_mhsn',
             '30_6pch' : pleiades + '/starIC/run11_30km/6pc_hsn',
             '30_3pcvh' : pleiades + '/starIC/run11_30km/3pc_vhsn',
             '40'     : pleiades + '/starIC/run11_40km/corrected_sndriving'}

color_dict['17'] = ps.orange
color_dict['25'] = ps.magenta
color_dict['r11 ibug']  = ps.magenta
color_dict['30']  = ps.purple
color_dict['30_3pc'] = ps.blue
color_dict['30_3pch'] = ps.orange
color_dict['30_3pcmh'] = ps.magenta
color_dict['30_6pch'] = 'C1'
color_dict['40']  = ps.black
ls_dict['17'] = '-'
ls_dict['25'] = '-'
ls_dict['r11 ibug']  = '--'
ls_dict['30'] = '-'
ls_dict['30_3pc'] = '-'
ls_dict['30_3pch'] = '-'
ls_dict['30_6pch'] = '-'
ls_dict['30_3pcmh'] = '-'
ls_dict['40'] = '-'

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

run13_feedback = {'Full Physics' : pleiades + '/starIC/run11_30km/final_sndriving',
                  'Ion Only' : stampede2 + '/run11_30km/ion_no-otrad-sn',
                  'PE+LW Only' : stampede2 + '/run11_30km/otrad_no-ion-sn',
                  'SN Only' : stampede2 + '/run11_30km/sn_only',
                  'Ion+Pe+LW - No SN' : stampede2 + '/run11_30km/otrad_ion-no-sn',
                  'SN+Ion' : stampede2 + '/run11_30km/sn_ion-no-otrad',
                  'Full - No RP' : stampede2 + '/run11_30km/sn_otrad_ion_noRP'}

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
                  'Hu sndriving' : pleiades +'/hu_metal/correct_feedback',
                  'Hu shortrad'  : pleiades +'/hu_metal/correct_shortrad',
#                  'Hu sndriv metal' : pleiades + '/hu_metal/correct_feedback',
#                  'Hu shortrad metal' : pleiades + '/hu_metal/correct_shortrad/',
#                  'run11 no radpressure' : pleiades + '/starIC/run11_largebox/no_radpressure',
                  'run11 original' : pleiades + '/starIC/run11_largebox/no_wind',
                  'run15 lbox' : pleiades + '/starIC/run15_largebox/better_IC'}

IC_comparisons = {'collapse' : stampede2 + '/run11_30km/collapse_IC',
                  'cool ramp' : stampede2 + '/run11_30km/cool_ramp',
                  'cool ramp2' : stampede2 + '/run11_30km/cool_ramp2',
                  'cool ramp3' : stampede2 + '/run11_30km/cool_ramp3',
                  'CR2 SN'     : stampede2 + '/run11_30km/cool_ramp2_sn',
                  'SNx5'       : stampede2 + '/run11_30km/normalx5',
                  'Fiducial'   : pleiades + '/starIC/run11_30km/final_sndriving'}
for k in IC_comparisons:
    ls_dict[k] = '-'
ls_dict['CR2 SN'] = '--'
ls_dict['SNx5']   = '--'
color_dict['collapse']   = ps.purple
color_dict['cool ramp']  = ps.blue
color_dict['cool ramp2'] = ps.orange
color_dict['cool ramp3'] = ps.green
color_dict['CR2 SN'] = ps.blue
color_dict['Fiducial'] = ps.black
color_dict['SNx5'] = ps.black

ls_dict['Hu et. al. 2017'] = '-'
color_dict['Hu et. al. 2017'] = ps.black
ls_dict['Hu sndriving'] = '-'
color_dict['Hu sndriving'] = ps.black
ls_dict['Hu shortrad'] = '-'
color_dict['Hu shortrad'] = ps.blue
color_dict['Hu sndriv metal'] = ps.purple
ls_dict['Hu sndriv metal'] = '-'

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


feedback_comparisons = {'ion_no-otrad-sn' : stampede2 + '/run11_30km/ion_no-otrad-sn',
                        'otrad_ion-no-sn' : stampede2 + '/run11_30km/otrad_ion-no-sn',
                        'otrad_no-ion-sn' : stampede2 + '/run11_30km/otrad_no-ion-sn',

                        'sn_ion_LW-no-PE' : stampede2 + '/run11_30km/sn_ion_LW-no-PE',
                        'sn_ion-no-otrad' : stampede2 + '/run11_30km/sn_ion-no-otrad',
                        'sn_ion_PE-no-LW' : stampede2 + '/run11_30km/sn_ion_PE-no-LW',

                        'sn_otrad_ion_noRP' : stampede2 + '/run11_30km/sn_otrad_ion_noRP',
                        'sn_otrad_ion_RPx2' : stampede2 + '/run11_30km/sn_otrad_ion_RPx2',
                        'sn_otrad_no-ion' : stampede2 + '/run11_30km/sn_otrad_no-ion',

                        'Full' : pleiades + '/starIC/run11_30km/final_sndriving',
                        'Shortrad' : pleiades + '/starIC/run11_30km/final_shortrad',

                        'SN only' : stampede2 + '/run11_30km/sn_only'}
feedback_plot_1 = ['Full']
x = 'Radiation Only - No SN'
feedback_comparisons[x] = feedback_comparisons['otrad_ion-no-sn']
ls_dict[x] = '-'; color_dict[x] = ps.purple
feedback_plot_1 += [x]
x = 'SN + PE + LW - No Ion'
feedback_comparisons[x] = feedback_comparisons['sn_otrad_no-ion']
ls_dict[x] = '-'; color_dict[x] = ps.orange
feedback_plot_1 += [x]
x= 'SN + Ion - No PE+LW'
feedback_comparisons[x] = feedback_comparisons['sn_ion-no-otrad']
ls_dict[x] = '-'; color_dict[x] = ps.blue
feedback_plot_1 += [x]
x = 'Full - No RP'
feedback_comparisons[x] = feedback_comparisons['sn_otrad_ion_noRP']
ls_dict[x] = '-'; color_dict[x] = 'green'
feedback_plot_1 += [x]
x = "SN + Ion + Pe - no LW"
feedback_comparisons[x] = feedback_comparisons['sn_ion_PE-no-LW']
ls_dict[x] = '--'; color_dict[x] = 'blue'
feedback_plot_1 += [x]
x = "SN + Ion + LW - no PE"
feedback_comparisons[x] = feedback_comparisons['sn_ion_LW-no-PE']
ls_dict[x] = '--'; color_dict[x] = 'red'
feedback_plot_1 += [x]

feedback_plot_2 = ['Full']
x = 'PE+LW Only'
feedback_comparisons[x] = feedback_comparisons['otrad_no-ion-sn']
ls_dict[x] = '-'; color_dict[x] = ps.purple
feedback_plot_2 +=[x]
x = 'Ionizing Only'
feedback_comparisons[x] = feedback_comparisons['ion_no-otrad-sn']
ls_dict[x] = '-'; color_dict[x] = ps.orange
feedback_plot_2 += [x]
x = 'SN Only'
feedback_comparisons[x] = feedback_comparisons['SN only']
ls_dict[x] = '-'; color_dict[x] = ps.blue
feedback_plot_2 += [x]


#for i,k in enumerate(feedback_comparisons.keys()):
#    color_dict[k] = plasma(i / (1.0 * len(feedback_comparisons.keys())))
#    ls_dict[k]    = '-'
ls_dict['Full'] = '-' ; color_dict['Full'] = 'black'
ls_dict['Shortrad'] = '--' ; color_dict['Shortrad'] = 'black'
ls_dict['SN only'] = '-'; color_dict['SN only'] = plasma(1.0)

color_dict['sn_ion-no-otrad'] = plasma(0.0); ls_dict['sn_ion-no-otrad'] = '-'
color_dict['sn_ion_LW-no-PE'] = plasma(0.0); ls_dict['sn_ion_LW-no-PE'] = '--'
color_dict['sn_ion_PE-no-LW'] = plasma(0.0); ls_dict['sn_ion_PE-no-LW'] = '-.'

color_dict['sn_otrad_no-ion'  ] = plasma(0.25); ls_dict['sn_otrad_no-ion'] = '-'
color_dict['sn_otrad_ion_RPx2'] = plasma(0.25); ls_dict['sn_otrad_ion_RPx2'] = '--'
color_dict['sn_otrad_ion_RPx10'] = plasma(0.25); ls_dict['sn_otrad_ion_RPx10'] = '-.'

color_dict['otrad_no-ion-sn'] = plasma(0.5); ls_dict['otrad_ion-no-sn'] = '-'
color_dict['otrad_no-ion-sn'] = plasma(0.5); ls_dict['otrad_no-ion-sn'] = '--'
color_dict['ion_no-otrad-sn'] = plasma(0.5); ls_dict['ion_no-otrad-sn'] = '-.'


ALL_DATA = {}
for s in DATA_PATHS.keys():
    ALL_DATA[s] = np.sort(glob.glob(DATA_PATHS[s] + '/DD*.h5'))
for s in feedback_comparisons.keys():
    ALL_DATA[s] = np.sort(glob.glob(feedback_comparisons[s] + '/DD*.h5'))

for s in IC_comparisons.keys():
    ALL_DATA[s] = np.sort(glob.glob(IC_comparisons[s] + '/DD*.h5'))

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
        print("If supplying no data set, must supply both time array and sfr array")
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
        sim_names = list(DATA_PATHS.keys())
    fig,ax = plt.subplots()

    for s in sim_names:
        # always load most recent file in every case
        print(s)
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




def plot_mass_loading(sim_names = None, species = 'total', z = 0.25, mass_loading = False, include_inflow = True):
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
        sim_names = list(DATA_PATHS.keys())

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
    ymin =  1.0E99
    ymax = -1.0E99
    for s in sim_names:
        # always load most recent file to check
        data = dd.io.load(ALL_DATA[s][-1])
        t_first = time_first_star(data)
        SFR     = data['time_data']['SFR_100'][-1] # Msun / yr

        # make arrays for plotting
        t  = np.zeros(np.size(ALL_DATA[s]))
        ml = np.zeros(np.size(ALL_DATA[s]))
        iml = np.zeros(np.size(ALL_DATA[s]))

        # now go through every data analysis dump
        for i,x in enumerate(ALL_DATA[s]):
            try:
                xdata = dd.io.load(ALL_DATA[s][i], '/gas_profiles/outflow/sphere')
                ixdata = dd.io.load(ALL_DATA[s][i], '/gas_profiles/inflow/sphere')
            except:
                print(("outflow rates load failed for " + ALL_DATA[s][i]))
                continue

            t[i]  = dd.io.load(ALL_DATA[s][i], '/meta_data/Time')

            # find the bin whose center is closest to the desired z
            zbin  = np.argmin( np.abs( xdata['centers_rvir'] - z ) )

            if species == 'total':
                try:
                    ml[i] = xdata[('gas','cell_mass')][zbin]
                    iml[i] = ixata[('gas','cell_mass')][zbin]
                except:
                    ml[i] = xdata[('gas','H_total_mass')][zbin] +\
                            xdata[('gas','He_total_mass')][zbin] +\
                            xdata[('gas','metal_mass')][zbin]
                    iml[i] = ixdata[('gas','H_total_mass')][zbin] +\
                             ixdata[('gas','He_total_mass')][zbin] +\
                             ixdata[('gas','metal_mass')][zbin]
            else:
                ml[i] = xdata[species][zbin]
                iml[i] = ixdata[species][zbin]

        # normalize t and plot only after SF occurs
        t         = t - t_first
        t_plot    = t[t >= 0.0]
        mass_plot = ml[t >= 0.0]
        imass_plot = iml[t >= 0.0]

        norm = 1.0
        if mass_loading:
            if SFR <= 0.0:
                norm = 0.0
            else:
                norm = 1.0 / SFR

        ax.plot(t_plot, mass_plot*norm, lw = ps.lw, ls = ls_dict[s], color = color_dict[s], label = s, drawstyle = 'steps-post')
        if include_inflow:
            ax.plot(t_plot, np.abs(imass_plot*norm), lw = ps.lw, ls = '--', color = color_dict[s], drawstyle = 'steps-post')
        ymax = np.max([np.max(mass_plot*norm),ymax])
        ymin = np.min([np.min(mass_plot*norm),ymin])

    sname = species
    if len(species) == 2:
        sname = species[1]

    ax.set_xlabel(r'Time (Myr)')
    if mass_loading:
        ax.set_ylabel(sname + r'Mass Loading Factor')
    else:
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

    ymin = np.min( [ymin, 1.0E-2*ymax]) # show at least 2 orders of mag
    ymin = np.max( [ymin, 1.0E-6*ymax]) # but no more than 6
    ax.set_ylim(ymin,ymax)

    outname = sname + '_mass_' + plottype + '_z%2f.png'%(z)

    fig.savefig(outname)
    plt.close()

    return            

def plot_mass(sim_names = None, species = 'HI'):


    if sim_names is None:
        sim_names = list(DATA_PATHS.keys())

    fig,ax = plt.subplots()

    ymin = 1.0E99
    ymax = -1.0E99

    for s in sim_names:
        # always load most recent file in every case
        print(s)
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
                    print(("Breaking for field ", fname, " in plot mass evolution"))
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

    ymin = np.max( [ymin, 1.0E-3*ymax])

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
        sim_names = list(DATA_PATHS.keys())

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
        sim_names = list(DATA_PATHS.keys())

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
#    all_s = ['Hu sndriving','Hu shortrad'] #, 'Hu sndriv metal' 'Hu shortrad metal']
#    all_s = ['r11 fiducial', 'r11 2xr_s', 'r11 30km', 'r11 40km']

#    all_s = ['Fiducial','collapse','cool_ramp','cool_ramp2','cool_ramp3','CR2 SN']
#    all_s = feedback_comparisons.keys()

    all_s = ['30','30_3pc','30_3pch','30_3pcmh','30_6pch']
#    all_s = ['30']
#    all_s = PERT_DATA_PATHS.keys()

#    all_s = feedback_plot_1
#    all_s = feedback_plot_2

    plot_sfr(sim_names = all_s)
    plot_sfr(sim_names = all_s, sampling = 1)
    plot_sfr(sim_names = all_s, sampling = 100)
    plot_snr(sim_names = all_s)
    print("done with sfr")

    if True:
        plot_mass(sim_names = all_s, species = 'HI')
        plot_mass(sim_names = all_s, species = 'HII')
        plot_mass(sim_names = all_s, species = 'H2')
#        plot_mass(sim_names = all_s, species = 'H2_total')
#        plot_mass(sim_names = all_s, species = 'He_total')
        plot_mass(sim_names = all_s, species = 'H_total')
        plot_mass(sim_names = all_s, species = 'stars')
        plot_mass(sim_names = all_s, species = 'total')
        plot_mass(sim_names = all_s, species = 'metals')
#        plot_mass(sim_names = all_s, species = 'Fe')
#        plot_mass(sim_names = all_s, species = 'O')
#        plot_mass(sim_names = all_s, species = 'C')
#        plot_mass(sim_names = all_s, species = 'Mg')
#        plot_mass(sim_names = all_s, species = 'N')
#        plot_mass(sim_names = all_s, species = 'Si')


        plot_sfr(sim_names = all_s)
        plot_sfr(sim_names = all_s, sampling = 1)
        plot_sfr(sim_names = all_s, sampling = 100)
        plot_snr(sim_names = all_s)

    if False:
        for species in ['total', 'metals', 'O', 'N', 'Fe', 'Sr']:
            print(('mass outflow for ' + species))
            plot_mass_loading(sim_names = all_s, species = species, z = 0.1)
            plot_mass_loading(sim_names = all_s, species = species, z = 0.25)
            plot_mass_loading(sim_names = all_s, species = species, z = 0.5)
            plot_mass_loading(sim_names = all_s, species = species, z = 1.0)

        species = 'total'
        print(('mass loading for ' + species))
        plot_mass_loading(sim_names = all_s, species = species, z = 0.1, mass_loading = True)
        plot_mass_loading(sim_names = all_s, species = species, z = 0.25, mass_loading = True)
        plot_mass_loading(sim_names = all_s, species = species, z = 0.5, mass_loading = True)
        plot_mass_loading(sim_names = all_s, species = species, z = 1.0, mass_loading = True)
        plot_mass_loading(sim_names = all_s, species = species, z = 500)

