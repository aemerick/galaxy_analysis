import os
import sys
#
#
# Use this script to run through all of the individual python scripts
#   / functions to generate the necessary plots for all of the things
#     in the method paper. Additional work will be requried after
#     this point to get to the final versions, however.
#
work_dir     = './'
output_dir   = './method_paper_plots/'

override_dsi = [119, 209, 299, 400] # plot these files!
t_o          = 119                # dataset of first star formation


if not os.path.exists(output_dir):
    os.mkdir(output_dir)

def get_dsi(x):
    if override_dsi is None:
        return x
    else:
        return override_dsi

# Fig. 1 and Fig. 2: Phase-on and Edge-on projections:
if not ('no_panel' in sys.argv):
    from galaxy_analysis.method_paper_plots import panel_plot

# need to pick which datasets to use
#    can optionally choose additional fields here using
#    field_list argument
    _dsi     = [46,196,346,546]
    ds_list = ["DD%0004i/DD%0004i"%(i,i) for i in get_dsi(_dsi)]
    panel_plot.plot_panel(ds_list = ds_list, axis = 'x', outdir = output_dir)
    panel_plot.plot_panel(ds_list = ds_list, axis = 'z', outdir = output_dir)


# Fig. 3: Scale height plot:
#    requires running
#    compute_scale_height.py to calculate the
#    gas scale height at each time step
#
#    this can be done as python compute_scale_height.py
#
if not ('no_scale' in sys.argv):
    from galaxy_analysis.method_paper_plots import compute_scale_height

    compute_scale_height.plot_all_data(t_o = t_o,
                                       dt = 5,
                                       t = [(x - t_o) for x in get_dsi(_dsi)] )
    compute_scale_height.plot_phase_comparison(t_o = t_o,
                                               t = 30, dt = 5,
                                               phases=['CNM','WNM','WIM','HIM'])

# Fig. 4:
#
from galaxy_analysis.method_paper_plots import sfr_snr
from galaxy_analysis.method_paper_plots import mass_plot

sfr_snr.plot(work_dir   = work_dir, outdir = output_dir)
mass_plot.plot_mass_evolution(work_dir = work_dir, outdir = output_dir)

#
# Fig. 5: Time average phase diagrams:
#
if not ( 'no_phase' in sys.argv):
    t_min = 390.0 # 300.0
    t_max = 399.0 # 350.0

    from galaxy_analysis.method_paper_plots import phase_diagrams
    phase_diagrams.plot_time_average_PD(work_dir, t_min = t_min, t_max = t_max, plots = ['nT'],
                                        outdir = output_dir)

#
# Fig. 6: mass and volume fraction evolution plots
#
from galaxy_analysis.plot import plot_mass_volume_fractions

plot_mass_volume_fractions.plot_fractions(outdir = output_dir)

#
# Fig. 7: G_o and Q_o 1D radial profiles
#
from galaxy_analysis.method_paper_plots import radiation_profiles
radiation_profiles.plot(t_min = t_min, t_max = t_max,
                        fields = ['G_o','Q0_flux'],
                        work_dir = work_dir, outdir = output_dir)


#
# Fig. 8: G_o and Q_o 2D phase diagrams
#
t_min = 408.1 # 200
t_max = 410.1 # 201
phase_diagrams.plot_time_average_PD(work_dir, t_min = t_min, t_max = t_max, 
                                    plots = ['G_o','Q_o'], outdir = output_dir)

#
# Fig. 9: Mass outflow, mass loading rates
#   and
# Fig. 11 (top panel only): Metal mass loading factor
#
from galaxy_analysis.method_paper_plots import mass_outflow
mass_outflow.plot_basic_outflow_and_loading(work_dir = work_dir, t_min = 0.0, t_max = 1000.0,
                                            outdir = output_dir)


#
# Fig. 10: time averaged radial velocity
#
t_min = 390.0
t_max = 420.1
from galaxy_analysis.method_paper_plots import time_average_velocity
time_average_velocity.plot(workdir = work_dir, outdir = output_dir,
                           t_min = t_min, t_max = t_max)

#
# Fig. 11: bottom panel
#
from galaxy_analysis.method_paper_plots import metal_retention
metal_retention.plot_metal_retention(workdir = work_dir, outdir = output_dir)

#
# Fig. 12:
#

#
# This needs to be ported from the ipython notebook
#

#
# Fig. 13:
#
from galaxy_analysis.method_paper_plots import schmidt_law

schmidt_law.schmidt_law(work_dir = work_dir, obs_method = True,
                        total_gas = False, tmin = 0.0) # maybe change tmin to after trans phase



#
#
# Appendix plots from simulation data
#
#

#
# Fig. E1
#
from galaxy_analysis.method_paper_plots import sfr_resolution
# from galaxy_analysis.method_paper_plots import 
comparison = {'3pcH2' : ('../3pc_H2', '3.6 pc', '--'),
              '6pcH2'  :  '../6pc_H2', '7.2 pc', '-.'),
              'Fiducial' : (work_dir, 'Fiducial', '-')}
sfr_resolution.sfr_resolution(work_dir = work_dir,
                              output_dir = output_dir, comparison=comparison)

mass_plot.plot_mass_resolution(work_dir = work_dir, output_dir = output_dir, comparison = comparison)
