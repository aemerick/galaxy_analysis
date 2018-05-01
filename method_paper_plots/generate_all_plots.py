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

override_dsi = [119, 129, 139, 149] # plot these files!
t_o          = 119                # dataset of first star formation


if not os.path.exists(output_dir):
    os.mkdir(output_dir)

def get_dsi(x):
    if override_dsi is None:
        return x
    else:
        return override_dsi

# Fig. 1 and Fig. 2: Phase-on and Edge-on projections:
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
from galaxy_analysis.method_paper_plots import compute_scale_height

compute_scale_height.plot_all_data(t_o = t_o, dt = 20, t = get_dsi(_dsi))
compute_scale_height.plot_phase_comparison(t_o=t_o,t=100,dt=20,phases=['CNM','WNM','WIM','HIM'])

# Fig. 4:
#
from galaxy_analysis.method_paper_plots import sfr_snr
from galaxy_analysis.method_paper_plots import mass_plot

sfr_snr.plot(work_dir   = work_dir, outdir = output_dir)
mass_plot.plot_mass_evolution(work_dir = work_dir, outdir = output_dir)

#
# Fig. 5: Time average phase diagrams:
#
t_min = 145.0 # 300.0
t_max = 151.0 # 350.0

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
t_min = 149.1 # 200
t_max = 150.1 # 201
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
# Fig. 13:
#
