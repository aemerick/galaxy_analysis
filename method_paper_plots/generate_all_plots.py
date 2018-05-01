#
#
# Use this script to run through all of the individual python scripts
#   / functions to generate the necessary plots for all of the things
#     in the method paper. Additional work will be requried after
#     this point to get to the final versions, however.
#
work_dir     = './'
output_dir   = './method_paper_plots/'

override_dsi = [129,139,140, 146] # plot these files!
t_o          = 119                # dataset of first star formation

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
panel_plot.plot_panel(ds_list = ds_list, axis = 'x')
panel_plot.plot_pabel(ds_lits = ds_list, axis = 'z')


# Fig. 3: Scale height plot:
#    requires running
#    compute_scale_height.py to calculate the
#    gas scale height at each time step
#
#    this can be done as python compute_scale_height.py
#
from galaxy_analysis.method_paper_plots import compute_scale_height

compute_scale_height.plot_all_data(t_o = t_o, dt = 20, t = get_dsi(_dsi))
compute_scale_height.plot_phase_comparison(t_o=t_o,dt=20,phases=['CNM','WNM','WIM','HIM'])

# Fig. 4:
#
from galaxy_analysis.method_paper_plots import sfr_snr
from galaxy_analysis.method_paper_plots import mass_plot

sfr_snr.plot(work_dir   = work_dir, outdir = output_dir)
mass_plot.plot_mass_evolution(work_dir = work_dir, outdir = output_dir)

#
# Fig. 5: Time average phase diagrams:
#
from galaxy_analysis.method_paper_plots import phase_diagram
phase_diagrams.plot_time_average_PD(work_dir, t_min = 300.0, t_max = 350.0, plots = ['nT'])

#
# Fig. 6: mass and volume fraction evolution plots
#
from galaxy_analysis.plot import plot_mass_volume_fractions

plot_mass_volume_fractions.plot_fractions()

#
# Fig. 8: G_o and Q_o radial profiles
#
phase_diagrams.plot_time_average_PD(workdir, t_min = 200.0, t_max = 201.0, plots = ['G_o','Q_o'])
