#
#
# Use this script to run through all of the individual python scripts
#   / functions to generate the necessary plots for all of the things
#     in the method paper. Additional work will be requried after
#     this point to get to the final versions, however.
#


# Fig. 1 and Fig. 2: Phase-on and Edge-on projections:
from galaxy_analysis.method_paper_plots import panel_plot

# need to pick which datasets to use
#    can optionally choose additional fields here using
#    field_list argument
dsi     = [46,196,346,546]
ds_list = ["DD%0004i/DD%0004i"%(i,i) for i in dsi]

panel_plot.plot_panel(ds_list = ds_list, axis = 'x')
panel_plot.plot_pabel(ds_lits = ds_list, axis = 'z')


# Fig. 3: Scale height plot:
#    requires running
#    compute_scale_height.py to calculate the
#    gas scale height at each time step
