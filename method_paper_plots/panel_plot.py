from galaxy_analysis.plot.plot_styles import *
rc('font', size=12) # make smaller for these plots
from galaxy_analysis.analysis import Galaxy

import yt
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

yt.enable_parallelism()

fns = ['DD0050/DD0050','DD0100/DD0100','DD0150/DD0150','DD0187/DD0187']
axis = 'z'
fields_list = ['number_density', 'temperature', 'H_p0_number_density', 'O_Number_Density']
nrow,ncol = 4,4

#nrow,ncol = 2,2
#fns = ['DD0100/DD0100','DD0150/DD0150']
#fields_list = ['number_density','temperature']

fig = plt.figure()

# See http://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html
# These choices of keyword arguments produce a four panel plot with a single
# shared narrow colorbar on the right hand side of the multipanel plot. Axes
# labels are drawn for all plots since we're slicing along different directions
# for each plot.
grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols   = (nrow, ncol),
                axes_pad      = 0.075,
                label_mode    = "L",
                share_all     = True,
                cbar_location = "right",
                cbar_mode     = "edge",
                cbar_size     = "7%",
                cbar_pad      = "2%")


# region = ds.disk([0.5,0.5,0.5], [0,0,1], 2.0*yt.units.kpc, 2.0*yt.units.kpc)
width  = (1.5,'kpc')
fig.set_size_inches(8*nrow,8*ncol)

for i, fn in enumerate(fns):
    for j, field in enumerate(fields_list):
        index = i + j * ncol

        # Load the data and create a single plot
        gal = Galaxy(fn.split('/')[0])

        region = gal.ds.disk([0.5,0.5,0.5], [0,0,1], 2.0*yt.units.kpc, 2.0*yt.units.kpc)

        if field == 'number_density':
            p = yt.ProjectionPlot(gal.ds, axis, field, weight_field = 'number_density',
                                  width = width, data_source = region)
            p.set_zlim(field, 1.0E-3, 1.0E3)
            p.set_cmap(field, 'viridis')
            if gal.ds.parameters['NumberOfParticles'] > 0:
                p.annotate_particles((1.0,'kpc'), p_size = 0.5, ptype='main_sequence_stars')

        elif field == 'temperature':
            p = yt.SlicePlot(gal.ds, axis, field, width = width, data_source = region)
            p.set_zlim(field, 1.0E2, 1.0E7)
            p.set_cmap(field, 'RdYlBu_r')
        elif field == 'H_p0_number_density':
            p = yt.ProjectionPlot(gal.ds, axis, field, width = width, data_source = region,
                                  weight_field = None)
            p.set_zlim(field, 1.0E16, 1.0E22)
            p.set_cmap(field, 'Greys')
        elif field == 'O_Number_Density':
            p = yt.ProjectionPlot(gal.ds, axis, field, width = width, data_source = region,
                                  weight_field=None)
            p.set_zlim(field, 1.0E12, 1.0E17)
            p.set_cmap(field, 'magma')

        p.set_axes_unit('kpc')
        p.set_buff_size(1600)

        # This forces the ProjectionPlot to redraw itself on the AxesGrid axes.
        plot = p.plots[field]
        plot.figure = fig
        plot.axes = grid[index].axes
#        plot.cax = grid.cbar_axes[index]

        if ((index+1)%(nrow)):
            # grid.cbar_axes[index/4].colorbar(plot.image)
        # try
            plot.cax = grid.cbar_axes[index/nrow]
        print np.size(grid.cbar_axes)

        # Finally, this actually redraws the plot.
        p._setup_plots()

#for cax in grid.cbar_axes:
#    cax.toggle_label(True)
#    cax.axis[cax.orientation].set_label('label')

plt.savefig('multiplot_4x4_' + axis+ '.png')
