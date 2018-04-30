from galaxy_analysis.plot.plot_styles import *
rc('font', size=12) # make smaller for these plots
import matplotlib
#matplotlib.use('Agg')

from galaxy_analysis.analysis import Galaxy

import yt
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import sys


yt.enable_parallelism()

default_list = ['DD0046/DD0046','DD0196/DD0196', 'DD0346/DD0346', 'DD0546/DD0546']
default_fields = ['number_density', 'temperature', 'H_p0_number_density', 'H2_p0_number_density']

def plot_panel(ds_list = default_list, fields_list = default_fields, axis = 'x',
               width = (1.55,'kpc'), fsize = 4)

    ncol = len(fields_list)
    nrow = len(ds_list)

    fig = plt.figure(figsize=(fsize*nrow, fsize*ncol))

    # See http://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html
    # These choices of keyword arguments produce a four panel plot with a single
    # shared narrow colorbar on the right hand side of the multipanel plot. Axes
    # labels are drawn for all plots since we're slicing along different directions
    # for each plot.
    grid = AxesGrid(fig, (0.075,0.075,0.875,0.875),
                    nrows_ncols   = (nrow, ncol),
                    axes_pad      = 0.1,
                    label_mode    = "L",
                    share_all     = True,
                    cbar_location = "right",
                    cbar_mode     = "edge",
                    cbar_size     = "8%",
                    cbar_pad      = "0%")



    for i, fn in enumerate(ds_list):
        for j, field in enumerate(fields_list):
            index = i + j * ncol

            # Load the data and create a single plot
            gal = Galaxy(fn.split('/')[0])
            region = gal.ds.disk([0.5,0.5,0.5], [0,0,1], 2.0*yt.units.kpc, 2.0*yt.units.kpc)

            if field == 'number_density':
                p = yt.ProjectionPlot(gal.ds, axis, field, weight_field = 'number_density',
                                      width = width, data_source = region)
#            p.set_zlim(field, 1.0E-3, 1.0E3)
                p.set_zlim(field, 0.9E-3, 0.9E3)
                p.set_cmap(field, 'viridis')
                p.set_colorbar_label(field, r'n (cm$^{-3}$)')
    
                if gal.ds.parameters['NumberOfParticles'] > 0:
                    p.annotate_particles((1.0,'kpc'), p_size = 0.5, ptype='main_sequence_stars')
 
            elif field == 'temperature':
                p = yt.SlicePlot(gal.ds, axis, field, width = width, data_source = region)
#            p.set_zlim(field, 1.0E2, 1.0E7)
                p.set_zlim(field, 1.1E2, 0.9E7)
                p.set_cmap(field, 'RdYlBu_r')
                p.set_colorbar_label(field, r'T (K)')
            elif field == 'H_p0_number_density':
                p = yt.ProjectionPlot(gal.ds, axis, field, width = width, data_source = region,
                                  weight_field = None)
#            p.set_zlim(field, 1.0E16, 1.0E22)
                p.set_zlim(field, 1.1E16, 0.9E22)
                p.set_cmap(field, 'Greys')
                p.set_colorbar_label(field, r'N$_{\rm HI}$ (cm$^{-2}$)')
            elif field == 'H2_p0_number_density':
                p = yt.ProjectionPlot(gal.ds, axis, field, width = width, data_source = region,
                                  weight_field = None)
#            p.set_zlim(field, 6.0E18, 1.0E21)
                p.set_zlim(field, 6.1E18, 0.9E21)
                p.set_cmap(field, 'plasma')
                p.set_colorbar_label(field, r'N$_{\rm H_2}$ (cm$^{-2}$)')
            elif field == 'O_Number_Density':
                p = yt.ProjectionPlot(gal.ds, axis, field, width = width, data_source = region,
                                  weight_field=None)
#            p.set_zlim(field, 1.0E12, 1.0E17)
                p.set_zlim(field, 1.1E12, 0.9E17)
                p.set_cmap(field, 'magma')
                p.set_colorbar_label(field, r'N$_{\rm O}$ (cm$^{-2}$)')

            p.set_axes_unit('kpc')
            p.set_buff_size(1600)

        # This forces the ProjectionPlot to redraw itself on the AxesGrid axes.
            plot = p.plots[field]
            plot.figure = fig
            plot.axes = grid[index].axes

            if ((index+1)%(nrow)):
                plot.cax = grid.cbar_axes[index/nrow]

        # Finally, this actually redraws the plot.
            p._setup_plots()

    #for cax in grid.cbar_axes:
    #    cax.toggle_label(True)
    #    cax.axis[cax.orientation].set_label('label')
    fig.set_size_inches(fsize*nrow,fsize*ncol)
    plt.savefig('multiplot_4x4_' + axis+ '_w%.2f.png'%(width[0]))

    del(fig)
    return


if __name__ == "__main__":

    if len(sys.argv) == 1:
        plot_panel(axis = 'x')
        plot_panel(axis = 'z')
    elif len(sys.argv) > 1:
        ds_list = [  "DD%0004i/DD%0004i"%(int(i),int(i)) for i in np.array(sys.argv)[1:]]
        plot_panel(ds_list = ds_list, axis = 'x')
        plot_panel(ds_list = ds_list, axis = 'z')
