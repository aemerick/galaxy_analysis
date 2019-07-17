from galaxy_analysis.plot.plot_styles import *
import yt
from galaxy_analysis.analysis import Galaxy
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
from yt.visualization.base_plot_types import get_multi_plot
import matplotlib.colorbar as cb
from matplotlib.colors import LogNorm
import sys, os

from multiprocessing import Pool
from contextlib import closing
import itertools

def panel_plot(dsname):

    if not os.path.isfile(dsname + '/' + dsname):
        print(dsname + " does not exist")
        return

    gal = Galaxy(dsname)

    ds = gal.ds # yt.load(dsname + '/' + dsname)
    data = ds.all_data()


    cmaps = {'number_density' : 'viridis', 'Temperature' : 'RdYlBu_r',
             'O_Number_Density' : 'magma'}

    unit  = {'number_density' : 'cm**(-3)', 'Temperature' : 'K',
             'O_Number_Density' : 'cm**(-2)'}

    zlim  = {'number_density' : (1.0E-3, 1.0E3), 'Temperature' : (100.0,1.0E7),
             'O_Number_Density' :  (1.0E12,1.0E16)}

    labels = {'number_density' : r'log(n [cm$^{-3}$])',
              'Temperature'     : r'log(T [K])',
              'O_Number_Density'    : r'log(O Column [cm$^{-2}$])'}

    img_axes = {'number_density' : 'z', 'Temperature' : 'z', 'O_Number_Density' : 'x'}
    width    = {'number_density' : (2,'kpc'), 'Temperature' : (2,'kpc'), 'O_Number_Density' : (24,'kpc')}
    log = {'number_density' : True, 'Temperature' : True, 'O_Number_Density' : False}
    weight = {'number_density':'Density','Temperature':'Density','O_Number_Density':None}

    fields = ['number_density','Temperature','O_Number_Density']

    fig, axes, colorbars = get_multi_plot(3,1,colorbar='horizontal',bw=4)

    for ax1 in axes:
        for ax2 in ax1:
            ax2.xaxis.set_visible(False)
            ax2.yaxis.set_visible(False)

    all_plots = {}
    all_frb   = {}

    buff = 1024

    for f in fields:
        if f == 'Temperature':
            all_plots[f] = yt.SlicePlot(ds, img_axes[f], f, width = width[f])
        else:
            all_plots[f] = yt.ProjectionPlot(ds, img_axes[f], f, width = width[f],
                                         weight_field = weight[f])
        all_plots[f].set_buff_size(buff)

        all_plots[f].set_cmap(f, cmaps[f])
        all_plots[f].set_unit(f, unit[f])
        all_plots[f].set_zlim(f, zlim[f][0], zlim[f][1])
        all_plots[f].set_log(f, log[f])

        #if img_axes[f] == 'z':
        #    all_plots[f].annotate_particles(0.9)
        all_frb[f] = all_plots[f].data_source.to_frb(width[f], buff)

    time = ds.current_time.to('Myr').value - 119.0

    l = 0.5
    scale = ds.domain_width[0].to('kpc').value / 2.0
    particle_x = scale*(data['particle_position_x'] - ds.domain_center[0]).value + l #.to('kpc').value
    particle_y = scale*(data['particle_position_y'] - ds.domain_center[1]).value + l #to('kpc').value

    pt = data['particle_type']

    px = particle_x[pt == 11]
    py = particle_y[pt == 11]


    axi = axj = 0
    plots = [ [None,None,None]]
    for f in fields:
        ax = axes[axi][axj]
        plot_data = np.array( all_frb[f][f] )

        #if f == 'O_Number_Density':
        #    plots[axi][axj] = ax.imshow(plot_data, origin='lower')
        #else:
        plots[axi][axj] = ax.imshow(plot_data, origin='lower',norm=LogNorm())
        plots[axi][axj].set_clim(zlim[f])

        plots[axi][axj].set_cmap(cmaps[f])
        if axj == 0:
#            xy = (0.07,0.07)
#            ax.text(xy[0], xy[1], k, color = 'white', fontsize = 28,
#                                     transform = ax.transAxes)
            if axi == 0:
              xy = (0.07,0.9)
              ax.text(xy[0],xy[1], "%.1f Myr"%(time), color = "white", fontsize = 26,
                     transform = ax.transAxes)

        color = 'white'
        if f == 'Temperature':
            color = 'black'
        ax.plot( [0.07 + 0.05, 0.17 + 0.05], [0.05,0.05], color = color,
                 lw = 3, transform = ax.transAxes)
        xy = (0.07, 0.075)
        ax.text(xy[0], xy[1], "%.1f kpc"%(width[f][0] / 10.0), color = color,
                fontsize = 16, transform = ax.transAxes)

        if img_axes[f] == 'z':
            ax.scatter(px, py, color = 'black', s = 0.1, alpha = 0.5, transform = ax.transAxes)


        if axi == 0:
            cbar = fig.colorbar(plots[axi][axj], cax=colorbars[axj], orientation='horizontal')
            cbar.set_label(labels[f])

        axj = axj + 1
        if axj >= 3:
            axj = 0
            axi = axi + 1

    outname = "./movie_panels/" + dsname + '_panel.png'

    fig.savefig(outname)
    return True

if __name__ == "__main__":

    dsmin = int(sys.argv[1])
    dsmax = int(sys.argv[2])
    dsi   = int(sys.argv[3])

    nproc = 1

    if len(sys.argv) > 4:
        nproc = int(sys.argv[4])

    if nproc == 1:
        for i in np.arange(dsmin,dsmax,dsi):
            panel_plot( "DD%0004i"%(i))
    else:
        ds_list = ["DD%0004i"%(i) for i in np.arange(dsmin,dsmax,dsi)]

        for sub_list in itertools.zip_longest(*(iter(ds_list),) * nproc):
            sub_list = list(sub_list)
            sub_list = [s for s in sub_list if s is not None]
            reduced_nproc = np.min( [len(sub_list), nproc] )

            pool = Pool(reduced_nproc)
            results = pool.map_async(panel_plot, sub_list)

            pool.close()
            pool.join()

            del(results)
