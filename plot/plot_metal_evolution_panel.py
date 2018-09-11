import yt
from yt import derived_field

import numpy as np
import matplotlib.pyplot as plt
from galaxy_analysis import Galaxy
import glob
from collections import OrderedDict
from yt.visualization.base_plot_types import get_multi_plot
import matplotlib.colorbar as cb
from matplotlib.colors import LogNorm
import sys, os

from multiprocessing import Pool
from contextlib import closing
import itertools



#
# As Copied from stack overflow with some minor renaming edits:
#
colors     = {'massive_star_winds' : 'white', 'AGB_winds' : 'C1', 'SN' : 'black', 'other_stars' : 'black'}
markers    = {'massive_star_winds' : '*',     'AGB_winds' : 'D', 'SN' : '*', 'other_stars' : '.'}
ps         = {'massive_star_winds' :  75, 'AGB_winds' : 300, 'SN' : 300, 'other_stars' : 15}
all_fields = ['number_density','Temperature','N_over_O_filtered','O_over_H_filtered', 'G_o', 'Q0_flux']

tol      = 1.0E-25
import matplotlib.colors as mpl_colors
class MidpointNormalize(mpl_colors.Normalize):
    """
    Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        mpl_colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

def make_filtered_field(ds, fieldname, filter_fields = [], tolerance = tol):
    """
    Creates a derived field that filters 'fieldname' if any field in
    'filter_fields' has a value below 'tolerance'. Filtered by setting
    cells where this is true to np.nan (i.e. they will show up as
    white / blank in projection / slice plots).

    Derived field will have the name  : ('gas',fieldname + "_filtered")
    Assumes that the original field is: ('gas',fieldname)
    """
    def _filtered_field(field, data):
        x = data[('gas',fieldname)]

        select = data[filter_fields[0]] < 0
        for f in filter_fields:
            select = select + (data[f] < tolerance)
        x[select] = np.nan

        return x

    ds.add_field(('gas',fieldname + '_filtered'), function = _filtered_field, units = "")

    return


def panel_plot(dsname, width = 1000.0, thickness = 20.0):

    cmaps = {'number_density' : 'viridis', 'Temperature' : 'RdYlBu_r',
             'N_over_O_filtered' : 'PRGn', 'O_over_H_filtered' : 'cubehelix'}

    unit  = {'number_density' : 'cm**(-3)', 'Temperature' : 'K'}

    zlim  = {'number_density' : (1.0E-3, 1.0E3), 'Temperature' : (100.0,1.0E7),
             'N_over_O_filtered' : (-2,2), 'O_over_H_filtered' : (-5,1)}


    labels = {'number_density' : r'log(n [cm$^{-3}$])',
              'Temperature'     : r'log(T [K])',
              'N_over_O_filtered'    : r'[N/O]', 'O_over_H_filtered' : r'[O/H]'}

    log = {'number_density' : True, 'Temperature' : True,
           'N_over_O_filtered' : False, 'O_over_H_filtered': False}


    if not os.path.isfile(dsname + '/' + dsname):
        print dsname + " does not exist"
        return

    gal = Galaxy(dsname)

    ds = gal.ds # yt.load(dsname + '/' + dsname)
    data = ds.all_data()

#    @derived_field(name="logNO", units="")
#    def _logNO(field, data):
#        return  np.log10(data['N_Abundance'] / data['O_Abundance'])
#    gal.ds.add_field(("gas", "logNO"), function=_logNO, units="")

#    make_filtered_field(gal.ds, 'logNO', ['O_Fraction','N_Fraction'])
    make_filtered_field(gal.ds, 'O_over_H', ['O_Fraction'])
    make_filtered_field(gal.ds, 'N_over_O', ['O_Fraction','N_Fraction'])

    dt = 5.0 * yt.units.Myr

    M            = data['birth_mass']
    t_o          = data['creation_time'].convert_to_units('Myr')
    MS_lifetime  = data[('io','particle_model_lifetime')].to('Myr')
    MS_death     = t_o + MS_lifetime
    px           = (data['particle_position_x'] - gal.ds.domain_center[0]).to('pc')
    py           = (data['particle_position_y'] - gal.ds.domain_center[1]).to('pc')
    pz           = (data['particle_position_z'] - gal.ds.domain_center[2]).to('pc')

    recent_death = (MS_death > gal.ds.current_time - dt) * (MS_death <= gal.ds.current_time + 0.001*yt.units.Myr)
    alive        = MS_death > gal.ds.current_time + 0.001*yt.units.Myr

    AGB           = M < 8.0
    massive_star  = (M > 8.0) * (M < 25.0)

    boxdim = np.array([width*1.25,width*1.25,thickness])*yt.units.pc
    region = gal.ds.box(gal.ds.domain_center - boxdim*0.5, gal.ds.domain_center + boxdim*0.5)


    fields = ["number_density","Temperature","N_over_O_filtered","O_over_H_filtered"]


    fig, axes, colorbars = get_multi_plot(4,1, colorbar="horizontal", bw=4)

    for ax1 in axes:
        for ax2 in ax1:
            ax2.xaxis.set_visible(False)
            ax2.yaxis.set_visible(False)

    all_plots = {}
    all_frb   = {}

    buff = 1024

    for f in fields:
        all_plots[f] = yt.ProjectionPlot(gal.ds, 'z', f, weight_field = 'number_density',
                                    data_source = region,
                                    width = (width,'pc'))

        all_plots[f].set_buff_size(buff)

        all_plots[f].set_cmap(f, cmaps[f])
        if f in unit.keys():
            all_plots[f].set_unit(f, unit[f])
        all_plots[f].set_zlim(f, zlim[f][0], zlim[f][1])
        all_plots[f].set_log(f, log[f])

        #if img_axes[f] == 'z':
        #    all_plots[f].annotate_particles(0.9)
        all_frb[f] = all_plots[f].data_source.to_frb((width,'pc'), buff)


    dt = 5.0 * yt.units.Myr
    # buffer around image. otherwise points plotted near edge of image my run a little outside
    # viewing area, causing weird shifts in plotting. Not sure how to control this otherwise
    buffer = 15.0 # in pc
    in_image     = (np.abs(pz) <= boxdim[2]*0.5) *\
                   (np.abs(px) <= (width*0.5 - buffer)) *\
                   (np.abs(py) <= (width*0.5 - buffer))

    pp = {}
    pp['massive_star_winds'] = in_image * alive * massive_star
    pp['AGB_winds']          = in_image * recent_death * AGB
    pp['SN']                 = in_image * recent_death * massive_star

    l = 0.5
    scale = ds.domain_width[0].to('pc').value / width
    particle_x = scale*(data['particle_position_x'] - ds.domain_center[0]).value + l #.to('kpc').value
    particle_y = scale*(data['particle_position_y'] - ds.domain_center[1]).value + l #to('kpc').value

    axi = axj = 0
    plots = [ [None,None,None,None]]
    time = gal.ds.current_time.to('Myr').value

    for f in fields:
        ax = axes[axi][axj]
        plot_data = np.array( all_frb[f][f] )

        if not log[f]:
            plots[axi][axj] = ax.imshow(plot_data, origin='lower')
        else:
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
        #if f == 'Temperature':
        #    color = 'black'
            ax.plot( [0.07 + 0.05, 0.17 + 0.05], [0.05,0.05], color = color,
                 lw = 3, transform = ax.transAxes)
            xy = (0.07, 0.075)
            ax.text(xy[0], xy[1], "%.1f pc"%(width / 10.0), color = color,
                    fontsize = 16, transform = ax.transAxes)

        #if img_axes[f] == 'z':
        for s in pp.keys():
            if np.size( particle_x[pp[s]]) > 0:
                ax.scatter(particle_x[pp[s]], particle_y[pp[s]], s = ps[s],
                           marker = markers[s], color = colors[s],
                           transform = ax.transAxes)

                print s, particle_x[pp[s]]

        if axi == 0:
            if f == 'N_over_O_filtered':
                cbar = fig.colorbar(plots[axi][axj], cax = colorbars[axj],
                                    ticks = [-2,-1,0,1,2], orientation = 'horizontal')
                cbar.ax.set_xticklabels(['-2','-1','0','1','2'])
                cbar.set_label(labels[f])

            else:
                cbar = fig.colorbar(plots[axi][axj], cax=colorbars[axj],
                                    orientation='horizontal')
                cbar.set_label(labels[f])

        axj = axj + 1
        if axj >= 4:
            axj = 0
            axi = axi + 1


    outname = "./metal_evolution/" + dsname + '_panel.png'
    fig.savefig(outname)
    plt.close()

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

        for sub_list in itertools.izip_longest(*(iter(ds_list),) * nproc):
            sub_list = list(sub_list)
            sub_list = [s for s in sub_list if s is not None]
            reduced_nproc = np.min( [len(sub_list), nproc] )

            pool = Pool(reduced_nproc)
            results = pool.map_async(panel_plot, sub_list)

            pool.close()
            pool.join()

            del(results)
