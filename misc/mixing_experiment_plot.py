import numpy as np
from galaxy_analysis.plot.plot_styles import *
import matplotlib.pyplot as plt
import yt
from yt.fields.api import ValidateParameter
from matplotlib.colors import LogNorm
from galaxy_analysis.analysis import Galaxy


from multiprocessing import Pool
from contextlib import closing
import itertools
import os
import glob
import sys

def plot_field(dsname, mass_field, paperstyle=False, show_colorbar=True):

    fontsize = 'large'

    field_to_plot   = mass_field + "_disk_fraction"
    field_parameter = "total_field_mass"

    gal   = Galaxy(dsname)

    dims = [700,700,56]
    pix  = np.array(dims) / 2.0
    dist = pix * np.min( gal.df['dx'].to('pc'))

    LE = gal.ds.domain_center - dist

    cg = gal.ds.smoothed_covering_grid(level     = int(np.max(gal.df['grid_level'])),
                              left_edge = LE,
                              dims      = dims)

    disk_mass = np.sum( (gal.disk[mass_field].to('Msun').value) )
    cg.set_field_parameter(field_parameter,
                            disk_mass * yt.units.Msun)

    def _disk_mass_fraction(field,data):
        if data.has_field_parameter(field_parameter):
            mtot = data.get_field_parameter(field_parameter)
            if hasattr(mtot, 'convert_to_units'):
                mtot = mtot.convert_to_units('Msun')
            else:
                mtot = mtot * yt.units.Msun
        else:
            raise ValueError
        mfield = mass_field
        return data[mfield].convert_to_units('Msun') / mtot

    gal.ds.add_field(field_to_plot, function = _disk_mass_fraction,
                     units = "", take_log = False, force_override=True,
                     validators=[ValidateParameter(field_parameter)])

    # sum the data over 2nd axes
    axis = 2
    #print np.shape(cg[field_to_plot])
    proj_data       = np.sum(cg[field_to_plot], axis=axis)
    n_proj_data     = np.max(cg['number_density'].value, axis=axis)
    #print np.shape(proj_data)

    #print np.min(proj_data), np.max(proj_data)

    if axis == 2:
        extent = [-dist[0],dist[0],-dist[1],dist[1]]
        fsize  = (8,8)
    elif axis == 0 or axis == 1:
        extent = [-dist[1],dist[1],-dist[2],dist[2]]
        proj_data = proj_data.T
        fsize  = (12.5,1)
        n_proj_data = n_proj_data.T

    pp = plt.imshow( proj_data, norm = LogNorm(), extent=extent)
    pp.set_clim(1.0E-10,1.0E-1)
    pp.set_cmap("magma")
    #plt.tight_layout()
    #pp.axes.set_xticks([-600,-400,-200,0,200,400,600])
    pp.axes.yaxis.set_visible(False)
    pp.axes.xaxis.set_visible(False)
    pp.axes.set_frame_on(False)
    pp.figure.set_size_inches(fsize)
    pp.figure.patch.set_facecolor('black')
    # annotate the size
    anncoord = (0.9,0.1)
    dim = np.shape(proj_data)
    x  = 0.5
    di = 250.0 / ( 1.8 )
    pp.axes.plot( [dim[0]*x, dim[0]*x + di],
                  [-dim[1]*0.825]*2, lw = 3, color = "white")

    pp.axes.annotate('250 pc',
                     xy=anncoord, xycoords='axes fraction',
                     xytext=anncoord, textcoords='axes fraction', color = "white",
                     # arrowprops=dict(facecolor='black', shrink=0.05),
                     horizontalalignment='right', verticalalignment='top', fontsize = fontsize)

    start = 220.0
    anncoord2 = (0.225,0.96)
    if not paperstyle:
        pp.axes.annotate("%.1f Myr"%(gal.ds.current_time.to('Myr').value - start),
                        xy=anncoord2, xycoords='axes fraction',
                        xytext=anncoord2, textcoords='axes fraction', color = "white",
                    # arrowprops=dict(facecolor='black', shrink=0.05),
                        horizontalalignment='right', verticalalignment='top', fontsize = fontsize)


#    pp.figure.savefig("./mixing_plots/" + dsname + "_" + mass_field + "_disk_fraction.png")
    outname = "./mixing_plots/" + dsname + "_" + mass_field + "_disk_fraction.png"
    if show_colorbar:
        cax = plt.colorbar(orientation="horizontal")
        cticks =  [1.0E-8, 1.0E-6, 1.0E-4, 1.0E-2, 1.0]
        cticks =  [1.0E-9, 1.0E-7, 1.0E-5, 1.0E-3, 1.0E-1]
        cax.set_ticks(cticks)
        cax.set_ticklabels( [r"10$^{%2i}$"%(np.log10(x)) for x in cticks] )
        cax.set_label("Tracer Fraction (Arbitrary Units)")


    plt.savefig(outname, bbox_inches="tight", pad_inches = 0.0)
    plt.close()

############
    pp = plt.imshow( n_proj_data, norm = LogNorm(),
                     extent=extent)
    pp.set_clim(1.0E-4,1.0E3)
    #plt.tight_layout()
    #pp.axes.set_xticks([-600,-400,-200,0,200,400,600])
    pp.axes.yaxis.set_visible(False)
    pp.axes.xaxis.set_visible(False)
    pp.axes.set_frame_on(False)
    pp.figure.set_size_inches(fsize)
    pp.figure.patch.set_facecolor('black')
    # annotate the size
    anncoord = (0.9,0.1)
    dim = np.shape(n_proj_data)
    x  = 0.5
    di = 250.0 / ( 1.8 )
    pp.axes.plot( [dim[0]*x, dim[0]*x + di],
                  [-dim[1]*0.825]*2, lw = 3, color = "white")


    pp.axes.annotate('250 pc',
                     xy=anncoord, xycoords='axes fraction',
                     xytext=anncoord, textcoords='axes fraction', color = "white",
                     # arrowprops=dict(facecolor='black', shrink=0.05),
                     horizontalalignment='right', verticalalignment='top', fontsize = fontsize)


    start = 220.0
    anncoord2 = (0.225,0.96)
    pp.axes.annotate("%.1f Myr"%(gal.ds.current_time.to('Myr').value - start),
                    xy=anncoord2, xycoords='axes fraction',
                    xytext=anncoord2, textcoords='axes fraction', color = "white",
                # arrowprops=dict(facecolor='black', shrink=0.05),
                    horizontalalignment='right', verticalalignment='top', fontsize = fontsize)


#    pp.figure.savefig("./mixing_plots/" + dsname + "_" + mass_field + "_disk_fraction.png")
    outname = "./mixing_plots/" + dsname + "_" + "ndens.png"
    if show_colorbar:
        cax = plt.colorbar(orientation = "horizontal")
        cticks =  [1.0E-4, 1.0E-2, 1.0, 100.0]
        cax.set_ticks(cticks)
        cax.set_ticklabels( [r"10$^{%2i}$"%(np.log10(x)) for x in cticks] )
        cax.set_label(r"n (cm$^{-3}$)")

    plt.savefig(outname, bbox_inches="tight", pad_inches = 0.0)

    return


if __name__ == "__main__":

    mass_field = "C_Mass"
    nproc      = 1
    ds_list    = None

    if len(sys.argv) >= 2:
        mass_field = str(sys.argv[1])

    if len(sys.argv) >= 4:
        imin = int(sys.argv[2])
        imax = int(sys.argv[3])
    else:
        ds_list = np.sort(glob.glob('DD????'))

    if len(sys.argv) >= 5:
        di   = int(sys.argv[4])

    if len(sys.argv) >= 6:
        nproc = int(sys.argv[5])

    if ds_list is None:
        ds_list = ["DD%0004i"%(i) for i in np.arange(imin,imax,di)]

    if nproc > 1:

        def _parallel_loop(dsname):
            plot_field(dsname, mass_field)
            return

        for sub_list in itertools.izip_longest(*(iter(ds_list),) * nproc):
            sub_list = list(sub_list)
            sub_list = [s for s in sub_list if s is not None]
            reduced_nproc = np.min( [len(sub_list), nproc] )

            pool    = Pool(reduced_nproc)
            results = pool.map_async(_parallel_loop, sub_list)
            pool.close()
            pool.join()

    else:
        for dsname in ds_list:
            plot_field(dsname, mass_field)
