from galaxy_analysis.analysis import time_average_phase_diagram as tapd
from galaxy_analysis.utilities import utilities
from galaxy_analysis.analysis import Galaxy
from galaxy_analysis.yt_fields import field_generators as fg
import yt
import numpy as np
import glob as glob
import os, sys

print(yt.enable_parallelism() , " ------------------------")


def plot_time_average_PD(wdir, t_min, t_max, nbin = 100,
                         plots = ['nT','G_o','Q_o'], outdir = None):

    x = utilities.select_data_by_time(wdir, tmin=0.0, tmax=np.inf)
    sim_files = {'files': x[0], 'Time' : x[1]}
    #                    'DD_files' : np.sort(glob.glob(wdir + '/DD????/DD????'))}
    #iremove = len(sim_files['files'])
    #sim_files['DD_files'] = sim_files['DD_files'][-iremove:]
    # --- get data sets associated with output files (this is gross - I am sorry)
    sim_files['DD_files'] =  np.array([ x.split('_galaxy_data')[0] + '/' + x.split('_galaxy_data')[0][-6:] for x in sim_files['files']])
    # --- make sure they exist
    # sim_files['DD_files'] =  [ x for x in sim_files['DD_files'] if os.path.isfile(x)]
    # --- select the ones in the time range we want
    sim_files['phase_files'] = sim_files['DD_files'][ (sim_files['Time'] >= t_min) *\
                                                      (sim_files['Time'] <= t_max)]

    gal = Galaxy(sim_files['DD_files'][0].split('/')[-1])

    if 'G_o' in plots:
        pd = tapd.time_average_phase_diagram(t_min, t_max,
                                      ds_list = sim_files['phase_files'], wdir = wdir,
                                         x_bins = np.linspace(0.0, 600.0, nbin)*yt.units.pc,
                                         xlog  = False, xfield = 'cylindrical_radius', ylabel = r'G$_{\rm o}$',
                                         y_bins = np.logspace(np.log10(0.00324), 3, nbin)*yt.units.pc/yt.units.pc, yfield = 'G_o',
                                         region_type = 'Disk', region_kwargs = {'normal':np.array([0.0,0.0,1.0]),
                                                                               'radius': 600*yt.units.pc,
                                                                               'height': (1.5*yt.units.pc*20),
                                                                               'center': np.array([0.5,0.5,0.5])},
                                         zlog=True, zlabel = r'Mass (M$_{\odot}$)',
#                                         region_type = 'FullBox', zlog = True,
                                         ylog=True,
                                         zlim = [5.0E-2, 1.0E5], zunit = 'Msun', cmap = 'cubehelix',
#                                         outname = 'G_o_r_cell_mass_2D_phase.png',
                                         outdir  = outdir)

    if ('Q_o' in plots) or ('Q0' in plots) or ('Q_0' in plots):
        pd = tapd.time_average_phase_diagram(t_min, t_max,
                                     ds_list = sim_files['phase_files'], wdir = wdir,
                                     x_bins = np.linspace(0.0, 600.0, nbin)*yt.units.pc,
                                     xlog  = False, xfield = 'cylindrical_radius', ylabel = r'HI Ionizing Flux (erg s$^{-1}$ cm$^{-2}$)',
                                     y_bins = np.logspace(-8, 1, nbin)*yt.units.erg/yt.units.s/yt.units.cm**2, yfield = 'Q0_flux',
                                     region_type = 'Disk', region_kwargs = {'normal':np.array([0.0,0.0,1.0]),
                                                                           'radius': 600*yt.units.pc,
                                                                           'height': (1.5*yt.units.pc)*20,
                                                                            'center': np.array([0.5,0.5,0.5])},
                                     zlog=True, zlabel = r'Mass (M$_{\odot}$)',
#                                         region_type = 'FullBox', zlog = True,
                                     ylog=True,
                                     zlim = [5.0E-2, 1.0E5], zunit = 'Msun', cmap = 'cubehelix',
#                                     outname = 'Q_o_r_cell_mass_2D_phase.png',
                                     outdir  = outdir)

    if 'nT' in plots:
        pd = tapd.time_average_phase_diagram(t_min, t_max,
                                     ds_list = sim_files['phase_files'], wdir = wdir,
                                     x_bins = np.logspace(-6, 3, 512)*yt.units.cm**(-3),
                                     xlog  = True, xfield = 'number_density',
                                     y_bins = np.logspace(0,7.5, 512) * yt.units.K, yfield = 'Temperature',
                                     region_type = 'sphere', region_kwargs = {'radius': 3.6*yt.units.kpc,
                                                                              'center': np.array([0.5,0.5,0.5])},
                                     zlog=True,
#                                      region_type = 'FullBox', zlog = True,
                                         ylog=True,
                                     zlim = [1.0E-2, 6.0E4], zunit = 'Msun', cmap = 'cubehelix',
#                                     outname = 'T_n_cell_mass_2D_phase.png',
                                     outdir  = outdir)


    return

if __name__ == "__main__":
    # do something
    t_min = float(sys.argv[1])
    t_max = float(sys.argv[2])

    if len(sys.argv) > 3:
        plots = [str(x) for x in sys.argv[3:]]
    else:
        plots = ['nT','G_o','Q_o']

    plot_time_average_PD('./', t_min, t_max, nbin = 100,
                         plots = plots)
