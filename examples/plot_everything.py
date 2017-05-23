import galaxy_analysis as ga
from galaxy_analysis.static_data\
     import IMAGE_COLORBAR_LIMITS as cbar_lim

from galaxy_analysis.static_data\
     import PLOT_LIMITS as plim

from galaxy_analysis.static_data\
     import FIELD_UNITS as field_units

from galaxy_analysis.static_data\
     import LABELS as axis_label

import yt
import glob
import os
import numpy as np
import sys

from joblib import Parallel, delayed
import multiprocessing

#
# move to yt field defines
# 
def main_sequence(pfilter,data):
    filter = data[(pfilter.filtered_type, "particle_type")] == 11
    return filter

yt.add_particle_filter("main_sequence", function=main_sequence, filtered_type="all", requires=["particle_type"])


def phase_plots(ds, to_plot = 'all', region = None):

    # construct a list of phase diagrams here. List contains list of tuples, where each
    # list is the following set of 4 items:  [ x-axis field, y-axis field, cbar field, weight_field ]
    #
    pp_def    = { 'T_n' :  [ ('gas','number_density'), ('enzo','Temperature'), ('gas','cell_mass'), None],
                  'T_n_Fe' :  [ ('gas','number_density'), ('enzo','Temperature'), ('gas','Fe_Mass'), None],
                  'T_n_O'  :  [ ('gas','number_density'), ('enzo','Temperature'), ('gas','O_Mass'), None],
                  'T_n_C'  :  [ ('gas','number_density'), ('enzo','Temperature'), ('gas','C_Mass'), None],
                  'P_T' :  [ ('enzo','Temperature'), ('gas','pressure'), ('gas', 'cell_mass'), None],
                  'Go_n' : [ ('gas','number_density'), ('gas','G_o'), ('gas','cell_mass'), None],
                  'Go_r' : [ ('index','magnitude_cylindrical_radius'), ('gas','G_o'), ('gas','cell_mass'), None]
                 }

    pdict = {}
    if to_plot == 'all':
        for n in ['T_n_Fe','T_n_O','T_n','P_T','Go_n','T_n_C']:
            pdict[n] = pp_def[n]

    if region is None:
        region = ds.all_data()

    print cbar_lim

    for name in pdict:
        pdef = pdict[name]

        pp = yt.PhasePlot( region, pdef[0], pdef[1], pdef[2],
                           weight_field = pdef[3])

        if pdef[0] == ('gas','magnitude_cylindrical_radius'):
            pp.set_log('magnitude_cylindrical_radius', False)

        # set units on each axis
        for f in pdef[:3]:
            pp.set_unit(f[1], field_units[f].units)

        # set plot limits for horizonatal and vertical axes
        if not plim[pdef[0]] is None:
            pp.set_xlim(plim[pdef[0]][0], plim[pdef[0]][1])
        if not plim[pdef[1]] is None:
            pp.set_ylim(plim[pdef[1]][0], plim[pdef[1]][1])

        # set color map limits
        print cbar_lim[pdef[2]], pdef[2]
        if not (cbar_lim[pdef[2]] is None):
            pp.set_zlim( pdef[2], cbar_lim[pdef[2]][0], cbar_lim[pdef[2]][1])

        pp.set_cmap( pdef[2], 'cubehelix')

        pp.save('./phase_plots/')

        del(pp)


    return

def projection_plots(ds, fields = None, axis=['x','z'], has_particles = None, thin = False,
                         width = 2.5*yt.units.kpc, ndx = 10):

    if has_particles is None:
        has_particles = ('io','particle_position_x') in ds.field_list 

    if fields is None:
        fields = [('gas','number_density'), ('enzo','Temperature')]

        if thin:
            fields += [('gas','a_rad_over_a_grav')]

    m = 0.0
    t = ds.current_time.convert_to_units('Myr').value
    data   = ds.all_data()
    dx_min = data['dx'].convert_to_units('pc')

    if has_particles:
        m = np.sum(data['particle_mass'][data['particle_type'] == 11].convert_to_units('Msun'))
        m = m.value

    for a in axis:

        #
        # set up thin projection definitions if needed
        #
        if thin:
            center = ds.domain_center.convert_to_units('kpc')
            depth  = ndx * dx_min

            if axis == 'x':
                dl = np.array([depth, width, width])
            elif axis == 'y':
                dl = np.array([width, depth, width])
            else:
                dl = np.array([width, width, depth])

            data_source = ds.box( center - dl, center + dl)
        else:
            data_source = data

        pp = yt.ProjectionPlot(ds, axis = a, fields = fields,
                               width = width, weight_field = 'density', data_source = data_source)
        pp.set_buff_size(2048)

        for f in [('gas','number_density'),('enzo','Temperature')]:
            pp.set_cmap(f, ga.static_data.CMAPS[f])
            pp.set_unit(f, field_units[f].units)
            pp.set_zlim(f, cbar_lim[f][0], cbar_lim[f][1])
            pp.set_colorbar_label(f, axis_label[f])

        if has_particles:
            pp.annotate_particles(0.9, p_size = 0.4, stride = 1, ptype = 'main_sequence')

        pp.annotate_title(" Time = %.1f Myr     M_star = %.3E Msun"%(t,m))

        if thin:
            outdir = './thin_proj/'
        else:
            outdir = './proj/'

        pp.save(outdir)
        del(pp)

    return

def slice_plots(ds, fields, axis = ['x','z'], has_particles = None):

    if has_particles is None:
        has_particles = ('io','particle_position_x') in ds.field_list

    m = 0.0
    t = ds.current_time.convert_to_units('Myr').value

    data = ds.all_data()
    if has_particles:
        m = np.sum(data['particle_mass'][data['particle_type'] == 11].convert_to_units('Msun'))
        m = m.value


    for a in axis:
        sp   = yt.SlicePlot(ds, axis = a, fields = fields,
                                width = (2.5,'kpc'))
        sp.set_buff_size(2048)

        for f in fields:
            sp.set_cmap(f, ga.static_data.CMAPS[f])
            sp.set_unit(f, field_units[f].units)
            sp.set_zlim(f, cbar_lim[f][0], cbar_lim[f][1])
            sp.set_colorbar_label(f, axis_label[f])


        if has_particles:
            sp.annotate_particles(0.9, p_size = 0.4, stride = 1, ptype = 'main_sequence')

        sp.annotate_title(" Time = %.1f Myr     M_star = %.3E Msun"%(t,m))
        sp.save('./slice/')

        del(sp)

    return



def _parallel_loop(dsname, fields, axis = ['x','z']):


    ds   = yt.load(dsname)
    yt.add_particle_filter("main_sequence", function=main_sequence, filtered_type="all", requires=["particle_type"])
    ds.add_particle_filter("main_sequence")
    ds = yt.load(dsname)
    ds.add_particle_filter("main_sequence")

    if ds.parameters['NumberOfParticles'] < 1:
        print "Functionality currently broken when no particles exist"
        return

    ga.yt_fields.field_generators.generate_derived_fields(ds)
    ds   = yt.load(dsname)

    data = ds.all_data()

    has_particles = ('io','particle_position_x') in ds.field_list


#    phase_plots(ds)

    slice_plots(ds, fields, axis, has_particles = has_particles)

    projection_plots(ds, fields = ['number_density','temperature'], axis = axis, has_particles = has_particles, 
                     thin = True)


    del(ds)

    return


def make_plots(ds_list, fields, axis = ['x', 'z'], n_jobs = None):

    if n_jobs is None:
        n_jobs = multiprocessing.cpu_count()

    print "Starting to make %i plots for each of %i datasets with %i processors"%(len(fields) * len(axis), len(ds_list), n_jobs)

    Parallel(n_jobs = n_jobs)(\
            delayed(_parallel_loop)(dsname,fields,axis) for dsname in ds_list)

    return


if __name__ == "__main__":

    fields = [('gas','number_density'), ('enzo','Temperature'),
              ('gas','velocity_magnitude'), ('gas','G_o'),
              ('gas','Pe_heating_rate_masked')]

    if not os.path.exists('./slice'):
        os.makedirs('./slice')

    if not os.path.exists('./proj'):
        os.makedirs('./proj')

    if not os.path.exists('./thin_proj'):
        os.makedirs('./thin_proj')

    all_ds_list = glob.glob('./DD????/DD????')
    all_ds_list = np.sort(all_ds_list)

    n_jobs = None

    if len(sys.argv) == 1:
        ds_list = all_ds_list

    elif len(sys.argv) == 2:
        ds_list = all_ds_list
        n_jobs  = int(sys.argv[1])

    elif len(sys.argv) == 3:
        imin, imax = [ int(x) for x in sys.argv[1:]]
        ds_list = all_ds_list[imin:imax]

    elif len(sys.argv) == 4:
        imin, imax, di = [int(x) for x in sys.argv[1:]]

        ds_list = all_ds_list[np.arange(imin,imax,di)]

    elif len(sys.argv) == 5:
        imin, imax, di, n_jobs = [int(x) for x in sys.argv[1:]]
        ds_list = all_ds_list[np.arange(imin,imax,di)]

    make_plots(ds_list, fields, n_jobs = n_jobs)
        
