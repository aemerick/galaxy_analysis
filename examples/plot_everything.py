import galaxy_analysis as ga
from galaxy_analysis.static_data\
     import IMAGE_COLORBAR_LIMITS as cbar_lim

from galaxy_analysis.static_data\
     import FIELD_UNITS as field_units

from galaxy_analysis.static_data\
     import LABELS as cbar_label

import yt
import glob
import os
import numpy as np
import sys

from joblib import Parallel, delayed
import multiprocessing

def _parallel_loop(dsname, fields, axis = ['x','z']):


    ds   = yt.load(dsname)
    ga.yt_fields.field_generators.generate_derived_fields(ds)
    ds   = yt.load(dsname)

    data = ds.all_data()

    for a in axis:
        sp   = yt.SlicePlot(ds, axis = a, fields = fields, 
                                width = (2.5,'kpc'))
        sp.set_buff_size(1664)

        for f in fields:
            sp.set_cmap(f, ga.static_data.CMAPS[f])
            sp.set_unit(f, field_units[f].units)
            sp.set_zlim(f, cbar_lim[f][0], cbar_lim[f][1])
            sp.set_colorbar_label(f, cbar_label[f])

        if ('io','particle_position_x') in ds.field_list:
            sp.annotate_particles(0.9, p_size = 0.25)

        sp.save('./slice/')

    del(sp)
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
        
