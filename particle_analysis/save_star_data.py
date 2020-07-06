import yt
import numpy as np
from galaxy_analysis.particle_analysis import particle_analysis
from galaxy_analysis.yt_fields import field_generators as fg
from galaxy_analysis.utilities import utilities

import glob
import sys
#
#
#


def save_dataset_particles(data,
                           particle_type = 'all_stars',
                           derived_fields = ['particle_model_lifetime'],
                           save_abundance_ratios = True,
                           abundance_denominators = ['H']):
    """
    Saves as dataset ALL star particle information to make 
    analysis with the star particles A LOT easier in datasets with
    many DM particles, since the way these are natively stored in Enzo
    requires loading ALL particle data for each field desired (at least once).

    Parameters
    -----------
    data   :    yt data object
    particle_type : yt-recognized particle type, string. Optional. Default 'all_stars'
    derived_fields : yt-recognized derived fields to also compute and save (for convenience).
                     Optional. Default = ['particle_model_lifetime']
    save_abundance_ratios : bool. also include derived abundance ratios for convenience.
                            Optional. Default True
    abundance_denominators : list of denominators in abundance ratios to compute.
                             Optional. Default : ['H']
    """


    particle_data = {}

    select = particle_analysis.particle_filter(particle_type, data = data)

    #
    # get native on-disk fields
    #
    fields = [x[1] for x in data.ds.field_list if x[0]=='all']

    for f in fields:
        particle_data[f] = data[('all', f)][select]

    for field in derived_fields:
        particle_data[field] = data[('all', field)][select]

    if save_abundance_ratios:
        elements = utilities.species_from_fields( data.ds.field_list )

        for denom in abundance_denominators:

            for e in elements:

                if e == denom:
                    continue

                field = 'particle_' + e + '_over_' + denom

                particle_data[field] = data[('all', field)][select]

    
    if yt.is_root():
        outname = str(ds) + "_star_particle_data.h5"
        yt.save_as_dataset(data.ds, outname, particle_data)

    return



if __name__ == "__main__":


    parallel = False
    if len(sys.argv) > 1:
        if "parallel" in sys.argv:
            parallel = True

    if parallel:
        yt.enable_parallelism()


    all_ds = np.sort(glob.glob("DD????/DD????"))

    if len(sys.argv) > 1:
        if 'datasets' in sys.argv:
            index = sys.argv.index('datasets')
            imin = int(sys.argv[index+1])
            imax = int(sys.argv[index+2])

            dslist = ["DD%0004i/DD%0004i"%(i,i) for i in np.arange(imin,imax+1,1)]

            dslist = [x for x in dslist if x in all_ds]
        else:
            dslist = all_ds



    for dsname in dslist:
        ds = fg.load(dsname)
        data = ds.all_data()

        save_dataset_particles(data)
