import numpy as np
import yt
from yt.extensions.astro_analysis.halo_analysis.api import HaloCatalog
from yt.data_objects.particle_filters import add_particle_filter
from yt.extensions.astro_analysis.halo_finding.rockstar.api import RockstarHaloFinder
from yt.fields.api import ValidateParameter

import sys

yt.enable_parallelism() # must be true for rockstar

@particle_filter("dark_matter", requires=["particle_type"])
def DarkMatter(pfilter,data):
    filter = data[('all','particle_type')] == 1 # DM in Enzo
    return filter

@particle_filter("max_res_dark_matter", requires=["particle_type"])
def MaxResDarkMatter(pfilter,dobj):
    #min_dm_mass = dobj.get_field_parameter('min_dm_mass')
    min_dm_mass = dobj.ds.parameters['min_dm_mass'] * yt.units.Msun
    return dobj['particle_mass'] <= 1.01 * min_dm_mass

def setup_ds(ds):
    ds.add_particle_filter("dark_matter")
    ds.add_particle_filter("max_res_dark_matter")
    return


assert(yt.communication_system.communicators[-1].size >= 3)

ROCKSTAR_OUTPUT_PREFIX = 'rockstar_halos/halos_'
HALOCATALOG_PREFIX     = 'halo_catalogs/catalog_'

def find_halos(dsname, wdir = './', *args, **kwargs):
    """
    Find the halos in a dataset. For simplicity
    this ONLY does the halo finding. We do other computations
    elsewhere.


    """
    ds = yt.load(wdir + dsname + '/' + dsname)
    setup_ds(ds)

    data = ds.all_data()
    min_dm_mass = data.quantities.extrema(('dark_matter','particle_mass'))[0].to('Msun').value

    ds.parameters['min_dm_mass'] = min_dm_mass

    ds.add_particle_filter('dark_matter')

    hc = HaloCatalog(data_ds = ds, finder_method = 'rockstar',
                     finder_kwargs = {'dm_only':True,
                                      'outbase' : ROCKSTAR_OUTPUT_PREFIX + str(ds),
                                      'particle_type':'dark_matter'},
                     output_dir = wdir + HALOCATALOG_PREFIX +str(ds))

    hc.create()

    #rhf = RockstarHaloFinder(ds, particle_type='max_res_dark_matter')
    #rhf.run()

    #halos = yt.load('rockstar_halos/halos_0.0.bin')
    #hc = HaloCatalog(halos_ds=halos, output_dir = widr + 'halo_catalogs/catalog_'+str(ds)))
    #hc.load()

    del(hc)
    del(ds)
    del(data)

    return

def compute_virial_quantities(dsname, wdir = './', *args, **kwargs):
    """
    Assuming halos are already generated, compute and store their
    virial quantities.
    """
    data_ds = yt.load(wdir + dsname + '/' + dsname)
    halos_ds = yt.load(wdir + ROCKSTAR_OUTPUT_PREFIX + dsname + '/halos_0.0.bin')

    hc = HaloCatalog(data_ds = data_ds, halos_ds = halos_ds)
    hc.add_filter('quantity_value', 'particle_mass', '>', 1E4, 'Msun')

    if ('enzo','Density') in data_ds.field_list:
        mass_field = 'matter_mass'
    else:
        # DM only simulation
        mass_field = 'particle_mass'
        
    hc.add_recipe("calculate_virial_quantities", ["radius", mass_field ])
    hc.create()

    return


if __name__ == '__main__':

    #dsname = 'DD0011'
    for dsname in [str(sys.argv[1])]:
        find_halos(dsname)
        compute_virial_quantities(dsname)
