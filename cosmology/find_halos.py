import numpy as np
import yt
yt.enable_parallelism()

from yt.extensions.astro_analysis.halo_analysis.api import HaloCatalog
from yt.data_objects.particle_filters import add_particle_filter
from yt.extensions.astro_analysis.halo_finding.rockstar.api import RockstarHaloFinder
from yt.data_objects.particle_filters import particle_filter
from yt.fields.api import ValidateParameter

from yt.extensions.astro_analysis.halo_analysis.halo_recipes import add_recipe

import sys, glob

def my_calculate_virial_quantities(hc, fields,
                                weight_field=None, accumulation=True,
                                radius_field="radius", factor=2.0,
                                overdensity_field=("gas", "overdensity"),
                                critical_overdensity=200):
    storage = "virial_quantities_profiles"
    pfields = [field for field in fields if field != radius_field]

    hc.add_callback("sphere", factor=factor)
    if pfields:
        hc.add_callback("profile", [radius_field], pfields,
                        weight_field=weight_field,
                        accumulation=accumulation,
                        storage=storage)
    hc.add_callback("profile", [radius_field], [overdensity_field],
                    weight_field="cell_volume", accumulation=True,
                    storage=storage)
    hc.add_callback("virial_quantities", fields,
                    overdensity_field=overdensity_field,
                    critical_overdensity=critical_overdensity,
                    profile_storage=storage)
    hc.add_callback("delete_attribute", storage)
    hc.add_callback("delete_attribute", "data_object")
    return

add_recipe("my_calculate_virial_quantities", my_calculate_virial_quantities)


def DarkMatter(pfilter,data):
    filter = data[('all','particle_type')] == 1 # DM in Enzo
    return filter

def MaxResDarkMatter(pfilter,dobj):
    min_dm_mass = dobj[('dark_matter','particle_mass')]
    #min_dm_mass = dobj.get_field_parameter('min_dm_mass')
    #min_dm_mass = dobj.ds.parameters['min_dm_mass'] * yt.units.Msun
    return dobj['particle_mass'] <= 1.01 * min_dm_mass

add_particle_filter("dark_matter", function=DarkMatter, \
                    filtered_type='all', requires=["particle_type"])

add_particle_filter("max_res_dark_matter", function=MaxResDarkMatter, \
                    filtered_type='dark_matter', requires=["particle_mass"])

def setup_ds(ds):
    ds.add_particle_filter("dark_matter")
    ds.add_particle_filter("max_res_dark_matter")
    return


#assert(yt.communication_system.communicators[-1].size >= 3)

ROCKSTAR_OUTPUT_PREFIX = 'rockstar_halos/halos_'
HALOCATALOG_PREFIX     = 'halo_catalogs/catalog_'

def find_halos(dsname, simfile = None, wdir = './', restart = False, *args, **kwargs):
    """
    Find the halos in a dataset. For simplicity
    this ONLY does the halo finding. We do other computations
    elsewhere.


    """
    es = yt.simulation(simfile, "Enzo")
    es.get_time_series(initial_redshift=30.0)
 
    print(es[0].derived_field_list)
    for ds in es:  
        setup_ds(ds)
    print(es[0].derived_field_list)

    #ds = yt.load(wdir + dsname + '/' + dsname)
    #setup_ds(ds)

    #data = ds.all_data()
    #min_dm_mass = data.quantities.extrema(('dark_matter','particle_mass'))[0].to('Msun').value

    #ds.parameters['min_dm_mass'] = min_dm_mass

    #hc = HaloCatalog(data_ds = ds, finder_method = 'rockstar',
    #                 finder_kwargs = {'dm_only':True,
    #                                  'outbase' : ROCKSTAR_OUTPUT_PREFIX + str(ds),
    #                                  'particle_type':'dark_matter'},
    #                 output_dir = wdir + HALOCATALOG_PREFIX +str(ds))

    #hc.create()

    rhf = RockstarHaloFinder(es, 
                             #num_readers=1, num_writers=1,
                             particle_type='max_res_dark_matter')
    rhf.run(restart=restart)

    #halos = yt.load('rockstar_halos/halos_0.0.bin')
    #hc = HaloCatalog(halos_ds=halos, output_dir = widr + 'halo_catalogs/catalog_'+str(ds)))
    #hc.load()

    #del(hc)
    #del(ds)
    #del(data)

    return

def compute_virial_quantities(dsname, wdir = './', *args, **kwargs):
    """
    Assuming halos are already generated, compute and store their
    virial quantities.
    """
    data_ds = yt.load(wdir + dsname + '/' + dsname)
    halos_ds = yt.load(wdir + ROCKSTAR_OUTPUT_PREFIX + dsname + '/halos_0.0.bin')

    hc = HaloCatalog(data_ds = data_ds, halos_ds = halos_ds,
                     output_dir = wdir + HALOCATALOG_PREFIX + str(data_ds))
    hc.add_filter('quantity_value', 'particle_mass', '>', 1E4, 'Msun')

    if ('enzo','Density') in data_ds.field_list:
        mass_field    = 'matter_mass'
        radius_field  = "radius"
    else:
        # DM only simulation
        mass_field = ('all',"particle_mass")
        radius_field = ('all','particle_radius')
        
    hc.add_recipe("my_calculate_virial_quantities", [radius_field, mass_field ], radius_field=radius_field)
    hc.create()

    return


if __name__ == '__main__':
    simfile = None

#    if len(sys.argv) > :
#        dsnames = ["RD%0004i"%(i) for i in np.arange(int(sys.argv[1]),int(sys.argv[2])+1,1)]
    if str(sys.argv[1]) == 'all':
        dsnames = np.sort(glob.glob('RD????'))
    elif "." in str(sys.argv[1]):
        dsnames = None
        simfile = str(sys.argv[1])
        restart = False
        if len(sys.argv) > 2:
            restart = bool( sys.argv[2] )
       
    else:
        dsnames = [str(x) for x in sys.argv[1:]]

    if dsnames is None:
        find_halos(None, simfile = simfile)
    else:

        for dsname in dsnames:
            find_halos(dsname, restart=restart)
#        compute_virial_quantities(dsname)
