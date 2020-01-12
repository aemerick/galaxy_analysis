import yt
import numpy as np
from collections import Iterable
from galaxy_analysis.utilities import utilities as util
#
#
#
from onezone.star import Star, StarList


__all__ = ['get_list_of_stars', 'get_star_property', 'Star', 'StarList']

def get_list_of_stars(ds, data, dummy_call = False, overload_type = None):

    types = {11 : 'star', 12 : 'WD', 13 : 'remnant', 1 : 'N/A'}

    M      = data['birth_mass'].value
    Z      = data['metallicity_fraction'].value
    t_form = data['creation_time'].convert_to_units('Myr').value
    id     = data['particle_index'].value

    if overload_type is None:
        PT     = np.abs( data['particle_type'].value )
    elif overload_type in [11,12,13]:
        PT = np.ones(np.size(M)) * overload_type
    else: # if True or rando value, assume MS
        PT = np.ones(np.size(M)) * 11

    yield_names = util.species_from_fields(ds.field_list, include_primordial=True)
    yield_names = ['m_tot','m_metal'] + yield_names

    # just set up a mostly empty dictionary with species names
    # so that the metal species properties get computed - don't need to actually
    # set up abundances correctly.
    yd = {}
    for k in yield_names:
        yd[k] = 0.0
    yd['m_tot'] = 1.0

    try:
        star_list = [Star(star_type = types[ptype], M=m, Z=z, tform = t, id = i, abundances=yd)\
                                    for m, z, t, i, ptype in zip(M, Z, t_form, id, PT)]
    except:
        if dummy_call:
            star_list = [Star(star_type = 'star', M = 1.0, Z = 0.01, tform = 0.0, id = 0)]
        else:
            # try again and actuall fail this time
            star_list = [Star(star_type = types[ptype], M=m, Z=z, tform = t, id = i, abundances=yd)\
                                    for m, z, t, i, ptype in zip(M, Z, t_form, id, PT)]


    M = data['particle_mass'].convert_to_units('Msun').value
    lifetime = data['dynamical_time'].convert_to_units('Myr').value

    for i, s in enumerate(star_list):

        if overload_type is None: # use actual lifetimes and masses from the simulation
            s.properties['lifetime'] = lifetime[i]
            s.M = M[i]

        s.set_SNII_properties()
        if ((s.M_o > ds.parameters['IndividualStarSNIaMinimumMass']) and\
            (s.M_o < ds.parameters['IndividualStarSNIaMaximumMass']) and\
            (s.M == 0.0)):
            s.set_SNIa_properties()

    AllStars = StarList(star_list)

    return AllStars

def get_model_yields(ds, data, AllStars = None, sum_only = True, overload_type = None):
    """
    Using a provided list of model stars ("AllStars"), or using a generated
    list of model stars, computes the model yields for each species and
    returns the cumulative amount of yields expected for each species.
    """

    if AllStars is None:
        AllStars = get_list_of_stars(ds, data, overload_type = overload_type)

    # now extract all of the yields into a dictionary of arrays
    total_wind_ejecta = {}
    model_sn_ejecta   = {}
    for k in AllStars[0].wind_ejecta_masses().keys():
        total_wind_ejecta[k] = np.array([x.wind_ejecta_masses()[k] for x in AllStars.stars_iterable])
        model_sn_ejecta[k]   = np.array([x.sn_ejecta_masses()[k] for x in AllStars.stars_iterable])

    bm = AllStars.property_asarray('M_o')
    lt = AllStars.property_asarray('lifetime')
    bt = AllStars.property_asarray('tform')
    pt = AllStars.property_asarray('star_type')

    current_time = ds.current_time.convert_to_units('Myr').value
    age          = current_time - bt

    # we need to adjust the wind yields for their time dependense
    AGB = (bm < ds.parameters['IndividualStarSNIIThreshold']) *\
           (pt == 11)
    other = (bm > ds.parameters['IndividualStarSNIIThreshold']) *\
           (pt == 11)

    factor = age / lt
    factor[factor > 1.0] = 1.0

    model_wind_ejecta = dict(total_wind_ejecta)
    for k in total_wind_ejecta.keys():
        model_wind_ejecta[k][AGB] = 0.0 # stars that will go AGB, but have not yet
        model_sn_ejecta[k][ pt == 11] = 0.0 # stars that have not died yet

        model_wind_ejecta[k][select] = model_wind_ejecta[k][select]*factor[select]

    all_ejecta = {}
    for k in model_wind_ejecta.keys():
        all_ejecta[k] = model_sn_ejecta[k] + model_wind_ejecta[k]

    # just return the total in each species and type, not arrays
    # of all stars
    if sum_only:
        for k in model_wind_ejecta.keys():
            all_ejecta[k] = np.sum(all_ejecta[k])
            model_wind_ejecta[k] = np.sum(model_wind_ejecta[k])
            model_sn_ejecta[k]   = np.sum(model_sn_ejecta[k])

    return all_ejecta, model_wind_ejecta, model_sn_ejecta

def get_star_property(ds, data, AllStars = None, property_names = None,
                      dummy_call = False, overload_type = None):
    """

    """

    if AllStars is None:
        AllStars = get_list_of_stars(ds, data, dummy_call = dummy_call, overload_type = overload_type)

    single_field = False

    if property_names is None:

        property_names = ['luminosity', 'L_FUV', 'L_LW', 'Q0', 'Q1', 'E0', 'E1',
                          'Teff', 'R', 'agb_phase_length', 'mechanical_luminosity']

    elif not isinstance(property_names, Iterable):
        single_field = True
        property_names = [property_names]

    elif np.size(property_names) < 2:
        single_field = True

    properties = {}
    for field in property_names:
        properties[field] = AllStars.property_asarray(field)

    if single_field:
        return properties[field]
    else:
        return properties

def generate_CMD():

    print("Placeholder / reminder to do this")

    return
