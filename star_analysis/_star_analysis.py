import yt
import numpy as np
from collections import Iterable
#
#
#
from onezone.star import Star, StarList


__all__ = ['get_list_of_stars', 'get_star_property', 'Star', 'StarList']

def get_list_of_stars(ds, data, dummy_call = False):

    types = {11 : 'star', 12 : 'WD', 13 : 'remnant', 1 : 'N/A'}

    M      = data['birth_mass'].value
    Z      = data['metallicity_fraction'].value
    t_form = data['creation_time'].convert_to_units('Myr').value
    id     = data['particle_index'].value
    PT     = np.abs( data['particle_type'].value )

    try:
        star_list = [Star(star_type = types[ptype], M=m, Z=z, tform = t, id = i)\
                                    for m, z, t, i, pt in zip(M, Z, t_form, id, PT)]
    except:
        if dummy_call:
            star_list = [Star(star_type = 'star', M = 1.0, Z = 0.01, tform = 0.0, id = 0)]
        else:
            # try again and actuall fail this time
            star_list = [Star(star_type = types[ptype], M=m, Z=z, tform = t, id = i)\
                                    for m,z,t,i,ptype in zip(M, Z, t_form, id, PT)]


    M = data['particle_mass'].convert_to_units('Msun').value
    lifetime = data['dynamical_time'].convert_to_units('Myr').value

    for i, s in enumerate(star_list):
        s.M = M[i]
        s.properties['lifetime'] = lifetime[i]

    AllStars = StarList(star_list)

    return AllStars

def get_star_property(ds, data, AllStars = None, property_names = None,
                      dummy_call = False):
    """

    """

    if AllStars is None:
        AllStars = get_list_of_stars(ds, data, dummy_call = dummy_call)

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
