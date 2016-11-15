import yt
import matplotlib.pyplot as plt
import numpy as np
import glob 

from galaxy_analysis.static_data import MAX_R

#
# doing spatial distributino stuff for the particles
#
# want to be able to compute particle radial profile for
# an abitrary field
#
def _generate_radial_bins(bins):

    if bins is None:
        dr   = 10.0 * yt.units.pc
        bins = np.arange(0.0, MAX_R.convert_to_units(dr.units)+dr, dr)
        bins = bins * dr.unit_quantity

    elif np.size(bins) == 1:

        if hasattr(bins, 'units'):
            dr = bins
            bins = np.arange(0.0, MAX_R.convert_to_units(dr.units)+dr, dr)
            bins = bins * dr.unit_quantity
        else:
            # take as number of bins and bin in parsecs
            nbins = bins
            bins = np.linspace(0.0, MAX_R.convert_to_units('pc'), nbins)
            bins = bins * yt.units.pc

    return bins

_particle_profile_fields = ['count', 'particle_mass', 'total-dispersion',
                            'r-dispersion','anisotropy','metallicity_fraction']

def radial_profile(ds, data, fields, weight_field = None, bins = None,
                                     mode = 'both'):
    """
    Mode is going to be either sum or average or both
    """

    r = data['particle_position_spherical_radius'].convert_to_units('pc')
    bins = _generate_radial_bins(bins)

    all_profiles = {}

    if np.size(mode) == 1:
        if mode == 'both':
            mode = ['sum','average']
        else:
            mode = [mode]

    for m in mode:

        profiles = {}

        for field in fields:

            p = np.zeros(np.size(bins) - 1)

            if field in ['number','count','N']:
                if m == 'average':
                    # average number in bin doesn't make sense
                    continue
                #
                # profile of particle count as function of radius
                #
                for i in np.arange(0, np.size(bins)-1):
                    p[i] = np.size(r[ (r <    bins[i+1])*
                                          (r >= bins[i])])
            elif 'dispersion' or 'anisotropy' in field:
                if m == 'sum':
                    continue

                if 'dispersion' in field:

                    disp_type = field.rsplit('-')[0]
                    if disp_type == 'total':
                        types = ['x','y','z']
                    else:
                        types = [disp_type]

                    val = compute_velocity_dispersion( data, types = [disp_type],
                                                        filter = (r<bins[i+1])*(r>=bins[i]))

                    if types == 'total':
                       p[i] = np.sqrt(val['x']*val['x'] + val['y']*val['y'] + val['z']*val['z'])
                    else:
                       p[i] = val[disp_type]
                else:
                    val = compute_velocity_dispersion(data, types =['r','theta'],
                                                      filter = (r<bins[i+1])*(r>=bins[i]))

                    p[i] = 1.0 - val['theta']**2 / val['r']**2
            else:
                for i in np.arange(0, np.size(bins)-1):
                    if m == 'sum':
                        p[i] = np.sum( data[field][ (r< bins[i+1])*
                                                    (r>=bins[i])])
                    elif m == 'average':
                        vals = data[field][ (r<bins[i+1])*(r>=bins[i])]

                        if np.size(vals) >= 1:
                            p[i] = np.average(vals)
                        else:
                            p[i] = 0.0



            profiles[field] = p

        all_profiles[m] = profiles

    centers = 0.5 * (bins[1:] + bins[:-1])

    return bins, centers, all_profiles

def cylindrical_profile(ds, data, field, weight_field = None):

    raise NotImplementedError

    return
