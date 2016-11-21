import numpy as np
import yt
import glob
import matplotlib.pyplot as plt

def vel_dispersion(v):
    return np.std(v)


def compute_velocity_dispersion(data, types = None, fields = None, filter = None):
    """
    Returns all possible velocity dispersons from all particles found in the
    data set. A particle filter can be passed using "filter" which is a list
    of booleans:
      e.x.
        > filter = data['particle_type'] == 11
        > particle_dispersion.compute_velocity_dispersion(data, filter = filter)
      --- 
        > filter = (data['particle_type'] == 11) * (data['particle_spherical_r'] < 100.0)
        > particle_dispersion.compute_velocity_dispersion(data, filter = filter)
    """

    types_to_fields = {'x': 'particle_velocity_x',
                       'y': 'particle_velocity_y',
                       'z': 'particle_velocity_z',
                       'r': 'particle_velocity_spherical_radius',
                       'theta': 'particle_velocity_spherical_theta',
                       'phi': 'particle_velocity_spherical_phi'}

    if types is None and fields is None:
        fields = types_to_fields.values()
        keys   = types_to_fields.keys()
    elif fields is None:
        fields = [ types_to_fields[x] for x in types ]
        keys   = types
    else:
        keys = fields

    dispersion = {}

    for i,x in enumerate(fields):

        if filter is not None:
            v = data[x][filter]
        else:
            v = data[x]

        if np.size(v) == 0:
            dispersion[keys[i]] = 0.0
        else:
            dispersion[keys[i]] = vel_dispersion( v.convert_to_units('km/s') )


    return dispersion

#
# function to do LOS velocity dispersion calculatoins
#     - take a list of positions
#     - take a list of velocities
#  then
#     - choose random line of sight
#     - project velocities onto LOS
#     - compute velocity disperions
#     - repeat for X random LOS and average together


if __name__ == '__main__':

    ds_list = np.sort( glob.glob('./DD????/DD????'))

    ds   = yt.load(ds_list[-1])
    data = ds.all_data()

    filter = data['particle_type'] == 11

    dispersions = compute_velocity_dispersion(data, filter = filter)

    dr = 50.0
    bins = np.arange(0.0, 750.0 + dr, dr) * yt.units.pc

    r = data['particle_position_spherical_radius'].convert_to_units('pc')

    beta = np.zeros(np.size(bins)-1)
    sigma = np.zeros(np.size(bins)-1)
    for i in np.arange(1, np.size(bins)):
        radial_filter = (r < bins[i]) * (r >= bins[i-1])

        dispersions = compute_velocity_dispersion(data, filter = filter*radial_filter)

        beta[i-1] = 1.0 - dispersions['theta'].value**2 / dispersions['r'].value**2
        sigma[i-1] = np.sqrt(dispersions['x']**2 + dispersions['y']**2 + dispersions['z']**2)

    centers = 0.5 * ( bins[1:] + bins[:-1])

    fig, ax = plt.subplots(figsize=(8,8))
    ax.plot(centers.value, beta, color = 'black', lw = 3)
    ax.set_xlabel('Radius (pc)')
    ax.set_ylabel(r'Anisotropy Parameter')
    ax.minorticks_on()
    plt.savefig('radial_anisotropy.png')
    plt.close()

    fig, ax = plt.subplots(figsize=(8,8))
    ax.plot(centers.value, sigma, color = 'black', lw =3)
    ax.set_xlabel('Radius (pc)')
    ax.set_ylabel(r'3D velocity dispersion (km/s)')
    ax.minorticks_on()
    plt.savefig('velocity_dispersion.png')
