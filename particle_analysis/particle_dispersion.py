import numpy as np
import yt


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

    for x,i in enumerate(fields):

        if filter is not None:
            v = data[x][filter]
        else:
            v = data[x]

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

    for x in dispersions:
        print x, dispersions[x]

#    fig, ax = plt.subplots(figsize=(8,8))
#    ax.plot(times/1.0E6, mass, color = 'black', lw = 3)
#    ax.set_xlabel('Time (Myr)')
#    ax.set_ylabel(r'Cumulative SFH (M$_{\odot}$)')
#    ax.set_ylim(1.0, np.max(mass)*5.0)
#    ax.semilogy()
#    ax.minorticks_on()
#    plt.savefig('sfh.png')

    
   
