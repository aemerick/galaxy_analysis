import yt.mods as yt
import matplotlib.pyplot as plt
import numpy as np

def sfrFromParticles(ds, data, selection = None, times = None):
    """
    Given a dataset, computes the star formation rate as a 
    function of time 
    """

    nstars = np.size( data['particle_mass'] )

    if selection is None:
        selection = [True]*nstars

    particle_mass = data['birth_mass'][selection].value * yt.units.Msun
    creation_time = data['creation_time'][selection].convert_to_units('Myr')
    currentTime   = ds.current_time.convert_to_units('Myr')

    if times is None:
        bin_spacing = 10.0 * yt.units.Myr
        times = np.linspace(np.min(creation_time), currentTime, bin_spacing)
    elif np.size(times) == 1:
        bin_spacing = times
        if not hasattr(bin_spacing, 'value'):
            bin_spacing = bin_spacing * yt.units.Myr

        times = np.linspace(np.min(creation_time), currentTime, bin_spacing)

    times = np.array(times)
    times = np.sort(times)
    sfr   = np.zeros(np.shape(times))

    for i,t in enumerate(times[1:]):
        dt = t - times[i-1]
        dm = np.sum(particle_mass[creation_time <= t]) -\
             np.sum(particle_mass[creation_time <= times[i-1]])

        sfr[i] = dm / dt

    return times, sfr


if __name__=='__main__':

    ds   = yt.load('./DD0100/DD0100')
    data = ds.all_data()

    times = np.arange(0.0*yt.units.Myr, ds.current_time.convert_to_units('Myr'), 10.0*yt.units.Myr)

    times, sfr = sfrFromParticles(ds, data, times = times)
    fig, ax = plt.subplots(figsize=(8,8))
    ax.plot(times.value, sfr, color = 'black', lw = 3)
    ax.set_xlabel('Time (Myr)')
    ax.set_ylabel('SFR (Msun/Myr)')
    ax.set_ylim(0, np.max(sfr)*1.1)
    plt.savefig('sfr.png')



