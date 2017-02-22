#import yt.mods as yt
import yt
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import glob

__all__ = ['sfhFromParticles']

def sfhFromParticles(ds, data, selection = None, times = None):
    """
    Given a dataset, computes the star formation history
    """

    nstars = np.size( data['particle_mass'] )

    if selection is None:
        selection = np.array([True]*nstars)

    particle_mass = data['birth_mass'][selection].value * yt.units.Msun
    creation_time = data['creation_time'][selection].convert_to_units('Myr')
    currentTime   = ds.current_time.convert_to_units('Myr')

    if times is None:
        bin_spacing = 2.0 * yt.units.Myr
        times = np.arange(np.min(creation_time) - bin_spacing*2.0, currentTime, bin_spacing)*yt.units.Myr
    elif np.size(times) == 1:
        bin_spacing = times
        if not hasattr(bin_spacing, 'value'):
            bin_spacing = bin_spacing * yt.units.Myr

        times = np.linspace(np.min(creation_time), currentTime, bin_spacing)
        times = times *yt.units.Myr

    mass  = np.zeros(np.shape(times))

    times = times.convert_to_units('yr')
    for i,t in enumerate(times):
        mass[i] = np.sum( particle_mass[creation_time <= t] )

    return times, mass

if __name__=='__main__':

    ds_list = np.sort( glob.glob('./DD????/DD????'))

    
    ds   = yt.load(ds_list[-1])
    data = ds.all_data()

    times = np.arange(0.0*yt.units.Myr, ds.current_time.convert_to_units('Myr'), 5.0*yt.units.Myr)
    times = times*yt.units.Myr

    times, mass = sfhFromParticles(ds, data, times = times)

    fig, ax = plt.subplots(figsize=(8,8))
    ax.plot(times/1.0E6, mass, color = 'black', lw = 3)
    ax.set_xlabel('Time (Myr)')
    ax.set_ylabel(r'Cumulative SFH (M$_{\odot}$)')
    ax.set_ylim(1.0, np.max(mass)*5.0)
    ax.semilogy()
    ax.minorticks_on()
    plt.savefig('sfh.png')
