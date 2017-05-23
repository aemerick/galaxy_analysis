#import yt.mods as yt
import yt
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import glob

__all__ = ['sfrFromParticles']

def sfrFromParticles(ds, data, selection = None, times = None, t_o = None):
    """
    Given a dataset, computes the star formation rate as a
    function of time
    """

    nstars = np.size( data['particle_mass'] )

    if selection is None:
        selection = np.array( [True]*nstars )

    particle_mass = data['birth_mass'][selection].value * yt.units.Msun
    creation_time = data['creation_time'][selection].convert_to_units('Myr')
    currentTime   = ds.current_time.convert_to_units('Myr')

    # set start time of array
    if t_o is None:
        t_o = np.min(creation_time)
    else:
        if not hasattr(t_o, 'value'):
            t_o = t_o * yt.units.Myr
        else:
            t_o = t_o.convert_to_units('Myr')

    if times is None:
        bin_spacing = 2.0 * yt.units.Myr
        times = np.arange(t_o - bin_spacing*2.0, currentTime, bin_spacing)*yt.units.Myr
    elif np.size(times) == 1:
        bin_spacing = times
        if not hasattr(bin_spacing, 'value'):
            bin_spacing = bin_spacing * yt.units.Myr
        else:
            bin_spacing = bin_spacing.convert_to_units('Myr')

        t_f = currentTime + bin_spacing

        times = np.arange(t_o, t_f, bin_spacing)
        if not hasattr(times, 'value'):
            times = times * yt.units.Myr

    sfr   = np.zeros(np.size(times)-1)

    times = times.convert_to_units('yr')

    for i in np.arange(np.size(times)-1):
        dt = times[i+1] - times[i]
        dm = np.sum(particle_mass[creation_time <= times[i+1]]) -\
             np.sum(particle_mass[creation_time <= times[i]])

        sfr[i] = dm / dt

    return times, sfr


if __name__=='__main__':
    log = True

    ds_list = np.sort( glob.glob('./DD????/DD????'))

    ds   = yt.load(ds_list[-1])
    data = ds.all_data()

    dt = 10.0*yt.units.Myr

    times = np.arange(0.0*yt.units.Myr, ds.current_time.convert_to_units('Myr')+dt, dt)
    times = times*yt.units.Myr

    times, sfr = sfrFromParticles(ds, data, times = times)
    fig, ax = plt.subplots(figsize=(8,8))

    center = 0.5 * (times[1:]+times[:-1])

    if log:
        ax.plot(center/1.0E6, sfr, color = 'black', lw = 3)
        ax.semilogy()
    else:
        ax.step(times[:-1]/1.0E6, sfr*1.0E4, color = 'black', lw = 3, where='pre')

    ax.set_xlabel('Time (Myr)')

    ax.set_ylabel(r'SFR (10$^{-4}$ M$_{\odot}$/yr)')

    if log:
        ax.set_ylim(1.0E-5, 1.0E-2)
    else:

        ax.set_ylim(0.0, np.max(sfr)*1.1*1.0E4)
    ax.minorticks_on()
    plt.savefig('sfr.png')



