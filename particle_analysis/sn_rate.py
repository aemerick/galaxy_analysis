import yt.mods as yt
import matplotlib.pyplot as plt
import numpy as np
import glob

def snr(ds, data, times = None, sn_type = 'II'):
    """
    Computes the supernova rate of the desired time for a given dataset
    as a function of time. The way the particle types and particle lifetimes
    are handled, this can be done for the entire galaxy history using a single
    snapshot, rather than having to sort through each dump.

    One can provide sample times using "times" argument, or leave it alone for
    a 10 Myr sample spacing from t = 0 to t = current_time. If a single value
    is provided, this is taken to be the sample spacing (dt), sampled over
    t = 0 to t = current_time. Units are assumed to be Myr if not provided.

    Accounts for direct collapse model in computing SNII rates using
    parameter file.
    """

    _core_collapse_labels = ["SNII", "II", "2", "SN_II", "TypeII", "Type 2",
                             "Type II", "type II", "typeII", 'core collapse']
    _snia_labels = ["SN1a", "SNIa", "Type1a", "TypeIa", "Type Ia", "Type 1a",
                     "type 1a", "type Ia", "type ia", "type1a", "typeIa"]

    #
    # Determine the particle type we are working with
    #
    if any([sn_type in x for x in _core_collapse_labels]):
        sn_particle_type = 13 # Remnant (core collapse + direct collapse)
    elif any([sn_type in x for x in _snia_labels]):
        sn_particle_type = 12 # WD and exploded WD have same type
    else:
        print "sn_type :" + sn_type + " not a valid option - check spelling"
        return -1

    if times is None:
        bin_spacing = 10.0 * yt.units.Myr
        times = np.linspace(np.min(creation_time), currentTime, bin_spacing)*yt.units.Myr
    elif np.size(times) == 1:
        bin_spacing = times
        if not hasattr(bin_spacing, 'value'):
            bin_spacing = bin_spacing * yt.units.Myr

        times = np.linspace(np.min(creation_time), currentTime, bin_spacing)
        times = times *yt.units.Myr


    # load particle properties
    birth_mass    = data['birth_mass'].value
    mass          = data['particle_mass'].convert_to_units("Msun").value
    creation_time = data['creation_time'].convert_to_units('Myr').value
    metallicity   = data['metallicity_fraction'].value
    lifetimes     = data['dynamical_time'].convert_to_units('Myr').value
    pt            = data['particle_type'].value
    current_time  = ds.current_time.convert_to_units('Myr').value

    # check to see if there are any SN candidates in the first place
    if not any([sn_particle_type == x for x in np.unique(pt)]):
        print "no supernova of type " + sn_type + " found"
        return times, np.zeros(np.size(times.value) - 1)

    # start the particle array slicer
    pcut = (pt == sn_particle_type)

    # looking for core collapse supernova rate
    if sn_particle_type == 13:

        # ignore stars that did not actually go supernova
        collapse_threshold = ds.parameters['IndividualStarDirectCollapseThreshold']
        if not any([x <= collapse_threshold for x in birth_mass[pcut]]):
            print "no core collapse supernova present, only direct collapse"
            return times, np.zeros(np.size(times.value) - 1)

        # slice!
        pcut *= (birth_mass <= collapse_threshold)

    elif sn_particle_type == 12:

        # SNIa are the ones that are just masless tracers, rest are WD
        if not any(mass[pcut] == 0.0):
            print "no Type Ia supernova, only white dwarfs"
            print "N_WD = %i -- Lowest mass = %.3f Msun"%(np.size(mass[pcut]), np.min(mass[pcut]))
            print "Current time = %.2E Myr - Next to explode at t = %.2E Myr"%(current_time, np.min(lifetimes[pcut] + creation_time[pcut]))
            return times, np.zeros(np.size(times.value) - 1)

        # slice!
        pcut *= (mass == 0.0)


    #
    # now get the explosion times for all supernova
    # when stars go SN, lifetime is set to be lifetime*huge_number
    # therefore, explosion time can be backed out as:
    #
    explosion_times = creation_time[pcut] + lifetimes[pcut]/ds.parameters['huge_number']
    explosion_times = explosion_times * yt.units.Myr

    times = times.convert_to_units('yr')
    snr   = np.zeros(np.size(times.value) - 1)

    # compute SNR
    for i in np.arange(np.size(times) - 1):
        dt = times[i+1] - times[i]
        dN = np.size( explosion_times[explosion_times <= times[i+1]]) -\
             np.size( explosion_times[explosion_times <= times[i]])

        snr[i] = dN / dt

    return times, snr


if __name__ == '__main__':
    # example usage - uses most recent data file

    log = False

    ds_list = np.sort( glob.glob('./DD????/DD????'))

    ds   = yt.load(ds_list[-1])
    data = ds.all_data()

    dt = 25.0
    times = np.arange(0.0, ds.current_time.convert_to_units('Myr').value + dt, dt)
    times = times*yt.units.Myr

    times, snrII = snr(ds, data, times = times, sn_type = 'TypeII')
    times, snrIa = snr(ds, data, times = times, sn_type = "TypeIa")

    center = 0.5 * (times[1:] + times[:-1])

    fig, ax = plt.subplots(figsize=(8,8))

    snialabel = 'Type Ia'
    sniilabel = 'Core Collapse'

    if log:
        ax.plot(center/1.0E6, snrII*1.0E6, color = 'black', lw = 3, ls = '-', label = sniilabel)
        ax.plot(center/1.0E6, snrIa*1.0E6, color = 'black', lw = 3, ls = '--', label = snialabel)
        x.semilogy()
    else:
        ax.step(times[:-1]/1.0E6, snrII*1.0E6, color ='black', lw = 3, ls = '-', label = sniilabel)
        ax.step(times[:-1]/1.0E6, snrIa*1.0E6, color ='orange', lw = 3, ls = '-', label = snialabel)

    ax.set_xlabel('Time (Myr)')

    ax.set_ylabel(r'SNR (10$^{-6}$ yr$^{-1}$)')
    ax.set_ylim( np.min( [np.min(snrIa), np.min(snrII)])*1.0E6,
                 np.max( [np.max(snrIa), np.max(snrII)])*1.5*1.0E6)
    ax.legend(loc ='best')
    plt.tight_layout()
    ax.minorticks_on()
    plt.savefig('snr.png')
