import matplotlib as mpl
mpl.use('Agg')

import yt
from yt import units as u
import numpy as np
import glob
import matplotlib.pyplot as plt

# only needed for plotting sampling error
ALLOW_SAMPLE_ERROR = True
try:
    from onezone import imf as onez_imf
except:
    print "Could not load onezone imf module - cannot do imf sampling error"
    ALLOW_SAMPLE_ERROR = False

__all__ = ['compute_IMF', 'IMF', 'scaled_IMF']

def compute_IMF(ds, data, mode = 'mass', tmin = 0.0, tmax = np.inf, **kwargs):
    """
    Wrapper around IMF function to compute IMF
    given Enzo data. Options to compute observed IMF using
    only "live" particles with current masses, or using
    birth mass to get the theoretical total IMF over the
    entire history of the simulation

    mode = 'mass' (default)
        Use current particle mass of live stars only
    mode = 'birth_mass'
        Use birth mass of ALL star particles (even
        ones that have moved on to greener pastures)
    """

    if mode == 'mass':
        M  = data['particle_mass'].convert_to_units('Msun')
        pt = data['particle_type']
        M  = M[pt == 11]
    elif mode == 'birth_mass':
        M = data['birth_mass'].value * u.Msun
    else:
        raise ValueError("Choose either 'mass' or 'birth_mass' for mode")

    t_o = data['creation_time'].to('Myr')

    M = M[ (t_o > tmin) * (t_o < tmax)]

    return IMF(M, m_min = ds.parameters['IndividualStarIMFLowerMassCutoff'],
                  m_max = ds.parameters['IndividualStarIMFUpperMassCutoff'],
                  **kwargs)

def IMF(M, m_min = 1.0, m_max = 100.0, bins = None, log = True):
    """
    Actually tabulates the IMF given a list
    of masses and binning properties. If log is true, assume
    whatever bin information is passes is log bins. If bins is a
    single value, it is taken as the number of bins.
    """

    if log:
        M     = np.log10(M)
        m_min = np.log10(m_min)
        m_max = np.log10(m_max)

    if bins is None:
        bins = np.linspace(m_min, m_max, 25)

    elif np.size(bins) == 1:
        bins = np.linspace(m_min, m_max, bins)

    hist, bins = np.histogram(M, bins = bins)
    cent       = 0.5 * (bins[1:] + bins[:-1])


    return hist, bins, cent


def scaled_IMF(centers, dNdm, alpha = None):
    """
    Returns scaled IMF to match experiemental
    """

    if alpha < 0:
        alpha = np.abs(alpha)

    y = (10**(centers))**(- alpha)
    A = dNdm[0] / y[0]

    return A*y

def determine_sample_error(Mstar, bins, nmodel = 100):
    """
    Runs `nmodel' samplings of the IMF to compute standard deviation
    from the samplings. This may take a while depending on the
    stellar mass sampled.
    """

    n_simulations = nmodel
    nbins         = np.size(bins)

    hist_total      = np.zeros(nbins-1)
    running_average = np.zeros(nbins-1)
    variance        = np.zeros(nbins-1)

    #mass_average = 0.0
    #mass_variance  = 0.0

    for sim_num in np.arange(1,n_simulations+1):

#        all_stars = [None] * n_cell
#        for i in np.arange():
#            all_stars[i] = imf.sample_IMF(imf.salpeter, M = m_cell, mass_mode = 'keep', alpha=1.35)
        imf =onez_imf.salpeter(alpha = 1.35, M_min = 1.0, M_max = 100.0)
        all_stars = imf.sample(M = Mstar)

#        all_stars = np.asarray(all_stars)
#        all_stars = np.hstack(all_stars)


        hist, bins = np.histogram(all_stars, bins=bins)

        hist_total += hist

        delta = hist - running_average

        running_average += delta / (1.0*sim_num)
        variance        += delta*(hist - running_average)

        M = np.sum(all_stars)
#        delta_mass       = M - mass_average
#        mass_average    += delta_mass / (1.0 * sim_num)
#        mass_variance   += delta_mass * ( M - mass_average)

#    mass_variance = mass_variance / ((n_simulations - 1)*1.0)
#    mass_std      = np.sqrt(mass_variance)

    std = np.sqrt(variance / ((n_simulations -1)*1.0))

    return std

def plot_IMF(ds, data = None, nbins=80, compute_std = False, tmin = 0.0, tmax = np.inf):
    """
    Plot the IMF. If data is not provided, assumes all data
    """
    if data is None:
        data = ds.all_data()

    Mstar = np.sum(data['birth_mass'].value)

    hist, bins, cent = compute_IMF(ds, data, mode = 'mass', bins = nbins, tmin=tmin, tmax = tmax)

    fig, ax = plt.subplots(1)
    ax.plot(bins[1:], hist, color = 'red', ls = '-', lw = 3, label = 'Current Masses of Main Sequence Stars',drawstyle='steps-post')

    hist, bins, cent = compute_IMF(ds, data, mode = 'birth_mass', bins = nbins, tmin= tmin, tmax = tmax)

    ax.plot(bins[1:], hist, color = 'blue', ls ='-', lw = 3, label = 'Initial Masses of all Stars',drawstyle='steps-post')

    theory = scaled_IMF(cent, hist, alpha = ds.parameters['IndividualStarSalpeterSlope'])

    ax.plot(cent, theory, color = 'black', ls = '--', lw = 3, label = "Salpeter")

    ax.semilogy()
    ax.set_ylim(ax.get_ylim())

    if compute_std:
        ymin = ax.get_ylim()[0]
        std = determine_sample_error(Mstar, bins = bins)
        low = theory - std
        x   = 0.1 * ymin * np.ones(np.size(std))
        low = np.max([[low],[x]],axis=0)[0]
        ax.fill_between(cent, low, theory + std, facecolor = 'grey', interpolate=True, alpha = 0.5)

    ax.set_xlabel(r'log(M [M$_{\odot}$])')
    ax.set_ylabel(r'dN/log(M)')
    fig.set_size_inches(8,8)
    ax.legend(loc='best')
    ax.minorticks_on()
    plt.tight_layout()
    fig.savefig('IMF.png')
    plt.close()

    return 

if __name__ == "__main__":
#
# Example usage - plotting IMF's for most recent
# dump in a directory - this example script assumes using
# Salpeter IMF
#
    ds_list = np.sort(glob.glob('./DD????/DD????'))
    ds      = yt.load(ds_list[-1])
    data    = ds.all_data()

    nbins = 80

    plot_IMF(ds, nbins = nbins, compute_std = True)
