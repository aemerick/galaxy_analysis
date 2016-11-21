import yt
from yt import units as u
import numpy as np
import glob
import matplotlib.pyplot as plt


def compute_IMF(ds, data, mode = 'mass', **kwargs):
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


if __name__ == "__main__":
#
# Example usage - plotting IMF's for most recent
# dump in a directory - this example script assumes using
# Salpeter IMF
#
    ds_list = np.sort(glob.glob('./DD????/DD????'))
    ds      = yt.load(ds_list[-1])
    data    = ds.all_data()

    nbins = 25

    hist, bins, cent = compute_IMF(ds, data, mode = 'mass', bins = nbins)

    fig, ax = plt.subplots(1)
    ax.plot(cent, hist, color = 'black', ls = '--', lw = 3, label = 'Mass')

    hist, bins, cent = compute_IMF(ds, data, mode = 'birth_mass', bins = nbins)

    ax.plot(cent, hist, color = 'black', ls ='-', lw = 3, label = 'Birth Mass')

    theory = scaled_IMF(cent, hist, alpha = ds.parameters['IndividualStarSalpeterSlope'])

    ax.plot(cent, theory, color = 'red', ls = ':', lw = 3, label = "Salpeter")

    ax.semilogy()
    ax.set_xlabel(r'log(M [M$_{\odot}$])')
    ax.set_ylabel(r'dN/log(M)')
    fig.set_size_inches(8,8)
    ax.legend(loc='best')
    ax.minorticks_on()
    plt.tight_layout()
    fig.savefig('IMF.png')
    plt.close()
