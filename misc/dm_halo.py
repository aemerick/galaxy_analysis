import numpy as np
from yt import physical_constants as const
import yt.units as u

#
# functions that define the background potential
# given parameters known to the data set already
#
#

def burkert_density(r, r_s, rho_o):
    """
    Burkert dark matter density profile
    """

    x = r / r_s

    density = rho_o / ( (x) * (1.0 + x)**2)

    return density.convert_to_cgs()

def burkert_potential(r, r_s, rho_o):
    """
    Burkert dark matter density potential
    """

    x = r / r_s

    G = const.G.convert_to_cgs()

    phi_o = 4.0 * np.pi * G * r_s**2 * rho_o

    phi   = -0.5 * phi_o * ( (1.0 / x) * (0.5 * np.log(1.0+x*x) + np.log(1.0+x) - np.arctan(x)) +\
                             (np.log(1.0+x) - 0.5 * np.log(1.0+x*x) - np.arctan(x) + np.pi/2.0))


    return phi.convert_to_cgs()


