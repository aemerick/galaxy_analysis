import numpy as np
from scipy import optimize as opt
from yt import physical_constants as const
import yt.units as u

#
# functions that define the background potential
# given parameters known to the data set already
#
#

def burkert_solve_virial(r_s, rho_o, rho_crit = 9.74E-30 * u.g/u.cm**3):
    """
    Returns M200 in Msun, R200 in same units as r_s.
    """


    f_M = lambda x : 1.5* (0.5*np.log(1.0 + x**2) +\
                     np.log(1.0 + x)  -\
                     np.arctan(x) )

    # Solve the profile for radius with average density equal to
    # 200 * rho_crit (i.e. solve for R200)
    if hasattr(rho_o, 'value'):
        rho_o    = (rho_o.to('g/cm**3')).value

    if hasattr(rho_crit, 'value'):
        rho_crit = (rho_crit.to('g/cm**3')).value

    eq_solve = lambda x : (rho_o)*f_M(x)/(x**3)-200.0*rho_crit

    R200 = r_s * opt.bisect(eq_solve, .1, 10000.0, xtol=1.0E-12)
    rho_crit = rho_crit * u.g / u.cm**3
    M200 = 4.0 * np.pi * R200**3 * (200.0 * rho_crit) / 3.0


    return M200.to('Msun') , R200

def burkert_density(r, r_s, rho_o):
    """
    Burkert dark matter density profile
    """

    x = r / r_s

    density = rho_o / ( (x) * (1.0 + x)**2)

    return density.to('g/cm**3')

def burkert_potential(r, r_s, rho_o):
    """
    Burkert dark matter density potential
    """

    x = r / r_s

    G = 1.0 * const.G
    G.convert_to_cgs()

    phi_o = 4.0 * np.pi * G * r_s**2 * rho_o

    phi   = -0.5 * phi_o * ( (1.0 / x) * (0.5 * np.log(1.0+x*x) + np.log(1.0+x) - np.arctan(x)) +\
                             (np.log(1.0+x) - 0.5 * np.log(1.0+x*x) - np.arctan(x) + np.pi/2.0))


    phi.convert_to_cgs()
    return phi
