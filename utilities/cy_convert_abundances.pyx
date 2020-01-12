"""

    Author : A. Emerick


    Cython versions of some of the routines in convert_abundances.py that can
    be used to speed things up if doing many conversions. Speed up increases
    with number of conversions, from factors of ~2 for <1E6 points to
    ~ for 1E8 points. Can get faster.

    Currently lacking setup.py to handle cython routines in this repo. Best way
    to deal with this for now is doing: 

    import pyximport; pyximport.install(setup_args={'include_dirs':[np.get_include()]},
                                        language_level=3)
"""


import numpy as np
cimport numpy as np


# --- These could be changed to pure cython variables / arrays to speed things
#     up even more, but I'm a little lazy
from galaxy_analysis.static_data import SOLAR_ABUNDANCE, MOLECULAR_WEIGHT, AMU



cpdef np.ndarray abundance_ratio_array(str x1e,
                           double[:] x1,
                           str x2e,
                           double[:] x2,
                           str input_type = 'abundance', str normalize='solar'):
    """
    Rrturns abundance ratios for elements given element symbol (x1e), array of values (x1),
    and second symbol (x2e) and array of values (x2). input_type is either "abundance"
    or "mass".
    """

    if x1.shape[0] != x2.shape[0]:
        print("Arrays are not of equal size")
        raise ValueError

    cdef double[:] result = np.zeros(x1.shape[0])
    cdef int i = 0

    for i in np.arange(x1.shape[0]):
        result[i] = cy_abundance_ratio( x1e, x1[i], x2e, x2[i],input_type=input_type)

    return np.array(result)

cdef double cy_elemental_abundance(str element, double mass):

    cdef double n = mass / (MOLECULAR_WEIGHT[element] * AMU)

    return n

cdef double cy_abundance_ratio(str x1e, double x1_abund, str x2e, double x2_abund, str input_type = 'abundance'):
    """
    Normalize abundance ratio to solar. x1 and x2 are tuples
    containing either element atomic number and abundance, or
    atomic symbol and abundance (abundance = number of particles).
    E.g.:
       normalize_abundance_ratio( ('Fe',1), ('H', 100) )
    or
       normalize_abundance_ratio( (26, 1), (1, 100) )

    Returns [x1/x2] where
       [x1/x2] = log10(x1/x2) - log10(x1_sun / x2_sun)

    Optionally, x1 and x2 can be mass (in cgs, unless astropy units
    are used) and
    """

    if input_type == 'mass':
        x1_abund = cy_elemental_abundance(x1e, x1_abund)
        x2_abund = cy_elemental_abundance(x2e, x2_abund)

    cdef double x1_solar = SOLAR_ABUNDANCE[x1e]
    cdef double x2_solar = SOLAR_ABUNDANCE[x2e]

    cdef double aratio = np.log10(x1_abund / x2_abund) - (x1_solar - x2_solar) # np.log10( x1_solar / x2_solar)

    return aratio
