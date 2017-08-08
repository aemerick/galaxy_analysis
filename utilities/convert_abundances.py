import numpy as np
#import yt.mods as yt
import yt

from galaxy_analysis.static_data import \
     MOLECULAR_WEIGHT, \
     AMU, \
     SOLAR_ABUNDANCE

def elemental_abundance(element, mass):
    n = mass / (MOLECULAR_WEIGHT[element] * AMU)

    return n


def abundance_ratio(x1, x2, input_type = 'abundance'):
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

    # convert to abundance first
    x1_abund = x1[1]
    x2_abund = x2[1]

    if input_type == 'mass' or hasattr(x1[1],'unit'):

        # if we have units, use them to convert to grams
        # and strip units
        # otherwise assume grams
        if hasattr(x1[1], 'unit') and hasattr(x2[1], 'unit'):
            x1_abund = x1[1].to(u.g).value
            x2_abund = x2[1].to(u.g).value


        x1_abund = elemental_abundance(x1[0], x1[1])
        x2_abund = elemental_abundance(x2[0], x2[1])

    #
    # solar dictionary gives log(e_x), where
    # log(e_x) = log(N_x / N_H) + 12.0
    #
    # so log(N_x / N_y) = log(e_x) - log(e_y)
    # 
    #
    #
    x1_solar = SOLAR_ABUNDANCE[x1[0]]
    x2_solar = SOLAR_ABUNDANCE[x2[0]]

    aratio = np.log10(x1_abund / x2_abund) - (x1_solar - x2_solar) # np.log10( x1_solar / x2_solar)

    return aratio

