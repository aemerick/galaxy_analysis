import numpy as np
#import yt.mods as yt
import yt

from galaxy_analysis.static_data import \
     MOLECULAR_WEIGHT, \
     AMU, \
     SOLAR_ABUNDANCE

from onezone import data_tables as DT

SN_YIELD_TABLE           = DT.StellarYieldsTable('SNII')
WIND_YIELD_TABLE         = DT.StellarYieldsTable('wind')
MASSIVE_STAR_YIELD_TABLE = DT.StellarYieldsTable('massive_star')

def return_yields(yield_type, M, Z, species):
    """
    Get the mass yields for a list of species for a star of mass M
    and metallicity Z as given in the adopted yield tables. This can
    be done for three yield_types: 
        'SNII' : Just yield from supernova associated with star (if any)
        'wind' : Yield from winds from star (AGB or otherwise)
        'all'  : sum of 'SNII' and 'wind'

    Stars above M > 25.0 are assumed to have NO CCSN yields
    """


    interp = {'SNII' : lambda x,y : SN_YIELD_TABLE.interpolate([x,y], species),
              'wind' : lambda x,y : WIND_YIELD_TABLE.interpolate([x,y], species),
              'massive_star' : lambda x,y : MASSIVE_STAR_YIELD_TABLE.interpolate([x,y], species)}

    interp['all'] = lambda x,y : interp['SNII'](x,y) + interp['wind'](x,y)

    if M > 25.0:
        yield_type = 'massive_star'

    return np.array(interp[yield_type](M,Z))

def get_yield_ratio(ele1, ele2, M, Z, yield_type = 'all'):
    """
    Get the abundance ratio between two elements for a star of mass M
    and metallicity Z as given in the adopted yield tables. This can
    be done for three yield_types:
        'SNII' : Just yield from supernova associated with star (if any)
        'wind' : Yield from winds from star (AGB or otherwise)
        'all'  : sum of 'SNII' and 'wind'

    Stars above M > 25.0 are assumed to have NO CCSN yields
    """

    yields = return_yields(yield_type, M, Z, [ele1, ele2])

    yields = (yields * yt.units.Msun).to('g').value

    return abundance_ratio((ele1,yields[0]), (ele2,yields[1]), input_type = 'mass')

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

def renormalize(aratio, e1, e2, to_solar = False):
    """
    Takes an abundance ratio (aratio) between elements e1 and e2, and reconverts
    the ratio to just straight abundance, instead of normalized to solar.
    """

    x1_solar = SOLAR_ABUNDANCE[e1]
    x2_solar = SOLAR_ABUNDANCE[e2]
    
    conversion = x1_solar - x2_solar
    
    if to_solar:
        conversion *= -1
        
    return aratio + conversion
