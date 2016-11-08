import yt.mods as yt
import numpy as np
from collections import Iterable, OrderedDict


# -- internal --
from galaxy_analysis.utilities import convert_abundances


def generate_aratio(ds, data, ratios, particle_type = 11):
    """
    Generate abundances for all particles of given particle type.
    If particle_type is None, does this for all particles
    """

    birth_mass = data['birth_mass'].value * yt.units.Msun
    ptype      = data['particle_type']

    if not isinstance(ratios, Iterable):
        ratios = [ratios]
    #
    # split string
    #
    aratios = OrderedDict()

    for ratio in ratios:
        if '/' in ratio:
            ele1, ele2 = ratio.rsplit('/')
        else:
            print "Must provide abundance ratio string of style: Fe/H"
            return

        enzo_name1 = ('io','particle_' + ele1 + '_fraction')
        enzo_name2 = ('io','particle_' + ele2 + '_fraction')

        mass1 = data[enzo_name1] * birth_mass
        mass2 = data[enzo_name2] * birth_mass

        aratios[ratio] = convert_abundances.abundance_ratio( (ele1, mass1), 
                                                             (ele2, mass2), 
                                                             'mass' )


        
    return aratios

