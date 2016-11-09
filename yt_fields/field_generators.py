import yt

from collections import Iterable

from galaxy_analysis.static_data import AMU,\
                 MOLECULAR_WEIGHT

def _mass_fraction_function_generator(asym):

    if not isinstance(asym, Iterable):
        asym = [asym]

    for a in asym:

        def _mass_fraction(field,data):
            ele_dens = data[('enzo', a + '_Density')].value
            ele_dens = ele_dens * data.ds.mass_unit / data.ds.length_unit**3
            ele_dens = ele_dens.convert_to_cgs()

            dens = data[('enzo','Density')].convert_to_cgs()

            mass_fraction = ele_dens / dens

            return mass_fraction

        yt.add_field(('gas', a + '_Fraction'), function = _mass_fraction, units='dimensionless')

    return

def _number_density_function_generator(asym):

    if not isinstance(asym, Iterable):
        asym = [asym]

    for a in asym:

        def _number_density(field,data):
            ele_dens = data[('enzo', a + '_Density')].value
            ele_dens = ele_dens * data.ds.mass_unit / data.ds.length_unit**3
            ele_dens = ele_dens.convert_to_cgs()

            n = ele_dens / (MOLECULAR_WEIGHT[a] * AMU * yt.units.g)

            return n.convert_to_cgs()

        yt.add_field(('gas', a + '_Number_Density'),
                     function = _number_density, units='cm**(-3)')


def _abundance_ratio_function_generator(ratios):

    if not isinstance(ratios, Iterable):
        ratios = [ratios]

    for r in ratios:

        ele1, ele2 = r.rsplit('/')

        def _abundance_ratio(field, data):

            mass1 = data[('gas', ele1 + '_Density')].value * data.ds.mass_unit / data.ds.length_unit**3
            mass1 = (mass1 * data['cell_volume']).convert_to_cgs()

            mass2 = data[('gas', ele2 + '_Density')].value * data.ds.mass_unit / data.ds.length_unit**3
            mass2 = (mass2 * data['cell_volume']).convert_to_cgs()

