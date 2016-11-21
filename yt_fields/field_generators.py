
__author__ = "Andrew Emerick"
__email__  = "aemerick11@gmail.com"

import yt
from yt.units import dimensions

from collections import Iterable

from galaxy_analysis.static_data import AMU,\
                 MOLECULAR_WEIGHT

from galaxy_analysis.utilities import convert_abundances
from galaxy_analysis.utilities import utilities as util
from galaxy_analysis import star_analysis


def _mass_function_generator(asym):

    if not isinstance(asym, Iterable):
        asym = [asym]

    def return_function(a):
        def _mass(field,data):
            ele_dens = data[('enzo', a + '_Density')].value
            ele_dens = ele_dens * data.ds.mass_unit / data.ds.length_unit**3
            ele_dens = ele_dens.convert_to_cgs()

            mass = (ele_dens * data['cell_volume']).convert_to_units('g')

            return mass

        return _mass

    nfields = 0
    for a in asym:

        yt.add_field(('gas', a + '_Mass'), function = return_function(a), units='g')

        nfields = nfields + 1

    return nfields

#
# Construct arbitrary mass fraction derived fields in yt
# using a loop to generate functions
#
def _mass_fraction_function_generator(asym):

    if not isinstance(asym, Iterable):
        asym = [asym]


    def return_function(a):
        def _mass_fraction(field,data):
            ele_dens = data[('enzo', a + '_Density')].value
            ele_dens = ele_dens * data.ds.mass_unit / data.ds.length_unit**3
            ele_dens = ele_dens.convert_to_cgs()

            dens = data[('enzo','Density')].convert_to_cgs()
            mass_fraction = ele_dens / dens
            return mass_fraction

        return _mass_fraction

    nfields = 0
    for a in asym:

        yt.add_field(('gas', a + '_Fraction'), function = return_function(a), units='auto',
                                               dimensions=dimensions.dimensionless)

        nfields = nfields + 1

    return nfields

#
# Construct arbitrary number density derived fields in yt
# using a loop to generate functions
#
def _number_density_function_generator(asym):

    if not isinstance(asym, Iterable):
        asym = [asym]


    def return_function(a):
        def _number_density(field,data):
            ele_dens = data[('enzo', a + '_Density')].value
            ele_dens = ele_dens * data.ds.mass_unit / data.ds.length_unit**3
            ele_dens = ele_dens.convert_to_cgs()
            n = ele_dens / (MOLECULAR_WEIGHT[a] * AMU * yt.units.g)

            return n.convert_to_cgs()

        return _number_density

    nfields = 0
    for a in asym:

        yt.add_field(('gas', a + '_Number_Density'),
                     function = return_function(a), units='cm**(-3)')
        nfields = nfields + 1

    return nfields

def _particle_abundance_function_generator(ratios):

    if not isinstance(ratios, Iterable):
        ratios = [ratios]

    def return_function(ele1, ele2, field1, field2):
        def _abundance_ratio(field, data):
            mass1 = data[field1].value
            mass1 = ((mass1 * data['birth_mass'].value) * yt.units.Msun).convert_to_units('g')

            mass2 = data[field2].value
            mass2 = ((mass2 * data['birth_mass'].value) * yt.units.Msun).convert_to_units('g')
            ratio = convert_abundances.abundance_ratio( (ele1, mass1.value), (ele2, mass2.value), 'mass')

            return ratio * yt.units.g / yt.units.g

        return _abundance_ratio

    nfields = 0
    for r in ratios:

        ele1, ele2 = r.rsplit('/')

        field1 = ('io','particle_' + ele1 + '_fraction')
        field2 = ('io','particle_' + ele2 + '_fraction')

        fieldname = 'particle_' + ele1 + '_over_' + ele2

        yt.add_field(('io', fieldname), function = return_function(ele1,ele2,field1,field2),
                              units = 'dimensionless', particle_type = True)
        nfields = nfields + 1

    return nfields

#
# Construct arbitrary abundance ratio fields in yt
# using a function generator
#
def _abundance_ratio_function_generator(ratios, H_mode = 'total'):

    if not isinstance(ratios, Iterable):
        ratios = [ratios]

    def _H_mass(data, mode):

        if mode == 'total':
            mass = data[('enzo','HI_Density')] + data[('enzo','HII_Density')]
            
            if ('enzo','H2I_Density') in data.ds.field_list:
                mass += data[('enzo','HM_Density')] + data[('enzo','H2I_Density')] +\
                        data[('enzo','H2II_Density')]

        elif mode == 'HI':
            mass = data[('enzo','HI_Density')]
        elif mode == 'HII':
            mass = data[('enzo','HII_Denisty')]

        return mass

    def return_function(ele1, ele2, field1, field2):
        def _abundance_ratio(field, data):

            if ele1 == 'H':
                mass1 = _H_mass(data, H_mode)
            else:
                mass1 = data[field1]

            if ele1 != 'H' and ele1 != 'He':
                mass1 = mass1.value * data.ds.mass_unit / data.ds.length_unit**3
            mass1 = (mass1 * data['cell_volume']).convert_to_units('g')

            if ele2 == 'H':
                mass2 = _H_mass(data, H_mode)
            else:
                mass2 = data[field2]

            if ele2 != 'H' and  ele2 != 'He':
                mass2 = mass2.value * data.ds.mass_unit / data.ds.length_unit**3
            mass2 = (mass2 * data['cell_volume']).convert_to_units('g')

            ratio = convert_abundances.abundance_ratio( (ele1, mass1.value), (ele2, mass2.value), 'mass')

            return ratio * yt.units.g / yt.units.g

        return _abundance_ratio

    nfields = 0
    for r in ratios:

        ele1, ele2 = r.rsplit('/')

        if ele1 != 'H' and ele2 != 'He':
            field1 = ('enzo', ele1 + '_Density')
        else:
            field1 = ('gas', ele1 + '_density')

        if ele2 != 'H' and ele2 != 'He':
            field2 = ('enzo', ele2 + '_Density')
        else:
            field2 = ('gas', ele2 + '_density')

        fieldname = ele1 + '_over_' + ele2

        yt.add_field(('gas', fieldname), function = return_function(ele1,ele2,field1,field2),
                              units = 'dimensionless')
        nfields = nfields + 1

    return nfields


def generate_stellar_model_fields(ds):

    #
    # luminosity, L_FUV, L_LW, Q0, Q1, E0, E1
    #

    field_names = ['luminosity','L_FUV','L_LW','Q0','Q1','E0','E1']

    units = {'luminosity' : yt.units.erg/yt.units.s,
             'L_FUV' : yt.units.erg/yt.units.s,
             'L_LW'  : yt.units.erg/yt.units.s,
             'Q0' : 1.0 /yt.units.s, 'Q1' : 1.0 / yt.units.s,
             'E0' : yt.units.eV, 'E1' : 1.0 / yt.units.eV}

    unit_label = {'luminosity': 'erg/s', 'L_FUV' : 'erg/s', 'L_LW' : 'erg/s',
                  'Q0' : '1/s', 'Q1' : '1/s', 'E0' : 'erg', 'E1': 'erg'}

    def _function_generator(field_name):
        def _function(field, data):

            p = star_analysis.get_star_property(ds, data, property_names = [field_name])
            p = p * units[field_name]
            return p

        return _function

    for field in field_names:
        yt.add_field(('io', 'particle_model_' + field),
                     function = _function_generator(field), units=unit_label[field],
                     particle_type = True)

    return

def _additional_helper_fields():

    nfields = 0


    def _H_total_mass(field, data):
        mass = data[('gas','H_p0_mass')] + data[('gas','H_p1_mass')]

        return mass

    def _He_total_mass(field, data):
        mass = data[('gas','He_p0_mass')] + data[('gas','He_p1_mass')] +\
               data[('gas','He_p2_mass')]

        return mass

    def _metal_total_mass(field, data):
        mass = data['Metal_Density'] * data['cell_volume']

        return mass.convert_to_units('g')

#    def _H2_total_mass(field, data):
#        mass = data[('gas',

    yt.add_field(('gas','H_total_mass'), function = _H_total_mass, units ='g')
    yt.add_field(('gas','He_total_mass'), function = _He_total_mass, units = 'g')
    yt.add_field(('gas','metal_mass'), function = _metal_total_mass, units = 'g')
#    yt.add_field(('gas','H2_total_mass'), function = _H2_total_mass, units = 'g')
#    yt.add_field(('gas','All_H_total_mass'), function = _all_H_total_mass, units = 'g')

    nfields = 3

    return nfields

def generate_derived_fields(ds):
    """
    Given a data set (to extract the on-disk field names), generate
    derived fields that will persist for the python session (i.e. not
    tied only to the passed data set).

    Right now, takes in all metal species tracer fields and constructs
    fields for their mass fraction, number density, and all possible
    interesting abundance ratios.

    NOTE: The derived fields will only exist for data sets loaded after
    this function call. If analysis is intended for passed data set,
    it will need to be reloaded for fields to exist.
    """

    fields = ds.field_list

    # lets figure out the metal tracers present
    metals = util.species_from_fields(fields)
    ratios = util.ratios_list(metals)

    print "tracer species present: ", metals
    nfields = _mass_function_generator(metals)
    print nfields, "mass fields defined"
    nfields = _mass_fraction_function_generator(metals)
    print nfields, "mass fraction fields defined"
    nfields = _number_density_function_generator(metals)
    print nfields, "number density fields defined"
    nfields = _abundance_ratio_function_generator(ratios, H_mode = 'total')
    print nfields, "abundance ratio fields defined"

    if ds.parameters['NumberOfParticles'] > 0:
        nfields =  _particle_abundance_function_generator(ratios)
        print nfields, "particle abundance ratio fields defined"

    print "no particle fields found"

    generate_stellar_model_fields(ds)

    nfields = _additional_helper_fields()
    print nfields, "additional helper fields defined"


    return
