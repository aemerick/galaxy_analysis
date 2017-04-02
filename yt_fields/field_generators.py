
__author__ = "Andrew Emerick"
__email__  = "aemerick11@gmail.com"

import yt
from yt.units import dimensions
import numpy as np

from collections import Iterable

from galaxy_analysis.static_data import AMU,\
                 MOLECULAR_WEIGHT

from galaxy_analysis.utilities import convert_abundances
from galaxy_analysis.utilities import utilities
from galaxy_analysis import star_analysis
from galaxy_analysis.misc import dm_halo

FIELDS_DEFINED = False

def _density_function_generator(asym):
    if not isinstance(asym, Iterable):
        asym = [asym]

    def return_function(a):
        def _density(field, data):
            dens = data[('enzo', a + '_Density')].value
            dens = dens * data.ds.mass_unit / data.ds.length_unit**3
            return dens.convert_to_units('g/cm**3')

        return _density

    for a in asym:
        yt.add_field(('gas', a + "_Density"), function = return_function(a), units = 'g/cm**3')

    return

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

        yt.add_field(('gas', a + '_Fraction'), function = return_function(a), units="")

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

    # make a metal number density field - make an assumption about the metal molecular weight
    def _metal_number_density(field,data):
        ele_dens = data[('enzo','Metal_Density')].convert_to_units('g/cm**3')
        n = ele_dens / (MOLECULAR_WEIGHT['metal'] * AMU * yt.units.g)

        return n.convert_to_cgs()

    yt.add_field(('gas', 'Metal_Number_Density'),
                 function = _metal_number_density, units = 'cm**(-3)')

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
                              units = "", particle_type = True)
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
                              units = "")
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
             'E0' : yt.units.eV, 'E1' : yt.units.eV}

    unit_label = {'luminosity': 'erg/s', 'L_FUV' : 'erg/s', 'L_LW' : 'erg/s',
                  'Q0' : '1/s', 'Q1' : '1/s', 'E0' : 'erg', 'E1': 'erg'}

    def _function_generator(field_name):
        def _function(field, data):
            if np.size(data['particle_mass']) == 1:
                # this is ugly, but a way to bypass yt's validation step
                # and the fact that the below will print errors during this
                with utilities.nostdout():
                    p = star_analysis.get_star_property(ds, data, property_names = [field_name],
                                                            dummy_call = True)
            else:
                p = star_analysis.get_star_property(ds, data, property_names = [field_name])
            p = p * units[field_name]

            return p

        return _function
 
    def _model_L0(field, data):
        Q0 = data[('io','particle_model_Q0')]
        E0 = data[('io','particle_model_E0')]

        return (E0 * Q0).convert_to_units('erg/s')

    def _model_L1(field, data):
        Q1 = data[('io','particle_model_Q1')]
        E1 = data[('io','particle_model_E1')]

        return (E1 * Q1).convert_to_units('erg/s')
   

    def _age(field, data):
        p = data[('io','creation_time')]
        t = data.ds.current_time

        return (t - p).convert_to_units('Myr')

    for field in field_names:
        yt.add_field(('io', 'particle_model_' + field),
                     function = _function_generator(field), units=unit_label[field],
                	     particle_type = True)

    yt.add_field(('io', 'particle_age'), function = _age, units = 'Myr',
                 particle_type = True)

    yt.add_field(('io','particle_model_L0'), function = _model_L0, units = 'erg/s',
                 particle_type = True)

    yt.add_field(('io','particle_model_L1'), function = _model_L1, units = 'erg/s',
                 particle_type = True)

    return

def _additional_helper_fields(fields):

    nfields = 0


    def _H_total_mass(field, data):
        mass = data[('gas','H_p0_mass')] + data[('gas','H_p1_mass')]

        return mass

    def _He_total_mass(field, data):
        mass = data[('gas','He_p0_mass')] + data[('gas','He_p1_mass')] +\
               data[('gas','He_p2_mass')]

        return mass

    def _pe_heating_cgs(field,data):
        pe = data[('enzo','Pe_heating_rate')].value

        energy_unit = data.ds.mass_unit * data.ds.velocity_unit**2
        pe = pe * energy_unit / data.ds.length_unit**3 / data.ds.time_unit

        return pe.convert_to_units('erg/s/cm**3')

    def _pe_heating_rate(field, data):
        pe = data[('gas','Pe_heating_rate')].convert_to_units('erg/s/cm**3')

        x = 1.0 * pe

        x[data['temperature'] > data.ds.parameters['IndividualStarFUVTemperatureCutoff']] = 0.0

        return x

    def _otlwcgs(field, data):
        lw = data[('enzo','OTLW_kdissH2I')].value / data.ds.time_unit
        return lw.convert_to_units('1/s') 

    def _G_o(field,data):

        pe  = data[('gas','Pe_heating_rate')].convert_to_units('erg/s/cm**3').value
        Z   = (data['Metal_Density'] / data['Density']).value
        n_H = (data['H_p0_number_density'] + data['H2_number_density'] + data['H_p1_number_density']).convert_to_units('cm**(-3)').value

        g_to_d = 0.68 - 3.08 * np.log10(Z / 0.014)
        d_to_g = 1.0 / (10.0**(g_to_d))

        D = d_to_g / 6.616595E-3

        epsilon = 0.01488637246 * (n_H)**(0.235269059)

        atten = np.exp( - 1.33E-21 * D * data['dx'].convert_to_units('cm').value * n_H)

        G_o = pe / (1.3E-24 * n_H * epsilon * D * atten)

        return G_o * (data['Density'] / data['Density'])


    def _metal_total_mass(field, data):
        mass = data['Metal_Density'] * data['cell_volume']

        return mass.convert_to_units('g')

    def _grav_pot(field,data):
        try:
            x = (data['PotentialField'] * -1.0).convert_to_units('erg/g')
        except:
            x = ( (data['GravPotential'].value * data.ds.velocity_unit**2)
                   * -1.0).convert_to_units('erg/g')

        return x

    def _tot_grav_pot(field,data):
        try:
            x = data['PotentialField']
        except:
            x = data['GravPotential'].value * data.ds.velocity_unit**2

        x = x + data[('index','DM_background_potential')]

        return x.convert_to_units('erg/g')

    def _pos_tot_grav_pot(field, data):

        return np.abs(data[('gas','total_gravitational_potential')])

    def _potential_energy(field,data):

        x = data[('gas','total_gravitational_potential')] * data['cell_mass']

        return x.convert_to_units('erg')

    def _grav_bound(field, data):
        PE = data[('gas','potential_energy')].convert_to_units('erg')
        TE = ( data[('gas','thermal_energy')] * data['cell_mass'].convert_to_units('g')).convert_to_units('erg')
        KE = ( data[('gas','kinetic_energy')] * data['cell_volume']).convert_to_units('erg')

        result = 1 * ((TE + KE) + PE < 0.0)

        return result


    def _mag_cyl_r(field,data):
        return np.abs( data[('index','cylindrical_radius')].convert_to_units('cm'))

    def _mag_cyl_z(field,data):
        return np.abs( data[('index','cylindrical_z')].convert_to_units('cm') )

    def _dm_density(field, data):
        r     = data[('index','spherical_r')].convert_to_units('cm')
        r_s   = (data.ds.parameters['DiskGravityDarkMatterR'] * yt.units.Mpc).convert_to_units('cm')
        rho_o = (data.ds.parameters['DiskGravityDarkMatterDensity'] * yt.units.g / yt.units.cm**3)

        rho = dm_halo.burkert_density(r, r_s, rho_o)

        return rho.convert_to_cgs()

    def _dm_potential(field, data):
        r     = data[('index','spherical_r')].convert_to_units('cm')
        r_s   = (data.ds.parameters['DiskGravityDarkMatterR'] * yt.units.Mpc).convert_to_units('cm')
        rho_o = (data.ds.parameters['DiskGravityDarkMatterDensity'] * yt.units.g / yt.units.cm**3)

        pot = dm_halo.burkert_potential(r, r_s, rho_o)

        return pot.convert_to_cgs()


    yt.add_field(('index','DM_background_density'), function = _dm_density, units = 'g/cm**3')
    yt.add_field(('index','DM_background_potential'), function = _dm_potential, units = 'erg/g')

    yt.add_field(('index','magnitude_cylindrical_radius'), function = _mag_cyl_r, units = 'cm')

    yt.add_field(('index','magnitude_cylindrical_z'), function = _mag_cyl_z, units = 'cm')
#    def _H2_total_mass(field, data):
#        mass = data[('gas',

    yt.add_field(('gas','Pe_heating_rate'), function = _pe_heating_cgs, units = 'erg/s/cm**3')
    yt.add_field(('gas','H_total_mass'), function = _H_total_mass, units ='g')
    yt.add_field(('gas','He_total_mass'), function = _He_total_mass, units = 'g')
    yt.add_field(('gas','metal_mass'), function = _metal_total_mass, units = 'g')
    yt.add_field(('gas','OTLW_kdissH2I'), function = _otlwcgs, units = '1/s')
    yt.add_field(('gas','Pe_heating_rate_masked'), function = _pe_heating_rate, units='erg/s/cm**3')
    yt.add_field(('gas','G_o'), function = _G_o, units = "")
#    yt.add_field(('gas','H2_total_mass'), function = _H2_total_mass, units = 'g')
#    yt.add_field(('gas','All_H_total_mass'), function = _all_H_total_mass, units = 'g')

    if ('enzo','PotentialField') in fields or ('enzo', 'GravPotential') in fields:
        yt.add_field(('gas','pos_gravitational_potential'), function=_grav_pot, units = 'erg/g')
        yt.add_field(('gas','total_gravitational_potential'), function=_tot_grav_pot, units = 'erg/g')
        yt.add_field(('gas','pos_total_gravitational_potential'), function=_pos_tot_grav_pot, units = 'erg/g')
        yt.add_field(('gas','potential_energy'), function=_potential_energy, units = 'erg')
        yt.add_field(('gas','gravitationally_bound'), function=_grav_bound, units = "")

    nfields = 5

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
    metals = utilities.species_from_fields(fields)
    ratios = utilities.ratios_list(metals)

    # make new functions to do correct units for species fields
    _density_function_generator(metals + ['Metal'])

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

    nfields = _additional_helper_fields(fields)
    print nfields, "additional helper fields defined"

    FIELDS_DEFINED = True
    return
