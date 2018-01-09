
__author__ = "Andrew Emerick"
__email__  = "aemerick11@gmail.com"

import yt
yt.funcs.mylog.setLevel(40)

from yt.units import dimensions
import numpy as np

from collections import Iterable

from galaxy_analysis.static_data import AMU,\
                 MOLECULAR_WEIGHT

from galaxy_analysis.utilities import convert_abundances
from galaxy_analysis.utilities import utilities
from galaxy_analysis import star_analysis
from galaxy_analysis.misc import dm_halo

from onezone import data_tables

SE_table = data_tables.StellarEvolutionData()

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

def _abundance_function_generator(asym):

    if not isinstance(asym, Iterable):
        asym = [asym]

    def return_function(a):
        def _abundance(field,data):
            mass = data[('gas', a + '_Mass')].convert_to_units('g').value
            abund = convert_abundances.elemental_abundance( a, mass)
            return abund

        return _abundance

    if not ('H' in asym):
        asym = asym + ['H']

    for a in asym:
        yt.add_field(('gas',a + '_Abundance'), function = return_function(a), units = "")

    if (('O' in asym) and ('Mg' in asym) and ('Si' in asym)):
        def _alpha_abundance(field, data):
            alpha = data[('gas','O_Abundance')] + data[('gas','Mg_Abundance')] + data[('gas','Si_Abundance')]
            return alpha / 3.0
        yt.add_field(('gas','alpha_Abundance'), function=_alpha_abundance, units="")

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

    if (('O' in asym) and ('Mg' in asym) and ('Si' in asym)):
        def _alpha_mass(field, data):
            alpha = data[('gas','O_Mass')] + data[('gas','Mg_Mass')] + data[('gas','Si_Mass')]
            return alpha

        yt.add_field(('gas','alpha_Mass'), function = _alpha_mass, units = "g") # mass of alpha elements

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

    if (('O' in asym) and ('Mg' in asym) and ('Si' in asym)):
        def _alpha_mass_fraction(field, data):
            alpha = data[('gas','O_Fraction')] + data[('gas','Mg_Fraction')] +\
                        data[('gas','Si_Fraction')]
            return alpha / 3.0

#        yt.add_field(('gas','alpha_Fraction'), function = _alpha_mass_fraction, units = "")

        if (('S' in asym) and ('Ca' in asym)):
            def _alpha_5(field, data):
                alpha = data[('gas','alpha_Fraction')]*3.0 + data[('gas','S_Fraction')] +\
                                                 data[('gas','Ca_Fraction')]
                return alpha / 5.0

#            yt.add_field( ('gas','alpha_5_Fraction'), function = _alpha_5, units = "")

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

    def _H_total_number_density(field,data):
        n_H = data['H_p0_number_density'] + data['H_p1_number_density']

        try:
            n_H += data['H_m1_number_density'] +\
             0.5 * (data['H2_p0_number_density'] + data['H2_p1_number_density'])

        except:
            n_H += np.zeros(np.shape(n_H))

        return n_H.convert_to_units('cm**(-3)')

    yt.add_field(('gas','H_total_number_density'),
                 function = _H_total_number_density, units = 'cm**(-3)')

    return nfields

def _particle_abundance_function_generator(asym, ds = None):

    if not (ds is None):
        if not (ds.parameters['NumberOfParticles'] > 0):
            return

    if not isinstance(asym, Iterable):
        asym = [asym]

    if not ('H' in asym):
        asym = asym + ['H']

    if not ('He' in asym):
        asym = asym + ['He']

    def return_function(element, fraction_field):
        def _abundance(field, data):
            mass = data[fraction_field].value * (data['birth_mass'].value *yt.units.Msun).convert_to_units('g').value
            abund = convert_abundances.elemental_abundance(element, mass)
            return abund
        return _abundance

    for a in asym:
        fraction_field = ('io','particle_' + a + '_fraction')
        yt.add_field(('io','particle_' + a + '_abundance'),
                      return_function(a, fraction_field), units = "", particle_type = True)

    if (('O' in asym) and ('Mg' in asym) and ('Si' in asym)):
        def _alpha_fraction(field, data):
            alpha = data[('io','particle_O_fraction')] +\
                    data[('io','particle_Mg_fraction')] +\
                    data[('io','particle_Si_fraction')]
            return alpha / 3.0

        def _alpha_abundance(field, data):
            alpha = data[('io','particle_O_abundance')] +\
                    data[('io','particle_Mg_abundance')] +\
                    data[('io','particle_Si_abundance')]

            return alpha / 3.0

#        yt.add_field(('io','particle_alpha_fraction'), function=_alpha_fraction, units = "", particle_type = True)
        yt.add_field(('io','particle_alpha_abundance'), function=_alpha_abundance, units="", particle_type=True)

        if ('S' in asym) and ('Ca' in asym):
            def _alpha_5_fraction(field,data):
                alpha = (data[('io','particle_alpha_fraction')]*3.0 + data[('io','particle_S_fraction')] +\
                                        data[('io','particle_Ca_fraction')]) / 5.0
                return alpha

            def _alpha_5(field,data):
                alpha = data[('io','particle_alpha_abundance')]*3.0 + data[('io','particle_S_abundance')]+\
                                        data[('io','particle_Ca_abundance')]
                return alpha / 5.0

#            yt.add_field( ('io','particle_alpha_5_fraction'), function = _alpha_5_fraction, units = "", particle_type =True)
            yt.add_field( ('io','particle_alpha_5_abundance'), function = _alpha_5, units = "", particle_type = True)

    return

def _particle_abundance_ratio_function_generator(ratios, ds = None):

    if not (ds is None):
        if not (ds.parameters['NumberOfParticles'] > 0):
            return

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

    def _alpha_return_function(base):
        def _alpha_over_x(field,data):
            alpha = data[('io','particle_alpha_abundance')]
            x     = data[('io','particle_' + base + '_abundance')]

            return convert_abundances.abundance_ratio(('alpha',alpha),(base,x), 'abundances')

        return _alpha_over_x

    def _alpha_5_return_function(base):
        def _alpha_5_over_x(field,data):
            alpha = data[('io','particle_alpha_5_abundance')]
            x     = data[('io','particle_' + base + '_abundance')]

            return convert_abundances.abundance_ratio(('alpha_5',alpha),(base,x), 'abundances')

        return _alpha_5_over_x

    denoms = [x.split('/')[1] for x in ratios]
    denoms = np.unique(denoms)
    for x in denoms:
        yt.add_field(('io','particle_alpha_over_' + x), function = _alpha_return_function(x), units = "", particle_type = True)
        yt.add_field(('io','particle_alpha_5_over_' + x), function = _alpha_5_return_function(x), units = "", particle_type = True)


#    def _alpha_over_Fe(field,data):
#        alpha = data[('io','particle_alpha_abundance')]
#        Fe    = data[('io','particle_Fe_abundance')]
#        return convert_abundances.abundance_ratio( ('alpha', alpha), ('Fe', Fe), 'abundances')

#    yt.add_field(('io','particle_alpha_over_Fe'), function = _alpha_over_Fe, units = "", particle_type = True)


#    def _alpha_5_over_Fe(field,data):
#        alpha = data[('io','particle_alpha_5_abundance')]
#        Fe    = data[('io','particle_Fe_abundance')]
#        return convert_abundances.abundance_ratio( ('alpha_5', alpha), ('Fe', Fe), 'abundances')

 #   yt.add_field(('io','particle_alpha_5_over_Fe'), function = _alpha_5_over_Fe, units = "", particle_type = True)


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

    def _return_alpha_over_x(element_name):
        def _alpha_over_x(field, data):
            alpha = data[('gas','alpha_Abundance')]
            x     = data[('gas',element_name + '_Abundance')]
            return convert_abundances.abundance_ratio(('alpha',alpha),(element_name,x),'abundances')

        return _alpha_over_x

    denoms = [x.split('/')[1] for x in ratios]
    denoms = np.unique(denoms)
    for x in denoms:
        yt.add_field(('gas','alpha_over_' + x), function = _return_alpha_over_x(x), units = "")

    return nfields


def generate_stellar_model_fields(ds):

    if not (ds.parameters['NumberOfParticles'] > 0):
        return

    #
    # luminosity, L_FUV, L_LW, Q0, Q1, E0, E1
    #

    field_names = ['luminosity','L_FUV','L_LW','Q0','Q1','E0','E1']

    units = {'luminosity' : yt.units.erg/yt.units.s,
             'L_FUV' : yt.units.erg/yt.units.s,
             'L_LW'  : yt.units.erg/yt.units.s,
             'Q0' : 1.0 /yt.units.s, 'Q1' : 1.0 / yt.units.s,
             'E0' : yt.units.eV, 'E1' : yt.units.eV, 'lifetime' : yt.units.s}

    unit_label = {'luminosity': 'erg/s', 'L_FUV' : 'erg/s', 'L_LW' : 'erg/s',
                  'Q0' : '1/s', 'Q1' : '1/s', 'E0' : 'erg', 'E1': 'erg', 'lifetime' : 's'}

    overload_type = {} # when generating stars, default is use type from simulation
    for k in units.keys():
        overload_type[k] = None # keep simulation type when making stars for all fields

    overload_type['lifetime'] = 11 # lifetime field will now be the lifetime of the original MS star
                                   # and NOT the lifetime from the simulation --- this is done by
                                   # overloading the particle type to 11, or main sequence. This is
                                   # also more useful, as the field would just be redundant
                                   # info if this isn't done. See star_analysis code.

    def _function_generator(field_name):
        def _function(field, data):
            if np.size(data['particle_mass']) == 1:
                # this is ugly, but a way to bypass yt's validation step
                # because throwing junk values into the routine below will cause problems
                with utilities.nostdout():
                    p = star_analysis.get_star_property(ds, data, property_names = [field_name],
                                                            dummy_call = True, overload_type = None)
            else:
                p = star_analysis.get_star_property(ds, data, property_names = [field_name],
                                                        overload_type = overload_type[field_name] )
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

    def _lifetime(field, data):
        m = data['birth_mass'].value
        z = data['metallicity_fraction'].value

        if np.size(m) == 1:  # get around yt's checking
            lt = np.shape(m) * yt.units.Myr
        else:
            lt = np.zeros(np.size(m))
            for i in np.arange(np.size(m)):
                lt[i] = SE_table.interpolate({'mass' : m[i], 'metallicity' : z[i]}, 'lifetime')
            lt = (lt * yt.units.s).convert_to_units('Myr')

        return lt

    yt.add_field(('io','particle_model_lifetime'), function = _lifetime, units = 'Myr',
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
        n_H = (data['H_p0_number_density'] + data['H_p1_number_density'] + data['H_m1_number_density'] +\
                0.5*(data['H2_p0_number_density'] + data['H2_p1_number_density'])).convert_to_units('cm**(-3)').value
        g_to_d = 0.68 - 3.08 * np.log10(Z / 0.014)
        d_to_g = 1.0 / (10.0**(g_to_d))
        D = d_to_g / 6.616595E-3
        epsilon = 0.01488637246 * (n_H)**(0.235269059)
        atten = np.exp( - 1.33E-21 * D * data['dx'].convert_to_units('cm').value * n_H)
        G_o = pe / (1.3E-24 * n_H * epsilon * D * atten)
        return G_o * (data['Density'] / data['Density'])

    def _G_eff(field,data):
        pe  = data[('gas','Pe_heating_rate')].convert_to_units('erg/s/cm**3').value
        Z   = (data['Metal_Density'] / data['Density']).value
        n_H = (data['H_p0_number_density'] + data['H_p1_number_density'] + data['H_m1_number_density'] +\
                0.5*(data['H2_p0_number_density'] + data['H2_p1_number_density'])).convert_to_units('cm**(-3)').value

        g_to_d = 0.68 - 3.08 * np.log10(Z / 0.014)

        d_to_g = 1.0 / (10.0**(g_to_d))

        D = d_to_g / 6.616595E-3

        epsilon = 0.01488637246 * (n_H)**(0.235269059)

        # atten = np.exp( - 1.33E-21 * D * data['dx'].convert_to_units('cm').value * n_H)

        G_eff = pe / (1.3E-24 * n_H * epsilon * D)

        return G_eff * (data['Density'] / data['Density'])


    def _FUV_flux(field, data):
        # 1.59E-3 converts from MW normalized flux density to flux dens in cgs
        G_o = data[('gas','G_o')] # relative to MW
        G   = (G_o.value * 1.59E-3) * yt.units.erg / yt.units.cm**2 /yt.units.s
        return G

    def _LW_flux(field, data):
        LW_energy = 12.8 * yt.units.eV
        H2Isigma  = 3.71E-18 * yt.units.cm**(2)

        kdissH2I = (data[('enzo','OTLW_kdissH2I')].value / data.ds.time_unit).convert_to_units('1/s')

        LW_flux = kdissH2I / H2Isigma * LW_energy

        return LW_flux.convert_to_units('erg/cm**2/s')

    def _Q0_flux(field, data):
        E_HI = 13.6 * yt.units.eV
        kph = data[('enzo','HI_kph')].convert_to_cgs()
        n   = data[('gas','H_p0_number_density')].convert_to_cgs()
        dt  = data.ds.parameters['dtPhoton']
        V   = data['cell_volume'].convert_to_cgs()
        dx  = data['dx'].convert_to_cgs()
        s   = 6.34629E-18 * yt.units.cm**(2) # cross section of HI at 13.6 eV

        tau   = s * n * dx
        denom = 1.0 - np.exp(-tau)

        Q = kph * n * V / denom  # this gives number of photons / s

        flux = Q * E_HI / dx**2

        return flux.convert_to_units('erg/cm**2/s')

    def _Q1_flux(ds,data):
        E_HeI = 24.6 * yt.units.eV
        kph = data[('enzo','HeI_kph')].convert_to_cgs()
        n   = data[('gas','H_p0_number_density')].convert_to_cgs()
        dt  = data.ds.parameters['dtPhoton']
        V   = data['cell_volume'].convert_to_cgs()
        dx  = data['dx'].convert_to_cgs()
        s   = 7.4300459E-18 * yt.units.cm**(2) # cross section of HeI at 24.6 eV

        tau   = s * n * dx
        denom = 1.0 - np.exp(-tau)

        Q = kph * n * V / denom  # this gives number of photons / s

        flux = Q * E_HeI / dx**2

        return flux.convert_to_units('erg/cm**2/s')


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

    def _gas_grav_pot(field,data):
        try:
            x = data['PotentialField']
        except:
            x = data['GravPotential'].value * data.ds.velocity_unit**2

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

    def _rad_accel(field, data):
        return np.sqrt(data['RadAccel1']**2 + data['RadAccel2']**2 + data['RadAccel3']**2).convert_to_units('cm/s**2')

    def _is_star_forming(field, data):
        n      = data[('gas','number_density')].convert_to_units('cm**(-3)')
        T      = data['Temperature']
        divv   = data[('gas','velocity_divergence')]
        l      = data['grid_level']

        answer = 1 * ((n > data.ds.parameters['StarMakerOverDensityThreshold']) *\
                       (T < data.ds.parameters['IndividualStarTemperatureThreshold']) *\
                        (divv < 0) *\
                          (l == data.ds.parameters['MaximumRefinementLevel']))

        return answer

    yt.add_field(('gas','is_star_forming'), function = _is_star_forming,
                         units = "")

    yt.add_field(("gas","a_rad"), function=_rad_accel, units="cm/s**2")

    yt.add_field(('index','DM_background_density'), function = _dm_density, units = 'g/cm**3')
    yt.add_field(('index','DM_background_potential'), function = _dm_potential, units = 'erg/g')

    yt.add_field(('index','magnitude_cylindrical_radius'), function = _mag_cyl_r, units = 'cm')

    yt.add_field(('index','magnitude_cylindrical_z'), function = _mag_cyl_z, units = 'cm')
#    def _H2_total_mass(field, data):
#        mass = data[('gas',

    yt.add_field(('gas','Pe_heating_rate'), function = _pe_heating_cgs, units = 'erg/s/cm**3')
    yt.add_field(('gas','H_total_mass'), function = _H_total_mass, units ='g')
    yt.add_field(('gas','H_Mass'), function = _H_total_mass, units = 'g') # define as same
    yt.add_field(('gas','He_total_mass'), function = _He_total_mass, units = 'g')
    yt.add_field(('gas','metal_mass'), function = _metal_total_mass, units = 'g')
    yt.add_field(('gas','OTLW_kdissH2I'), function = _otlwcgs, units = '1/s')
    yt.add_field(('gas','Pe_heating_rate_masked'), function = _pe_heating_rate, units='erg/s/cm**3')
    yt.add_field(('gas','G_o'), function = _G_o, units = "")
    yt.add_field(('gas','G_eff'), function = _G_eff, units = "")
    yt.add_field(('gas','FUV_flux'), function = _FUV_flux, units = "erg/s/cm**2")
    yt.add_field(('gas','LW_flux'), function = _LW_flux, units = "erg/s/cm**2")
    yt.add_field(('gas','Q0_flux'), function = _Q0_flux, units = "erg/s/cm**2")
    yt.add_field(('gas','Q1_flux'), function = _Q1_flux, units = "erg/s/cm**2")
#    yt.add_field(('gas','H2_total_mass'), function = _H2_total_mass, units = 'g')
#    yt.add_field(('gas','All_H_total_mass'), function = _all_H_total_mass, units = 'g')

    if ('enzo','PotentialField') in fields or ('enzo', 'GravPotential') in fields:
        yt.add_field(('gas','pos_gravitational_potential'), function=_grav_pot, units = 'erg/g')
        yt.add_field(('gas','gas_gravitational_potential'), function=_gas_grav_pot, units = 'erg/g')
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
    nfields = _abundance_function_generator(metals)

    if ds.parameters['NumberOfParticles'] > 0:
        nfields =  _particle_abundance_ratio_function_generator(ratios, ds)
        print nfields, "particle abundance ratio fields defined"
        _particle_abundance_function_generator(metals, ds)


    generate_stellar_model_fields(ds)

    nfields = _additional_helper_fields(fields)
    print nfields, "additional helper fields defined"

    FIELDS_DEFINED = True
    return


def load_and_define(name):
    """
    Wrapper around yt to load a data set and define gradient
    fields and particle filters which must be defined for each 
    simulation file separately (unlike the above)
    """

    ds = yt.load(name)

    generate_gradient_fields(ds)

    def _grav_accel_x(field,data):
        return data[('gas','gas_gravitational_potential_gradient_x')].convert_to_units('cm/s**2')
    def _grav_accel_y(field,data):
        return data[('gas','gas_gravitational_potential_gradient_y')].convert_to_units('cm/s**2')
    def _grav_accel_z(field,data):
        return data[('gas','gas_gravitational_potential_gradient_z')].convert_to_units('cm/s**2')
    def _grav_accel(field,data):
        return np.sqrt(data[('gas','a_grav_x')]**2 + data[('gas','a_grav_y')]**2 + data[('gas','a_grav_z')]**2)

    def _a_rad_a_grav(field,data):
        a = data[('gas','a_rad')] / data[('gas','a_grav')]
    
        a[data[('gas','a_grav')] == 0.0] = 0.0

        return a

    ds.add_field(('gas','a_grav_x'), function = _grav_accel_x, units = 'cm/s**2', sampling_type='cell')
    ds.add_field(('gas','a_grav_y'), function = _grav_accel_y, units = 'cm/s**2', sampling_type='cell')
    ds.add_field(('gas','a_grav_z'), function = _grav_accel_z, units = 'cm/s**2', sampling_type='cell')
    ds.add_field(('gas','a_grav'),   function = _grav_accel,   units = 'cm/s**2', sampling_type='cell')
    ds.add_field(('gas','a_rad_over_a_grav'), function = _a_rad_a_grav, units = '', sampling_type = 'cell')

    generate_particle_filters(ds)
    
    return ds

def generate_gradient_fields(ds):
    """
    generate gas self gravity gradient fields and rename them to
    something sensible
    """

    ds.add_gradient_fields(("gas","gas_gravitational_potential"))

    return

def generate_particle_filters(ds):
    """
    Make filter definitions for the various particle types:
        Main Sequence :
        White Dwarf   :
        SNIa remnant  :
        SNII remnant  :
        AGB phase (likely very few or none since short lived) :
    """ 

    @yt.particle_filter(requires=["particle_type"], filtered_type='all')
    def main_sequence_stars(pfilter, data):
        filter = data[(pfilter.filtered_type, "particle_type")] == 11
        return filter

    ds.add_particle_filter('main_sequence_stars')


    return
