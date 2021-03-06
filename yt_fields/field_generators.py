
__author__ = "Andrew Emerick"
__email__  = "aemerick11@gmail.com"

import yt
yt.funcs.mylog.setLevel(40)
from yt.fields.api import ValidateDataField, ValidateParameter
from yt.units import dimensions
import numpy as np

from collections import Iterable

from galaxy_analysis.static_data import AMU,\
                 MOLECULAR_WEIGHT

from galaxy_analysis.utilities import convert_abundances
from galaxy_analysis.utilities import utilities
from galaxy_analysis import star_analysis
from galaxy_analysis.misc import dm_halo
from galaxy_analysis import physics

from galaxy_analysis.yt_fields import ionization

from onezone import data_tables, radiation



GRACKLE_IMPORTED = True
try:
    import pygrackle
    from pygrackle.grackle_wrapper import \
         calculate_cooling_time
except:
    GRACKLE_IMPORTED = False


SE_table = data_tables.StellarEvolutionData()

FIELDS_DEFINED = False

def _density_function_generator(asym):
    if not isinstance(asym, Iterable):
        asym = [asym]

    def return_function(a):
        def _density(field, data):
            dens = data[('enzo', a + '_Density')].value
            dens = dens * data.ds.mass_unit / data.ds.length_unit**3
            return dens.to('g/cm**3')

        return _density

    for a in asym:
        yt.add_field(('gas', a + "_Density"), sampling_type = 'cell',
                     function = return_function(a), units = 'g/cm**3')

    return

def _abundance_function_generator(asym):

    if not isinstance(asym, Iterable):
        asym = [asym]

    def return_function(a):
        def _abundance(field,data):
            mass = data[('gas', a + '_Mass')].to('g').value
            abund = convert_abundances.elemental_abundance( a, mass)
            return abund

        return _abundance

    if not ('H' in asym):
        asym = asym + ['H']

    for a in asym:
        yt.add_field(('gas',a + '_Abundance'), sampling_type = 'cell',
                      function = return_function(a), units = "")

    if (('O' in asym) and ('Mg' in asym) and ('Si' in asym)):
        def _alpha_abundance(field, data):
            alpha = data[('gas','O_Abundance')] + data[('gas','Mg_Abundance')] + data[('gas','Si_Abundance')]
            return alpha / 3.0
        yt.add_field(('gas','alpha_Abundance'), sampling_type = 'cell', function=_alpha_abundance, units="")

    return


def _mass_function_generator(asym):

    if not isinstance(asym, Iterable):
        asym = [asym]

    def return_function(a):
        def _mass(field,data):
            ele_dens = data[('enzo', a + '_Density')].value
            ele_dens = ele_dens * data.ds.mass_unit / data.ds.length_unit**3
            ele_dens.convert_to_cgs()

            mass = (ele_dens * data['cell_volume']).to('g')

            return mass

        return _mass

    nfields = 0
    for a in asym:

        yt.add_field(('gas', a + '_Mass'), sampling_type = 'cell', function = return_function(a), units='g')
        nfields = nfields + 1

    if (('O' in asym) and ('Mg' in asym) and ('Si' in asym)):
        def _alpha_mass(field, data):
            alpha = data[('gas','O_Mass')] + data[('gas','Mg_Mass')] + data[('gas','Si_Mass')]
            return alpha

        yt.add_field(('gas','alpha_Mass'), sampling_type = 'cell', function = _alpha_mass, units = "g") # mass of alpha elements

        nfields = nfields + 1

    return nfields

#
# Construct arbitrary mass fraction derived fields in yt
# using a loop to generate functions
#
def _mass_fraction_function_generator(ds, asym):

    if not isinstance(asym, Iterable):
        asym = [asym]


    def return_function(a):
        def _mass_fraction(field,data):
            ele_dens = data[('enzo', a + '_Density')].value
            ele_dens = ele_dens * data.ds.mass_unit / data.ds.length_unit**3
            ele_dens.convert_to_cgs()

            dens = data[('enzo','Density')].to('g/cm**3')
            mass_fraction = ele_dens / dens
            return mass_fraction

        return _mass_fraction

    def return_function_2(a):
        def _mass_fraction_2(field,data):
            ele_dens = data[('enzo', a + '_Density_2')].value
            ele_dens = ele_dens * data.ds.mass_unit / data.ds.length_unit**3
            ele_dens.convert_to_cgs()

            dens = data[('enzo','Density')].to('g/cm**3')
            mass_fraction = ele_dens / dens
            return mass_fraction

        return _mass_fraction_2



    nfields = 0
    for a in asym:
        yt.add_field(('gas', a + '_Fraction'),  sampling_type = 'cell', function = return_function(a), units="")
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


    if 'IndividualStarTrackAGBMetalDensity' in ds.parameters:
        if ds.parameters['IndividualStarTrackAGBMetalDensity']:
            def _AGB_mass_fraction(field,data):
                ele_dens = data[('enzo', 'AGB_Metal_Density')].value
                ele_dens = ele_dens * data.ds.mass_unit / data.ds.length_unit**3
                ele_dens.convert_to_cgs()

                dens = data[('enzo','Density')].to('g/cm**3')
                mass_fraction = ele_dens / dens
                return mass_fraction

            yt.add_field(('gas', 'AGB_Mass_Fraction'), sampling_type = 'cell', function = _AGB_mass_fraction, units="")

        if ds.parameters['IndividualStarTrackWindDensity']:
            def _IWIND_mass_fraction(field,data):
                ele_dens = data[('enzo', 'Intermediate_Wind_Metal_Density')].value
                ele_dens = ele_dens * data.ds.mass_unit / data.ds.length_unit**3
                ele_dens.convert_to_cgs()

                dens = data[('enzo','Density')].to('g/cm**3')
                mass_fraction = ele_dens / dens
                return mass_fraction
            yt.add_field(('gas', 'Intermediate_Wind_Mass_Fraction'), sampling_type = 'cell', function = _IWIND_mass_fraction, units="")
            def _MWIND_mass_fraction(field,data):
                ele_dens = data[('enzo', 'Massive_Wind_Metal_Density')].value
                ele_dens = ele_dens * data.ds.mass_unit / data.ds.length_unit**3
                ele_dens.convert_to_cgs()

                dens = data[('enzo','Density')].to('g/cm**3')
                mass_fraction = ele_dens / dens
                return mass_fraction
            yt.add_field(('gas', 'Massive_Wind_Mass_Fraction'), sampling_type = 'cell', function = _MWIND_mass_fraction, units="")

        if ds.parameters['IndividualStarTrackSNMetalDensity']:
            def _SNII_mass_fraction(field,data):
                ele_dens = data[('enzo', 'SNII_Metal_Density')].value
                ele_dens = ele_dens * data.ds.mass_unit / data.ds.length_unit**3
                ele_dens.convert_to_cgs()

                dens = data[('enzo','Density')].to('g/cm**3')
                mass_fraction = ele_dens / dens
                return mass_fraction

            yt.add_field(('gas', 'SNII_Mass_Fraction'), sampling_type = 'cell', function = _SNII_mass_fraction, units="")

            def _SNIa_mass_fraction(field,data):
                ele_dens = data[('enzo', 'SNIa_Metal_Density')].value
                ele_dens = ele_dens * data.ds.mass_unit / data.ds.length_unit**3
                ele_dens.convert_to_cgs()

                dens = data[('enzo','Density')].to('g/cm**3')
                mass_fraction = ele_dens / dens
                return mass_fraction

            yt.add_field(('gas', 'SNIa_Mass_Fraction'), sampling_type = 'cell', function = _SNIa_mass_fraction, units="")

            if ds.parameters['IndividualStarSNIaModel'] == 2:
                def _SNIa_sCh_mass_fraction(field,data):
                    ele_dens = data[('enzo', 'SNIa_sCh_Metal_Density')].value
                    ele_dens = ele_dens * data.ds.mass_unit / data.ds.length_unit**3
                    ele_dens.convert_to_cgs()

                    dens = data[('enzo','Density')].to('g/cm**3')
                    mass_fraction = ele_dens / dens
                    return mass_fraction

                yt.add_field(('gas', 'SNIa_sCh_Mass_Fraction'), sampling_type = 'cell', function = _SNIa_sCh_mass_fraction, units="")
                def _SNIa_SDS_mass_fraction(field,data):
                    ele_dens = data[('enzo', 'SNIa_SDS_Metal_Density')].value
                    ele_dens = ele_dens * data.ds.mass_unit / data.ds.length_unit**3
                    ele_dens.convert_to_cgs()

                    dens = data[('enzo','Density')].to('g/cm**3')
                    mass_fraction = ele_dens / dens
                    return mass_fraction

                yt.add_field(('gas', 'SNIa_SDS_Mass_Fraction'), sampling_type = 'cell', function = _SNIa_SDS_mass_fraction, units="")
                def _SNIa_HeRS_mass_fraction(field,data):
                    ele_dens = data[('enzo', 'SNIa_HeRS_Metal_Density')].value
                    ele_dens = ele_dens * data.ds.mass_unit / data.ds.length_unit**3
                    ele_dens.convert_to_cgs()

                    dens = data[('enzo','Density')].to('g/cm**3')
                    mass_fraction = ele_dens / dens
                    return mass_fraction

                yt.add_field(('gas', 'SNIa_HeRS_Mass_Fraction'), sampling_type = 'cell', function = _SNIa_HeRS_mass_fraction, units="")

            if ds.parameters['IndividualStarPopIIIFormation']:
                def _PopIII_mass_fraction(field,data):
                    ele_dens = data[('enzo', 'PopIII_Metal_Density')].value
                    ele_dens = ele_dens * data.ds.mass_unit / data.ds.length_unit**3
                    ele_dens.convert_to_cgs()

                    dens = data[('enzo','Density')].to('g/cm**3')
                    mass_fraction = ele_dens / dens
                    return mass_fraction

                def _PopIII_PISNe_mass_fraction(field,data):
                    ele_dens = data[('enzo', 'PopIII_PISNe_Metal_Density')].value
                    ele_dens = ele_dens * data.ds.mass_unit / data.ds.length_unit**3
                    ele_dens.convert_to_cgs()

                    dens = data[('enzo','Density')].to('g/cm**3')
                    mass_fraction = ele_dens / dens
                    return mass_fraction

                yt.add_field(('gas', 'PopIII_Mass_Fraction'), sampling_type = 'cell', function = _PopIII_mass_fraction, units="")
                yt.add_field(('gas', 'PopIII_PISNe_Mass_Fraction'), sampling_type = 'cell', function = _PopIII_PISNe_mass_fraction, units="")


                for a in asym:
                    yt.add_field(('gas', a + '_PopIII_Fraction'),  sampling_type = 'cell', function = return_function_2(a), units="")
                    nfields = nfields + 1

            if ds.parameters['IndividualStarRProcessModel']:
                def _RProcess_mass_fraction(field,data):
                    ele_dens = data[('enzo', 'RProcess_Metal_Density')].value
                    ele_dens = ele_dens * data.ds.mass_unit / data.ds.length_unit**3
                    ele_dens.convert_to_cgs()

                    dens = data[('enzo','Density')].to('g/cm**3')
                    mass_fraction = ele_dens / dens
                    return mass_fraction

                yt.add_field(('gas', 'RProcess_Mass_Fraction'), sampling_type = 'cell', function = _RProcess_mass_fraction, units="")


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
            ele_dens.convert_to_cgs()
            n = ele_dens / (MOLECULAR_WEIGHT[a] * AMU * yt.units.g)

            n.convert_to_cgs()
            return n

        return _number_density

    nfields = 0
    for a in asym:

        yt.add_field(('gas', a + '_Number_Density'), sampling_type = 'cell',
                     function = return_function(a), units='cm**(-3)')
        nfields = nfields + 1

    # make a metal number density field - make an assumption about the metal molecular weight
    def _metal_number_density(field,data):
        ele_dens = data[('enzo','Metal_Density')].to('g/cm**3')
        n = ele_dens / (MOLECULAR_WEIGHT['metal'] * AMU * yt.units.g)

        n.convert_to_cgs()
        return n

    yt.add_field(('gas', 'Metal_Number_Density'),sampling_type = 'cell',
                 function = _metal_number_density, units = 'cm**(-3)')

    def _H_total_number_density(field,data):
        n_H = data['H_p0_number_density'] + data['H_p1_number_density']

        try:
            n_H += data['H_m1_number_density'] +\
             0.5 * (data['H2_p0_number_density'] + data['H2_p1_number_density'])

        except:
            n_H += np.zeros(np.shape(n_H))

        return n_H.to('cm**(-3)')

    yt.add_field(('gas','H_total_number_density'), sampling_type = 'cell',
                 function = _H_total_number_density, units = 'cm**(-3)')

    return nfields

def _ionization_state_generator(metals):

    temp_ions     = ionization.get_ions()

    all_ions      = [x for x in temp_ions if ionization.get_elements(x) in metals]


    def return_function(ion):
        def _ion_density(field,data):
            ele    = ionization.get_elements(ion)
            n_ele  = data[('gas', ele + '_Number_Density')].value

            n      = data['H_total_number_density'].to('cm**(-3)').value
            n      = np.log10(n)
            T      = data['Temperature'].to('K').value
            T      = np.log10(T)

            f_ion  = 10.0**(ionization.get_ion_fraction(n, T, ion))

            return (n_ele * f_ion) * yt.units.cm**(-3)

        return _ion_density

    nfields = 0
    for ion in all_ions:

        yt.add_field(('gas', ion + '_Number_Density'), sampling_type = 'cell',
                     function = return_function(ion), units='cm**(-3)')
        nfields = nfields + 1

    return nfields

def _generate_rates(ds):
    """
    Generate reaction rate equations
    """

    kunit = 1.0 #

    def k8(T):
        T = T.to('K').value
        k8 = 1.35E-9 * (T**(9.8493E-2) + 3.2852E-1 *\
                          T**(5.561E-1) + 2.881E-7 * T**2.1826) /\
            (1.0 + 6.191E-3 * T**1.0461 + 8.9712E-11*T**3.0424 +\
               3.2576E-14 * T**3.7741)
        return k8  # now in cgs

    def k10(T): # have arg T to look same as other functions
        k10 = 6.0E-10

        return k10

    def k19(T):
        T = T.to('K').value
        k19 = 5.0E-7 * np.sqrt(100.0 / T)
        return k19

    def k22(T):
        T = T.to('K').value
        # for GLover 2008 three body rate ONLY
        k22 = 7.7E-31 / T**0.464
        return k22

    def k13(T):
        T = T.to('K').value
        k13 = 10.0**(-178.4239 - 68.42243 * np.log10(T)
                        + 43.20243 * np.log10(T)**2
                        - 4.633167 * np.log10(T)**3
                        + 69.70086 * np.log10(1.0 + 40870.38 / T)
                        - (23705.7 / T))
        ###############
#        above  is for use with Glover 2008 three body rate
#
#        T_eV = (T / yt.physical_constants.k_b).to('eV').value
#        T_lim = 0.3
#
#        k13 = np.ones(np.shape(T)) * 1.0E-20
#
#        k13[ T > T_lim] = 1.0670825E-10*T_eV**(2.012) /\
#               (np.exp(4.463/T_eV) * (1.0 + 0.2472 * T_eV)**3.512)

        return k13

    def k11(T):
        T_eV = (T / yt.physical_constants.k_b).to('eV').value
        T_lim = 0.3
        k11 = np.ones(np.shape(T)) * 1.0E-20

        log_T = np.log(T.to('K').value)

        k11[ T_eV > T_lim] = (np.exp(-21237.15/T) *\
                (- 3.3232183E-7
                 + 3.3735382E-7  * log_T
                 - 1.4491368E-7  * log_T**2
                 + 3.4172805E-8  * log_T**3
                 - 4.7813720E-9  * log_T**4
                 + 3.9731542E-10 * log_T**5
                 - 1.8171411E-11 * log_T**6
                 + 3.5311932E-13 * log_T**7))

        return k11

    def k12(T):
        T_eV = (T / yt.physical_constants.k_b).to('eV').value
        T_lim = 0.3
        k12  = np.ones(np.shape(T)) * 1.0E-20

        k12[T>T_lim] = 4.4886E-9*T**(0.109127)*np.exp(-101858.0/T)

        return k12

#    def k29(T)
#
#        return k29

    reaction_units = 1.0 / yt.units.cm**3 / yt.units.s
    ru_label = '1/s/cm**3'
    def _k8_reaction_rate(field, data):
        rr = k8(data['Temperature'].to('K'))
        rr = 2.0 * rr * data[('gas','H_m1_number_density')].to('cm**(-3)').value *\
                        data[('gas','H_p0_number_density')].to('cm**(-3)').value
        return rr * reaction_units

    def _k10_reaction_rate(field,data):
        rr = k10(data['Temperature'].to('K'))
        rr = rr * data[('gas','H2_p1_number_density')].to('cm**(-3)').value *\
                  data[('gas','H_p0_number_density')].to('cm**(-3)').value
        return rr * reaction_units

    def _k19_reaction_rate(field,data):
        rr = k19(data['Temperature'].to('K'))
        rr = rr * data[('gas','H2_p1_number_density')].to('cm**(-3)').value *\
                  data[('gas','H_m1_number_density')].to('cm**(-3)').value
        return rr * reaction_units

    def _k22_reaction_rate(field, data):
        rr = k22(data['Temperature'].to('K'))
        rr = rr * (data[('gas','H_p0_number_density')].to('cm**(-3)').value)**3
        return rr * reaction_units

    yt.add_field(('gas','k8_rr'), sampling_type = 'cell',
                 function = _k8_reaction_rate, units = ru_label)
    yt.add_field(('gas','k10_rr'), sampling_type = 'cell',
                 function = _k10_reaction_rate, units = ru_label)
    yt.add_field(('gas','k19_rr'), sampling_type = 'cell',
                 function = _k19_reaction_rate, units = ru_label)
    yt.add_field(('gas','k22_rr'), sampling_type = 'cell',
                 function = _k22_reaction_rate, units = ru_label)


    # rates
    # scoef
    #  2.0 * (    k8  *   HM * HI         - set with interp
    #             k10 * H2II * HI * 0.5   - set with interp
    #             k19 * H2II * HM * 0.5   - set with interp
    #             k22 * HI   * (HI*HI)    - set with interp
    #
    # acoef
    #    k13*HI + k11*HII + k12*de + k29 + k31shield

    # idust
    #   + 2 * H2dust * HI * rhoH
    #
    # H2I = (scoef*dtit + H2I) / (1.0 + acoef*dtit)
    #
    #    passes density field


    return

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
            mass = data[fraction_field].value * (data['birth_mass'].value *yt.units.Msun).to('g').value
            abund = convert_abundances.elemental_abundance(element, mass)
            return abund
        return _abundance

    for a in asym:
        fraction_field = ('all','particle_' + a + '_fraction')
        yt.add_field(('all','particle_' + a + '_abundance'),
                      return_function(a, fraction_field), units = "", particle_type = True)

    if (('O' in asym) and ('Mg' in asym) and ('Si' in asym)):
        def _alpha_fraction(field, data):
            alpha = data[('all','particle_O_fraction')] +\
                    data[('all','particle_Mg_fraction')] +\
                    data[('all','particle_Si_fraction')]
            return alpha / 3.0

        def _alpha_abundance(field, data):
            alpha = data[('all','particle_O_abundance')] +\
                    data[('all','particle_Mg_abundance')] +\
                    data[('all','particle_Si_abundance')]

            return alpha / 3.0

#        yt.add_field(('all','particle_alpha_fraction'), function=_alpha_fraction, units = "", particle_type = True)
        yt.add_field(('all','particle_alpha_abundance'), function=_alpha_abundance, units="", particle_type=True)

        if ('S' in asym) and ('Ca' in asym):
            def _alpha_5_fraction(field,data):
                alpha = (data[('all','particle_alpha_fraction')]*3.0 + data[('all','particle_S_fraction')] +\
                                        data[('all','particle_Ca_fraction')]) / 5.0
                return alpha

            def _alpha_5(field,data):
                alpha = data[('all','particle_alpha_abundance')]*3.0 + data[('all','particle_S_abundance')]+\
                                        data[('all','particle_Ca_abundance')]
                return alpha / 5.0

#            yt.add_field( ('all','particle_alpha_5_fraction'), function = _alpha_5_fraction, units = "", particle_type =True)
            yt.add_field( ('all','particle_alpha_5_abundance'), function = _alpha_5, units = "", particle_type = True)


    def _particle_above_chiaki_threshold(field,data):
        C_f  = data[('all','particle_C_fraction')].value
        Fe_f = data[('all','particle_Fe_fraction')].value
        H_f  = data[('all','particle_H_fraction')].value

        return physics.chiaki_threshold(C_f, Fe_f, H_f, return_value=False).astype(np.float64)

    def _particle_chiaki_value(field,data):
        C_f  = data[('all','particle_C_fraction')].value
        Fe_f = data[('all','particle_Fe_fraction')].value
        H_f  = data[('all','particle_H_fraction')].value

        return physics.chiaki_threshold(C_f, Fe_f, H_f, return_value=True)

    yt.add_field( ('all','particle_above_chiaki_threshold'), function = _particle_above_chiaki_threshold, units = "", particle_type = True)
    yt.add_field( ('all','particle_chiaki_value'), function = _particle_chiaki_value, units = "", particle_type = True)

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
            mass1 = ((mass1 * data['birth_mass'].value) * yt.units.Msun).to('g')

            mass2 = data[field2].value
            mass2 = ((mass2 * data['birth_mass'].value) * yt.units.Msun).to('g')
            ratio = convert_abundances.abundance_ratio( (ele1, mass1.value), (ele2, mass2.value), 'mass')

            return ratio * yt.units.g / yt.units.g

        return _abundance_ratio

    nfields = 0
    for r in ratios:
        ele1, ele2 = r.rsplit('/')

        field1 = ('all','particle_' + ele1 + '_fraction')
        field2 = ('all','particle_' + ele2 + '_fraction')

        fieldname = 'particle_' + ele1 + '_over_' + ele2

        yt.add_field(('all', fieldname), function = return_function(ele1,ele2,field1,field2),
                              units = "", particle_type = True)
        nfields = nfields + 1

    def _alpha_return_function(base):
        def _alpha_over_x(field,data):
            alpha = data[('all','particle_alpha_abundance')]
            x     = data[('all','particle_' + base + '_abundance')]

            return convert_abundances.abundance_ratio(('alpha',alpha),(base,x), 'abundances')

        return _alpha_over_x

    def _alpha_5_return_function(base):
        def _alpha_5_over_x(field,data):
            alpha = data[('all','particle_alpha_5_abundance')]
            x     = data[('all','particle_' + base + '_abundance')]

            return convert_abundances.abundance_ratio(('alpha_5',alpha),(base,x), 'abundances')

        return _alpha_5_over_x

    denoms = [x.split('/')[1] for x in ratios]
    denoms = np.unique(denoms)
    if ('all','particle_alpha_abundance') in ds.derived_field_list:
        for x in denoms:
            yt.add_field(('all','particle_alpha_over_' + x), function = _alpha_return_function(x), units = "", particle_type = True)
            yt.add_field(('all','particle_alpha_5_over_' + x), function = _alpha_5_return_function(x), units = "", particle_type = True)


#    def _alpha_over_Fe(field,data):
#        alpha = data[('all','particle_alpha_abundance')]
#        Fe    = data[('all','particle_Fe_abundance')]
#        return convert_abundances.abundance_ratio( ('alpha', alpha), ('Fe', Fe), 'abundances')

#    yt.add_field(('all','particle_alpha_over_Fe'), function = _alpha_over_Fe, units = "", particle_type = True)


#    def _alpha_5_over_Fe(field,data):
#        alpha = data[('all','particle_alpha_5_abundance')]
#        Fe    = data[('all','particle_Fe_abundance')]
#        return convert_abundances.abundance_ratio( ('alpha_5', alpha), ('Fe', Fe), 'abundances')

 #   yt.add_field(('all','particle_alpha_5_over_Fe'), function = _alpha_5_over_Fe, units = "", particle_type = True)


    return nfields

#
# Construct arbitrary abundance ratio fields in yt
# using a function generator
#
def _abundance_ratio_function_generator(ratios, metals, H_mode = 'total'):

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
            mass1 = (mass1 * data['cell_volume']).to('g')

            if ele2 == 'H':
                mass2 = _H_mass(data, H_mode)
            else:
                mass2 = data[field2]

            if ele2 != 'H' and  ele2 != 'He':
                mass2 = mass2.value * data.ds.mass_unit / data.ds.length_unit**3
            mass2 = (mass2 * data['cell_volume']).to('g')

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
                              units = "", sampling_type = 'cell')
        nfields = nfields + 1


    if ('O' in metals) and ('Mg' in metals) and ('Si' in metals):
        def _return_alpha_over_x(element_name):
            def _alpha_over_x(field, data):
                alpha = data[('gas','alpha_Abundance')]
                x     = data[('gas',element_name + '_Abundance')]
                return convert_abundances.abundance_ratio(('alpha',alpha),(element_name,x),'abundances')

            return _alpha_over_x

        denoms = [x.split('/')[1] for x in ratios]
        denoms = np.unique(denoms)
        for x in denoms:
            yt.add_field(('gas','alpha_over_' + x), sampling_type = 'cell', function = _return_alpha_over_x(x), units = "")



    return nfields


def generate_stellar_model_fields(ds):

    if not (ds.parameters['NumberOfParticles'] > 0):
        return

    #
    # luminosity, L_FUV, L_LW, Q0, Q1, E0, E1
    #

    field_names = ['luminosity','L_FUV','L_LW','Q0','Q1','E0','E1', 'Teff','R']

    units = {'luminosity' : yt.units.erg/yt.units.s,
             'L_FUV' : yt.units.erg/yt.units.s,
             'L_LW'  : yt.units.erg/yt.units.s,
             'Q0' : 1.0 /yt.units.s, 'Q1' : 1.0 / yt.units.s,
             'E0' : yt.units.erg, 'E1' : yt.units.erg, 'lifetime' : yt.units.s,
             'Teff' : yt.units.K, 'R' : yt.units.cm}

    unit_label = {'luminosity': 'erg/s', 'L_FUV' : 'erg/s', 'L_LW' : 'erg/s',
                  'Q0' : '1/s', 'Q1' : '1/s', 'E0' : 'erg', 'E1': 'erg', 'lifetime' : 's', 'Teff' : 'K',
                  'R' : 'cm'}

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
            if np.size(data[(field.name[0],'particle_mass')]) == 1:
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
        Q0 = data[(field.name[0],'particle_model_Q0')]
        E0 = data[(field.name[0],'particle_model_E0')]
        return (E0 * Q0).to('erg/s')

    def _model_L1(field, data):
        Q1 = data[(field.name[0],'particle_model_Q1')]
        E1 = data[(field.name[0],'particle_model_E1')]
        return (E1 * Q1).to('erg/s')

    def _age(field, data):
        p = data[(field.name[0],'creation_time')]
        t = data.ds.current_time
        return (t - p).to('Myr')

    def _model_L_1_3eV(field, data):
        Teff = data[(field.name[0],'particle_model_Teff')].value
        flux = np.array([radiation.BB_flux(1.0, 3.0, T) for T in Teff])

        flux = flux * yt.units.erg / yt.units.s / yt.units.cm**2

        SA = 4.0 * np.pi * data[(field.name[0],'particle_model_R')]**2

        return flux * SA

    for field in field_names:
        yt.add_field(('all', 'particle_model_' + field),
                     function = _function_generator(field), units=unit_label[field],
                	     particle_type = True)

    def _lifetime(field, data):
        ptname = field.name[0]
        m = data[(ptname,'birth_mass')].value
        z = data[(ptname,'metallicity_fraction')].value
        ispopiii = data[(ptname,'particle_is_popiii')]

        if np.size(m) == 1:  # get around yt's checking
            lt = np.shape(m)
        else:
            lt = np.zeros(np.size(m))
            for i in np.arange(np.size(m)):
                if m[i] < 0 and z[i] < 0:
                    lt[i] = 0.0
                elif ispopiii[i]:
                    lt[i] = (physics.popIII_lifetime(m[i]) * yt.units.yr).to('Myr').value
                else:
                    lt[i] = SE_table.interpolate({'mass' : m[i], 'metallicity' : z[i]}, 'lifetime')
                    lt[i] = (lt[i] * yt.units.s).to('Myr').value


        return lt * yt.units.Myr


    def _is_popIII(field,data):

        pt = data[(field.name[0],'particle_type')]

        am_i_popiii = ((pt == 14)).astype(np.float64)

        if data.ds.parameters['IndividualStarPopIIIFormation'] > 0:

            if data.ds.parameters['PopIIIMetalCriticalFraction'] > 0:
                am_i_popiii = am_i_popiii + (data[(field.name[0],'metallicity_fraction')] > data.ds.parameters['PopIIIMetalCriticalFraction']).astype(np.float64)
            else:
                am_i_popiii = am_i_popiii + np.logical_not(data[(field.name[0],'particle_above_chiaki_threshold')])


        return am_i_popiii.astype(np.float64)

    yt.add_field(('all','particle_is_popiii'), function = _is_popIII, units ='', particle_type=True)

    yt.add_field(('all','particle_model_lifetime'), function = _lifetime, units = 'Myr',
                 particle_type = True)

    yt.add_field(('all', 'particle_age'), function = _age, units = 'Myr',
                 particle_type = True)

    yt.add_field(('all','particle_model_L0'), function = _model_L0, units = 'erg/s',
                 particle_type = True)

    yt.add_field(('all','particle_model_L1'), function = _model_L1, units = 'erg/s',
                 particle_type = True)

    yt.add_field(('all','particle_model_L_1_3eV'), function = _model_L_1_3eV, units = 'erg/s',
                 particle_type = True)

    return


def _grackle_fields(ds):
    """
    Fields that require use of pygrackle
    """

    cdata                       = pygrackle.chemistry_data()
    cdata.use_grackle           = 1

    enzo_to_grackle             = { 'MultiSpecies' : 'primordial_chemistry',
                                    'MetalCooling' : 'metal_cooling',
                                    'self_shielding_method' : 'self_shielding_method',
#
#                                   The below is broken in pygrackle - not sure of prob
#                                    'H2_self_shielding' : 'H2_self_shielding',
                                    'grackle_data_file' : 'grackle_data_file',
                                    'DensityUnits' : 'density_units',
                                    'LengthUnits' : 'length_units',
                                    'TimeUnits'   : 'time_units',
                                    'ComovingCoordinates': 'comoving_coordinates',
                                    'with_radiative_cooling' : 'with_radiative_cooling',
                                    'UVbackground' : 'UVbackground'}

    for k in enzo_to_grackle:
        setattr(cdata, enzo_to_grackle[k], ds.parameters[k])

    cdata.a_units = 1.0
    cdata.a_value = 1.0
    cdata.velocity_units = cdata.length_units / cdata.time_units
    #cdata.energy_units = (cdata.length_units / cdata.time_units)**2.0
    cdata.initialize()

    def _H2_self_shielding_length(field, data):
        return data['dx'].to('cm')
    ds.add_field(('gas','H2_self_shielding_length'), function = _H2_self_shielding_length, units='cm')

    def _cooling_time(field, data):

        #
        # This checks if yt is doing its fake-data error checking and
        # gives dummy result... a bit of a hack....
        #
        # compute field in grackle grid-by-grid
        field_list      = data.ds.derived_field_list + data.ds.field_list
        flat_fields = {}

        fc = pygrackle.FluidContainer(cdata,10) # dummy container for now
        for f1, f2, conv in pygrackle.fluid_container._needed_fields(fc):
            if f2 not in field_list:
                raise pygrackle.fluid_container.FieldNotFound(f2)
            else:
                flat_fields[f2[1]] = np.zeros(np.size( data['Density']))

        for f1, f2,conv in pygrackle.fluid_container._needed_fields(fc):
            flat_fields[f2[1]] = ((1.0*data[f2]).value).flatten() / conv

        flat_fields[f2[1]] = np.zeros(np.size( data['Density'].flatten()))
        flat_fields['cooling_time'] = np.zeros(np.size(flat_fields[f2[1]]))

        # compute a new FC every 8192 zones
        imin = 0
        imax = 0
        di   = 8192
        ncells = np.size(flat_fields['cooling_time'])
        while imax < ncells:
            imin = 1*imax
            imax = np.min( [imax + di, ncells] )
            fc   = pygrackle.FluidContainer(cdata, imax - imin)

            for f1, f2, conv in pygrackle.fluid_container._needed_fields(fc):
                fc[f1][:] = flat_fields[f2[1]][imin:imax]
            fc.calculate_cooling_time()


            flat_fields['cooling_time'][imin:imax] = fc['cooling_time']

        flat_fields['cooling_time'] = flat_fields['cooling_time'] * cdata.time_units

        return flat_fields['cooling_time'].reshape( np.shape(data['Density'].value)) * yt.units.s

    def _neg_cooling_time(field,data):
        return -1.0 * data['cooling_time']

    ds.add_field(('gas','cooling_time'), function = _cooling_time, units = 's')
    ds.add_field(('gas','neg_cooling_time'), function = _neg_cooling_time, units = 's')

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

        return pe.to('erg/s/cm**3')

    def _pe_heating_rate_masked(field, data):
        pe = data[('gas','Pe_heating_rate')].to('erg/s/cm**3')

        x = 1.0 * pe

        x[data['temperature'] > data.ds.parameters['IndividualStarFUVTemperatureCutoff']] = 0.0

        return x

    def _otlwcgs(field, data):
        if ('enzo','OTLW_kdissH2I') in data.ds.field_list:
            lw = data[('enzo','OTLW_kdissH2I')].value / data.ds.time_unit
        else:
            lw = np.zeros(np.shape(data['Density'])) / data.ds.time_unit

        return lw.to('1/s')

    def _G_o(field,data):
        pe  = data[('gas','Pe_heating_rate')].to('erg/s/cm**3').value
        Z   = (data['Metal_Density'] / data['Density']).value
        n_H = (data['H_p0_number_density'] + data['H_p1_number_density'] + data['H_m1_number_density'] +\
                0.5*(data['H2_p0_number_density'] + data['H2_p1_number_density'])).to('cm**(-3)').value

        logZ   = np.log10(Z / 0.014)
        g_to_d = np.zeros(np.shape(logZ))
        g_to_d[logZ <= -0.73] = 0.68 - 3.08*logZ[logZ <= -0.73]
        g_to_d[logZ  > -0.73] = 2.21 - 1.00*logZ[logZ  > -0.73]
        d_to_g = 1.0 / (10.0**(g_to_d))
        D = d_to_g / 6.616595E-3
        epsilon = 0.01488637246 * (n_H)**(0.235269059)
        atten = np.exp( - 1.33E-21 * D * data['dx'].to('cm').value * n_H)
        G_o = pe / (1.3E-24 * n_H * epsilon * D * atten)
        return G_o * (data['Density'] / data['Density'])

    def _G_eff(field,data):
        pe  = data[('gas','Pe_heating_rate')].to('erg/s/cm**3').value
        Z   = (data['Metal_Density'] / data['Density']).value
        n_H = (data['H_p0_number_density'] + data['H_p1_number_density'] + data['H_m1_number_density'] +\
                0.5*(data['H2_p0_number_density'] + data['H2_p1_number_density'])).to('cm**(-3)').value

        logZ   = np.log10(Z / 0.014)
        g_to_d = np.zeros(np.shape(logZ))
        g_to_d[ logZ <= -0.73] = 0.68 - 3.08*logZ[logZ <= -0.73]
        g_to_d[ logZ  > -0.73] = 2.21 - 1.00*logZ[logZ  > -0.73]
        d_to_g = 1.0 / (10.0**(g_to_d))
        D = d_to_g / 6.616595E-3

        epsilon = 0.01488637246 * (n_H)**(0.235269059)

        # atten = np.exp( - 1.33E-21 * D * data['dx'].to('cm').value * n_H)

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

        if ('enzo','OTLW_kdissH2I') in data.ds.field_list:
            kdissH2I = (data[('enzo','OTLW_kdissH2I')].value / data.ds.time_unit).to('1/s')
        else:
            kdissH2I = (np.zeros(np.shape(data['Density'])) / data.ds.time_unit).to('1/s')

        LW_flux = kdissH2I / H2Isigma * LW_energy

        return LW_flux.to('erg/cm**2/s')

    def _Q0_flux(field, data):
        E_HI = 13.6 * yt.units.eV
        kph = data[('enzo','HI_kph')]
        kph.convert_to_cgs()
        n   = data[('gas','H_p0_number_density')]
        n.convert_to_cgs()
        dt  = data.ds.parameters['dtPhoton']
        V   = data['cell_volume']
        V.convert_to_cgs()
        dx  = data['dx']
        dx.convert_to_cgs()
        s   = 6.34629E-18 * yt.units.cm**(2) # cross section of HI at 13.6 eV

        tau   = s * n * dx
        denom = 1.0 - np.exp(-tau)

        Q = kph * n * V / denom  # this gives number of photons / s

        flux = Q * E_HI / dx**2

        return flux.to('erg/cm**2/s')

    def _Q1_flux(ds,data):
        E_HeI = 24.6 * yt.units.eV
        kph = data[('enzo','HeI_kph')]
        kph.convert_to_cgs()
        n   = data[('gas','H_p0_number_density')]
        n.convert_to_cgs()
        dt  = data.ds.parameters['dtPhoton']
        V   = data['cell_volume']
        V.convert_to_cgs()
        dx  = data['dx']
        dx.convert_to_cgs()
        s   = 7.4300459E-18 * yt.units.cm**(2) # cross section of HeI at 24.6 eV

        tau   = s * n * dx
        denom = 1.0 - np.exp(-tau)

        Q = kph * n * V / denom  # this gives number of photons / s

        flux = Q * E_HeI / dx**2

        return flux.to('erg/cm**2/s')


    def _metal_total_mass(field, data):
        mass = data['Metal_Density'] * data['cell_volume']

        return mass.to('g')

    def _grav_pot(field,data):
        try:
            x = (data['PotentialField'] * -1.0).to('erg/g')
        except:
            try:
                x = ( (data['GravPotential'].value * data.ds.velocity_unit**2)
                   * -1.0).to('erg/g')
            except:
                x = ( (data['Grav_Potential'].value * data.ds.velocity_unit**2)
                   * -1.0).to('erg/g')


        return x

    def _tot_grav_pot(field,data):
        try:
            x = data['PotentialField']
        except:
            try:
                x = data['GravPotential'].value * data.ds.velocity_unit**2
            except:
                x = data['Grav_Potential'].value * data.ds.velocity_unit**2

        x = x + data[('index','DM_background_potential')]

        return x.to('erg/g')

    def _gas_grav_pot(field,data):
        try:
            x = data['PotentialField']
        except:
            try:
                x = data['GravPotential'].value * data.ds.velocity_unit**2
            except:
                x = data['Grav_Potential'].value * data.ds.velocity_unit**2

        return x.to('erg/g')


    def _pos_tot_grav_pot(field, data):

        return np.abs(data[('gas','total_gravitational_potential')])

    def _potential_energy(field,data):

        x = data[('gas','total_gravitational_potential')] * data['cell_mass']

        return x.to('erg')

    def _grav_bound(field, data):
        PE = data[('gas','potential_energy')].to('erg')
        TE = ( data[('gas','thermal_energy')] * data['cell_mass'].to('g')).to('erg')
        KE = ( data[('gas','kinetic_energy')] * data['cell_volume']).to('erg')

        result = 1 * ((TE + KE) + PE < 0.0)

        return result*1.0


    def _mag_cyl_r(field,data):
        return np.abs( data[('index','cylindrical_radius')].to('cm'))

    def _mag_cyl_z(field,data):
        return np.abs( data[('index','cylindrical_z')].to('cm') )

    def _dm_density(field, data):
        r     = data[('index','spherical_r')].to('cm')
        r_s   = (data.ds.parameters['DiskGravityDarkMatterR'] * yt.units.Mpc).to('cm')
        rho_o = (data.ds.parameters['DiskGravityDarkMatterDensity'] * yt.units.g / yt.units.cm**3)

        rho = dm_halo.burkert_density(r, r_s, rho_o)

        rho.convert_to_cgs()
        return rho

    def _dm_potential(field, data):
        r     = data[('index','spherical_r')].to('cm')
        r_s   = (data.ds.parameters['DiskGravityDarkMatterR'] * yt.units.Mpc).to('cm')
        rho_o = (data.ds.parameters['DiskGravityDarkMatterDensity'] * yt.units.g / yt.units.cm**3)

        pot = dm_halo.burkert_potential(r, r_s, rho_o)

        pot.convert_to_cgs()
        return pot

    def _rad_accel(field, data):
        return np.sqrt(data['RadAccel1']**2 + data['RadAccel2']**2 + data['RadAccel3']**2).to('cm/s**2')

    def _is_star_forming(field, data):
        n      = data[('gas','number_density')].to('cm**(-3)')
        T      = data['Temperature']
        divv   = data[('gas','velocity_divergence')]
        l      = data['grid_level']

        answer = 1 * ((n > data.ds.parameters['StarMakerOverDensityThreshold']) *\
                       (T < data.ds.parameters['IndividualStarTemperatureThreshold']) *\
                        (divv < 0) *\
                          (l == data.ds.parameters['MaximumRefinementLevel']))

        return answer

    def _above_chiaki_threshold(field,data):
        C_f  = data[('gas','C_Fraction')].value
        Fe_f = data[('gas','Fe_Fraction')].value
        H_f  = data[('gas','H_fraction')].value

        return physics.chiaki_threshold(C_f, Fe_f, H_f, return_value=False).astype(np.float64)

    def _chiaki_value(field,data):
        C_f  = data[('gas','C_Fraction')].value
        Fe_f = data[('gas','Fe_Fraction')].value
        H_f  = data[('gas','H_fraction')].value

        return physics.chiaki_threshold(C_f, Fe_f, H_f, return_value=True)


    yt.add_field(("gas","a_rad"), sampling_type = 'cell',function=_rad_accel, units="cm/s**2")

    yt.add_field(('index','DM_background_density'), sampling_type = 'cell',function = _dm_density, units = 'g/cm**3')
    yt.add_field(('index','DM_background_potential'), sampling_type = 'cell',function = _dm_potential, units = 'erg/g')

    yt.add_field(('index','magnitude_cylindrical_radius'), sampling_type = 'cell',function = _mag_cyl_r, units = 'cm')

    yt.add_field(('index','magnitude_cylindrical_z'),sampling_type = 'cell', function = _mag_cyl_z, units = 'cm')
#    def _H2_total_mass(field, data):
#        mass = data[('gas',

    yt.add_field(('gas','Pe_heating_rate'),sampling_type = 'cell', function = _pe_heating_cgs, units = 'erg/s/cm**3')
    yt.add_field(('gas','H_total_mass'), sampling_type = 'cell',function = _H_total_mass, units ='g')
    yt.add_field(('gas','H_Mass'),sampling_type = 'cell', function = _H_total_mass, units = 'g') # define as same
    yt.add_field(('gas','He_total_mass'), sampling_type = 'cell',function = _He_total_mass, units = 'g')
    yt.add_field(('gas','metal_mass'), sampling_type = 'cell',function = _metal_total_mass, units = 'g')

    yt.add_field(('gas','OTLW_kdissH2I'), sampling_type = 'cell',function = _otlwcgs, units = '1/s',
                 validators=ValidateDataField(('enzo','OTLW_kdissH2I')))
    yt.add_field(('gas','LW_flux'), sampling_type = 'cell',function = _LW_flux, units = "erg/s/cm**2",
                 validators=ValidateDataField(('enzo','OTLW_kdissH2I')))

    yt.add_field(('gas','above_chiaki_threshold'), sampling_type = 'cell',function = _above_chiaki_threshold,
                 units="")
    yt.add_field(('gas','chiaki_value'),sampling_type = 'cell', function = _chiaki_value,
                 units="")

    yt.add_field(('gas','is_star_forming'),sampling_type = 'cell', function = _is_star_forming,
                         units = "")

    yt.add_field(('gas','Pe_heating_rate_masked'),sampling_type = 'cell', function = _pe_heating_rate_masked, units='erg/s/cm**3')
    yt.add_field(('gas','G_o'), sampling_type = 'cell',function = _G_o, units = "")
    yt.add_field(('gas','G_eff'), sampling_type = 'cell',function = _G_eff, units = "")
    yt.add_field(('gas','FUV_flux'), sampling_type = 'cell',function = _FUV_flux, units = "erg/s/cm**2")
    yt.add_field(('gas','Q0_flux'),sampling_type = 'cell', function = _Q0_flux, units = "erg/s/cm**2")
    yt.add_field(('gas','Q1_flux'),sampling_type = 'cell', function = _Q1_flux, units = "erg/s/cm**2")
#    yt.add_field(('gas','H2_total_mass'), function = _H2_total_mass, units = 'g')
#    yt.add_field(('gas','All_H_total_mass'), function = _all_H_total_mass, units = 'g')

    if (('enzo','PotentialField') in fields) or (('enzo', 'GravPotential') in fields) or (('enzo','Grav_Potential') in fields):
        yt.add_field(('gas','pos_gravitational_potential'), sampling_type = 'cell',function=_grav_pot, units = 'erg/g')
        yt.add_field(('gas','gas_gravitational_potential'), sampling_type = 'cell',function=_gas_grav_pot, units = 'erg/g')
        yt.add_field(('gas','total_gravitational_potential'),sampling_type = 'cell', function=_tot_grav_pot, units = 'erg/g')
        yt.add_field(('gas','pos_total_gravitational_potential'), sampling_type = 'cell',function=_pos_tot_grav_pot, units = 'erg/g')
        yt.add_field(('gas','potential_energy'),sampling_type = 'cell', function=_potential_energy, units = 'erg')
        yt.add_field(('gas','gravitationally_bound'), sampling_type = 'cell',function=_grav_bound, units = "")

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
    print("defining for the following metals ", metals)
    # make new functions to do correct units for species fields
    _density_function_generator(metals + ['Metal'])

    print("tracer species present: ", metals)
    nfields = _mass_function_generator(metals)
    print(nfields, "mass fields defined")
    nfields = _mass_fraction_function_generator(ds, metals)
    print(nfields, "mass fraction fields defined")
    nfields = _number_density_function_generator(metals)
    print(nfields, "number density fields defined")

    if not (ionization._ion_table is None):
        nfields = _ionization_state_generator(metals)
        print(nfields, "ionization state fields defined")

    nfields = _abundance_ratio_function_generator(ratios, metals, H_mode = 'total')
    print(nfields, "abundance ratio fields defined")
    nfields = _abundance_function_generator(metals)

    if ds.parameters['NumberOfParticles'] > 0:
        if ('all','particle_' + metals[0] + '_fraction') in ds.field_list:

            nfields =  _particle_abundance_ratio_function_generator(ratios, ds)
            print(nfields, "particle abundance ratio fields defined")
            _particle_abundance_function_generator(metals, ds)


    generate_stellar_model_fields(ds)


    nfields = _additional_helper_fields(fields)
    print(nfields, "additional helper fields defined")

    #generate_grackle_fields(ds)

    FIELDS_DEFINED = True
    return



def load_and_define(name):
    """
    Wrapper around yt to load a data set and define gradient
    fields and particle filters which must be defined for each
    simulation file separately (unlike the above)
    """


    ds = yt.load(name)

    if not FIELDS_DEFINED:
        generate_derived_fields(ds)
        ds = yt.load(name)
        generate_derived_fields(ds)

    gradient_available = generate_gradient_fields(ds)

    if gradient_available:
        def _grav_accel_x(field,data):
            return data[('gas','gas_gravitational_potential_gradient_x')].to('cm/s**2')
        def _grav_accel_y(field,data):
            return data[('gas','gas_gravitational_potential_gradient_y')].to('cm/s**2')
        def _grav_accel_z(field,data):
            return data[('gas','gas_gravitational_potential_gradient_z')].to('cm/s**2')
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

    #generate_particle_filters(ds)
    #generate_grackle_fields(ds)

    return ds


def load(name):
    return load_and_define(name)

def generate_grackle_fields(ds):

    if not GRACKLE_IMPORTED:
        print("Grackle's python wrapper (pygrackle) was not imported successfully")

    if ds.parameters['use_grackle']:
        _grackle_fields(ds)

    return

def generate_gradient_fields(ds):
    """
    generate gas self gravity gradient fields and rename them to
    something sensible
    """

    if ("gas","gas_gravitational_potential") in ds.derived_field_list:
        ds.add_gradient_fields(("gas","gas_gravitational_potential"))
        gradient_available = True
    else:
        gradient_available = False

    return gradient_available

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
    def all_stars(pfilter, data):
        filter = data[(pfilter.filtered_type, "particle_type")] >= 11

        return filter

    @yt.particle_filter(requires=["particle_type"], filtered_type='all_stars')
    def all_popIII_stars(pfilter, data):
        filter = data[(pfilter.filtered_type,"particle_is_popiii")].astype(np.bool)

        return filter

    @yt.particle_filter(requires=["particle_type"], filtered_type='all_stars')
    def all_popII_stars(pfilter, data):
        filter = np.logical_not(data[(pfilter.filtered_type,"particle_is_popiii")])

        return filter

    @yt.particle_filter(requires=["particle_type"], filtered_type='all_stars')
    def main_sequence_stars(pfilter, data):
        filter = (data[(pfilter.filtered_type, "particle_type")] == 11) +\
                 (data[(pfilter.filtered_type, "particle_type")] == 15)

        return filter

    @yt.particle_filter(requires=["particle_type"], filtered_type='all_stars')
    def main_sequence_popIII_stars(pfilter, data):
        filter = data[(pfilter.filtered_type, "particle_type")] == 14
        return filter

    @yt.particle_filter(requires=["particle_type"], filtered_type='all_stars')
    def remnant_stars(pfilter, data):
        filter = data[(pfilter.filtered_type, "particle_type")] == 13
        return filter

    @yt.particle_filter(requires=["particle_type",'birth_mass'], filtered_type='all_stars')
    def low_mass_stars(pfilter, data):
        filter = data[(pfilter.filtered_type, "particle_type")] == 11
        filter = filter * (data[(pfilter.filtered_type,"birth_mass")] > 2.0) * (data[(pfilter.filtered_type,"birth_mass")] < 8.0)
        return filter

    @yt.particle_filter(requires=["particle_type",'birth_mass'], filtered_type='all_stars')
    def low_mass_unresolved_stars(pfilter, data):
        filter = data[(pfilter.filtered_type, "particle_type")] == 15
        return filter

    @yt.particle_filter(requires=["particle_type",'birth_mass'], filtered_type='all_stars')
    def white_dwarf(pfilter, data):
        filter = data[(pfilter.filtered_type, "particle_type")] == 12
        return filter



    ds.add_particle_filter('main_sequence_stars')
    ds.add_particle_filter('remnant_stars')
    ds.add_particle_filter('low_mass_stars')
    ds.add_particle_filter('low_mass_unresolved_stars')
    ds.add_particle_filter("white_dwarf")
    ds.add_particle_filter('main_sequence_popIII_stars')

    ds.add_particle_filter("all_popII_stars")
    ds.add_particle_filter("all_popIII_stars")
    ds.add_particle_filter("all_stars")

    #
    #
    # End of life filteres for non-snia
    #
    #

    @yt.particle_filter(requires=['particle_type','birth_mass'], filtered_type='all_stars')
    def all_remnants(pfilter, data):
        filter = data[(pfilter.filtered_type,"particle_type")] == 13
        return filter

    @yt.particle_filter(requires=['particle_type','birth_mass'], filtered_type='all_remnants')
    def popIII_remnant(pfilter, data):

        if ('IndividualStarPopIIIFormation' in data.ds.parameters) and\
           ('PopIIIMetalCriticalFraction' in data.ds.parameters):
            if data.ds.parameters['IndividualStarPopIIIFormation'] > 0:

                if data.ds.parameters['PopIIIMetalCriticalFraction'] > 0:
                    filter = data[(pfilter.filtered_type,'metallicity_fraction')] < data.ds.parameters['PopIIIMetalCriticalFraction']
                else: # use the Chiaki threshold:

                    filter = np.logical_not( data[(pfilter.filtered_type,'particle_above_chiaki_threshold')] )

        else:
            filter = np.logical_not(data[(pfilter.filtered_type, "birth_mass")] == data[(pfilter.filtered_type, "birth_mass")])


        return filter

    @yt.particle_filter(requires=['particle_type','birth_mass'], filtered_type='all_remnants')
    def popIII_ccsne_remnant(pfilter, data):

        if ('IndividualStarPopIIIFormation' in data.ds.parameters) and\
           ('PopIIIMetalCriticalFraction' in data.ds.parameters):
            if data.ds.parameters['IndividualStarPopIIIFormation'] > 0:

                if data.ds.parameters['PopIIIMetalCriticalFraction'] > 0:
                    filter = data[(pfilter.filtered_type,'metallicity_fraction')] < data.ds.parameters['PopIIIMetalCriticalFraction']
                else: # use the Chiaki threshold:

                    filter = np.logical_not( data[(pfilter.filtered_type,'particle_above_chiaki_threshold')] )

                filter = filter * ((data[(pfilter.filtered_type,'birth_mass')] >= ds.parameters['TypeIILowerMass']) *\
                               (data[(pfilter.filtered_type,'birth_mass')] <= ds.parameters['TypeIIUpperMass']))
        else:
            filter = np.logical_not(data[(pfilter.filtered_type, "birth_mass")] == data[(pfilter.filtered_type, "birth_mass")])


        return filter


    @yt.particle_filter(requires=['particle_type','birth_mass'], filtered_type='all_remnants')
    def popIII_pisne_remnant(pfilter, data):

        if ('IndividualStarPopIIIFormation' in data.ds.parameters) and\
           ('PopIIIMetalCriticalFraction' in data.ds.parameters):
            if data.ds.parameters['IndividualStarPopIIIFormation'] > 0:

                if data.ds.parameters['PopIIIMetalCriticalFraction'] > 0:
                    filter = data[(pfilter.filtered_type,'metallicity_fraction')] < data.ds.parameters['PopIIIMetalCriticalFraction']
                else: # use the Chiaki threshold:

                    filter = np.logical_not( data[(pfilter.filtered_type,'particle_above_chiaki_threshold')] )

                filter = filter * ((data[(pfilter.filtered_type,'birth_mass')] >= ds.parameters['PISNLowerMass']) *\
                               (data[(pfilter.filtered_type,'birth_mass')] <= ds.parameters['PISNUpperMass']))
        else:
            filter = np.logical_not(data[(pfilter.filtered_type, "birth_mass")] == data[(pfilter.filtered_type, "birth_mass")])


        return filter

    @yt.particle_filter(requires=['particle_type','birth_mass'], filtered_type='all_remnants')
    def popIII_direct_collapse_remnant(pfilter, data):

        if ('IndividualStarPopIIIFormation' in data.ds.parameters) and\
           ('PopIIIMetalCriticalFraction' in data.ds.parameters):
            if data.ds.parameters['IndividualStarPopIIIFormation'] > 0:

                if data.ds.parameters['PopIIIMetalCriticalFraction'] > 0:
                    filter = data[(pfilter.filtered_type,'metallicity_fraction')] < data.ds.parameters['PopIIIMetalCriticalFraction']
                else: # use the Chiaki threshold:

                    filter = np.logical_not( data[(pfilter.filtered_type,'particle_above_chiaki_threshold')] )

                filter = filter *  (np.logical_not((data[(pfilter.filtered_type,'birth_mass')] >= ds.parameters['PISNLowerMass']) *\
                                                  (data[(pfilter.filtered_type,'birth_mass')] <= ds.parameters['PISNUpperMass'])) *\
                                    np.logical_not((data[(pfilter.filtered_type,'birth_mass')] >= ds.parameters['TypeIILowerMass']) *\
                                                  (data[(pfilter.filtered_type,'birth_mass')] <= ds.parameters['TypeIIUpperMass'])))
        else:
            filter = np.logical_not(data[(pfilter.filtered_type, "birth_mass")] == data[(pfilter.filtered_type, "birth_mass")])


        return filter

    @yt.particle_filter(requires=["particle_type",'birth_mass'], filtered_type='all_remnants')
    def ccsne_remnant(pfilter, data):

        filter = ((data[(pfilter.filtered_type, "birth_mass")] <= data.ds.parameters['IndividualStarDirectCollapseThreshold']) *\
                  (data[(pfilter.filtered_type, "birth_mass")] >= data.ds.parameters['IndividualStarAGBThreshold']))

        if ('IndividualStarPopIIIFormation' in data.ds.parameters) and\
           ('PopIIIMetalCriticalFraction' in data.ds.parameters):
            if data.ds.parameters['IndividualStarPopIIIFormation'] > 0:

                if data.ds.parameters['PopIIIMetalCriticalFraction'] > 0:
                    filter = filter * data[(pfilter.filtered_type,'metallicity_fraction')] < data.ds.parameters['PopIIIMetalCriticalFraction']
                else: # use the Chiaki threshold:
                    filter = filter * data[(pfilter.filtered_type,'particle_above_chiaki_threshold')].astype(np.bool)

        return filter

    @yt.particle_filter(requires=["particle_type",'birth_mass'], filtered_type='all_remnants')
    def direct_collapse_remnant(pfilter, data):

        filter = data[(pfilter.filtered_type, "birth_mass")] > data.ds.parameters['IndividualStarDirectCollapseThreshold']

        if ('IndividualStarPopIIIFormation' in data.ds.parameters) and\
           ('PopIIIMetalCriticalFraction' in data.ds.parameters):
            if data.ds.parameters['IndividualStarPopIIIFormation'] > 0:

                if data.ds.parameters['PopIIIMetalCriticalFraction'] > 0:
                    filter = filter * data[(pfilter.filtered_type,'metallicity_fraction')] < data.ds.parameters['PopIIIMetalCriticalFraction']
                else: # use the Chiaki threshold:
                    filter = filter * data[(pfilter.filtered_type,'particle_above_chiaki_threshold')].astype(np.bool)

        return filter


    ds.add_particle_filter('all_remnants')
    ds.add_particle_filter('popIII_ccsne_remnant')
    ds.add_particle_filter('popIII_pisne_remnant')
    ds.add_particle_filter('ccsne_remnant')
    ds.add_particle_filter('direct_collapse_remnant')
    ds.add_particle_filter('popIII_direct_collapse_remnant')


    #
    # this is the same as selecting for WD's
    #
    #@yt.particle_filter(requires=["particle_type",'birth_mass'], filtered_type='all_stars')
    #def agb_remnant(pfilter, data):
    #    filter = data[(pfilter.filtered_type, "particle_type")] == 12
    #
    #    filter = filter * ( (data[(pfilter.filtered_type,'birth_mass')] >= data.ds.parameters['IndividualStarWDMinimumMass']) *
    #                        (data[(pfilter.filtered_type,'birth_mass')] <= data.ds.parameters['IndividualStarWDMaximumMass']))
    #    return filter
    #ds.add_particle_filter("agb_remnant")


    #
    #
    # SNIa particle filters
    #
    #

    @yt.particle_filter(requires=["particle_type",'birth_mass','snia_sch_metal_fraction','snia_sds_metal_fraction',
                                  'snia_hers_metal_fraction',
                                  'snia_metal_fraction'], filtered_type='all_stars')
    def snia_progenitor(pfilter,data):
        filter = ((data[(pfilter.filtered_type, "birth_mass")] >= data.ds.parameters['IndividualStarSNIaMinimumMass']) *\
                  (data[(pfilter.filtered_type, "birth_mass")] <= data.ds.parameters['IndividualStarSNIaMaximumMass']))

        filter = filter * ( (data[(pfilter.filtered_type, 'snia_sch_metal_fraction')] < 0) +\
                            (data[(pfilter.filtered_type, 'snia_sds_metal_fraction')] < 0) +\
                            (data[(pfilter.filtered_type, 'snia_hers_metal_fraction')] < 0) +\
                            (data[(pfilter.filtered_type, 'snia_metal_fraction')] < 0) )
        return filter

    @yt.particle_filter(requires=["particle_type",'birth_mass','snia_sch_metal_fraction','snia_sds_metal_fraction',
                                  'snia_hers_metal_fraction',
                                  'snia_metal_fraction'], filtered_type='all_stars')
    def snia_dds_progenitor(pfilter,data):
        filter = ((data[(pfilter.filtered_type, "birth_mass")] >= data.ds.parameters['IndividualStarSNIaMinimumMass']) *\
                   (data[(pfilter.filtered_type, "birth_mass")] <= data.ds.parameters['IndividualStarSNIaMaximumMass']))

        filter = filter * (data[(pfilter.filtered_type, 'snia_metal_fraction')] < 0)

        return filter

    @yt.particle_filter(requires=["particle_type",'birth_mass','snia_sch_metal_fraction','snia_sds_metal_fraction',
                                  'snia_hers_metal_fraction',
                                  'snia_metal_fraction'], filtered_type='all_stars')
    def snia_sch_progenitor(pfilter,data):
        filter = ((data[(pfilter.filtered_type, "birth_mass")] >= data.ds.parameters['IndividualStarSNIaMinimumMass']) *\
                  (data[(pfilter.filtered_type, "birth_mass")] <= data.ds.parameters['IndividualStarSNIaMaximumMass']))

        filter = filter * (data[(pfilter.filtered_type, 'snia_sch_metal_fraction')] < 0)

        return filter

    @yt.particle_filter(requires=["particle_type",'birth_mass','snia_sch_metal_fraction','snia_sds_metal_fraction',
                                  'snia_hers_metal_fraction',
                                  'snia_metal_fraction'], filtered_type='all_stars')
    def snia_hers_progenitor(pfilter,data):
        filter = ((data[(pfilter.filtered_type, "birth_mass")] >= data.ds.parameters['IndividualStarSNIaMinimumMass']) *\
                  (data[(pfilter.filtered_type, "birth_mass")] <= data.ds.parameters['IndividualStarSNIaMaximumMass']))

        filter = filter * (data[(pfilter.filtered_type, 'snia_hers_metal_fraction')] < 0)

        return filter

    @yt.particle_filter(requires=["particle_type",'birth_mass','snia_sch_metal_fraction','snia_sds_metal_fraction',
                                  'snia_hers_metal_fraction',
                                  'snia_metal_fraction'], filtered_type='all_stars')
    def snia_sds_progenitor(pfilter,data):
        filter = ((data[(pfilter.filtered_type, "birth_mass")] >= data.ds.parameters['IndividualStarSNIaMinimumMass']) *\
                  (data[(pfilter.filtered_type, "birth_mass")] <= data.ds.parameters['IndividualStarSNIaMaximumMass']))

        filter = filter * (data[(pfilter.filtered_type, 'snia_sds_metal_fraction')] < 0)

        return filter

    #
    # And select those that have actually exploded by
    # which ones have zero mass
    #
    if ('all','snia_sch_metal_fraction') in ds.field_list:
        ds.add_particle_filter("snia_progenitor")
        ds.add_particle_filter("snia_dds_progenitor")
        ds.add_particle_filter("snia_hers_progenitor")
        ds.add_particle_filter("snia_sds_progenitor")
        ds.add_particle_filter("snia_sch_progenitor")



    return
