import yt
import numpy as np

def load_and_define_fields(dsname):
    """
    Works just like yt.load, but also defines fields for you. Returns the 
    dataset object that would normally be returned by yt.load. Example:

       ds = load_and_define_fields("DD0000/DD0000")

    Defines 7 different fields:

       ('gas','Pe_heating_rate') : Local photoelectric heating rate (erg/s/cm**3)
                                   uncorrected for temperature. Above a threshold T,
                                   in the simulation, this field is zeroed. This
                                   is usually at about 10^4 K.
       ('gas','Pe_heating_rate_masked') : The correct heating rate as applied in sim
                                          corrected for temperature (erg/s/cm**3)
       ('gas','G_o')      : ISRF normalized to the Milky Way value (unitless)
       ('gas','FUV_flux') : FUV band flux (in erg/s/cm**2)
       ('gas','LW_flux')  : LW band flux (in erg/s/cm**2)
       ('gas','Q0_flux')  : HI ionizing photon flux (in erg/s/cm**2)
       ('gas','Q1_flux')  : HeI ionizing photon flux (in erg/s/cm**2)

    Parameters:
        dsname : string pointing to filepath of data set
    Returns:
        ds : yt data set object
    """

    def _pe_heating_rate(field, data):
        """
        Correct units for PE heating rate from dataset
        """
        eu = data.ds.mass_unit * data.ds.velocity_unit**2
        lu = data.ds.length_unit
        tu = data.ds.time_unit

        pe = data[('enzo','Pe_heating_rate')] * eu / lu**3 / tu
        pe = pe.convert_to_units('erg/s/cm**3')

        return pe

    def _pe_heating_rate_masked(field,data):
        """
        Mask PE heating rate for temperatures above threshold, like it
        is actually done in code.
        """

        x = data[('gas','Pe_heating_rate')]
        x[data['temperature'] > data.ds.parameters['IndividualStarFUVTemperatureCutoff']] = 0.0

        return x

    def _G_o(field,data):

        pe  = data[('gas','Pe_heating_rate')].convert_to_units('erg/s/cm**3').value
        Z   = (data['Metal_Density'] / data['Density']).value
        n_H = (data['H_p0_number_density'] + data['H_p1_number_density'] + data['H_m1_number_density'] +\
               data['H2_p0_number_density'] + data['H2_p1_number_density']).convert_to_units('cm**(-3)').value

        g_to_d = 0.68 - 3.08 * np.log10(Z / 0.014)
        d_to_g = 1.0 / (10.0**(g_to_d))

        D = d_to_g / 6.616595E-3

        epsilon = 0.01488637246 * (n_H)**(0.235269059)

        atten = np.exp( - 1.33E-21 * D * data['dx'].convert_to_units('cm').value * n_H)

        G_o = pe / (1.3E-24 * n_H * epsilon * D * atten)

        return G_o * (data['Density'] / data['Density'])

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


    ds = yt.load(dsname)
    yt.add_field(('gas','Pe_heating_rate'), function = _pe_heating_rate, units = 'erg/s/cm**3')
    yt.add_field(('gas','Pe_heating_rate_masked'), function = _pe_heating_rate_masked, units = 'erg/s/cm**3')
    yt.add_field(('gas','G_o'), function = _G_o, units = "")
    yt.add_field(('gas','FUV_flux'), function = _FUV_flux, units = "erg/s/cm**2")
    yt.add_field(('gas','LW_flux'), function = _LW_flux, units = "erg/s/cm**2")
    yt.add_field(('gas','Q0_flux'), function = _Q0_flux, units = "erg/s/cm**2")
    yt.add_field(('gas','Q1_flux'), function = _Q1_flux, units = "erg/s/cm**2")
    ds = yt.load(dsname)

    return ds
