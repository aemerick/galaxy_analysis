"""
   Author: Britton Smith
   Modified: Andrew Emerick


   Generate modified high-z Cloudy tables for Grackle using
   an updated model for the LW background (and scale d rates
   for H- and H2 dissociation) spliced to merge with
   the lower-z UVB for consistency.
"""

import h5py
import numpy as np

from h5_in_memory import H5InMemory

#
# Lyman Werner model to use
#

#LW_model = 'RFT14'  # original used here and Enzo
#LW_model = 'JW2012' 
LW_model = 'Qin2020' # most up to date model 

def cos_ramp(z, off, on):
    my_z = np.clip(z, on, off)
    x = np.pi*(my_z - off)/(on - off)
    return -0.5*np.cos(x)+0.5

def interpvals(xtab, ytab, x, log=True):
    i = np.digitize(x, xtab) - 1
    i = np.clip(i, 0, xtab.size-2)
    if log:
        m = np.log10(ytab[i+1] / ytab[i]) / np.log10(xtab[i+1] / xtab[i])
        return np.power(10, m * np.log10(x / xtab[i]) + np.log10(ytab[i]))
    else:
        m = (ytab[i+1] - ytab[i]) / (xtab[i+1] - xtab[i])
        return m * (x - xtab[i]) + ytab[i]

def k31_RFT14(z):
    k31 = np.zeros_like(z)
    k31[z > 50] = 0.
    logJ = -23.56688 + 4.56213e-1 * (1.0+z) - \
            2.67982e-2 * np.power(1.0+z, 2.0) + \
            5.88234e-4 * np.power(1.0+z, 3.0) - \
            5.05576e-6 * np.power(1.0+z, 4.0)
    fil = (z <= 50) & (z > 6)
    k31[fil] = 1.13e8 * np.power(10.0, logJ[fil])
    k31[z <= 6] = 1.13e8 * 1.0e-21
    return k31

def k31_JW2012(z):
    """
    Power-law fit to Qin+2020 - reference model
    taken from digitizing Figure 5
    """
    
    k31 = np.zeros_like(z)
    k31[ z > 50] = 0.

    logJ21 = -2.356 + 0.45620 * z +\
             -0.02680 * z * z +\
              5.882E-4 * z * z * z+\
             -5.056E-6 * z * z * z * z

    fil = (z <= 50) & (z > 6)

    k31[fil] = 1.13e8 * 1.0E-21 * np.power(10.0, logJ21[fil])
    k31[z <= 6] = 1.13e8 * 1.0E-21

    return k31

def k31_Qin2020(z):
    """
    Power-law fit to Qin+2020 - reference model
    taken from digitizing Figure 5
    """
    
    k31 = np.zeros_like(z)
    k31[ z > 50] = 0.

    logJ21 = 1.6685 - (1.5068E-1)*z + (2.3357E-3)*z*z +\
             (-6.8115E-5) * z * z * z

    fil = (z <= 50) & (z > 6)

    k31[fil] = 1.13e8 * 1.0E-21 * np.power(10.0, logJ21[fil])
    k31[z <= 6] = 1.13e8 * 1.0E-21

    return k31

if __name__ == "__main__":
    # filename = "CloudyData_UVB=HM2012.h5"
    # output_filename = "CloudyData_HM2012_highz.h5"
    filename = "CloudyData_UVB=HM2012_shielded.h5"
    output_filename = "CloudyData_HM2012_highz_shielded.h5"

    all_filenames = [ ("CloudyData_UVB=HM2012.h5","CloudyData_HM2012_highz.h5"),
                      ("CloudyData_UVB=HM2012_shielded.h5","CloudyData_HM2012_highz_shielded.h5")]

    for filename, output_filename in all_filenames:
        fd = H5InMemory(filename)

        ztable = fd['UVBRates']['z'].value
        Xtable = ztable + 1
        zi = 50.
        lXi = np.log10(zi + 1)
        dlX = np.diff(np.log10(ztable + 1)).mean()

        # z+1 values for new table
        Xnew = np.logspace(0, lXi, int(lXi / dlX) + 1)
        znew = Xnew - 1

        # interpolate HM2012 k27, k28, and k31 values to new z table
        k27_HM2012 = fd['UVBRates']['Chemistry']['k27'].value
        k27_int = interpvals(Xtable, k27_HM2012, Xnew)
        k28_HM2012 = fd['UVBRates']['Chemistry']['k28'].value
        k28_int = interpvals(Xtable, k28_HM2012, Xnew)
        k31_HM2012 = fd['UVBRates']['Chemistry']['k31'].value
        k31_int = interpvals(Xtable, k31_HM2012, Xnew)

        if LW_model == 'RFT14':
            # John's k31 model
            k31_model = k31_RFT14(znew)
        elif LW_model == 'JW2012':
            k31_model = k31_JW2012(znew)
        elif LW_model == 'Qin2020':
            k31_model = k31_Qin2020(znew)
        else:
            print(LW_model + " not understood")
            raise ValueError

        chem_rates = {}

        # Combine with a cosine ramp to make new model
        # z_off = 50 and z_on = 7 gives a nice smooth transition
        highz_ramp = cos_ramp(np.log10(Xnew), np.log10(50+1), np.log10(7+1))
        chem_rates['k31'] = np.power(
             10, (np.log10(k31_int) *      highz_ramp +
                 np.log10(k31_model)  * (1 - highz_ramp)))

        # scale k27 and k28 to John's k31 model according to a 30,000 K blackbody
        k27_scale = 3.862260e+01 / 1.553719e+00
        k28_scale = 5.400897e+00 / 1.553719e+00
        k27_model = k31_model * k27_scale
        k28_model = k31_model * k28_scale
        chem_rates['k27'] = np.power(
            10, (np.log10(k27_int) *      highz_ramp +
                 np.log10(k27_model)  * (1 - highz_ramp)))
        chem_rates['k28'] = np.power(
            10, (np.log10(k28_int) *      highz_ramp +
                 np.log10(k28_model)  * (1 - highz_ramp)))

        # Set the rest of the high-z photo rates to zero.
        chem_ramp = cos_ramp(np.log10(Xnew), np.log10(15+1), np.log10(10+1))
        for rate in ['k24', 'k25', 'k26', 'k29', 'k30']:
            rate_int = interpvals(
                Xtable, fd['UVBRates']['Chemistry'][rate].value, Xnew)
            chem_rates[rate] = rate_int * chem_ramp

        # Set photo-heating rates similarly
        photo_rates = {}
        for rate in ['piHI', 'piHeI', 'piHeII']:
            rate_int = interpvals(
                Xtable, fd['UVBRates']['Photoheating'][rate].value, Xnew)
            photo_rates[rate] = rate_int * chem_ramp

        # Update values in h5 object and save new file
        for rate in chem_rates:
            del fd['UVBRates']['Chemistry'][rate]
            fd['UVBRates']['Chemistry'].create_dataset(
                rate, data=chem_rates[rate])
        for rate in photo_rates:
            del fd['UVBRates']['Photoheating'][rate]
            fd['UVBRates']['Photoheating'].create_dataset(
                rate, data=photo_rates[rate])

        del fd['UVBRates']['z']
        fd['UVBRates'].create_dataset('z', data=znew)

        fd.save(output_filename)
