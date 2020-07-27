"""
  Gather yields from PopIII Star models

"""


import numpy as np
from galaxy_analysis.static_data import asym_to_anum, anum_to_asym
import re as _re


def _remove_num(val):
    pattern='[0-9]'
    return _re.sub(pattern,'',val)

def _remove_num_from_list(list):
    """
    Removes all numbers within a string for a list of strings.
    """
    list = [_remove_num(i) for i in list] 
    return list

def generate_heger_woosley_2010(filepath = './heger_woosley_2010_yields.dat',
                                explosion_energy = 1.2,
                                mixing=0.001,
                                piston="S4"):
    """
    Takes the full, machine readable table from Heger + Woosley 2010 and computes
    the total yields for each element (summing over isotopes).

    Heger and Woosley 2010 provides many different yield models, varying
    explosion energy, mixing, and piston location. They provide a
    'standard' model with explosion energy of 1.2E51, mixing value of 0.1,
    and piston postion (S4).

    Parameters
    ----------
    filepath : string, optional
        Path to the Heger and Woosley 2010 full, machine-readable yield table.

    explosion_energy : float, optional
        Explosion energy of model in units of 1.0E51 erg. Default is
        the 'standard' model. Available values are:
        0.3,  0.6,  0.9,  1.2,  1.5,  1.8,  2.4,  3. ,  5. , 10.
        Default : 1.2

    mixing : float, optional
        Mixing parameter. Default is 'standard'. Available values:
        0.     , 0.001  , 0.00158, 0.00251
        Default : 0.001

    piston : string, optional
        Piston location. Available options are 'S4' and 'Ye'.
        Default: 'S4'

    Returns
    -------

    final_yields : 2D array
        2D numpy array containing the yields for ALL elements from H
        to Bi for all sampled masses. Yields not included in table are
        zeroed.

    masses :
        masses contained in the model
    """

    raw_heger_yields = np.genfromtxt(filepath, skip_header=17,
                                     dtype=[('mass',"f8"),('energy',"f8"),
                                            ("piston","|U2"),("mixing","f8"),
                                            ("isotope","|U5"),("yield","f8")])

    if not (explosion_energy in np.unique(raw_heger_yields["energy"])):
        print("Explosion energy %4.2E not available"%(explosion_energy))
        print("Available values:  ", np.unique(raw_heger_yields["energy"]))
        raise ValueError

    if not (piston in np.unique(raw_heger_yields["piston"])):
        print("Piston model %5S not available"%(piston))
        print("Available values:  ", np.unique(raw_heger_yields["piston"]))
        raise ValueError

    if not (mixing in np.unique(raw_heger_yields["mixing"])):
        print("Mixing parameter %4.2E not available"%(mixing))
        print("Available values:  ", np.unique(raw_heger_yields["mixing"]))
        raise ValueError

    model_select =   (raw_heger_yields["energy"] == explosion_energy) *\
                     (raw_heger_yields["piston"] == piston) *\
                     (raw_heger_yields["mixing"] == mixing)

    selected_yields = raw_heger_yields[model_select]

    all_isotopes    = np.unique(selected_yields["isotope"])
    all_elements    = np.unique(_remove_num_from_list( all_isotopes ))

    isotope_element = np.array(_remove_num_from_list( selected_yields["isotope"])) # big list of element name for each val in table

    masses          = np.unique(selected_yields["mass"])


    # now make a uniform table that sums over isotopes and has a slot for
    # each element in the heger table at each mass
    processed_yields     = np.zeros((np.size(all_elements),np.size(masses)))

    for i in np.arange( np.size(all_elements) ):
        for mi in np.arange( np.size(masses)):
            select = (selected_yields['mass'] == masses[mi])*(isotope_element==all_elements[i])
            yields_to_sum = selected_yields['yield'][select]

            if np.size(yields_to_sum)>0:
                processed_yields[i][mi] = np.sum(yields_to_sum)
            else:
                processed_yields[i][mi] = 0.0

    # get atomic numbers that exist
    table_anum    = np.array([asym_to_anum[str(x).strip()] for x in all_elements])
    all_anum      = np.arange(1,84,1) # H to Bi


    # final yield table is : Total, Total metals, H, He, ..... , Bi
    final_yields = np.zeros( (np.size(all_anum)+2, np.size(masses)))

    for i in np.arange(np.size(all_anum)):
        if all_anum[i] in table_anum:
            final_yields[i + 2] = processed_yields[table_anum == all_anum[i]]

    final_yields[0] = np.sum(final_yields, axis = 0)     # total yields
    final_yields[1] = np.sum(final_yields[4:], axis = 0) # total metal yields


    return final_yields, masses


    return

def generate_nomoto(filepath = './nomoto_latex_table.dat'):
    """
    Takes the Z=0 portion of Table 2 of Nomoto et. al. 2006
    and sums over all elemental isotopes to combine a table
    that gives yields for all elements from H to Bi and includes
    a total yield and metal yield value.

    This spans metal free star masses from 13 to 40 Msun
    """

    raw_nomoto_yields = np.genfromtxt(filepath,
                                  delimiter = '&', names = True,
                                  dtype = "|U2,f8,f8,f8,f8,f8,f8,f8")

    nomoto_e          = np.unique(raw_nomoto_yields['element'])
    masses            = np.array([float(y) for y in raw_nomoto_yields.dtype.names[1:]])
    nomoto_yields     = np.zeros(
                           (np.size(nomoto_e),7)
                         )

    j = 0
    for i in np.arange( np.size(nomoto_e)):

        for j in np.arange(np.size(raw_nomoto_yields['element'])):

            if raw_nomoto_yields['element'][j] == nomoto_e[i]:
                nomoto_yields[i] +=  np.array(list(raw_nomoto_yields[j])[1:])


    # get atomic numbers that exist
    nomoto_anum    = np.array([asym_to_anum[str(x).strip()] for x in nomoto_e])
    all_anum       = np.arange(1,84,1) # H to Bi

    final_yields = np.zeros( (np.size(all_anum)+2, 7))

    for i in np.arange(np.size(all_anum)):
        if all_anum[i] in nomoto_anum:
            final_yields[i + 2] = nomoto_yields[nomoto_anum == all_anum[i]]

    final_yields[0] = np.sum(final_yields, axis = 0)
    final_yields[1] = np.sum(final_yields[4:], axis = 0)


    return final_yields, masses


def generate_heger_woosley_2002(filepath = './heger_latex_table.dat'):
    """
    Heger & Woosley 2002 - Table 3

    Table gives yields in He core masses. Using Eq 1 in HW2002 to
    convert this to progenitor masses. This spans
    progenitor masses from 140 to 260 Msun.

    """
# anum element 65 70 75 80 85 90 95 100 105 110 115 120 125 130
    raw_heger_yields = np.genfromtxt(filepath,
                                     delimiter = '&', names = True,
                                     dtype = "i2,|U4,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8")

    heger_e      = [str(x).strip() for x in np.unique(raw_heger_yields['element'])]
    He_masses    = np.array([float(y) for y in raw_heger_yields.dtype.names[2:]])

    # table gives He core masses - Eq 1 in Heger&Woosley2002
    masses       = He_masses*(24.0/13.0) + 20.0
    heger_yields = np.zeros( (np.size(heger_e),np.size(masses)))

    j = 0
    for i in np.arange(np.size(heger_e)):
        for j in np.arange(np.size( raw_heger_yields['element'])):
            if str(raw_heger_yields['element'][j]).strip() == heger_e[i]:
                heger_yields[i] += np.array(list(raw_heger_yields[j])[2:])

    # get atomic numbers
    heger_anum = np.array([asym_to_anum[str(x).strip()] for x in heger_e])
    all_anum   = np.arange(1,84,1) # H to Bi

    # final yield table is : Total, Total metals, H, He, ..... , Bi
    final_yields = np.zeros( (np.size(all_anum)+2,np.size(masses)))

    for i in np.arange(np.size(all_anum)):
        if all_anum[i] in heger_anum:
            final_yields[i+2] = heger_yields[ heger_anum == all_anum[i] ]

    final_yields[0] = np.sum(final_yields, axis=0)       # all elements
    final_yields[1] = np.sum(final_yields[4:], axis = 0) # all metals

    return final_yields, masses


def generate_table(low_mass_model = 'heger_woosley_2010',
                   pisne_model    = 'heger_woosley_2002'):
    """
    Combine yield tables together
    """

    if low_mass_model == 'heger_woosley_2010':
        low_mass_yields, low_mass_masses = generate_heger_woosley_2010()
    elif low_mass_model == 'nomoto':
        low_mass_yields, low_mass_masses = generate_nomoto()
    else:
        print("Low mass model type not recognized")
        raise ValueError

    if pisne_model == 'heger_woosley_2002':
        pisne_yields,  pisne_masses = generate_heger_woosley_2002()


    # combine into a single table

    final_table  = np.concatenate( [low_mass_yields, pisne_yields], axis=1)
    masses_final = np.concatenate( [low_mass_masses, pisne_masses], axis=0)

    return final_table, masses_final

def generate_ascii_table(outname = './popIII_yields.dat',
                   low_mass_model = 'heger_woosley_2010',
                   pisne_model    = 'heger_woosley_2002'):

    final_table, masses_final = generate_table(low_mass_model, pisne_model)

    header = "# M Z m_tot m_metal " + ' '.join( [anum_to_asym[x] for x in np.arange(1,84,1)]) + '\n'

    output = open(outname, 'w')

    output.write(header)

    metallicity = 0.0
    for j in np.arange(np.size(masses_final)):

        output.write("%.3f %.3f "%(masses_final[j], metallicity))

        for i in np.arange(np.size(final_table[:,0])):
            output.write("%8.8E "%( final_table[i,j]))

        output.write("\n")

    output.close()

    return final_table, masses_final


if __name__ == '__main__':

    table, m = generate_ascii_table(outname='temp.dat')

    print(np.shape(table))
    print(np.shape(table.reshape(np.shape(table)[0],1,np.shape(table)[1])))
    print(np.shape(table.reshape(np.shape(table)[1],1,np.shape(table)[0])))
    print(np.shape(m))
