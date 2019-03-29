import numpy as np
from scipy.interpolate import griddata

try:
    _ion_table = np.genfromtxt("/home/aemerick/code/galaxy_analysis/data/lt00HM12_h1_he2_c3_c4_o4_o6_ne8_mg2_si4", names = True)
else:
    _ion_table = None

def get_ion_fraction(n, T, ion):

    f = get_function(ion)

    _n = n.flatten()
    _T = T.flatten()

    fraction = f(_n,_T)

    return fraction.reshape(np.shape(n))

def get_function(ion):

    if not ('f' in ion):
        ion = 'f' + ion

    return lambda n,T : griddata( (_ion_table['nh'], _ion_table['T']),
                     _ion_table[ion], (n,T), fill_value = -99)

#    return  interp2d(_ion_table['nh'],
#                     _ion_table['T'] ,
#                     _ion_table[ion]  )

def get_ions():
    fields = _ion_table.dtype.names
    return [x[1:] for x in fields if ( not (x == 'nh')) and ( not (x=='T'))]

def list_ions():
    print get_ions()
    return

def get_elements(ions):
    """
    Strips roman numerals and returns elements
    """

    strip_strings = ['III','II','I','IV','V','VII','VIII']
    if np.size(ions) > 1:
        ele = [None] * len(ions)
        for i,x in enumerate(ions):
            ele[i] = [ x.strip(y) for y in strip_strings][-1]
    else:
        ele = [ions.strip(y) for y in strip_strings][-1]

    return ele
