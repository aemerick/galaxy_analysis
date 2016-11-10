import numpy as np

def species_from_fields(fields):
    """
    Returns a list of the individual species tracer fields
    given a list of enzo fields
    """

    element_fields = [x[1] for x in fields if ('_Density' in x[1] and\
                                               len(x[1]) <=10)]
    elements = np.sort([x.rsplit('_')[0] for x in element_fields])
    metals   = [x for x in elements if (x != 'H' and x != 'He' and x != 'HI')]


    return metals

def ratios_list(species):
    """
    Returns a list of abundance ratio names given a list
    of metal tracer species. Names are really all potentially
    interesting abundance ratios that can be derived from the given
    element names. This includes:

        1) all combinations of x/H , x/Mg, x/Fe, x/Na, x/Ni, if
           denominators are present [H is always assumed present]
        2) Y/Ba, Ba/Eu (if present)
        3) Fe/x, Ba/x (if present)

    where x is any element name from the passed in list
    """

    # remove H and He if they are there
    species = [x for x in species if x != 'H' and x != 'He']


    denominators = ['H', 'Mg', 'Fe', 'Na', 'Ni']

    ratios = []

    for d in denominators:

        if not d in species and d != 'H':
            continue

        for s in species:
            if s == d:
                continue

            ratios += [s + '/' + d]

    special_ratios = [('Y','Ba'), ('Ba','Eu')]

    for sr in special_ratios:
        if sr[0] in species and sr[1] in species:
            ratios += [sr[0] + '/' + sr[1]]

    # a couple more special ones:
    special_numerator = ['Fe', 'Ba']

    for n in special_numerator:

        if not n in species:
            continue

        for s in species:

            if s == n:
                continue

            ratios += [n + '/' + s]

    # remove any duplicates
    ratios = list(np.sort(np.unique(ratios)))

    return ratios
