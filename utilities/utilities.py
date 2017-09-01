import numpy as np
import sys
import contextlib

from galaxy_analysis.static_data import asym_to_anum

# to compute weighted statistics on a distribution
import statsmodels.stats.weightstats as stats
import statsmodels.tsa.stattools as stattools

rowcoldict = {2 : (1,1), 3: (1,3), 4:(2,2),
              5 : (2,3), 6: (2,3), 7:(2,4),
              8 : (2,4), 9: (2,5), 10:(2,5),
              11: (3,4), 12: (3,4), 13: (4,4),
              14: (4,4), 15: (4,4), 16: (4,4),
              17: (5,4), 18: (5,4), 19: (5,4), 20 : (5,4)}

def extract_nested_dict(dict, key_list):
    """
    Given a list of kwargs, extracts the information requested
    from a nested dictionary
    """
    x = dict
    for k in key_list:
        x = x[k]

    return x


def nested_haskey(x, keys):
    """
    For a nested dictionary 'x' and list of keys, checks
    if all keys exist in the nested key path.
    """

    if len(keys) == 1:
        return (keys[0] in x)

    if keys[0] in x:
        return nested_haskey( x[keys[0]], keys[1:])
    else:
        return False

def filter_dict(field, dictionary, level = 2):
    """
    Filter through nested dictionarys like one would do with
    arrays. Only works if sub-dictionaries all have same 
    kwargs as desired by 'field'
    """
    return [(x,dictionary[x][field]) for x in dictionary.keys]

class _DummyFile(object):
    """
    A do nothing write for use with below
    """
    def write(self, x): pass

@contextlib.contextmanager
def nooutput():
    """
    Silce stdout and stderr
    """
    save_stdout = sys.stdout
    save_stderr = sys.stderr
    sys.stdout  = _DummyFile()
    sys.stderr  = _DummyFile()
    yield
    sys.stdout = save_stdout
    sys.stderr = save_stderr
    return

@contextlib.contextmanager
def nostdout():
    """
    Silence a called function's stdout using:
        with nostdout():
            foo()
    """

    save_stdout = sys.stdout
    sys.stdout  = _DummyFile()
    yield
    sys.stdout = save_stdout
    return

#
# could probably do the check for scalar input stuff
# with a with statement in a similar way to the above
# but would be a little trickier
#

def sort_by_anum(names):
    """
    For a list of either element symbols or 
    arbitrary strings beginning with element symbols (with a
    consistent pattern), sort by atomic number. Returns
    the sorted list.
    """

    # find common part of string if it exists...
    if any( np.array([len(x) for x in names]) > 2):
        # assume all have the same common string
        index = np.argmax([len(x) for x in names])
        print index
        sample = names[index]
        print sample
        s_to_remove = sample[2:] # assume longest string is 2 letter sym
        print s_to_remove
        i = -1 * len(s_to_remove)
        filtered = [x[:i] for x in names]
        print filtered
    else:
        filtered = [x for x in names]

    # get list of atomic numbers
    anum = [ asym_to_anum[x] for x in filtered ]

    # get sorting array
    sort = np.argsort(anum)

    return list(np.array(names)[sort])

def compute_stats(x, return_dict = False, acf = False):

    if not return_dict:
        return np.min(x), np.max(x), np.average(x), np.median(x), np.std(x)

    d = {}
    d['min'] = np.min(x)
    d['max'] = np.max(x)
    d['mean'] = np.average(x)
    d['median'] = np.median(x)
    d['std'] = np.std(x)
    d['Q1']  = np.percentile(x, 25.0)
    d['Q3']  = np.percentile(x, 75.0)
    d['inner_quartile_range'] = d['Q3'] - d['Q1']
    d['variance'] = np.var(x)

    if (len(x) > 0) and acf:
        d['acf'] = np.array(stattools.acf(x))
    else:
        d['acf'] = 0.0
    return d

def compute_weighted_stats(x, w, return_dict = True):
    """
    Returns a dictionary containing stats information to be easiy chunked
    into a data set. Does this by constructing a weighted stats object
    using the statsmodels 
    """

    s = stats.DescrStatsW(x, weights = w)

    if not return_dict:
        return s

    d = {}
    d['mean']     = s.mean()
    d['std' ]     = s.std()
    q             = s.quantile( np.array([0.25, 0.75]))
    d['Q1']       = q[0]
    d['Q3']       = q[1]
    d['inner_quartile_range'] = q[1] - q[0]
    d['variance'] = s.var()
    d['min']      = np.min(x)
    d['max']      = np.max(x)

    return d

def chemistry_species_from_fields(fields):
    """
    Returns a list of the individual chemical species fields
    used in the non-equillibrium chemistry solver

    Example, returns ["HI", "HII", "HeI", "HeII", "HeIII"] if
    primordial chemistry is set to 1
    """

    species = [x[1] for x in fields if ('_Density' in x[1] and\
                                        ('I_' in x[1] or 'HM' in x[1]))]

    species = np.sort([x.rsplit('_')[0] for x in species])

    return sort_by_anum(species)

def species_from_fields(fields, include_primordial = False):
    """
    Returns a list of the individual species tracer fields
    given a list of enzo fields
    """

    element_fields = [x[1] for x in fields if ('_Density' in x[1] and\
                                               len(x[1]) <=10)]
    elements = np.sort([x.rsplit('_')[0] for x in element_fields])

    ignore = ['H', 'He', 'HI', 'HM', 'HD']

    metals   = [x for x in elements if x not in ignore]

    if include_primordial:
        metals = ['H','He'] + metals

    return sort_by_anum(metals)

def abundance_ratios_from_fields(fields):
    """
    Returns all abundance ratio fields that are defined (does not define
    any)

    Assume abundance ratio fields are all:
      ('gas','X_over_Y')
    where X and Y are atomic symbols. Returns None if there are no
    fields
    """

    fields = [x[1] for x in fields\
            if (      ('_over_' in x[1])\
                 and (not 'particle' in x[1])\
                 and (len(x[1]) < 11))]

    if len(fields) == 0:
        fields = None

    return fields

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
    species = sort_by_anum(species)
    species = [x for x in species if x != 'H' and x != 'He']


    denominators = sort_by_anum(['H', 'Mg', 'Fe', 'Na', 'Ni', 'O', 'C'])

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


