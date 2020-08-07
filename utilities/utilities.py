import numpy as np
import sys
import contextlib
import glob
import deepdish as dd

from galaxy_analysis.static_data import asym_to_anum

from astroML.time_series import ACF_EK

def map_to_pixels(x0,x1,y0=None, y1=None):
    """
    Map a single axis (or two axes with optional
    additional kwargs) in coordinates. Really this can be
    used to do any general linear mapping, but written to be
    used for mapping data coordinates to pixel coordinates for overplotting
    data on already-existing plot images.

    Returns a callable function
    """
    def _map_to_pixels(coord0, coord1):
        m = (coord1[1] - coord0[1]) / (coord1[0] - coord0[0])
        return (lambda x : m * (x - coord0[0]) + coord0[1])

    p_x = _map_to_pixels(x0,x1)

    if (y0 is None) and (y1 is None):
        return p_x
    elif ( (not (y0 is None)) and (not (y1 is None))):
        p_y = _map_to_pixels(y0,y1)
        return p_x, p_y
    else:
        print("Must supply two additional points for y-axis to generate additional function")
        raise ValueError

    return

def bin_edges(x):
    xnew = np.zeros(len(x) + 1)
    dx   = np.zeros(len(x) + 1)
    dx[1:-1]   = x[1:] - x[:-1]
    dx[0] = dx[1]
    dx[-1] = dx[-2]
    xnew[:-1] = x - 0.5*dx[:-1]
    xnew[-1] = x[-1] + 0.5*dx[-1]
    return xnew


def window_average(x, y, dx, where = 'pre', edge_case = 'trim',
                   sample_points = None, window=None,
                   method = 'average'):
    """
    Given input x and y values (x and y must be the same size),
    computes the window average using a width of dx at each x value.

    Parameters:
    -------
    x : np.ndarray
        Array of x values

    y : np.ndarray
        Array of y values

    dx : int,float
        Width (in x) to sample and average over

    where : string, optional
        Where to place the window. Default: 'pre':

            'pre'  : average over [x-dx,x)
            'mid'  : average over [x-dx/2, x+dx/2)
            'post' : average over [x,x+dx)

    edge_case : string, optional
           How to handle when the window is near the min / max
       of the x values. Default : trim

           'trim'   : reduce window size appropriately at edges.
                      Example, if 'pre' is used, for all x where
                      x - dx < xmin, the window is
                      [xmin, x).
           'fix'    : fix window size. This will flatten result
                      at edges. Example if 'pre' is used, for all
                      x where x - dx < xmin, the window
                      [xmin, xmin + dx) is used.

    sample_points : np.ndarray, optional
        Use custom sample points to do the averaging. Defaults to
        the input x array.

    window : string, optional
        Not yet implemented. Option to choose different window
        types.


    method : string, optional
        What to do with the points in the window. Default : average.
        Options include, average, sum, median, min, and max
    """

    if np.size(x) != np.size(y):
        print("x and y must be the same length: ",np.size(x),np.size(y))

        if np.size(x) == (np.size(y)+1):
            print("But it looks like x is just greater by one. Assuming bins and using centers")
            xvals = 0.5 * (x[1:] + x[-1])
        else:
            raise(ValueError)
    else:
        xvals = 1.0*x

    if sample_points is None:
        sample_points = 1.0*x

    output = np.zeros(np.size(sample_points))

    if edge_case == 'fix':
        print("fix edge case is not yet implemented... sorrry")
        raise(RuntimeError)

    if where == 'pre':

        select_function = lambda array, z : (array >= z-dx) * (array < z)

    elif where == 'mid':

        select_function = lambda array, z : (array >= z-dx*0.5) * (array < z+dx*0.5)

    elif where == 'post':

        select_function = lambda array, z : (array >= z) * (array < z+dx)


    if method == 'average':
        yfunction = lambda z : np.average(z)
    elif method == 'sum':
        yfunction = lambda z : np.sum(z)
    elif method == 'median':
        yfunction = lambda z : np.median(z)
    elif method == 'max':
        yfunction = lambda z : np.max(z)
    elif method == 'min':
        yfunction = lambda z : np.min(z)

    for i in np.arange(np.size(output)):

        select    = select_function(xvals, sample_points[i])
        output[i] = yfunction(y[select])


    return output


def simple_rebin(bins, y, new_bins, method = "sum", force=False):
    """
    Re-bin data (y) placed in evenly spaced bins (bins)
    where len(bins) = len(y) + 1. It is assumed the newly
    binned data should sum over the old bins, (method = 'sum'),
    but this can be changed to average over the values in the
    bins with method = 'average'.

    At the moment, new_bins must be such that there is even overlap
    with the previous bins, AND that it spans the full domain
    of the previous bins. new_bins can either be a float, representing
    a desired bin spacing, or an array representing the
    new bins.

    force computes new data even if bins are not evenly spaced. prints
    just the error message instead of failing.

    returns new bins and new values
    """

    if not (len(bins) == len(y) + 1):
        print("Length of bins must be greater than y by 1")
        raise ValueError

    dx = np.unique(np.diff(bins))
    if np.size(dx) > 1:
        dx = np.unique(np.round(np.diff(bins),10))
        if np.size(dx) > 1:
            print("Original bins must be evenly spaced")
            print(dx)
            print(np.diff(bins))
            print(bins)

            if not force:
                raise ValueError

    dx = dx[0] # unique returns an array even if only 1 value

    if np.size(new_bins) <= 1:
        new_dx   = new_bins
        new_bins = np.arange(np.min(bins), np.max(bins)+0.5*new_dx, new_dx)
    else:
        new_dx = np.unique(np.diff(new_bins))[0]

    new_y = np.zeros(np.size(new_bins)-1)

    if method == 'sum':
        func = np.sum
    elif method == 'average' or method == 'avg':
        func = np.average
    else:
        print("method must be either 'sum' or 'average'")
        raise ValueError

    tempy = np.zeros(np.size(y)+1) # bit of a hack unfortunately
    tempy[:-1] = y
    tempy[-1]  = y[-1]
    for i in np.arange(np.size(new_y)):
        new_y[i] = func( tempy[ (bins >= new_bins[i]) * (bins < new_bins[i+1]) ])


    return new_bins, new_y

def acf(x, y, dy = None, bins = None):

    # just pick small percentage error
    if dy is None:
        dy = 0.001 * y

    if bins is None:
        bins = np.linspace(np.min(x), np.max(x), 51)

    acf, acf_err, bins = ACF_EK(x, y, dy, bins)

    return acf, acf_err, bins

rowcoldict = {2 : (1,1), 3: (1,3), 4:(2,2),
              5 : (2,3), 6: (2,3), 7:(2,4),
              8 : (2,4), 9: (2,5), 10:(2,5),
              11: (3,4), 12: (3,4), 13: (4,4),
              14: (4,4), 15: (4,4), 16: (4,4),
              17: (5,4), 18: (5,4), 19: (5,4), 20 : (5,4),
              21: (4,6), 22: (4,6), 23: (4,6), 24 : (4,6),
              25: (5,5)}

def get_property(field_path, file_list=None, dir = '.', tmin = None, tmax = None, data_list = None,
                 times = None, self_contained = False):
    """
    field_path can either be a list of strings corresponding to the kwargs in the nested dictionary
    to pull from, or just a single string with kwargs separated by '/'. Must have a leading '/'.

       example: either is fine:
             field_path = ['meta_data','M_HI']
        or
             firld_path = '/meta_data/M_HI'
    """

    return_times = False
    if file_list is None:
        return_times = True

    if file_list is None:
        file_list, t = select_data_by_time( dir = dir, tmin = tmin, tmax = tmax,
                                            data_list = None, times = None, self_contained = False)

    if isinstance(field_path, str):
        if field_path[0] != '/':
            field_path = '/' + field_path
    else:
        field_path = "/" + "/".join(field_path)

    if self_contained:
        if not isinstance(file_list, str):
            print("If loading self-contained data, file_list MUST be the name of the file")
            raise ValueError

        if data_list is None:
            print("If loading self-contained data, data_list must be kwargs of 'files' in top level to loop over")
            raise ValueError

        load_data = lambda dname : dd.io.load(file_list, '/' + str(dname) + field_path)

        x = np.array( list(map( load_data, data_list)))
    else:
        load_data = lambda fname : dd.io.load(fname, field_path)
        x = np.array( list(map( load_data, file_list )) )

    if return_times:
        return x, t
    else:
        return x

def select_data_by_time(dir = '.', tmin = None, tmax = None,
                        data_list = None, times = None, self_contained = False):
    """
    Get list of data files that lie within the desired time range. Can pass a limited
    set of files or precomputed time array (if available). If self_contained is true,
    then rather than globbing for all files in a directory, assumes the data
    file is self contained.
    """

    if (tmax is None) and (tmin is None):
        print("Must set a tmin OR tmax, or both")
        raise ValueError

    if tmin is None: tmin = np.inf
    if tmax is None: tmax = np.inf

    if data_list is None:
        data_list = np.sort(glob.glob(dir + '/DD????_galaxy_data.h5'))

    times = np.zeros(len(data_list))

    if self_contained:
        for i, d in enumerate(data_list):
            times[i] = sc_data[d]['general']['Time']
    else:
        for i, d in enumerate(data_list):
            times[i] = dd.io.load(d, '/meta_data/Time')

    # now select which data sets to average
    if tmin is None and tmax is None:
        data_list = data_list
    elif tmin is None:
        data_list = data_list[ times < tmax ]
        times     = times[times<tmax]
    elif tmax is None:
        data_list = data_list[  times >= tmin ]
        times     = times[times>=tmin]
    else:
        data_list = data_list[(times<tmax)*(times>=tmin)]
        times     = times[(times<tmax)*(times>=tmin)]

    return data_list, times



def extract_nested_dict(dict, key_list):
    """
    Given a list of kwargs, extracts the information requested
    from a nested dictionary
    """
    x = dict

    if isinstance(key_list, str) or isinstance(key_list, tuple):
        x = x[key_list]
    elif isinstance(key_list, list):
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

def filter_dict(field, dict, level = 2):
    """
    Filter through nested dictionarys like one would do with
    arrays. Only works if sub-dictionaries all have same
    kwargs as desired by 'field'
    """
    return [(x,dict[x][field]) for x in dict.keys()]

def extract_nested_dict_aslist(dict, key_list, loop_keys = None, self_contained = True):
    """
    If one wants to extract a value deep in a nested dictionary
    across many "rows". Example, your dictionary looks like this maybe:

                                "time" - single_value
                               /
                        "DD0000"
                       /       \
                      /         "HI_Mass" - single_value
                     /
                    /           "time" - single_value
                   /           /
    dictionary ------- "DD0001"
                   \          \
                    \           "HI_Mass"      - single_value
                     \
                      \         "time" - single_value
                       \       /
                        "DD0002"
                               \
                                "HI_mass" - single value

    And you want to be able to plot the HI_mass for each data set as
    a function of, say, time. You would want each of these as a 1D array.
    You can do this as:

        times   = extract_nested_dict_aslist(dictionary, ["time"])
        HI_mass = extract_nested_dict_aslist(dictionary, ["HI_mass"])

    and if you only wanted it for a certain selection of the top level keys:

        times = extract_nested_dict_aslist(dictionary, ["time"], loop_keys = ['DD0000','DD0001'])

    and finally, if the desired key was one or more levels down (arbitrary),
    say something like "abundances/CNM/Fe_over_H", which would be the Fe_over_H
    abundance in the CNM, then one can do:

        CNM_Fe_over_H = extract_nested_dict_aslist(dictionary, ['abundances','CNM','Fe_over_H'])

    """
    if self_contained: # data sets are all contained in overarching dictionary
        if loop_keys is None:
            loop_keys = dict.keys()

        return [extract_nested_dict( dict[k], key_list) for k in loop_keys]
    else:
        # data sets are not contained, and are separate HDF5 files
        # --- this is a bit ugly of a process, but hopefully we can limit doing
        # --- this too much in the future
        if loop_keys is None:
            print("Currently must set loop_keys to file paths of HDF5 files to loop over")
            raise ValueError

        if not all(isinstance(x, str) for x in key_list):
            print("Currently cannot handle keys that are not all strings")
            print("doing this would require some more heavy lifting. It is on my to-do list")
            print("to get rid of non-string dictionary keys in the analysis")
            raise ValueError

        # must prepend first key with / to load using deepdish in below
        if not '/' in key_list[0]:
            key_list[0] = '/' + key_list[0]

        key_list_string = '/'.join(key_list)
        all_data = [dd.io.load(fpath, key_list_string) for fpath in loop_keys]

        return all_data


def extract_nested_dict_asarray(dict, key_list, loop_keys = None, self_contained = True):
    """
    Wrapper on extract_nested_dict_aslist to return a numpy array instead.
    """

    return np.array(extract_nested_dict_aslist(dict, key_list,
                    loop_keys = loop_keys, self_contained = self_contained))

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
        #print index
        sample = names[index]
        #print sample
        s_to_remove = sample[2:] # assume longest string is 2 letter sym
        #print s_to_remove
        i = -1 * len(s_to_remove)
        filtered = [x[:i] for x in names]
        #print filtered
    else:
        filtered = [x for x in names]

    # get list of atomic numbers
    anum = [ asym_to_anum[x] for x in filtered ]

    # get sorting array
    sort = np.argsort(anum)

    return list(np.array(names)[sort])

def compute_stats(x, return_dict = False, acf = False):
    """

    OLD FUNCTION. HERE FOR BACKWARDS COMPATABILITY.

    Please use `compute_statistics`
    """

    if not return_dict: # for bakcwards compatability - will remove eventually
        return np.min(x), np.max(x), np.average(x), np.median(x), np.std(x)

    d = {}
    d['min'] = np.min(x)
    d['max'] = np.max(x)
    d['mean'] = np.average(x)
    d['median'] = np.median(x)
    d['std'] = np.std(x)
    d['decile_1'], d['Q1'], d['Q3'], d['decile_9'] = np.percentile(x, [0.1,0.25,0.75,0.9])
    d['inner_quartile_range'] = d['Q3'] - d['Q1']
    d['d9_d1_range']          = d['decile_9'] - d['decile_1']
    d['variance'] = np.var(x)

    if (len(x) > 0) and acf:
        d['acf'] = np.array(stattools.acf(x))
    else:
        d['acf'] = 0.0

    return d


def weighted_quantile(values, quantiles, weight=None, values_sorted=False):
    """
        As taken from:
        https://stackoverflow.com/questions/21844024/weighted-percentile-using-numpy

        Very close to np.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: np.array with data
    :param quantiles: array-like with many quantiles needed
    :param weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of initial array
    :return: np.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if weight is None:
        weight = np.ones(len(values))
    weight = np.array(weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), 'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        weight = weight[sorter]

    weighted_quantiles = np.cumsum(weight) - 0.5 * weight
    weighted_quantiles /= np.sum(weight)

    return np.interp(quantiles, weighted_quantiles, values)

def weighted_percentile(values, percentiles, weight=None, values_sorted=False):
    """
    See `weighted_quantile` for more information.
    """
    return weighted_quantile(values, percentiles/100.0, weight=weight, values_sorted=values_sorted)


def compute_statistics(values, weights=None, limits = [-np.inf, np.inf]):
    """
    Compute a set of statistics for a given set of values with optional
    weights and limits.

    Provide an empty set of values or "None" to see which statistics
    are computed

    Parameters
    ----------
    values    : 1D-array
                Array of values to compute. Supply an empty array or
                'None' to see what statistics are computed.
    weights   : 1D array (optional)
                Optional weights for values. Must be same length as values.
                Default : None
    limits    : list, tuple (optional)
                Minimum and maximum values to consider (2D list, tuple).
                Default considers all values. Default: [-np.inf, np.inf].

    Returns:
    --------
    stats     : dictionary
                Dictionary of computed statistics. Provide an empty 'values'
                to see what statistics are computed.
    """

    #
    # construct default dictionary. This needs to be updated if
    # new stats are computed below.
    #
    stats = {'limits' : [], 'number' : 0, 'median' : 0, 'mean' : 0,
             '1-sigma' : 0, '2-sigma' : 0, '3-sigma' : 0,
             '1-sigma_range' : 0, '2-sigma_range' : 0, '3-sigma_range' : 0,
             'IQR' : 0, 'inter_quartile_range' : 0, 'inter_decile_range' : 0,
             'mean_median_diff' : 0,
             'variance' : 0,
             'std': 0, 'min' : 0, 'max' : 0}
    #
    # percentiles to compute, with string formatter
    # to eliminate trailing zeros and decimals in percentile
    #
    # e.g. percentile_0.135  and percentile_10
    #
    percentiles_list = [0.135, 2.275, 10, 15.865, 25, 50, 75, 84.135, 90, 97.725, 99.865]
    perc_string = lambda x : ("%.3f"%(x)).rstrip("0").rstrip('.')
    for pval in percentiles_list:
        stats['percent_' + perc_string(pval)] = 0


    # don't throw an error if values is empty
    if values is None or len(values) == 0:
        return stats

    # select out based on limits
    mask = (values >= limits[0] ) * (values < limits[-1])
    _values = values[mask]

    nvalues = np.size(_values)
    if weights is None:
        _weights = np.ones(nvalues)
    else:

        if len(weights) != len(mask):
            print("Weights must have same size (%i) as `values` array (%i)"%(np.size(weights),np.size(values)))
            raise ValueError

        _weights = weights[mask]

    stats['limits'] = limits
    stats['number'] = nvalues


    stats['mean']      = np.average(_values, weights=_weights)
    stats['variance']  = np.sum(weights/np.sum(weights) * (_values - stats['mean'])**2)
    stats['std']       = np.sqrt(stats['variance'])

    for pval in percentiles_list:
        stats['percent_' + perc_string(pval)] =\
                                  weighted_percentile(_values,pval,_weights)

    stats['median'] = stats['percent_50']

    # copy relevant percentiles into 1, 2, and 3 sigma interval pairs
    stats['1-sigma'] = [stats['percent_15.865'], stats['percent_84.135']]
    stats['2-sigma'] = [stats['percent_2.275'] , stats['percent_97.725']]
    stats['3-sigma'] = [stats['percent_0.135'] , stats['percent_99.865']]

    #
    # derived statistics giving ranges
    #
    for k in ['1-sigma','2-sigma','3-sigma']:
        stats[k + '_range'] = stats[k][1] - stats[k][0]
    stats['IQR']                  = stats['percent_75'] - stats['percent_25']
    stats['inter_quartile_range'] = stats['IQR']
    stats['inter_decile_range']   = stats['percent_90'] - stats['percent_10']
    stats['mean_median_diff']     = stats['mean'] - stats['median']

    # min, max
    stats['min'] = np.min(_values)
    stats['max'] = np.max(_values)


    return stats

def binned_statistics(values_x, values_y, limits = None, nbins = 100,
                      bins = None, weights = None,
                      y_limits = [-np.inf, np.inf],
                      empty_value = np.nan,
                      return_binned_values = False):
    """
    Computes binned statistics of `values_y` as a function of `values_x` with optional
    weights.

    Parameters :
    -------------
    values_x, values_y : 1D array
            1D arrays of same length to compute statistics on. Bins `values_y` as
            a function of `values_x`.

    limits : list, tuple, optional
            Min and max values in x dimension to bin along. If not provided,
            defaults to min and max values of `values_x`.
            Ignored if `bins` provided. Default : None

    nbins   : int, optional
            Number of bins. Ignored if `bins` provided. Default : 100

    bins    : 1D array, optional
            Bins in `values_x` to use. Overrides `limits` and `nbins` arguements.
            Default : None

    weights : 1D array, optional
            Weights to use when computing statistics. Must be same size as
            `values_x`. Default : None

    y_limits : list, tuple, optional
            Limits for `values_y` below / above which these values will be
            ignored. Default is to accept all `values_y`. Default : [-np.inf,np.inf]

    empty_value : optional
             If a bin is empty, stats for that bin (except `number`)
             are set to this value. Default is `np.nan` to mask out
             when plotting.

    return_binned_values : bool, optional
            Return a list of the `values_y` in each bin as well in the
            output dictionary. Default : False.

    Returns:
    ---------
    """

    if bins is None:
        if limits is None:
            limits = [np.min(values_x), np.max(values_x)]
        bins = np.linspace(limits[0], limits[-1], nbins+1)
    else:
        limits = [bins[0], bins[-1]]
        nbins  = np.size(bins)


    mask = (values_x >= limits[0]) * (values_x < limits[-1])
    _values_x = values_x[mask]
    _values_y = values_y[mask]

    if weights is None:
        weights = np.ones(np.size(_values_x))
    else:
        _weights = weights[mask]

    # get kwargs of statistics that will be computed
    # to build dictionary of lists to house the binned stats
    empty_stats = compute_statistics(None)

    binned_stats = {'bins'        : bins,      # bin edges
                    'lbins'       : bins[:-1], # left bin edges
                    'rbins'       : bins[1:],  # right bin edges
                    'bins_mid'    : 0.5*(bins[:-1]+bins[1:]), # bin centers
                    'bins_diff'   : (bins[:-1]-bins[1:]),     # db
                    'bin_limits'  : limits,
                    'y_limits'    : y_limits,
                    'total_number': np.size(_values_x),
                    'binned_y'    : [],    # empty for now, optional to generate
                    'nbins'       : nbins}

    #
    # initialize the statistics (better than appending to a list)
    #

    # stats to skip and not save in binned statistic.
    # these are all lists / tuples and redundant with other
    # information
    skip_stats = ['limits','1-sigma','2-sigma','3-sigma']

    for k in empty_stats.keys():
        if k in skip_stats:
            continue
        binned_stats[k] = np.zeros(binned_stats['nbins']-1)

    if return_binned_values:
        binned_stats['binned_y'] = [None] * binned_stats['nbins']

    #
    # generate the masks for all of the points
    #
    bin_indexes = np.digitize(_values_x, bins, right = False)
    print(binned_stats['nbins'],np.min(bin_indexes),np.max(bin_indexes))
    for i in np.arange(binned_stats['nbins']-1):
        #
        # Generate mask and get y values in this bin
        #
        mask     = bin_indexes == i
        masked_y = _values_y[mask]

        #
        # compute and save all stats for this bin.
        # this function returns
        #
        stats = compute_statistics(masked_y, weights = weights[mask],
                                   limits = y_limits)
        if stats['number'] == 0:
            #
            # Set the stats to the empty value
            # (by default, np.nan to mask out when plotting)
            #
            for k in stats.keys():
                if k in skip_stats:
                    continue
                binned_stats[k][i] = empty_value
            stats['number'] = 0

        else:
            for k in stats.keys():
                if k in skip_stats:
                    continue
                # print(i,k,stats[k])
                binned_stats[k][i] = stats[k]

        # save masked y values if requested
        if return_binned_values:
            binned_stats['binned_y'][i] = masked_y


    return binned_stats

def compute_weighted_stats(x, w, return_dict = True):
    """

    OLD FUNCTION. HERE FOR BACKWARDS COMPATABILITY.

    Please use `compute_statistics`

    Returns a dictionary containing statistics information that can
    be easily tossed into a data set. Uses exclusively numpy-based routines
    to compute the relevant weighted statistics.
    """

    _x = 1.0 * x
    if hasattr(x,'value'):
        _x = x.value

    _w = 1.0 * w
    if hasattr(w,'value'):
        _w = w.value

    d = {}
    d['w-avg_mean'] = np.average(                _x, weights=_w) # not necessarily the actual mean we might be interested in
    d['mean']       = 1.0*d['w-avg_mean'] # for backwards compatability (but this is bad)
    d['variance']   = np.average( (_x-d['mean'])**2, weights=_w)
    d['std']        = np.sqrt(d['variance'])

    # no weighted quantiles in numpy - need to use defined function
    q             = weighted_quantile(_x, [0.1, 0.25, 0.5, 0.75, 0.9], weight=_w)
    d['decile_1']  = q[0] # decile 1
    d['Q1']        = q[1] # quartile 1
    d['median']    = q[2]
    d['Q3']        = q[3]
    d['decile_9']  = q[4]
    d['inner_quartile_range'] = d['Q3'] - d['Q1']
    d['d9_d1_range']        = d['decile_9'] - d['decile_1']

    d['min']      = np.min(_x)
    d['max']      = np.max(_x)

    if hasattr(x, 'value'):
        for k in d.keys():
            d[k] = d[k] * x.unit_quantity

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

def fraction_fields(ds):
    """
    Which non-species fraction fields are available to be computed?
    """

    list = ['AGB','Intermediate_Wind','Massive_Wind','SNII','SNIa','SNIa_SDS',
            'SNIa_HeRS','SNIa_sCh','PopIII','PopIII_PISNe','RProcess']

    available = []
    if ds.parameters['IndividualStarTrackAGBMetalDensity']:
        available.append('AGB')
    if ds.parameters['IndividualStarTrackWindDensity']:
        available.append('Intermediate_Wind')
        available.append('Massive_Wind')
    if ds.parameters['IndividualStarTrackSNMetalDensity']:
        available.append("SNII")
        available.append("SNIa")
        if ds.parameters['IndividualStarSNIaModel'] == 2:
            available.append("SNIa_SDS")
            available.append("SNIa_HeRS")
            available.append("SNIa_sCh")
    if ds.parameters["IndividualStarPopIIIFormation"]:
        available.append("PopIII")
        available.append("PopIII_PISNe")

    if ds.parameters["IndividualStarRProcessModel"]:
        available.append("RProcess")

    return available

def abundance_ratios_from_fields(fields, select_denom = None, select_num = None):
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

    # only select a few of these fields with the chosen denominators
    if not (select_denom is None):
        select_denom = ['_over_' + x for x in select_denom]
        fields = [x for x in fields if any([y in x for y in select_denom])]

    # only select a few of these fields with the chosen numerators
    if not (select_num is None):
        select_num = [x + '_over_' for x in select_num]
        fields = [x for x in fields if any([y in x for y in select_num])]

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
#    if 'C' in species and 'O' in species and 'Mg' in species:
#        species = species + ['alpha']


    denominators = sort_by_anum(['H', 'Mg', 'Fe', 'Na', 'Ni','N', 'O', 'C', 'Ba'])
    # denominators = denominators + ['alpha']

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
