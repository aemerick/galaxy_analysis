import numpy as np
from scipy.optimize import curve_fit, minimize
from scipy.special import erf

def CDF_distance(cdf1, cdf2):
    """
    Very simple function to return distance between two CDFs.
    """ 
    return np.max( np.abs(cdf2 - cdf1))

def compute_cdf(x, y):
    """
    Given x and y data points, compute the CDF
    """
    
    
    
    return CDF


class fit_function():

    def __init__(self, name, function = None):
        self.name  = name

        if function is not None:
            self._f    = function
        self.p0    = None
        self.popt  = None
        self.pcov  = None

        return

    def fit_function(self, xdata, ydata, method = 'curve_fit', data_cdf = None, *args, **kwargs):
        """
        Fit function to data. By default, this uses the scipy method 'curve_fit', but the 
        'KS' can be provided as method to minimize the distance between the CDFs of the two 
        functions. If 'KS' is used, the data CDF is integrated using numpy trapz, but it would
        possibly be better to provide the CDF computed separately, using `data_cdf' argument.
        """
        if 'p0' in kwargs.keys():
            self.p0 = kwargs['p0']

        if method == 'curve_fit': # use scipy curve fitting
            self.popt, self.pcov = curve_fit(self._f, xdata, ydata,
                                                     *args, **kwargs)
        elif method == 'KS':
            # optimize the fit by minimizing distance between CDF
            if not 'p0' in kwargs.keys():
                print "Must supply initial guess if using KS method"
                raise ValueError
            else:
                del kwargs['p0']

            # compute the CDF from the data to fit
            if data_cdf is None:
                data_cdf = compute_cdf(xdata, ydata)
                        
            fit_function = lambda x0 : CDF_distance(self._CDF(xdata, *x0), data_cdf)
            result       = minimize(fit_function, self.p0, *args, **kwargs)

            self.popt = result.x
            self.pcov = None

        return self.popt, self.pcov

    def assign_function(self, function):
        self._f = function
        return




class lognormal(fit_function):

    def __init__(self, fix_mean = None):
        fit_function.__init__(self, 'lognormal')
        
        self._mu = fix_mean
        
        return

    def _CDF(self, x, mu, sigma):
        
        if not (self._mu is None):
            mu = self._mu
        
        CDF = 0.5 * (1.0 + erf( (np.log(x) - mu)/(sigma*np.sqrt(2.0))))
        
        
        return CDF
    
    def _f(self, x, mu, sigma):
        """
        Actual Function
        """
        
        if not (self._mu is None):
            mu = self.mu
        
        fx  = (1.0 / (x * sigma * np.sqrt(2.0*np.pi)))
        fx *= np.exp( -1.0 * (np.log(x) - mu)**2 / (2.0 * sigma * sigma))

        return fx

    def print_parameters(self):
        print "mu (mean of logged data) and sigma (standard deviaion of logged data)"
