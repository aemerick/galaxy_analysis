import numpy as np
from scipy.optimize import curve_fit, minimize
from scipy.special import erf
from copy import copy

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


class general_functions():

    def __init__(self, name, function = None):
        self.name  = name

        if function is not None:
            self._f    = function
        self.p0    = None
        self.popt  = None
        self.pcov  = None

        self.fit_in_logspace = False

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
            if self.fit_in_logspace:
                self.popt, self.pcov = curve_fit(self._logf, xdata, np.log10(ydata),
                                                 *args, **kwargs)

                #print self.name, self.popt, self._logf(xdata,*self.popt), np.log10(ydata)

            else:

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

            func = lambda x0 : CDF_distance(self._CDF(xdata, *x0), data_cdf)
            result       = minimize(func, self.p0, *args, **kwargs)

            self.popt = result.x
            self.pcov = None


        return self.popt, self.pcov

#    def assign_function(self, function):
#        self._f = function
#        return


class power_law(general_functions):

    def __init__(self):
        general_functions.__init__(self,'powerlaw')
        self.fit_in_logspace = True

        return

    def _f(self, x, k, a):
        return 10.0**(self._logf(x,k,a))

    def _logf(self, x, k, a):
        return -1.0 *k * np.log10(x) + np.log10(a)


    def _CDF(self, x, k, a, norm = False):

        CDF = np.zeros(np.size(x))

        CDF = a /((-k) + 1.0) * x**(-k + 1)

        if norm:
            CDF = CDF / np.max(CDF)
        return CDF

class truncated_powerlaw(general_functions):

    def __init__(self):
        general_functions.__init__(self, 'truncated_powerlaw')
        self.fit_in_logspace = True
        return

    def _f(self, x, k, a, x_t):
        fx = np.zeros(np.size(x))

        fx[ x < x_t] = 0.0
        fx[ x > x_t] = 10.0**(self._logf(x[x>x_t], k, a, x_t))

        return fx

    def _logf(self, x, k, a, x_t):
        fx = np.zeros(np.size(x))

        fx[x < x_t] = -99999
        fx[x > x_t] = -1.0 * k * np.log10(x[x>x_t]) + np.log10(a)
        return fx

    def _CDF(self, x, k, a, x_t, norm = False):

        CDF = np.zeros(np.size(x))
        CDF[x < x_t] = 0.0
        CDF[x > x_t] = power_law._CDF(x[x > x_t], k, a)

        return CDF

class gaussian(general_functions):
    def __init__(self, fix_mean = None):
        general_functions.__init__(self, 'gaussian')
        self._mu = fix_mean
        return


    def _CDF(self, x, mu, sigma):
        if not (self._mu is None):
            mu = self._mu

        CDF = 0.5 * (1.0 + erf( x / np.sqrt(2.0)))
        return CDF

    def _f(self, x, mu, sigma):
        if not (self._mu is None):
            mu = self.mu

        fx = (1.0 / (np.sqrt(2.0 * np.pi)*sigma) *\
               np.exp(- (x - sigma)*(x - sigma) / (2.0 * sigma * sigma)))

        return fx

    def print_parameters(self):
        print "mu (mean of data) and sigma (standard deviation of data)"

class lognormal_powerlaw(general_functions):
    """
    Following Chen, Burkhart, Goodman, and Collins 2018
    """

    def __init__(self):
        """
        Fit a lognormal + power law tail PDF
        """

        general_functions.__init__(self,'lognormal_powerlaw')
        self.fit_in_logspace = True
        return

#    def _f(self, x, mu, alpha, p_o, N):
    def _f(self, x, mu, alpha, sigma):
        N = 1.0
        fx = np.zeros(np.size(x))
        xlog  = np.log(x)
        self.xt    = self.full_mean * np.exp(self.st)

#        sigma = np.sqrt(-0.5*mu) # try this

        s     = np.log(x / self.full_mean)  #/ self.full_mean)
        mu    = mu - np.log( self.full_mean)

#        sigma = np.sqrt(-0.5 * mu)

        p_o = 1.0/(np.sqrt(2.0*np.pi)*sigma*self.xt) * np.exp(-1.0*(self.st-mu)**2 / (2.0*sigma*sigma) + alpha*self.st)

        fx[ s < self.st] = N / (np.sqrt(2.0*np.pi)*sigma*x[s<self.st]) * np.exp(-1.0 * (s[s<self.st] - mu)**2 / (2.0*sigma*sigma))
        fx[ s > self.st] = N * p_o * np.exp(-alpha * s[s>self.st])

        self.N     = N
        self.sigma = sigma
        self.p_o   = p_o

        return fx

    def _logf(self, x, *args):
        fvals = np.log10(self._f(x, *args))

        return fvals

    def fit_function(self, xdata, ydata, method = 'curve_fit', data_cdf = None, *args, **kwargs):
        """
        Fit function to data. By default, this uses the scipy method 'curve_fit', but the 
        'KS' can be provided as method to minimize the distance between the CDFs of the two 
        functions. If 'KS' is used, the data CDF is integrated using numpy trapz, but it would
        possibly be better to provide the CDF computed separately, using `data_cdf' argument.
        """
        if 'p0' in kwargs.keys():
            self.p0 = kwargs['p0']

        min_error = np.inf
        all_xt = np.logspace( np.log10(xdata[np.argmax(ydata)]), np.log10(np.max(xdata)), np.size(xdata)*2)
        for xt in all_xt:

            self.st = np.log( xt / self.full_mean)
            try:
                if self.fit_in_logspace:
                    self.popt, self.pcov = curve_fit(self._logf, xdata, np.log10(ydata),
                                                 *args, **kwargs)

                    #print self.name, self.popt, self._logf(xdata,*self.popt), np.log10(ydata)
                else:
                    self.popt, self.pcov = curve_fit(self._f, xdata, ydata,
                                                         *args, **kwargs)
            except:
                continue

            y_fit = self._f(xdata, *self.popt)
            error = np.sum(  (y_fit - ydata)**2 / ydata )

            if error < min_error:
                optimal_st    = 1.0*self.st
                optimal_popt  = copy(self.popt)
                optimal_pcov  = copy(self.pcov)
                min_error  = error

        self.popt = optimal_popt
        self.pcov = optimal_pcov
        self.st   = optimal_st

        return self.popt, self.pcov

class lognormal(general_functions):

    def __init__(self, fix_mean = None):
        general_functions.__init__(self, 'lognormal')
        self._mu = fix_mean
        self.fit_in_logspace = True
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

    def _logf(self, x, *args):
        return np.log10( self._f(x, *args))

    def print_parameters(self):
        print "mu (mean of logged data) and sigma (standard deviaion of logged data)"


by_name = {'log-normal' : lognormal,
           'powerlaw'   : power_law,
           'truncated_powerlaw' : truncated_powerlaw,
           'gaussian'   : gaussian,
           'lognormal_powerlaw' : lognormal_powerlaw}

