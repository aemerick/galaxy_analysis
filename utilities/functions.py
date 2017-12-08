import numpy as np
from scipy.optimize import curve_fit


class fit_function():

    def __init__(self, name, function = None):
        self.name  = name

        if function is not None:
            self._f    = function
        self.p0    = None
        self.popt  = None
        self.pcov  = None

        return

    def fit_function(self, xdata, ydata, *args, **kwargs):
        if 'p0' in kwargs.keys():
            self.p0 = p0

        self.popt, self.pcov = curve_fit(self._f, xdata, ydata,
                                                 *args, **kwargs)

        return self.popt, self.pcov

    def assign_function(self, function):
        self._f = function
        return




class lognormal(fit_function):

    def __init__(self):
        fit_function.__init__(self, 'lognormal')
        return

    def _f(self, x, mu, sigma):
        """
        Actual Function
        """
        fx  = (1.0 / (x * sigma * np.sqrt(2.0*np.pi)))
        fx *= np.exp( -1.0 * (np.log(x) - mu)**2 / (2.0 * sigma * sigma))

        return fx

    def print_parameters(self):
        print "mu (mean of logged data) and sigma (standard deviaion of logged data)"
