from galaxy_analysis.plot.plot_styles import *
import yt
import numpy as np
from copy import copy
import deepdish as dd
from scipy.optimize import brentq
from galaxy_analysis.utilities import utilities

from galaxy_analysis.utilities import functions

import matplotlib.pyplot as plt

from scipy import integrate
from scipy.stats import distributions  # may not need
from scipy import stats

def load_distribution_data(ldata, dsname, phase, field, centers = None):
    """
    Extracts some information from the dataset
    """
    rdata  = ldata[dsname][phase]['mass_fraction'][field]
    y      = ldata[dsname][phase]['mass_fraction'][field]['hist']
    mean   = ldata[dsname][phase]['mass_fraction'][field]['mean']
    std    = ldata[dsname][phase]['mass_fraction'][field]['std']
    q1     = ldata[dsname][phase]['mass_fraction'][field]['Q1']
    q3     = ldata[dsname][phase]['mass_fraction'][field]['Q3']
    q10    = ldata[dsname][phase]['mass_fraction'][field]['decile_1']
    q90    = ldata[dsname][phase]['mass_fraction'][field]['decile_9']

    if rdata['Q3'] is None or rdata['Q1'] is None:
        rdata['iqr']           = None
        rdata['q90_q10_range'] = None
    else:
        rdata['iqr']    = np.log10(q3) - np.log10(q1) # ldata[dsname][phase]['mass_fraction'][field]['inner_quartile_range']
        rdata['q90_q10_range'] = np.log10(q90) - np.log10(q10)

    rdata['label']  = ldata[dsname]['general']['Time']
    rdata['median'] = np.log10(np.interp(0.5, np.cumsum(y)/(1.0*np.sum(y)), centers))

    return rdata


def fit_multifunction_PDF(bins, y, data):
    """
    Fit lognormal + power law
    """

    success = {}
    rdict   = {}

    centers = 0.5 * (bins[1:] + bins[:-1])

    def _error(yfitvals, yvals):
#        return np.sum( np.abs(yvals-yfitvals)*yvals)/np.sum(yvals)
        return np.sum( np.abs(yvals[yvals>0] - yfitvals[yvals>0])**2/yvals[yvals>0] )
#        return np.sum( np.abs(yvals[yvals>0]-yfitvals[yvals>0])   ) / (np.sum(1.0/yvals[yvals>0]))

    try:
    #if True:
        lognormal_alone = fit_PDF(bins*1.0, y*1.0, data = data, function_to_fit = 'log-normal')
        success['lognormal'] = True
        lognormal_alone['error'] = _error( lognormal_alone['fit_function']._f(centers, *lognormal_alone['popt']) , lognormal_alone['norm_y'])
        rdict['lognormal']   = lognormal_alone
    except:
        success['lognormal'] = False

    try:
#    if True:
        ln_pl = fit_PDF(bins*1.0, y*1.0, data=data, function_to_fit = 'lognormal_powerlaw')
        success['lognormal_powerlaw'] = True
        ln_pl['error'] = _error(ln_pl['fit_function']._f(centers, *ln_pl['popt']) , ln_pl['norm_y'])
        rdict['lognormal_powerlaw'] = ln_pl
    except:
        success['lognormal_powerlaw'] = False

    powerlaw_alone  = fit_PDF(bins, y, data = data, function_to_fit = 'powerlaw')
    success['powerlaw'] = True
    powerlaw_alone['error'] = _error(powerlaw_alone['fit_function']._f( centers[centers>centers[np.argmax(powerlaw_alone['norm_y'])]] , *powerlaw_alone['popt']) ,
                                         powerlaw_alone['norm_y'][ centers > centers[np.argmax(powerlaw_alone['norm_y'])]]  )

    rdict['powerlaw'] = powerlaw_alone

#    try:
    truncated_powerlaw  = fit_PDF(bins, y, data = data, function_to_fit = 'truncated_powerlaw')
    success['truncated_powerlaw'] = True
    truncated_powerlaw['error'] = _error(truncated_powerlaw['fit_function']._f( centers , *truncated_powerlaw['popt']) ,
                                         truncated_powerlaw['norm_y'])

    rdict['truncated_powerlaw'] = truncated_powerlaw

    if all([not success[k] for k in success.keys()]):
        print "Cannot find a fit"
        raise ValueError

    min_error = np.inf
    for k in success.keys():

        if success[k]:
            rdict[k]['name'] = k
            if rdict[k]['error'] < min_error:
                min_error = rdict[k]['error']
                min_key   = k


#    if 'lognormal_powerlaw' in rdict.keys():
#        min_key = 'lognormal_powerlaw'

    return rdict[min_key]

def fit_PDF(bins, y, data = None, function_to_fit = None, p0 = None, bounds = None,
                     fit_method = 'curve_fit', remove_residual = None):
    """
    fit_function is either a function object from utilities.functions OR
    the name of one of these functions
    """

    centers = 0.5 * (bins[1:] + bins[:-1])
    binsize = (bins[1:] - bins[:-1])

    norm_y = y / binsize
    if not (remove_residual is None):
        norm_y = norm_y - remove_residual

    # if data is provided we can do some cool things
    if p0 is None:
        if not (data is None):
            if function_to_fit == 'log-normal' or function_to_fit == 'lognormal_powerlaw':
                u_guess   = np.log( data['mean'] / (np.sqrt(1.0 + data['std']**2 / (data['mean']**2))))
                std_guess = np.sqrt( np.log(1.0 + data['std']**2 / data['mean']**2))

                p0 = [u_guess, std_guess]
                bounds = ( [p0[0]*1.0 - 10, 0.0], [p0[0]*1.0+10, 10.0])


            if function_to_fit == 'lognormal_powerlaw':
#                p0 = [p0[0],p0[1],2,1]
#                bounds = ( [bounds[0][0], bounds[0][1], 1, 0.0], [bounds[1][0],bounds[1][1],np.inf,np.inf])
#                p0 = [u_guess - 2, 2.1, 1.0, 1.0]
#                bounds = ([u_guess - 20, 1.0, 0.0, 0.0],[ np.min([u_guess+20,0]), 10, np.inf, np.inf])
                 p0 = [u_guess - 3, 2.0, 1.0, 0.1]
                 bounds = ( [p0[0] - 5, 0.1, 0.0, 0.0], [u_guess+3, 10.0, 10.0, np.inf])

#                bounds = (  [p0[0] - np.log(100), 0.0], [p0[0] + np.log(100), 8] )
#                 p0 = [p0[0], p0[1], 1.0]
#                 bounds = ( [bounds[0][0],bounds[0][1], 0.00000001,], [bounds[1][0],bounds[1][1], 10.0])

        if function_to_fit == 'powerlaw' or function_to_fit == 'truncated_powerlaw':
            p0     = [2, 1.0E-5]
            bounds = ( [1,np.inf], [0.0, np.inf] )
            bounds = ( [1,0], [np.inf,1.0])

        if function_to_fit == 'truncated_powerlaw':
            p0      = [p0[0]*1,p0[1]*1, centers[np.argmax(norm_y)-1] ]
            bounds  = ([bounds[0][0]*1, 1*bounds[0][1], centers[np.argmax(norm_y)-10]], [1*bounds[1][0],1*bounds[1][1],centers[np.argmax(norm_y)+20]])



    # choose the region to fit over
    #   for all, this will be where there is data

    selection = (norm_y > 0)
    if function_to_fit == 'log-normal':
        selection = selection * (centers > 10.0**(np.log10(centers[np.argmax(norm_y)]) - 1.0)  )*\
                                (centers < 10.0**(np.log10(centers[np.argmax(norm_y)]) + 1.0)  )
    elif function_to_fit == 'powerlaw':
        selection = selection * ( centers > centers[np.argmax(norm_y)] )  # Fit everything to the right of the peak


    x_to_fit = centers[selection]
    y_to_fit = norm_y[selection]

    # set fitting function
    if isinstance(function_to_fit, str):
        function_to_fit= functions.by_name[function_to_fit]()


    if function_to_fit.name == 'lognormal_powerlaw':
        function_to_fit.full_mean = data['mean']

    # If we are using KS test to fit, then need to compute the CDF from the data
    if fit_method == 'KS':
        data_cdf  = np.cumsum(y[selection])
    else:
        data_cdf  = None

    # Fit -
    popt, pcov = function_to_fit.fit_function(x_to_fit, y_to_fit, p0 = p0,
                                              bounds = bounds, method = fit_method,
                                              data_cdf = data_cdf)

    if function_to_fit.name == 'lognormal_powerlaw':
        print "----", popt


    fit_dictionary = {'popt' : copy(popt), 'pcov': copy(pcov),
                      'fit_function' : function_to_fit, 'norm_y' : norm_y*1.0,
                      'fit_x' : x_to_fit*1.0, 'fit_y' : y_to_fit*1.0,
                      'fit_result' : lambda xx : function_to_fit._f(xx, *popt) }

    return fit_dictionary

colors = ['C' + str(i) for i in np.arange(9)]
lss     = ['-','--']

def plot_all_elements(data, dsname, phase, elements = None, **kwargs):

    if elements is None:
        elements = ['Ba','Y','As','Sr','Mn','Na','Ca','N','Ni','Mg','S','Si','Fe','C','O']

    fig, ax = plt.subplots(2)
    fig.set_size_inches(24,16)

    bins    = data[dsname][phase]['mass_fraction']['bins']
    centers = 0.5 * (bins[1:]+ bins[:-1])


    ci = li = 0
    for i, e in enumerate(elements):
        field = e + '_Fraction'

        ds_data = load_distribution_data(data, dsname, phase, field, centers = centers) # subselect

        fit_dict = fit_multifunction_PDF(1.0*bins, 1.0*ds_data['hist'], ds_data)
                   #fit_PDF(bins, ds_data['hist'], data = ds_data, **kwargs)


        plot_histogram(ax[0], np.log10(bins), fit_dict['norm_y']/np.max(fit_dict['norm_y']),
                                    color = colors[ci], ls = lss[li], lw = line_width)

        ax[0].plot(np.log10(centers),
                   #np.log10(fit_dict['fit_x']), 
                   fit_dict['fit_result'](centers) / np.max(fit_dict['norm_y']),
                                      lw = line_width, ls = lss[li], color = colors[ci])


        plot_histogram(ax[1],np.log10(bins),
                   #np.log10(fit_dict['fit_x']), 
                   fit_dict['norm_y']/np.max(fit_dict['norm_y']) - fit_dict['fit_result'](centers)/np.max(fit_dict['norm_y']),
                                      lw = line_width, ls = lss[li], color = colors[ci])



        print e, fit_dict['name'], fit_dict['popt']
        ci = ci + 1
        if ci >= np.size(colors):
            ci = 0
            li = li + 1

    for i in [0,1]:
        ax[i].set_xlim(-14, -1.5)
        ax[i].set_ylim(1.0E-5, 9.0)
        ax[i].semilogy()

    ax[0].set_ylabel('Peak Normalized PDF from Data')
    ax[1].set_ylabel('Peak Normalized PDF from Fit')

    fig.subplots_adjust(hspace = 0)
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible = False)
    plt.minorticks_on()

    fig.savefig(dsname + '_' + phase + '_all_elements_fit.png')
    plt.close()

    return


if __name__ == '__main__':

    data = {"DD0400" : dd.io.load('gas_abundances_5Myr.h5', "/DD0400") }
    plot_all_elements(data, 'DD0400', 'Molecular', function_to_fit = 'lognormal_powerlaw') #'lognormal_powerlaw')
    plot_all_elements(data, 'DD0400', 'CNM', function_to_fit = 'lognormal_powerlaw') #'lognormal_powerlaw')
    plot_all_elements(data, 'DD0400', 'WNM', function_to_fit = 'lognormal_powerlaw') #'lognormal_powerlaw')
    plot_all_elements(data, 'DD0400', 'WIM', function_to_fit = 'lognormal_powerlaw') #'lognormal_powerlaw')
    plot_all_elements(data, 'DD0400', 'HIM', function_to_fit = 'lognormal_powerlaw') #'lognormal_powerlaw')
    plot_all_elements(data, 'DD0400', 'Disk', function_to_fit = 'lognormal_powerlaw') #'lognormal_powerlaw')





# load_distribution_data(gas_data, "DD0500", "HIM", "O_Fraction", centers = centers)


