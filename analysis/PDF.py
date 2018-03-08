from galaxy_analysis.plot.plot_styles import *
import yt
import numpy as np
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


    try:
        lognormal_alone = fit_PDF(bins, y, data = data, function_to_fit = 'log-normal')
        success['lognormal'] = True
        lognormal_alone['error'] = np.sum(np.abs(lognormal_alone['fit_result'](centers[lognormal_alone['norm_y']>0]) - lognormal_alone['norm_y'][lognormal_alone['norm_y']>0]))
        rdict['lognormal']   = lognormal_alone
    except:
        success['lognormal'] = False

    try:
        powerlaw_alone  = fit_PDF(bins, y, data = data, function_to_fit = 'powerlaw')
        success['powerlaw'] = True
        powerlaw_alone['error'] = np.sum(np.abs(\
                                  powerlaw_alone['fit_result']( centers[centers>centers[np.argmax(powerlaw_alone['norm_y'])]] ) -\
                                  powerlaw_alone['norm_y'][centers>centers[np.argmax(powerlaw_alone['norm_y'])]]))

        rdict['powerlaw'] = powerlaw_alone
    except:
        success['powerlaw'] = False

    try:
        x_power_law  = centers[np.argmax(lognormal_alone['norm_y'])]
        residual = lognormal_alone['norm_y']*1.0
        residual[centers < x_power_law]     = 0.0
        remove_residual = lognormal_alone['norm_y'] - residual
        powerlaw_fit2 = fit_PDF(bins , y, data = data, function_to_fit = 'powerlaw', remove_residual = remove_residual)

        full_result = np.zeros(np.size(centers))
        full_result[centers>x_power_law] = powerlaw_fit2['fit_function']._f(centers[centers>x_power_law], *powerlaw_fit2['popt'])
        residual = lognormal_alone['norm_y'] - full_result
        residual[ residual < 0] = 0.0
        
        remove_residual = lognormal_alone['norm_y'] - residual

        ### Try a double fit
        lognormal_fit2 = fit_PDF(bins, y, data=data, function_to_fit = 'log-normal', remove_residual = remove_residual)

        def final_function2(x, x_power_law):
            result = np.zeros(np.size(x))

            result = lognormal_fit2['fit_function']._f(x,*lognormal_fit2['popt'])
            result[x>x_power_law] += powerlaw_fit2['fit_function']._f(x[x>x_power_law],*powerlaw_fit2['popt'])

            return result


        fit_dict = {'popt' : [lognormal_fit2['popt'], powerlaw_fit2['popt']],
                    'pcov' : [lognormal_fit2['pcov'], powerlaw_fit2['pcov']],
                    'norm_y' : lognormal_alone['norm_y'], 'norm_y_resid' : powerlaw_fit2['norm_y'],
                    'fit_x' : lognormal_fit2['fit_x'], 'fit_y' : lognormal_fit2['fit_y'],
                    'fit_result' : lambda x : final_function2(x, x_power_law)}

        fit_dict['error'] = np.sum(np.abs(fit_dict['fit_result'](centers[fit_dict['norm_y']>0]) - fit_dict['norm_y'][fit_dict['norm_y']>0]))


        rdict['lognormal_powerlaw2'] = fit_dict
        success['lognormal_powerlaw2'] = True
    except:
        print "Double function fitting failing ------ try 2"
        success['lognormal_powerlaw2'] = False


    try:
        ### Try a double fit
        lognormal_fit = fit_PDF(bins, y, data=data, function_to_fit = 'log-normal', remove_residual = None)

        full_result = lognormal_fit['fit_function']._f(centers, *lognormal_fit['popt'])
        residual = lognormal_fit['norm_y'] - full_result
        residual[ residual < 0] = 0.0
        residual[ centers < centers[np.argmax(full_result)] ] = 0.0

        remove_residual = (full_result - residual)


        x_power_law     = centers[np.argmax(residual)-5]

        for i in np.arange(np.size(residual)):
            if all( residual[i:i+5] > 0 ):
                x_power_law = centers[i]
                break
        if i == np.size(residual) -1:
            print '-------------------------------------------------'
        # compute the residual
        powerlaw_fit  = fit_PDF(bins, y, data=data, function_to_fit = 'powerlaw',
                                  remove_residual = remove_residual)

        def final_function(x, x_power_law):
            result = np.zeros(np.size(x))

            result = lognormal_fit['fit_function']._f(x,*lognormal_fit['popt'])
            result[x>x_power_law] += powerlaw_fit['fit_function']._f(x[x>x_power_law],*powerlaw_fit['popt'])

            return result


        fit_dict = {'popt' : [lognormal_fit['popt'], powerlaw_fit['popt']],
                    'pcov' : [lognormal_fit['pcov'], powerlaw_fit['pcov']],
                    'norm_y' : lognormal_fit['norm_y'], 'norm_y_resid' : powerlaw_fit['norm_y'],
                    'fit_x' : lognormal_fit['fit_x'], 'fit_y' : lognormal_fit['fit_y'],
                    'fit_result' : lambda x : final_function(x, x_power_law)}

        fit_dict['error'] = np.sum(np.abs(fit_dict['fit_result'](centers[fit_dict['norm_y']>0]) - fit_dict['norm_y'][fit_dict['norm_y']>0]))

        rdict['lognormal_powerlaw'] = fit_dict
        success['lognormal_powerlaw'] = True
    except:
        print "Double function fitting failing"
        success['lognormal_powerlaw'] = False

    if all([not success[k] for k in success.keys()]):
        print "Cannot find a fit"
        raise ValueError

    min_error = np.inf
    for k in success.keys():
        if success[k]:

            if rdict[k]['error'] < min_error:
                min_error = rdict[k]['error']
                min_key   = k

    return rdict[min_key]

def fit_PDF(bins, y, data = None, function_to_fit = None, p0 = None, bounds = None,
                     fit_method = 'curve_fit', remove_residual = None):
    """
    fit_function is either a function object from utilities.functions OR
    the name of one of these functions
    """

    # if data is provided we can do some cool things
    if p0 is None:
        if not (data is None):
            if function_to_fit == 'log-normal' or function_to_fit == 'lognormal_powerlaw':
                u_guess   = np.log( data['mean'] / (np.sqrt(1.0 + data['std']**2 / (data['mean']**2))))
                std_guess = np.sqrt( np.log(1.0 + data['std']**2 / data['mean']**2))
                p0 = [u_guess, std_guess]

                bounds = ( [p0[0] - 20, p0[0] + 20], [-np.inf,np.inf])
                bounds = ( [p0[0] - 20, 0.0], [p0[0]+20, 10.0])


            if function_to_fit == 'lognormal_powerlaw':
                p0 = [p0[0],p0[1],2,1]
                bounds = ( [bounds[0][0], bounds[0][1], 1, 0.0], [bounds[1][0],bounds[1][1],np.inf,np.inf])

        if function_to_fit == 'powerlaw':
            p0     = [2, 1.0E-5]
            bounds = ( [1,np.inf], [0.0, np.inf] )
            bounds = ( [1,0], [np.inf,1.0])


    centers = 0.5 * (bins[1:] + bins[:-1])
    binsize = (bins[1:] - bins[:-1])


    # choose the region to fit over
    #   for all, this will be where there is data

    norm_y = y / binsize
    if not (remove_residual is None):
        norm_y = norm_y - remove_residual

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
        if function_to_fit == 'lognormal_powerlaw':
            function_to_fit= functions.by_name[function_to_fit]( centers[np.argmax(norm_y)] )  # initialize function by name
        else:
            function_to_fit= functions.by_name[function_to_fit]()

    # If we are using KS test to fit, then need to compute the CDF from the data
    if fit_method == 'KS':
        data_cdf  = np.cumsum(y[selection])
    else:
        data_cdf  = None

    # Fit -
    popt, pcov = function_to_fit.fit_function(x_to_fit, y_to_fit, p0 = p0,
                                              bounds = bounds, method = fit_method,
                                              data_cdf = data_cdf)

    fit_dictionary = {'popt' : popt, 'pcov': pcov,
                      'fit_function' : function_to_fit, 'norm_y' : norm_y,
                      'fit_x' : x_to_fit, 'fit_y' : y_to_fit, 
                      'fit_result' : lambda x : function_to_fit._f(x, *popt) }

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

        fit_dict = fit_multifunction_PDF(bins, ds_data['hist'], ds_data)
                   #fit_PDF(bins, ds_data['hist'], data = ds_data, **kwargs)


        plot_histogram(ax[0], np.log10(bins), fit_dict['norm_y']/np.max(fit_dict['norm_y']),
                                    color = colors[ci], ls = lss[li], lw = line_width)

        ax[0].plot(np.log10(centers),
                   #np.log10(fit_dict['fit_x']), 
                   fit_dict['fit_result'](centers)/np.max(fit_dict['norm_y']),
                                      lw = line_width, ls = lss[li], color = colors[ci])
        plot_histogram(ax[1],np.log10(bins),
                   #np.log10(fit_dict['fit_x']), 
                   fit_dict['norm_y']/np.max(fit_dict['norm_y']) - fit_dict['fit_result'](centers)/np.max(fit_dict['norm_y']),
                                      lw = line_width, ls = lss[li], color = colors[ci])



        print e, fit_dict['popt']
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

