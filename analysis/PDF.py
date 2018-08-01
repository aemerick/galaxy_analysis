from galaxy_analysis.plot.plot_styles import *
rc('font', size=22)#, ftype=42)

import yt
import numpy as np
from copy import copy
import sys
import deepdish as dd
from scipy.optimize import brentq
from galaxy_analysis.utilities import utilities
import h5py
from galaxy_analysis.utilities import functions

import matplotlib.pyplot as plt

from scipy import integrate
from scipy.stats import distributions  # may not need
from scipy import stats

def gather_time_series(datafile, dsarray, phase, field, centers = None):

    if 'over' in field:
        cname = 'abins'
    else:
        cname = 'bins'

    i = 0
#    ldata = {dsarray[0]: dd.io.load(datafile, "/" + dsarray[0])}
    data  = load_distribution_data(datafile, dsarray[0], phase, field, centers = cname)
    time_series = {}
    for k in data.keys():
        if k == 'hist':
            continue
        time_series[k] = np.zeros(np.size(dsarray))

    for dsname in dsarray:
#        ldata = {dsname : dd.io.load(gas_file, "/" + dsname)}

        data  = load_distribution_data(datafile, dsname, phase, field, centers = cname)

        for k in data.keys():

            if np.size(data[k]) > 1: # cannot handle > 1D data for now
                continue

            time_series[k][i] = data[k]

        i = i + 1

    return time_series

def plot_normalized(ax, masses, field, xnorm = 'median'):
    return

def resolution_study(abundance_filename,
                     work_dir = './', output_dir = None,
                     comparison = None):

    phases = ['CNM','WNM','WIM','HIM']

    if output_dir is None:
        output_dir = work_dir

    if comparison is None:
        labels = {'3pcH2' : '3.6 pc', '6pcH2' : '7.2 pc' , 'Fiducial' : 'Fiducial'}
        lstyle = {'3pcH2' : '--', '6pcH2' : '-.' , 'Fiducial' : '-'}
        dirs   = {'3pcH2' : '../3pc_H2/abundances/', '6pcH2' : '../6pc_H2/abundances/', 'Fiducial' : ''}
        for k in dirs.keys():
            dirs[k] = work_dir + dirs[k]
    else:
        dirs = {}
        labels = {}
        lstyle = {}
        for k in comparison.keys():
            dirs[k] = work_dir + comparison[0]
            labels[k] = comparison[1]
            lstyle[k] = comparison[2]

    fig, ax = plt.subplots(2,2,sharex=True,sharey=True)
    fig.set_size_inches(12,12)
    fig.subplots_adjust(hspace=0.0,wspace=0.0)
    #
    # Plot individual panels showing:
    #   a time evolution of
    #
    time_data = {}
    for k in labels.keys():
        f = h5py.File(dirs[k] + abundance_filename,'r')
        dsarray = np.sort([str(x) for x in f.keys() if 'DD' in x])
        f.close()
        time_data[k] = {}
        time_data[k]['times'] = np.array([float(x.strip('DD')) for x in dsarray])
        time_data[k]['times'] = time_data[k]['times'] - time_data[k]['times'][0]
        for field in ['O_Fraction','Ba_Fraction']:
            time_data[k][field] = {}
            for phase in phases:
                time_data[k][field][phase] = gather_time_series(dirs[k] + abundance_filename,
                                                                dsarray, phase, field)

    for phase in phases:
        ls = '-'
        if phase == 'HIM':
            ls = ':'
        ax[(0,0)].plot(time_data['3pcH2']['times'],
                       time_data['3pcH2']['O_Fraction'][phase]['q90_q10_range'],
                        lw = 3, color = color_dict[phase], ls = ls)
        ax[(1,0)].plot(time_data['3pcH2']['times'],
                       time_data['3pcH2']['Ba_Fraction'][phase]['q90_q10_range'],
                       lw = 3, color = color_dict[phase], ls = ls)
        ax[(0,1)].plot(time_data['6pcH2']['times'],
                       time_data['6pcH2']['O_Fraction'][phase]['q90_q10_range'],
                         lw = 3, color = color_dict[phase], ls = ls)
        ax[(1,1)].plot(time_data['6pcH2']['times'],
                       time_data['6pcH2']['Ba_Fraction'][phase]['q90_q10_range'],
                         lw = 3, color = color_dict[phase], ls = ls)

    for a1 in ax:
        for a2 in a1:
            a2.set_xlim(0, 600)
            a2.set_ylim(0, 3)
            a2.minorticks_on()

    for i in [0,1]:
        ax[(i,0)].set_ylabel(r'log(90%) - log(10%) [dex]')
        ax[(1,i)].set_xlabel(r'Time (Myr)')

    x = 400
    y = 2.7
    size = 20
    ax[(0,0)].text(x, y, r'O  - 3.6 pc', color = 'black', size = size)
    ax[(1,0)].text(x, y, r'Ba - 3.6 pc', color = 'black', size = size)
    ax[(0,1)].text(x, y, r'O  - 7.2 pc', color = 'black', size = size)
    ax[(1,1)].text(x, y, r'Ba - 7.2 pc', color = 'black', size = size)

    #plt.tight_layout()
    fig.savefig('O_Ba_spread_resolution-comparison.png')
    plt.close()

    return

def element_by_element_panel(datafile, galaxy_file, dsname, show_fit = True,
                             elements = '3x3', outname = None):
    """
    Constructs a panel plot showing all phases for a single element
    on a single panel (so 4x4 plot: 15 elements + total metallicity field).
    Each plot is normalized to the median mass fraction (i.e. log(median) = 0)
    and centered on this value (so shareaxis is on). Each plot is normalized
    such that the integral over the disk PDF / max(disk PDF) = 1, and so that
    integral(phase)/integral(disk) gives mass fraction of that phase.

    Requires galaxy analysis output file (For mass fractions)
    """

    if outname is None:
        outname = galaxy_file.split('_galaxy')[0] + '_element_by_element.png'

    all_phases = ['CNM', 'WNM', 'WIM', 'HIM']

    if elements is None:
        elements   = get_element_list(datafile, dsname)
        xp,yp = 4,4
    elif elements == '3x3':
        elements   = ['O','C','Fe','Mg','N','Mn','Y','Sr','Ba']
        xp,yp = 3,3

    xsize, ysize = xp*6, yp*6

    fig, ax = plt.subplots(xp,yp, sharex=True, sharey=True)
    fig.set_size_inches(xsize,ysize)
    fig.subplots_adjust(hspace = 0.0, wspace = 0.0)


    gasdata = dd.io.load(galaxy_file, '/gas_meta_data/masses')
    masses = {}
    for k in gasdata.keys():
        try:
            masses[k] = gasdata[k]['Total']
        except:
          print k

#    disk_mass = masses['Disk']
#    for k in masses.keys():
#        masses[k] = masses[k] / disk_mass
# ----
    axi = 0
    axj = 0

    for element in elements:
        centers = 'bins'
        index = (axi,axj)

        field = element + '_Fraction'
        print element, axi, axj
        disk_data = load_distribution_data(datafile, dsname, 'Disk', field, centers = centers)
        disk_data['norm_y'] = disk_data['hist'] / disk_data['binsize']
        xbins = disk_data['bins']
        xcent = disk_data['centers']
        xnorm = disk_data['median']
        xplot  = np.log10(xbins) - xnorm
        xplotc = np.log10(xcent) - xnorm

        mass_norm = masses['Disk'] #- masses['Molecular']
        yplot     = disk_data['norm_y'] * masses['Disk'] / mass_norm
        ynorm     = np.max(yplot)
        yplot     = yplot / ynorm

        # plot disk here
        plot_histogram(ax[index], xplot, yplot, lw = line_width, color = 'black', label = 'Disk')
        if show_fit:
            fit = fit_multifunction_PDF(xbins, disk_data['hist'], disk_data)
            ax[index].plot(xplotc, fit['fit_result'](xcent)/ynorm*masses['Disk']/mass_norm,
                           lw = line_width, color = 'black', ls = '--')

        for phase in all_phases:
            data  = load_distribution_data(datafile, dsname, phase, field, centers = 'bins')
            data['norm_y'] = data['hist'] / disk_data['binsize']
            yplot = data['norm_y'] / ynorm * masses[phase] / mass_norm
            plot_histogram(ax[index], xplot, yplot, lw = line_width, color = color_dict[phase],
                           label=phase)

            if show_fit:
                fit = fit_multifunction_PDF(xbins, data['hist'], data)
                ax[index].plot( xplotc, fit['fit_result'](xcent) / ynorm * masses[phase]/mass_norm,
                                lw = line_width, color = color_dict[phase], ls = '--')

        ax[index].set_xlim(-3.8,3.3)
        ax[index].set_ylim(5.0E-7, 1.0)
        ax[index].semilogy()
        ax[index].plot([0,0],[1.0E-7,2.0], lw = 0.75 * line_width, ls = '-.', color = 'black')
        ax[index].plot([-1,-1],[1.0E-7,2.0], lw = 0.5 * line_width, ls = ':', color = 'black')
        ax[index].plot([1,1],[1.0E-7,2.0], lw = 0.5 * line_width, ls = ':', color = 'black')

        # ax[index].annotate()
        ax[index].minorticks_on()
        xytext = (-3.4,0.2)
        ax[index].text(xytext[0], xytext[1], element,  color = 'black', size = 40)
        axj = axj + 1
        if axj >= yp:
            axj = 0
            axi = axi + 1

    if (xp,yp) == (3,3):
        ax[(2,1)].set_xlabel(r'Median Normalized Mass Fraction')
        ax[(1,0)].set_ylabel(r'Peak-Normalized PDF')

#    for i in [0,1,2]:
#        ax[(2,i)].set_xticks([-3,-2,-1,0,1,2,3], [-3,-2,-1,0,1,2,3])

    ax[(1,0)].legend(loc = 'lower left', ncol = 1, fancybox=True, framealpha = 0.4)
    fig.savefig(outname)
    plt.close()

    return

def plot_time_evolution(datafile, dsarray, denominator = None):
    """
    Make a panel plot showing the evolution of:
        Median abundance           : top row
        log(median) - log(average) : middle row
        q90 - q10                  : bottom row
    For the following elements as columns (7):
        Ba, Sr, Fe, Mn, Mg, O, N

    Either as mass fractions (denominator = None) or as abundance
    ratios (e.g. denominator = "Fe")
    """

    all_phases = ['CNM','WNM','WIM']
    elements   = ['Ba','Sr','Fe','Mn','Ni','Mg','N','O']

    all_time_data = {}
    for e in elements:
        all_time_data[e] = {}


        if denominator is None:
            field = e + '_Fraction'
        elif e == denominator:
            continue # skip
        else:
            field = e + '_over_' + denominator

        for phase in all_phases:
            all_time_data[e][phase] =\
              gather_time_series(datafile, dsarray, phase, field)


    nrow = 3
    ncol = 8

    fig, ax = plt.subplots(nrow,ncol)
    fig.set_size_inches(ncol*6,nrow*6)

    axi = 0
    axj = 0

    times = np.arange(0, np.size(dsarray), 1) * 5

    color = {'Molecular' : 'C1', 'CNM' : 'C0', 'WNM' : 'C2', 'WIM':'C4','HIM' : 'C3'}

    for e in elements:

        if not (e == denominator):
            for phase in all_phases:

                median = all_time_data[e][phase]['median']
                mean   = all_time_data[e][phase]['mean']
                d9d1   = all_time_data[e][phase]['q90_q10_range']

                ax[(0,axj)].plot(times, median, lw = 3, color = color[phase], ls = '-')
                ax[(1,axj)].plot(times, median - mean, lw = 3, color = color[phase], ls = '-')
                ax[(2,axj)].plot(times, d9d1, lw = 3, color = color[phase], ls = '-')

        if not (denominator is None):
            if denominator == 'H':
                ax[(0,axj)].set_ylim(-8, 0)
                ax[(1,axj)].set_ylim(0, 1)
                ax[(2,axj)].set_ylim(0, 3)
            else:
                ax[(0,axj)].set_ylim(-4, 2)
                ax[(1,axj)].set_ylim(-0.2, 0.2)
                ax[(2,axj)].set_ylim(0, 2)

        for i in [0,1,2]:
            ax[(i,axj)].set_xlim(0, 500.0)
            xy = (0.9,0.1)
            ax[(i,axj)].annotate(e, xy,xy, textcoords='axes fraction')

        axj = axj + 1

    ax[(0,0)].set_ylabel('Median (dex)')
    ax[(1,0)].set_ylabel('Med - Mean (dex)')
    ax[(2,0)].set_ylabel('d9 - d1 (dex)')

    fig.subplots_adjust(hspace = 0)
    outname = 'time_evolution_'
    if denominator is None:
        outname += 'fraction'
    else:
        outname += denominator

    fig.savefig(outname + '.png')

    plt.close()
    
    return

def get_element_list(file_name, dsname):
    keys = dd.io.load(file_name, '/' + '/'.join([dsname,'CNM','mass_fraction'])).keys()

    x = [y.split('_over_')[0] for y in keys if '_over_Fe' in y]
    x = x + ['Fe']

    return utilities.sort_by_anum(x)

def load_distribution_data(file_name, dsname, phase, field, centers = None):
    """
    Extracts some information from the dataset
    """

    # load the data from file
    ldata = dd.io.load(file_name, '/' + '/'.join([dsname,phase,'mass_fraction',field]))

    rdata  = ldata
    y      = ldata['hist']
    mean   = ldata['mean']
    std    = ldata['std']
    q1     = ldata['Q1']
    q3     = ldata['Q3']
    q10    = ldata['decile_1']
    q90    = ldata['decile_9']

    if not (centers is None):

        if not (np.size(centers) == np.size(y)):
            if centers == 'bins' or centers == 'abins':
                bins = dd.io.load(file_name, '/' + '/'.join([dsname,phase,'mass_fraction',centers]))

            centers = 0.5 * (bins[1:] + bins[:-1])

    if rdata['Q3'] is None or rdata['Q1'] is None:
        rdata['iqr']           = None
        rdata['q90_q10_range'] = None
    else:

        if (q1 > 0 and q3 > 0 and q10 > 0 and q90 > 0) and (not 'over' in field):
            # safe to assume these are not logged data if all positive
            rdata['iqr']    = np.log10(q3) - np.log10(q1) # ldata['inner_quartile_range']
            rdata['q90_q10_range'] = np.log10(q90) - np.log10(q10)
        else:
            rdata['iqr'] = q3 - q1
            rdata['q90_q10_range'] = q90 - q10

    rdata['label']  = dd.io.load(file_name, '/' + '/'.join([dsname,'general','Time']))
    rdata['median'] = np.interp(0.5, np.cumsum(y)/(1.0*np.sum(y)), centers)
    if not 'over' in field:
        rdata['median'] = np.log10(rdata['median'])

    rdata['centers'] = centers
    rdata['binsize'] = bins[1:] - bins[:-1]
    rdata['bins'] = bins
    return rdata


def fit_multifunction_PDF(bins, y, data):
    """
    Fit lognormal + power law
    """

    success = {}
    rdict   = {}

    centers = 0.5 * (bins[1:] + bins[:-1])

    def _error(yfitvals, yvals):
        return np.sum( np.abs(yvals[yvals>0] - yfitvals[yvals>0])**2/yvals[yvals>0] )

    # lognormal
    try:
        lognormal_alone = fit_PDF(bins*1.0, y*1.0, data = data, function_to_fit = 'log-normal')
        success['lognormal'] = True
        lognormal_alone['error'] = _error( lognormal_alone['fit_function']._f(centers, *lognormal_alone['popt']) , lognormal_alone['norm_y'])
        rdict['lognormal']   = lognormal_alone
    except RuntimeError:
        success['lognormal'] = False

    # lognormal powerlaw
    ln_pl = fit_PDF(bins*1.0, y*1.0, data=data, function_to_fit = 'lognormal_powerlaw')
    success['lognormal_powerlaw'] = True
    ln_pl['error'] = _error(ln_pl['fit_function']._f(centers, *ln_pl['popt']) , ln_pl['norm_y'])
    rdict['lognormal_powerlaw'] = ln_pl

    if False: #try:
        gau_pl = fit_PDF(bins*1.0, y*1.0, data=data, function_to_fit = 'gaussian_powerlaw')
        success['gaussian_powerlaw'] = True
        gau_pl['error'] = _error(gau_pl['fit_function']._f(centers, *gau_pl['popt']), gau_pl['norm_y'])
        rdict['gaussian_powerlaw'] = gau_pl	
    #except RuntimeError:
    #    success['gaussian_powerlaw'] = False

    # powerlaw
    try:
        powerlaw_alone  = fit_PDF(bins, y, data = data, function_to_fit = 'powerlaw')
        success['powerlaw'] = True
#    powerlaw_alone['error'] = _error(powerlaw_alone['fit_function']._f( centers[centers>centers[np.argmax(powerlaw_alone['norm_y'])]] , *powerlaw_alone['popt']) ,
#                                         powerlaw_alone['norm_y'][ centers > centers[np.argmax(powerlaw_alone['norm_y'])]]  )
        powerlaw_alone['error'] = _error(powerlaw_alone['fit_function']._f( centers, *powerlaw_alone['popt']), powerlaw_alone['norm_y'])
        rdict['powerlaw'] = powerlaw_alone
    except RuntimeError:
        success['powerlaw'] = False

    # truncated powerlaw
    #truncated_powerlaw  = fit_PDF(bins, y, data = data, function_to_fit = 'truncated_powerlaw')
    #success['truncated_powerlaw'] = True
    #truncated_powerlaw['error'] = _error(truncated_powerlaw['fit_function']._f( centers , *truncated_powerlaw['popt']) ,
    #                                     truncated_powerlaw['norm_y'])
    #rdict['truncated_powerlaw'] = truncated_powerlaw

#    if all([not success[k] for k in success.keys()]):
#        print "Cannot find a fit"
#        raise ValueError

    min_error = np.inf
    for k in success.keys():

        if success[k]:
            rdict[k]['name'] = k
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
                bounds = ( [p0[0]*1.0 - 10, 0.0], [p0[0]*1.0+10, 30.0])


            if function_to_fit == 'lognormal_powerlaw':
                 p0 = [u_guess - 3, 2.0, np.sqrt(-0.5 * (u_guess - 3 - np.log(data['mean'])))  ]
                 bounds = ( [p0[0] - 5, 1.0, 0.01], [u_guess+3, 20.0, 10.0])

        if function_to_fit == 'powerlaw' or function_to_fit == 'truncated_powerlaw':
            p0     = [2, 1.0E-5]
            bounds = ( [1.0,0], [np.inf,1.0])

        if function_to_fit == 'truncated_powerlaw':
            p0      = [p0[0]*1,p0[1]*1, centers[np.max([np.argmax(norm_y)-1,0])] ]
            bounds  = ([bounds[0][0]*1, 1*bounds[0][1], centers[np.max([np.argmax(norm_y)-10,0])]],
                                                        [1*bounds[1][0],1*bounds[1][1],centers[np.min([np.argmax(norm_y)+20,np.size(centers)-1])]])
        if function_to_fit == 'gaussian_powerlaw':
            p0     = [data['mean'], 2.0, data['std']]
            bounds = [ (p0[0]/100.0, 0.01, p0[2] / 1000.0), (p0[0]*100.0, 10.0, p0[2]*100.0)]

    # choose the region to fit over
    #   for all, this will be where there is data

    selection = (norm_y > 0)
#    if function_to_fit == 'log-normal':
#        selection = selection * (centers > 10.0**(np.log10(centers[np.argmax(norm_y)]) - 1.0)  )*\
#                                (centers < 10.0**(np.log10(centers[np.argmax(norm_y)]) + 1.0)  )
#    elif function_to_fit == 'powerlaw':
#        selection = selection * ( centers > centers[np.argmax(norm_y)] )  # Fit everything to the right of the peak


    x_to_fit = centers[selection]
    y_to_fit = norm_y[selection]

    # set fitting function
    if isinstance(function_to_fit, str):
        function_to_fit= functions.by_name[function_to_fit]()


    if function_to_fit.name == 'lognormal_powerlaw' or function_to_fit.name == 'gaussian_powerlaw':
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

    fit_dictionary = {'popt' : copy(popt), 'pcov': copy(pcov),
                      'fit_function' : function_to_fit, 'norm_y' : norm_y*1.0,
                      'fit_x' : x_to_fit*1.0, 'fit_y' : y_to_fit*1.0,
                      'fit_result' : lambda xx : function_to_fit._f(xx, *popt) }

    return fit_dictionary

colors = ['C' + str(i) for i in [0,1,2,4,5,6,7,8,9]]
lss     = ['-','--']

def plot_all_elements(file_name, dsname, phase, elements = None, **kwargs):

    if elements is None:
        elements = ['Ba','Y','As','Sr','Mn','Na','Ca','N','Ni','Mg','S','Si','Fe','C','O']

    elements = ['Ba', 'Sr', 'Fe', 'Mn', 'Mg', 'O', 'N']

    fig, ax = plt.subplots(2)
    fig.set_size_inches(24,16)

    bins    = dd.io.load(file_name, '/' + '/'.join([dsname,phase,'mass_fraction','bins']))
    centers = 0.5 * (bins[1:]+ bins[:-1])


#    sort   = np.argsort(ymax_order)
#    sort_e = np.array(elements)[sort]

    label_pos = np.empty((np.size(elements),))
    label_pos[::2] = 1
    label_pos[1::2] = -1

#    label_pos = label_pos[sort]

    ci = li = 0
    for i, e in enumerate(elements):
        field = e + '_Fraction'

        ds_data = load_distribution_data(file_name, dsname, phase, field, centers = centers) # subselect

        fit_dict = fit_multifunction_PDF(1.0*bins, 1.0*ds_data['hist'], ds_data)
                   #fit_PDF(bins, ds_data['hist'], data = ds_data, **kwargs)


        plot_histogram(ax[0], np.log10(bins), fit_dict['norm_y']/np.max(fit_dict['norm_y']),
                                    color = colors[ci], ls = lss[li], lw = line_width)

        ax[0].plot(np.log10(centers),
                   #np.log10(fit_dict['fit_x']), 
                   fit_dict['fit_result'](centers) / np.max(fit_dict['norm_y']),
                                      lw = line_width, ls = lss[li], color = colors[ci])


        select = fit_dict['norm_y'] > 0
        chisqr = (fit_dict['fit_result'](centers[select]) - fit_dict['norm_y'][select])**2 / fit_dict['norm_y'][select]
        error  = np.abs(fit_dict['fit_result'](centers[select]) - fit_dict['norm_y'][select]) / fit_dict['norm_y'][select]

        ax[1].plot(np.log10(centers[select]), error, lw = line_width, ls = lss[li], color = colors[ci])


####
        xtext = np.log10(centers[np.argmax(fit_dict['norm_y'])]) - 0.1 - 0.05
        xa    = np.log10(centers[np.argmax(fit_dict['norm_y'])])
        ya    = 1.0
        pos = label_pos[i]
        if pos > 0:
            xtext = xtext + 0.05
            ytext = 2.0
        xy = (xtext, ytext)
        xya = (xa,ya)
        ax[0].annotate(e, xy = xya, xytext=xy, color = colors[ci],
                       arrowprops=dict(arrowstyle="-", connectionstyle="arc3"))
####

        print e, fit_dict['name'], fit_dict['popt'], np.sqrt(-0.5 * fit_dict['popt'][0]), np.sum(chisqr) / (1.0*np.size(fit_dict['norm_y'][select]))
        ci = ci + 1
        if ci >= np.size(colors):
            ci = 0
            li = li + 1

    for i in [0,1]:
        ax[i].set_xlim(-14, -1.5)
        ax[i].set_ylim(1.0E-5, 9.0)
        ax[i].semilogy()
    ax[1].set_ylim(0.01, 10.0)

    ax[0].set_ylabel('Peak Normalized PDF from Data')
    ax[1].set_ylabel('Peak Normalized PDF from Fit')

    fig.subplots_adjust(hspace = 0)
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible = False)
    plt.minorticks_on()

    fig.savefig(dsname + '_' + phase + '_all_elements_fit.png')
    plt.close()

    return

def plot_phase_panel(file_name, dsname, elements = None, plot_fit = True, **kwargs):

    phases = ['CNM','WNM','WIM','HIM','Disk']

    if elements is None:
        elements = ['Ba','Y','As','Sr','Mn','Na','Ca','N','Ni','Mg','S','Si','Fe','C','O']

    elements = ['Ba', 'Sr', 'Fe', 'Mn', 'Mg', 'O', 'N']
    #elements = ['Mg','O']
    fig, ax = plt.subplots(3,2)
    fig.set_size_inches(36,6*3)

    bins    = dd.io.load(file_name,'/'+'/'.join([dsname,phases[0],'mass_fraction','bins']))
    centers = 0.5 * (bins[1:]+ bins[:-1])


#    sort   = np.argsort(ymax_order)
#    sort_e = np.array(elements)[sort]

    label_pos = np.empty((np.size(elements),))
    label_pos[::2] = 1
    label_pos[1::2] = -1

#    label_pos = label_pos[sort]

    ax_indexes = [(0,0),(0,1), (1,0), (1,1), (2,0), (2,1)]
    for axi, phase in enumerate(phases):
        ci = li = 0

        axi = ax_indexes[axi]
        for i, e in enumerate(elements):
            field = e + '_Fraction'

            ds_data = load_distribution_data(file_name, dsname, phase, field, centers = centers) # subselect

            if plot_fit:
                fit_dict = fit_multifunction_PDF(1.0*bins, 1.0*ds_data['hist'], ds_data)
                   #fit_PDF(bins, ds_data['hist'], data = ds_data, **kwargs)


                y = fit_dict['norm_y'] / np.max(fit_dict['norm_y'])
                plot_histogram(ax[axi], np.log10(bins), y,
                                            color = colors[ci], ls = lss[li], lw = line_width)

                ax[axi].plot(np.log10(centers),
                           #np.log10(fit_dict['fit_x']), 
                           fit_dict['fit_result'](centers) / np.max(fit_dict['norm_y']),
                                              lw = line_width, ls = '--', color = colors[ci])
            else:
                y = ds_data['hist'] / (bins[1:] - bins[:-1])
                plot_histogram(ax[axi], np.log10(bins), y / np.max(y), color = colors[ci],
                                       lw = line_width, ls = '--')
####
            xtext = np.log10(centers[np.argmax(y) ]) - 0.1 - 0.05
            xa    = np.log10(centers[np.argmax(y)])
            ya    = 1.0
            pos = label_pos[i]
            if pos > 0:
                xtext = xtext + 0.05
                ytext = 2.0
            xy = (xtext, ytext)
            xya = (xa,ya)
            ax[axi].annotate(e, xy = xya, xytext=xy, color = colors[ci],
                           arrowprops=dict(arrowstyle="-", connectionstyle="arc3"))
####

            if plot_fit:
                if fit_dict['name'] == 'lognormal_powerlaw':
                    N = fit_dict['fit_function'].N
                else:
                    N = 0
                print phase, e, fit_dict['name'], fit_dict['popt'], "%5.5E"%(N)

            ci = ci + 1
            if ci >= np.size(colors):
                ci = 0
                li = li + 1

    for i, axi in enumerate(ax_indexes):
        ax[axi].set_xlim(-25, -1.5)
        ax[axi].set_ylim(1.0E-5, 9.0)
        ax[axi].semilogy()
        xy = (-3.0,1.0)
        ax[axi].annotate(phases[i], xy = xy, xytext=xy)

    for i in [0,1,2]:
        ax[(i,0)].set_ylabel('Peak Normalized PDF')       
        plt.setp( ax[(i,1)].get_yticklabels(), visible = False)
    #ax[1].set_ylim(0.01, 10.0)
    for i in [0,1]:
        ax[(2,i)].set_xlabel('log(Z)')         

    #ax[1].set_ylabel('Peak Normalized PDF')

    fig.subplots_adjust(hspace = 0, wspace = 0)
   # plt.setp([a.get_yticklabels() for a in fig.axes[1:]], visible = False)
    plt.minorticks_on()

    fig.savefig(dsname + '_multiphase_all_elements_fit.png')
    plt.close()

    return

def plot_phase_abundance_panel(file_name, dsname, elements = None, denominator = 'H', **kwargs):
    
    
    phases = ['CNM','WNM','WIM','HIM','Disk']

    if elements is None:
        elements = ['Ba','Y','As','Sr','Mn','Na','Ca','N','Ni','Mg','S','Si','Fe','C','O']

    elements = ['Ba', 'Sr', 'Fe', 'Mn', 'Mg', 'O', 'N']
    #elements = ['Mg','O']
    fig, ax = plt.subplots(3,2)
    fig.set_size_inches(36,6*3)

    bins    = dd.io.load(file_name,'/'+'/'.join([dsname,phases[0],'mass_fraction','abins']))
    centers = 0.5 * (bins[1:]+ bins[:-1])


#    sort   = np.argsort(ymax_order)
#    sort_e = np.array(elements)[sort]

    label_pos = np.empty((np.size(elements),))
    label_pos[::2] = 1
    label_pos[1::2] = -1

#    label_pos = label_pos[sort]

    ax_indexes = [(0,0),(0,1), (1,0), (1,1), (2,0), (2,1)]
    for axi, phase in enumerate(phases):
        ci = li = 0

        axi = ax_indexes[axi]
        for i, e in enumerate(elements):
            if e == denominator:
                continue
            field = e + '_over_' + denominator

            ds_data = load_distribution_data(file_name, dsname, phase, field, centers = centers) # subselect

            #fit_dict = fit_multifunction_PDF(1.0*bins, 1.0*ds_data['hist'], ds_data)
                   #fit_PDF(bins, ds_data['hist'], data = ds_data, **kwargs)

            yplot = np.cumsum(ds_data['hist'])  / (1.0*np.sum(ds_data['hist']))
            binsize = (10**(bins[1:]) - 10**(bins[:-1]))
            yplot = ds_data['hist'] / binsize

#            select = yplot > 0
#            xplot = bins[select]
#            yplot = yplot[select]

            plot_histogram(ax[axi], bins, yplot / np.max(yplot), #fit_dict['norm_y']/np.max(fit_dict['norm_y']),
                                        color = colors[ci], ls = lss[li], lw = line_width, label = e)

            #ax[axi].plot(np.log10(centers),
                       #np.log10(fit_dict['fit_x']), 
            #           fit_dict['fit_result'](centers) / np.max(fit_dict['norm_y']),
            #                              lw = line_width, ls = lss[li], color = colors[ci])

####
            xtext = np.log10(centers[np.argmax(ds_data['hist'])]) - 0.1 - 0.05
            xa    = np.log10(centers[np.argmax(ds_data['hist'])])
            ya    = 1.0
            pos = label_pos[i]
            if pos > 0:
                xtext = xtext + 0.05
                ytext = 2.0
            xy = (xtext, ytext)
            xya = (xa,ya)
#            ax[axi].annotate(e, xy = xya, xytext=xy, color = colors[ci],
#                           arrowprops=dict(arrowstyle="-", connectionstyle="arc3"))
####

            print phase, e, # fit_dict['name'], fit_dict['popt']
            ci = ci + 1
            if ci >= np.size(colors):
                ci = 0
                li = li + 1

    for i, axi in enumerate(ax_indexes):

        if denominator == 'H':
            ax[axi].set_xlim(-8,0)
        else:
            ax[axi].set_xlim(-3,3)
        ax[axi].set_ylim(1.0E-5, 5.0)
        ax[axi].semilogy()
        xy = (-3.0,1.0)
        ax[axi].annotate(phases[i], xy = xy, xytext=xy)

    ax[(0,0)].legend(loc = 'upper right')
    for i in [0,1,2]:
        ax[(i,0)].set_ylabel('Peak Normalized PDF')       
        plt.setp( ax[(i,1)].get_yticklabels(), visible = False)
    #ax[1].set_ylim(0.01, 10.0)
    for i in [0,1]:
        ax[(2,i)].set_xlabel('[X/' + denominator + ']')         

    #ax[1].set_ylabel('Peak Normalized PDF')

    fig.subplots_adjust(hspace = 0, wspace = 0)
   # plt.setp([a.get_yticklabels() for a in fig.axes[1:]], visible = False)
    plt.minorticks_on()

    fig.savefig(dsname + '_multiphase_elements_over_' + denominator + '.png')
    plt.close()

    return



if __name__ == '__main__':
    """
    Runs the analysis and plots a full panel of all phases and
    the full disk for selected elements. Additional arguments
    to run are:
       1) ds_number
       2) ds_number_start ds_number_end  (assume go up by 1)
       3) ds_number_start ds_number_end di

    e.x. To run starting with DD0400 to DD0500 going up by 5:
        $python ./PDF.py 400 500 5
    """


    if len(sys.argv) == 1:
        all_ds = ["DD0400"]
    elif len(sys.argv) == 2:
        all_ds = ["DD%0004i"%(int(sys.argv[1]))]
    elif len(sys.argv) == 3 or len(sys.argv) == 4:

        if len(sys.argv) == 4:
            di = int(sys.argv[3])
        else:
            di = 1

        all_ds = ["DD%0004i"%(x) for x in np.arange( int(sys.argv[1]),
                                                     int(sys.argv[2])+di/2.0, di)]

    individual_fail = False
    gas_file = 'gas_abundances.h5'
#    try:
#        data = {all_ds[0] : dd.io.load(gas_file, "/" + all_ds[0]) }
#    except:
#        individual_fail = True
#        all_data = dd.io.load(gas_file)


#    dsarray = ["DD%0004i"%(x) for x in np.arange(119, 202, 1)] #np.arange(50, 555, 250)]
#    plot_time_evolution(gas_file, dsarray, denominator = None)
#    plot_time_evolution(gas_file, dsarray, denominator = 'O')
#    plot_time_evolution(gas_file, dsarray, denominator = 'Fe')
#    plot_time_evolution(gas_file, dsarray, denominator = 'N')




    for dsname in all_ds:
        print "Beginning on " + dsname

#        if individual_fail:
#            data = {dsname : all_data[dsname]}
#        else:
#            data = {dsname : dd.io.load(gas_file, "/" + dsname) }

        plot_phase_panel(gas_file, dsname, plot_fit = False)


#        if True:
#             plot_phase_abundance_panel(data, dsname, denominator = 'Fe')
#            plot_phase_abundance_panel(data, dsname, denominator = 'Mg')
#            plot_phase_abundance_panel(data, dsname, denominator = 'O')
#
#
        if False:
            plot_all_elements(data, dsname, 'Molecular', function_to_fit = 'lognormal_powerlaw') #'lognormal_powerlaw')
            plot_all_elements(data, dsname, 'CNM', function_to_fit = 'lognormal_powerlaw') #'lognormal_powerlaw')
            plot_all_elements(data, dsname, 'WNM', function_to_fit = 'lognormal_powerlaw') #'lognormal_powerlaw')
            plot_all_elements(data, dsname, 'WIM', function_to_fit = 'lognormal_powerlaw') #'lognormal_powerlaw')
            plot_all_elements(data, dsname, 'HIM', function_to_fit = 'lognormal_powerlaw') #'lognormal_powerlaw')
            plot_all_elements(data, dsname, 'Disk', function_to_fit = 'lognormal_powerlaw') #'lognormal_powerlaw')


