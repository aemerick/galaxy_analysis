import matplotlib
matplotlib.use('agg')

import matplotlib.pyplot as plt
import numpy as np
import gizmo_analysis as ga
import utilities as ga_ut

import sys

FIRE_elements      = ['h','he','c','n','o','ne','mg','si','s','ca','fe']
FIRE_metals = ['c','n','o','ne','mg','si','s','ca','fe']

#
# wrapper to load data set and set FIRE
# abundance tables
#
def load_with_FIRE(sim_index, wdir, return_model = False, model_Z = 1.0):
    """
    Convenience wrapper to load a dataset and set the appropriate initial
    abundances and FIRE yield model tables to do age tracer postprocessing
    with default FIRE2 yields
    """

    initial_part = ga.io.Read.read_snapshots(['gas'],'index',0,simulation_directory=wdir)
    initial_abundances = {}

    # -- note: np.average isn't necessary since they all should have the same value...
    for e in FIRE_elements:
        initial_abundances[e]    = np.average(initial_part['gas'].prop('massfraction.' + e))
    initial_abundances['metals'] = np.average(initial_part['gas'].prop('massfraction.metals'))


    part = ga.io.Read.read_snapshots(['gas','star'],
                              'index',
                              sim_index,
                              simulation_directory = wdir)


    FIRE_yield_model = ga.agetracers.FIRE2_yields(model_Z = model_Z # yield table metallicity in solar units,
                                            )
    FIRE_yield_table = ga.agetracers.construct_yield_table(FIRE_yield_model, # object with a `yields` function
                                                              part.ageprop.age_bins/1000.0) # in Gyr

    part.set_yield_table(FIRE_yield_table,          # the 2D table
                         FIRE_yield_model.elements  # list of elements we want to be able to post process
                                                # these can be anything included in the yield model
                        )

    # finally, set the initial abundances:
    #      As generated above, this is a dictionary corresponding to the initial
    #      mass fractions of each element (and all metals). If an element is missing,
    #      it is assumed to have an initial mass fraction of 0
    part.set_initial_abundances(initial_abundances)

    if return_model:
        return part, FIRE_yield_model, FIRE_yield_table
    else:
        return part


def compute_error(part, element, particle_type = 'star', filter=True):
    """
    For a single element, compute and bin the error
    """

    if filter:
        select = part[particle_type].prop('metallicity.o') > -3.8
    else:
        select = part[particle_type].prop('metallicity.o') == part[particle_type].prop('metallicity.o')

    error = np.abs(part[particle_type].prop('metallicity.agetracer.'+element)[select] - #+np.log10(0.68) -
                   part[particle_type].prop('metallicity.' + element)[select])

    return error

def bin_error(part, element, particle_type='star',
              amin = 0.00, amax = 3.0, da = 0.001, logbins=False):
    """
    Return a distribution of the error
    """

    error = compute_error(part,element,particle_type)
    print(element, np.min(error), np.max(error))

    if logbins:
        bins = np.arange(np.log10(amin), np.log10(amax), da)
        _error = np.log10(error)
    else:
        bins = np.arange(amin, amax+0.5*da, da)
        _error = error

    hist, bins = np.histogram(_error, bins=bins)

    stats = {'bins' : bins,
             'hist' : hist,
             'cbins': 0.5*(bins[1:]+bins[:-1]),
             'median' : np.percentile(error,50.0),
             'onesigma' : np.percentile(error,68.0),
             'twosigma' : np.percentile(error,95.0),
             'threesigma' : np.percentile(error,99.7)}

    return  stats


def compute_all_errors(part,particle_type='star',logbins=False,amin=0.0,amax=3,da=0.001):
    """

    """

    all_hist  = {}
    all_stats = {}

    for e in FIRE_metals:
        all_stats[e] = bin_error(part,e,particle_type,amin=amin,amax=amax,da=da,logbins=logbins)

    return all_stats


def generate_analysis(runs):

    all_part = {}
    all_data = {}
    for runname in runs.keys():
        all_part[runname] = load_with_FIRE(600,"./" + runname, model_Z = 0.1)
        all_data[runname] = {}
        all_data[runname] = compute_all_errors(all_part[runname],'star')

    fs = 5
    fig,ax = plt.subplots(3,3, sharex=True,sharey=True)
    fig.set_size_inches(3*fs,3*fs)
    fig.subplots_adjust(hspace=0,wspace=0)


    axi, axj = 0,0
    for e in FIRE_metals:
        axindex = (axi,axj)

        for runname in runs.keys():
        
            ploty = np.cumsum(all_data[runname][e]['hist']*1.0)
            ploty = ploty/ploty[-1]

            ax[axindex].plot(all_data[runname][e]['cbins'], ploty,
                    lw = 3, label = runs[runname])
        ax[axindex].set_ylim(0.0,1.0)
        ax[axindex].set_xlim(0.001, 1.0)
        ax[axindex].semilogx()
 

        axj = axj + 1
        if axj >= 3:
            axj = 0
            axi = axi + 1

    for i in np.arange(3):
        ax[(i,0)].set_ylabel("Cumulative Histogram")
        ax[(2,0)].set_xlabel("Error [dex]")

    fig.savefig("rate_error_panel.png")
    #
    # now compute the metallicities for each and the difference
    #

    def average_stats(runname, elements):
        one = two = three = 0.0
        for e in elements:
            one = one + all_data[runname][e]['onesigma']
            two = two + all_data[runname][e]['twosigma']
            three = three  + all_data[runname][e]['threesigma']

        n = 1.0 * len(elements)

        return one / n, two / n, three / n

    f = open('rates_results.dat','w')

    f.write("name alpha_one alpha_two alpha_three wind_one wind_two wind_three ia_one ia_two ia_three\n")

    for run in runs.keys():

        f.write(run)

        for elist in [ ['o','mg','si','ca'],  ['c','n'], ['fe']]:
            one, two, three = average_stats(run, elist)
            f.write(" %.4f %.4f %.4f"%(one,two,three))
        f.write("\n")
        
    f.close()


def compare_distributions(runs):

    all_part = {}
    all_data = {}

    amin = 0.01
    amax = 2.0
    da = 0.01

    for runname in runs.keys():
        all_part[runname] = load_with_FIRE(600,"./" + runname, model_Z = 0.1)
        all_data[runname] = {}
        all_data[runname] = compute_all_errors(all_part[runname],'star',
                                               amin = amin, amax = amax, da = da,
                                               logbins=True)


    fs = 6
    fig,ax = plt.subplots(1,3, sharex=True,sharey=True)
    fig.set_size_inches(3*fs,1*fs)
    fig.subplots_adjust(wspace=0)

    for axi, element_list in enumerate([ ['o','mg','si','ca'], ['c','n'], ['fe']]):


        for runname in runs.keys():
            # compute average
            avg_hist = np.zeros(np.size(all_data[runname]['c']['hist']))
            count = 0

            for e in element_list:
                avg_hist += all_data[runname][e]['hist']
                count = count + 1
            avg_hist = avg_hist / (1.0*count)
        
            dz = all_data[runname][e]['bins'][1:] - all_data[runname][e]['bins'][:-1]

            ax[axi].plot(all_data[runname][e]['cbins'],
                         avg_hist / (1.0*np.sum(avg_hist)) / dz,
                         lw = 3, label = runs[runname])


        ax[axi].set_ylim(0,7)
        ax[axi].set_xlim(np.log10(amin),np.log10(amax))
#        ax[axi].semilogx()
        ax[axi].set_xlabel(r"log$_{10}$(Error [dex])")

    ax[0].set_ylabel("dN/d(log(Error))")


    xy = (0.8,0.1)
    ax[0].annotate(r'$\alpha$', xy, xy, xycoords='axes fraction')
    ax[1].annotate("C+N", xy, xy, xycoords='axes fraction')
    ax[2].annotate("Fe", xy, xy, xycoords='axes fraction')

    ax[0].legend(loc='best')


    fig.savefig('distributions.png', bbox_inches='tight', pad_inches=0.0)    


    return


def compare_runtime(runs):

    return

def plot_sigma(runs):

    data = np.genfromtxt("rates_results.dat",names=True)

    xvals = np.arange(6)

    fig, ax = plt.subplots(1,3, sharey=True)
    fig.set_size_inches(18,7)
    fig.subplots_adjust(wspace=0)

    colors = {'one': 'C0', 'two' : 'C1', 'three' :'C2'}
    labels = {'one': r'1 $\sigma$', 'two': r'2 $\sigma$', 'three' : r'3 $\sigma$'}

    annotation = {'alpha': r"$\alpha$ (CCSNe)", 'wind' : "C + N (Winds)", 'ia' : "Fe (Ia's)"}


    for i, metal in enumerate(['alpha','wind','ia']):
        for sigma in ['one','two','three']:
            ax[i].plot(xvals, data[metal + '_' + sigma], color = colors[sigma], lw = '2', ls = ':')
            ax[i].scatter(xvals, data[metal + '_' + sigma], s = 30, color = colors[sigma], label = labels[sigma])

    
        ax[i].set_ylim(0.0,1.0)
        ax[i].set_xlim(-0.5, 5.5)
        ax[i].set_xticks(xvals)
        ax[i].set_xticklabels( list(runs.values()), rotation = 'vertical', fontsize = 12)


        xy = (0.8,0.05)
        ax[i].annotate(annotation[metal], xy,xy,xycoords='axes fraction')

    ax[0].set_ylabel("Simulation - AgeTracer [dex]")
    ax[0].legend(loc='best')



   
    fig.savefig("sigma_comparison.png", bbox_inches='tight', pad_inches=0.05)


    return

if __name__ == "__main__":

    runs = {"test0" : "R = 1.0", 'test1' : "R = 0.5",
            "test2" : "R = 0.1", "test3" : "R = -1",
            'test4' : "R = -10", "test5" : "R = -100"}


    generate = False
    compare_runtime = False
    compare_sigma = False
    plot_distributions = False
    if len(sys.argv) > 1:
        if 'generate' in sys.argv:
            generate = True

        if 'runtime' in sys.argv:
            compare_runtime = True

        if 'compare_sigma' in sys.argv:
            compare_sigma = True

        if 'plot_distributions' in sys.argv:
            plot_distributions = True


    else:
        generate = True


    if generate:
        generate_analysis(runs)

    if compare_runtime:
        compare_runtime(runs)


    if compare_sigma:
        plot_sigma(runs)

    if plot_distributions:
        compare_distributions(runs)
