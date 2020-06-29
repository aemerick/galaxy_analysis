import matplotlib
matplotlib.use('agg')

import matplotlib.pyplot as plt
import numpy as np
import gizmo_analysis as ga
import utilities as ga_ut


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


def compute_error(part, element, particle_type = 'star'):
    """
    For a single element, compute and bin the error
    """

    error = np.abs(part[particle_type].prop('metallicity.agetracer.'+element) -
                   part[particle_type].prop('metallicity.' + element))

    return error

def bin_error(part, element, particle_type='star',
              amin = 0.001, amax = 1.0, da = 0.001):
    """
    Return a distribution of the error
    """

    error = compute_error(part,element,particle_type)
    bins = np.arange(amin, amax+0.5*da, da)
    hist, bins = np.histogram(error, bins=bins)

    stats = {'bins' : bins,
             'hist' : hist,
             'cbins': 0.5*(bins[1:]+bins[:-1]),
             'median' : np.percentile(error,50.0),
             'onesigma' : np.percentile(error,68.0),
             'twosigma' : np.percentile(error,95.0),
             'threesigma' : np.percentile(error,99.7)}

    return  stats


def compute_all_errors(part,particle_type='star'):
    """

    """

    all_hist  = {}
    all_stats = {}

    for e in FIRE_metals:
        all_stats[e] = bin_error(part,e,particle_type)

    return all_stats



if __name__ == "__main__":

    runs = {"test0" : "R = 1.0", 'test1' : "R = 0.5",
            "test2" : "R = 0.1"}

    all_part = {}
    all_data = {}
    for runname in runs.keys():
        all_part[runname] = load_with_FIRE(600,"./" + runname)
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
