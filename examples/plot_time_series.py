import numpy as np
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
import glob
import deepdish as dd
import sys
import os

# some general plot styles for consistency
from galaxy_analysis.plot import plot_styles as ps

colors = {'Disk' : ps.black,
              'CNM'  : ps.purple,
              'WNM'  : ps.purple,
              'HIM'  : ps.purple,
              'Molecular' : 'red',
          'FullBox'  : ps.magenta,
          'stars'    : ps.orange,
          'GravBound' : ps.blue,
          'OutsideBox' : ps.green }

ls = {'Disk' : '-',
         'CNM' : '-',
         'WNM' : '--',
         'HIM' : ':',
         'Molecular' : '-',
      'FullBox' : '-',
      'stars'   : '-',
      'GravBound' : '-',
      'OutsideBox' : '-'}


def plot_sequestering(directory = './', fields = None, elements = None,
                      fraction = None):
    """
    Given a directory, goes through all data outputs in that
    directory and plots the time evolution of the mass contained
    in specified gas phases for a given species. The default behavior
    is to plot all fields for each element on a different plot for each
    element.

    fields : list, default is to use all of: Disk, CNM, WNM, HIM, FullBox,
             stars, Molecular, OutsideBox, GravBound

    elements : list, default is to loop through all elements: H, He, C, N,
               O, Mg, Ca, Si, Mn, S, Fe, Y, Eu, and Metals

    fraction : optional, string. Plot the mass fraction. Normalize all
               lines by one of the fields listed above.
    """

    output_dir = directory + '/sequestering/'

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if not fraction is None:
        output_dir += fraction + '/'

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    sfields = fields
    if sfields is None:
        sfields = ['Disk', 'CNM', 'WNM', 'HIM', 'FullBox',
                   'stars', 'Molecular', 'OutsideBox', 'GravBound']
    #
    # get all data
    #
    all_output = np.sort(glob.glob(directory + '/DD*.h5'))

    all_elements = elements
    if all_elements is None:
        # this is a hacky way of getting the metal fields - should just save this to file
        metal_fields = dd.io.load(all_output[0], '/gas_meta_data/masses/FullBox')
        exclude      = ['H','HI','HII','H2','Metals','Total','He','HeI','HeII','HeIII']
        metal_fields = utilities.sort_by_anum([x for x in metal_fields if (not any([y in x for y in exclude]))])

        individual_metals
        all_elements = ['H','He','Metals'] + individual_metals

    for element in all_elements:
        fig, ax = plt.subplots()

        # construct dictionary and time array
        plot_data = {}
        t         = np.zeros(len(all_output))

        for s in sfields:
            plot_data[s] = np.zeros(len(all_output))

        # loop through all output, gathering time and mass values
        for i in np.arange(len(all_output)):
            t[i] = dd.io.load(all_output[i], '/meta_data/Time')
            x    = dd.io.load(all_output[i], '/gas_meta_data/masses')

            # now loop through all fields
            for s in sfields:
                if element == 'Metals' and s == 'stars':
                    plot_data[s][i] = x[s]['metals']
                else:
                    plot_data[s][i] = x[s][element]

        # normalize if necessary
        norm = np.ones(np.size(plot_data[sfields[0]]))
        if not fraction is None:
            norm = 1.0 / plot_data[fraction]

        ymax = np.max(plot_data['FullBox'] * norm)
        if fraction is None:
            ymax = 5.0 * ymax

        ymin = 1.0E-8 * ymax

        # try and set something reasonable for minimum if it exists
        ymin = np.max( [ymin, np.min(plot_data['stars'] * norm)*10.0] )

        for s in sfields:

            # ignore fields that are too tiny to show up on plot
            if np.max(plot_data[s] * norm) < ymin:
                continue

            ax.plot(t, plot_data[s] * norm, lw = ps.lw, label=s, ls = ls[s], color = colors[s])


        ax.set_xlabel(r'Time (Myr')

        if not fraction is None:
            ax.set_ylabel('Mass Fraction of ' + element + ' relative to ' + fraction)
        else:
            ax.set_ylabel(element + r' Mass (M$_{\odot}$)')

        ax.set_ylim(ymin, ymax)
        ax.semilogy()
        plt.minorticks_on()
        fig.set_size_inches(8,8)
        plt.tight_layout()
        ax.legend(loc='best', ncol = 2)

        outname = output_dir + '/' + element + '_sequestering.png'
        if not fraction is None:
            outname = output_dir + '/' + element + '_' + fraction + '_fraction_sequestering.png'

        fig.savefig(outname)
        plt.close(fig)
        del(plot_data)

    return


def plot_surface_density_profiles(directory = './', normalize = False):

    output_dir = directory + '/radial_profiles/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_dir = output_dir + 'surface_density/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    all_output = np.sort(glob.glob(directory + '/DD*.h5'))

    da = dd.io.load(all_output[-1])
    fields = [k[1] for k in da['gas_profiles']['surface_density']['disk'].keys() if k != 'xbins']

    ls = ['-','--',':']
    color = ['black',ps.purple,ps.blue,ps.orange]

    for ele in fields:
        lsi = 0
        ci  = 0

        norm = 1.0
        ymin = 1.0E99
        ymax = -1.0E99


        k1 = 'gas'
        try:
            temp = da['gas_profiles']['surface_density']['disk'][(k1,ele)]
        except:
            k1 = 'enzo'

        if normalize:
            da = dd.io.load(all_output[0])
            norm = da['gas_profiles']['surface_density']['disk'][(k1,ele)]

        fig, ax = plt.subplots()

        for i in np.floor(np.linspace(0, len(all_output) - 1, 12)):
            i = int(i)
            da = dd.io.load( all_output[i])

            t = da['meta_data']['Time']

            y = da['gas_profiles']['surface_density']['disk'][(k1,ele)]
            x = da['gas_profiles']['surface_density']['disk']['xbins']
            x = (x[1:] + x[:-1])*0.5

            ax.plot(x, y * norm, lw = 3, ls = ls[lsi], color = color[ci], label = "t = %0.1f Myr"%(t))

            ymax = np.max( [ymax, np.max(y * norm)])
            ymin = np.min( [ymin, np.min(y * norm)])

            lsi = lsi + 1
            if lsi >= len(ls):
                lsi = 0
                ci = ci + 1

        ax.set_xlabel(r'Radius (pc)')
        ax.set_ylabel(r'$\Sigma_{\rm ' + ele + '}$ (M$_{\odot}$ pc$^{-2}$)')
        ax.semilogy()
#        ax.semilogx()
        ax.set_xlim(0.0, 2000.0)

        ymin = np.max( [ymin, 1.0E-5*ymax])
        ax.set_ylim(ymin, ymax)

        plt.minorticks_on()
        fig.set_size_inches(8,8)
        plt.tight_layout()
        ax.legend(loc='best', ncol = 2)

        outname = output_dir + ele
        if normalize:
            outname = outname + '_norm'
        outname = outname + '_radial_profile.png'
        fig.savefig(outname)
        plt.close(fig)

    return




def plot_mass_profiles(directory = './', normalize = False):
    """
    Plot cumulative radial profiles at 6 different times (evenly spaced)
    throughout the simulation.
    """

    output_dir = directory + '/radial_profiles/'

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    all_output = np.sort(glob.glob(directory + '/DD*.h5'))

    da = dd.io.load(all_output[-1])
    fields = [k[1] for k in da['gas_profiles']['accumulation']['sphere'].keys() if k != 'xbins']

    ls = ['-','--',':']
    color = ['black',ps.purple,ps.blue,ps.orange]

    for ele in fields:
        lsi = 0
        ci  = 0

        norm = 1.0

        if normalize:
            da = dd.io.load(all_output[0])
            norm = da['gas_profiles']['accumulation']['sphere'][('gas',ele)]

        fig, ax = plt.subplots()

        for i in np.floor(np.linspace(0, len(all_output) - 1, 12)):
            i = int(i)
            da = dd.io.load( all_output[i])

            t = da['meta_data']['Time']

            y = da['gas_profiles']['accumulation']['sphere'][('gas',ele)]
            x = da['gas_profiles']['accumulation']['sphere']['xbins']
            x = (x[1:] + x[:-1])*0.5

            ax.plot(x, np.cumsum(y) * norm, lw = 3, ls = ls[lsi], color = color[ci], label = "t = %0.1f Myr"%(t))

            lsi = lsi + 1
            if lsi >= len(ls):
                lsi = 0
                ci = ci + 1

        ax.set_xlabel(r'Radius (kpc)')
        ax.set_ylabel(r'Mass (M$_{\odot}$)')
        ax.semilogy()
#        ax.semilogx()
        ax.set_xlim(0.0, 15.0)

        plt.minorticks_on()
        fig.set_size_inches(8,8)
        plt.tight_layout()
        ax.legend(loc='best', ncol = 2)

        outname = output_dir + ele
        if normalize:
            outname = outname + '_norm_'
        outname = outname + 'radial_profile.png'
        fig.savefig(outname)
        plt.close(fig)

    return



if __name__ == "__main__":

    directory = './'
    if len(sys.argv) == 2:
        directory = sys.argv[1]


    plot_sequestering(directory = directory)
    print "completed total sequestering"
    plot_surface_density_profiles(directory = directory)
    print "completed surface density profiles"
    plot_mass_profiles(directory = directory)
    print "completed mass profiles"

    # Disk and FullBox are likely the only two fractions that make
    # sense given the way the data is constructed. Other fractions are
    # certainly possible but will be misleading without recomputing
    # analysis (i.e. GravBound isn't useful unless I'm doing
    # fraction of a given field that is GravBound, which would
    # require re-computing, not just dividing by GravBound)
    plot_sequestering(directory = directory, fraction = 'Disk')
    plot_sequestering(directory = directory, fraction = 'FullBox')
