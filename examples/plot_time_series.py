import numpy as np
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
import glob
import deepdish as dd
import sys
import os

purple  = '#7d03a8'
magenta = '#cb4679'
blue    = '#0c0887'
orange  = '#fdc328'
black   = 'black'


colors = {'Disk' : black,
              'CNM'  : purple,
              'WNM'  : purple,
              'HIM'  : purple,
              'Molecular' : purple,
          'FullBox'  : magenta,
          'stars'    : orange,
          'GravBound' : blue,
          'OutsideBox' : 'green' }

ls = {'Disk' : '-',
         'CNM' : '--',
         'WNM' : '-.',
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
    
    all_elements = elements
    if all_elements is None:   
        all_elements = ['H','He','Fe','Metals','O','C','N','Eu','Mg','S','Y','Ca','Si','Mn']

    #
    # get all data 
    #
    all_output = np.sort(glob.glob(directory + '/DD*.h5'))

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
        ymin = 1.0E-10 * ymax

        for s in sfields:

            # ignore fields that are too tiny to show up on plot
            if np.max(plot_data[s] * norm) < ymin:
                continue

            ax.plot(t, plot_data[s] * norm, lw = 4, label=s, ls = ls[s], color = colors[s])


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
        ax.legend(loc='best', )

        outname = output_dir + '/' + element + '_sequestering.png'
        if not fraction is None:
            outname = output_dir + '/' + element + '_' + fraction + '_fraction_sequestering.png')
        
        fig.savefig(outname)
        plt.close(fig)
        del(plot_data)

    return


if __name__ == "__main__":

    directory = './'
    if len(sys.argv) == 2:
        directory = sys.argv[1]

    plot_sequestering(directory = directory)

    # Disk and FullBox are likely the only two fractions that make
    # sense given the way the data is constructed. Other fractions are
    # certainly possible but will be misleading without recomputing
    # analysis (i.e. GravBound isn't useful unless I'm doing
    # fraction of a given field that is GravBound, which would
    # require re-computing, not just dividing by GravBound)
    plot_sequestering(directory = directory, fraction = 'Disk')
    plot_sequestering(directory = directory, fraction = 'FullBox')
