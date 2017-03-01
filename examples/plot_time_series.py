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
          'GravBound' : blue}

ls = {'Disk' : '-',
         'CNM' : '--',
         'WNM' : '-.',
         'HIM' : ':',
         'Molecular' : '-',
      'FullBox' : '-',
      'stars'   : '-',
      'GravBound' : '-'}


def plot_sequestering(directory = './'):
    output_dir = directory + '/sequestering/'

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    sfields = ['Disk', 'CNM', 'WNM', 'HIM', 'FullBox',
               'stars', 'Molecular', 'OutsideBox', 'GravBound']
    all_elements = ['H','He','Fe','Metals','O','C','N','Eu','Mg','S','Y','Ca','Si','Mn']

    #
    # get all data
    #
    all_output = np.sort(glob.glob(directory + '/DD*.h5'))
    
    for element in all_elements:
        fig, ax = plt.subplots()

        plot_data = {}
        t         = np.zeros(len(all_output))

        for s in sfields:
            plot_data[s] = np.zeros(len(all_output))

        for i in np.arange(len(all_output)):
            t[i] = dd.io.load(all_output[i], '/meta_data/Time')
            x    = dd.io.load(all_output[i], '/gas_meta_data/masses')

            for s in sfields:
                if element == 'Metals' and s == 'stars':
                    plot_data[s][i] = x[s]['metals']
                else:
                    plot_data[s][i] = x[s][element]

        for s in sfields:
            ax.plot(t, plot_data[s], lw = 4, label=s, ls = ls[s], color = colors[s])

        ax.set_xlabel(r'Time (Myr')
        ax.set_ylabel(element + r' Mass (M$_{\odot}$)')
        ax.semilogy()
        plt.minorticks_on()
        fig.set_size_inches(8,8)
        plt.tight_layout()
        ax.legend(loc='best')

        fig.savefig(output_dir + '/' + element + '_sequestering.png')
        plt.close(fig)
        del(plot_data)

    return


if __name__ == "__main__":

    directory = './'
    if len(sys.argv) == 2:
        directory = sys.argv[1]

    plot_sequestering(directory = directory)
