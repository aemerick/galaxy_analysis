import numpy as np
import matplotlib.pyplot as plt
import glob
import deepdish as dd
import sys
import os

def plot_sequestering(directory = './'):
    output_dir = directory + '/sequestering/'

    if not os.path.exists(output_dir)
        os.makedirs(output_dir)

    sfields = ['Disk', 'CNM', 'WNM', 'HIM', 'FullBox',
               'stars', 'Molecular', 'OutsideBox']
    all_elements = ['H','He','Fe','metals','O','C','N','Eu','Mg','S','Y','Ca','Si','Mn']

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
                plot_data[s][i] = x[s][element]


        ax.set_xlabel(r'Time (Myr')
        ax.set_ylabel(element + r' Mass (M$_{\odot}$)')
        ax.semilogy()
        plt.minorticks_on()
        fig.set_size_inches(8,8)
        plt.tight_layout()

        fig.savefig(output_dir + '/' + element + '_sequestering.png')
        plt.close(fig)
        del(plot_data)

    return


if __name__ = "__main__":

    directory = './'
    if len(sys.argv) == 2:
        directory = sys.argv[1]

    plot_mass_loading(directory = directory)
