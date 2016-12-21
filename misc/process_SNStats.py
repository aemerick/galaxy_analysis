"""
    process_SNstats

    Author: A. Emerick

    Notes: script and functions to process output 
           SN statistics from IndividualStarPrintSNStats 
           parameter in Enzo. These stats are meant to understand
           the typical densities in which SN energy is injected,
           and to post-process check if SN are resolved at the given
           resolution
"""

import numpy as np
import matplotlib.pyplot as plt
import glob

import subprocess

_col_names = ['grid_id', 'pid', 'time', 'mass', 'birth_mass', 'metallicity',
             'm_eject', 'cells', 'volume', 'fraction', 'injected_mass',
             'energy_injected', 'ISM_mass', 'max_rho', 'avg_rho',
             'total_metal_mass', 'z_avg']

_dtypes = [int, int, float, float, float, float,
          float, int, float, float, float,
          float, float, float, float,
          float, float]

_ndtype = []
for i in np.arange(len(_col_names)):
    _ndtype += [(_col_names[i], _dtypes[i])]


def filter_data(data):
    """
    Filter data, returning only the non-repeating values
    """

    unique_sn = np.unique(data['pid'])
    num_unique = np.size(unique_sn)

    filtered_data = np.zeros( num_unique, dtype = _ndtype)

    for i, pid in enumerate(unique_sn):

        selection = data['pid'] == pid

        filtered_data[i] = (-1, pid, np.average(data['time'][selection]), np.average(data['mass'][selection]),
                        np.average(data['birth_mass'][selection]), np.average(data['metallicity'][selection]),
                        np.average(data['m_eject'][selection]), np.sum(data['cells'][selection]),
                        np.sum(data['volume'][selection]), np.sum(data['fraction'][selection]),
                        np.sum(data['injected_mass'][selection]), np.sum(data['energy_injected'][selection]),
                        np.sum(data['ISM_mass'][selection]), np.max(data['max_rho'][selection]),
                        np.average(data['avg_rho'][selection]), np.sum(data['total_metal_mass'][selection]),
                        np.average(data['z_avg'][selection]) ) #, dtype=_ndtype)
                        


    return filtered_data

def read_data(filename, filter = True):
    """
    Read in a single file and filter to individual SN events if
    filter is true.
    """



    return

def read_all_data(directory, filter = True):
    """
    Go through all files in directory and combine
    """

#    outputs = glob.glob(directory + '/*.o*')
    # merge the output files into one with just the data we want
    bash_command = "grep -h 'IndividualStarSNStats' " + directory + "/*.o* > combined_output.txt"
    subprocess.call(bash_command, shell=True)

    data = np.genfromtxt(directory + '/combined_output.txt', dtype = _ndtype)

    if filter:
        filtered_data = filter_data(data)
        return filtered_data
    else:
        return data

def save_data(directory, data):
    np.savetxt(directory + '/filtered_data.txt', data)
    return


def plot_density(data, estimate_n = True):
    fig, ax = plt.subplots(figsize=(8,8))

    avg_dens = data['avg_rho']
    max_dens = data['max_rho']

    if estimate_n:
        m_p = 1.6737236E-24
        mu  = 1.3

        avg_dens = data['avg_rho'] / (m_p*mu)
        max_dens = data['max_rho'] / (m_p*mu)

        bins = np.linspace(-3,3,30)
    else:
        bins = 10


    avg_hist, bins = np.histogram(np.log10(avg_dens), bins = bins)
    max_hist, bins = np.histogram(np.log10(max_dens), bins = bins)
    cent = 0.5 * (bins[1:] + bins[:-1])

    ax.step(cent, avg_hist, lw = 3, ls = '-', color = 'black', label = 'average density', where = 'pre')
    ax.step(cent, max_hist, lw = 3 , ls = '-', color = 'red', label = 'max density', where = 'pre')

    ax.set_xlabel(r'log[ $\rho$ (g cm$^{-3}$)]')
    ax.set_ylabel(r'Count')
    
    ax.legend(loc='best')
    plt.tight_layout()
    ax.minorticks_on()
    plt.savefig('density_distribution.png')

    plt.close()
    
    return


if __name__ == "__main__":

    directory = '.'

    data = read_all_data(directory)
    save_data(directory, data)
    plot_density(data)
