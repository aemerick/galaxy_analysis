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
from galaxy_analysis.plot.plot_styles import *
import numpy as np
import matplotlib.pyplot as plt
import glob

import subprocess


_col_names = ['name', 'grid_id', 'pid', 'time', 'mass', 'birth_mass', 'metallicity',
             'm_eject', 'cells', 'volume', 'fraction', 'injected_mass',
             'energy_injected', 'ISM_mass', 'max_rho', 'avg_rho',
             'total_metal_mass', 'z_avg']

_dtypes = [str, int, int, float, float, float, float,
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

    filtered_data = np.zeros( num_unique, dtype = _ndtype[1:])

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

def read_all_data(directory = '.', filter = True):
    """
    Go through all files in directory and combine
    """

    # grep for IndividualStarSNStats flag in output files, combine to one
    bash_commands = ["grep --no-filename -e '^IndividualStarSNStats' " + directory + "/*.o* > temp_combined_output_o.txt",
                     "grep --no-filename -e '^IndividualStarSNStats' " + directory + "/*.e* > temp_combined_output_e.txt",
                     "cat " + directory + "/temp_combined_output_e.txt " + directory + "/temp_combined_output_o.txt > " + directory + "/temp_combined_output.txt",
                     "sed '/P/d' " + directory + "/temp_combined_output.txt > " + directory + "/temp_combined_output_2.txt",
                     "sed '/Load/d' " + directory + "/temp_combined_output_2.txt > " + directory + "/combined_output.txt",
                     "rm " + directory + "/temp_combined_output*.txt"]

    # execute commands
    for bc in bash_commands:
        subprocess.call(bc, shell=True)

    data = np.genfromtxt(directory + '/combined_output.txt', dtype = _ndtype)


    if filter:
        filtered_data = filter_data(data)
        return filtered_data
    else:
        return data

def save_data(data, directory = '.'):
    np.savetxt(directory + '/filtered_data.txt', data)
    return


def R_PDS(n, z, E_51 = 1.0):

    z         = z / 0.017     # in solar

    n, z = np.asarray(n), np.asarray(z)
    scalar_input = False
    if n.ndim == 0:
        n, z = n[None], z[None]
        scalar_input = True

    R_PDS = np.zeros(np.size(n))

    R_PDS[ z  < 0.01] = 49.3 * E_51**(0.25) * n[z < 0.01]**(-0.5)
    R_PDS[ z >= 0.01] = 18.5 * E_51**(2.0/7.0) * n[z >= 0.01]**(-3.0/7.0) * z[z>=0.01]**(-1.0/7.0)

    if scalar_input:
        return np.squeeze(R_PDS)

    return R_PDS


def rho_to_n(rho, mu= 1.3):
    m_p = 1.6737236E-24

    return rho / (m_p * mu)

def plot_rpds(data, dx = None, ncell = 3, norm = False):
    fig, ax = plt.subplots(figsize=(8,8))

    avg_dens = data['avg_rho']
    max_dens = data['max_rho']

    avg_dens = rho_to_n(avg_dens)
    max_dens = rho_to_n(max_dens)

    ncell_vol = 4.0 * np.pi * ncell**3 # well... not exactly

    R_avg = R_PDS(avg_dens, data['metallicity'], E_51 = 1 / (ncell_vol))
    R_max = R_PDS(max_dens, data['metallicity'], E_51 = 1 / (ncell_vol))

    avg_hist, bins = np.histogram(np.log10(R_avg), bins = np.linspace(0.0, 3.0, 25))
    max_hist, bins = np.histogram(np.log10(R_max), bins = bins)
    cent = 0.5 * (bins[1:] + bins[:-1])

    ax.set_xlabel(r'log[Supernova PDS Radius (pc)]')
    if norm:
        ax.set_ylabel(r'Fraction of Total')
    else:
        ax.set_ylabel(r'Count')

    normalization = 1.0
    if norm:
        normalization = 1.0 / (1.0 * np.sum(avg_hist))
    #ax.step(cent, avg_hist * normalization , lw = 3, ls = '-', color = 'black', label = r'<$n$>', where = 'post')
    plot_histogram(ax, bins, avg_hist*normalization, lw = line_width, ls = '-', color = 'black', label = r'<$n$>')

    if norm:
        normalization = 1.0 / (1.0 * np.sum(max_hist))

    #ax.step(cent, max_hist * normalization, lw = 3, ls = '-', color = 'orange',   label = r'$n_{\rm max}$', where='post')
    plot_histogram(ax, bins, max_hist*normalization, lw = line_width, ls = '-', color = 'orange', label = r'$n_{\rm max}$')



    ax.set_ylim(0, np.max([np.max(avg_hist*normalization),np.max(max_hist*normalization)])*1.5)
    ax.set_xlim(np.floor(np.min(bins)),np.ceil(np.max(bins)))

    if dx is not None:
        logdx = np.log10(dx)
        ax.plot( [logdx,logdx], ax.get_ylim(), ls = '--', color = 'black', lw =3) 
        logdx = np.log10(4.5*dx)
        ax.plot( [logdx,logdx], ax.get_ylim(), ls = '-', color = 'black', lw =3)

    # find fraction that are unresolved
    if True:
        avg_unres = np.size(R_avg[R_avg < 4.5 *dx]) / (1.0 *np.size(R_avg))
        max_unres = np.size(R_max[R_max < 4.5 *dx]) / (1.0 *np.size(R_max))

        txy_1 = (0.75*logdx, ax.get_ylim()[1]*0.8)
        txy_2 = (txy_1[0]  , txy_1[1] - 0.05)

#        ax.annotate('%0.2f %%'%(100.0*avg_unres), xy = txy_1, xytext = txy_1,
#                        color = 'black')
#        ax.annotate('%0.2f %%'%(100.0*max_unres), xy = txy_2, xytext = txy_2,
#                        color = 'orange')

        print "avg    max   ", 100*avg_unres, 100 * max_unres
        ax.annotate(r"Resolved (4.5$\times$ dx)", xy=(logdx-0.15,txy_1[1]),xytext=(logdx-0.15,txy_1[1]),
                      color='black', rotation=90)
        ax.annotate(r"Cell Width", xy=(np.log10(dx)-0.15,txy_1[1]),
                                   xytext=(np.log10(dx)-0.15,txy_1[1]),
                                   color = 'black', rotation=90)



    ax.legend(loc='best')
    plt.tight_layout()
    ax.minorticks_on()

    if norm:
        plt.savefig('sn_radius_dist_fraction.png')
    else:
        plt.savefig('sn_radius_distribution.png')

    plt.close()

    return

def plot_density(data, estimate_n = True):
    fig, ax = plt.subplots(figsize=(8,8))

    avg_dens = data['avg_rho']
    max_dens = data['max_rho']

    if estimate_n:
        avg_dens = rho_to_n(avg_dens)
        max_dens = rho_to_n(max_dens)

        bins = np.linspace(-5,3,40)
    else:
        bins = 10


    avg_hist, bins = np.histogram(np.log10(avg_dens), bins = bins)
    max_hist, bins = np.histogram(np.log10(max_dens), bins = bins)
    cent = 0.5 * (bins[1:] + bins[:-1])

#    ax.step(cent, avg_hist, lw = 3, ls = '-', color = 'black', label = r'<$n$>', where = 'post')
    plot_histogram(ax, bins, avg_hist, lw = line_width, ls = '-', color = 'black', label = r'<$n$>')
#    ax.step(cent, max_hist, lw = 3 , ls = '-', color = 'orange', label = r'$n_{\rm max}$', where = 'post')
    plot_histogram(ax, bins, max_hist, lw = line_width, ls = '-', color = 'orange', label = r'$n_{\rm max}$')

    ax.set_xlabel(r'log[ n (cm$^{-3}$)]')
    ax.set_ylabel(r'Count')

    ax.set_ylim(0, np.max([np.max(avg_hist),np.max(max_hist)])*1.5)
    ax.set_xlim(np.floor(np.min(bins)), np.ceil(np.max(bins)))
    ax.legend(loc='best')
    plt.tight_layout()
    ax.minorticks_on()
    plt.savefig('density_distribution.png')

    plt.close()

    return


if __name__ == "__main__":

    directory = '.'

    data = read_all_data(directory)
    save_data(data, directory)
    plot_density(data)
    plot_rpds(data, dx = 1.8, ncell = 1)
    plot_rpds(data, dx = 1.8, norm = True, ncell = 1)
