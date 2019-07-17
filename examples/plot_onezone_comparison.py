import galaxy_analysis as gal
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import glob as glob
import os
import deepdish as dd
import sys
from galaxy_analysis.utilities import utilities

def sort_labels(l):

    elements = [x for x in l if len(x) <=2]
    other    = [x for x in l if len(x) > 2]

    return other + utilities.sort_by_anum(elements)

def load_galaxy_data(data_path = './../'):
    """
    Grab and load the most recent galaxy analysis data file
    """

    dfiles = glob.glob(data_path + 'DD????_galaxy_data.h5')
    if len(dfiles) < 1 and data_path == './../':
        dfiles = glob.glob(data_path + '../DD????_galaxy_data.h5')

    if len(dfiles) < 1:
        print("No galaxy analysis files found")
        return

    dfiles = np.sort(dfiles)

    gal    = dd.io.load(dfiles[-1])

    return gal

def load_onezone_data(data_path = './'):
    """
    Grab and load all of the summary outputs that exist
    """

    dfiles = glob.glob(data_path + 'run????_summary_output.txt')
    dfiles = np.sort(dfiles)

    all_data = [None]*len(dfiles)

    for i in np.arange(len(dfiles)):
        all_data[i] = np.genfromtxt(dfiles[i], names = True)

    return all_data

def compute_onezone_stats(all_data, element):
    """
    Go through all summary data and compute min, max, average, and std
    of all the chemical evolution data for a given field
    """

    if element == 'metals' or element == 'Metals' or element == 'Z':
        field_name = 'm_metal_mass'
    else:
        field_name = element + '_mass'

    mass = np.zeros(len(all_data))
    for i in np.arange(len(all_data)):
        mass[i] = all_data[i][field_name][-1]

    return utilities.compute_stats(mass, return_dict = True)

def save_onezone_stats(stats_dict):

    f = open('onezone_element_stats.dat','w')

    f.write('#element avg min max std frac_min frac_max frac_std\n')
    for k in list(stats_dict.keys()):
        x = stats_dict[k]
        f.write("%-8s  %5.5E %5.5E %5.5E %5.5E %5.5E %5.5E %5.5E\n"%(k,x['mean'],x['min'],x['max'],x['std'], x['min']/x['mean'], x['max']/x['mean'],x['std']/x['mean']))

    f.close()

    return

def compute_all_data():
    """
    Construct two dictionaries containing all the relevant data
    """

    gal = load_galaxy_data()

    simulation_data = gal['gas_meta_data']['masses']['FullBox']

    onezone_data = {}
    onez_data_files = load_onezone_data()

    # construct dictionary from onez model statistics
    elements = [x for x in list(simulation_data.keys()) if (not any([x in ['HI','HII','HeI','HeII','HeIII','Total Tracked Metals','H2','H2I']]))]
    for k in elements:
        if k == 'Total':
            continue
        onezone_data[k] = compute_onezone_stats( onez_data_files, k)

    save_onezone_stats(onezone_data)

    return simulation_data, onezone_data

def plot_variance(sim_data, onez_data):

    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    labels = sort_labels(list(onez_data.keys()))
    x = np.arange(1,len(labels)+1)

    yz = np.array([ onez_data[l]['std']/onez_data[l]['mean'] for l in labels])
    ax.scatter(x, yz, c = 'C0', s = 20, label = r'Onezone')
    ax.set_xticks(x)
    ax.set_ylim(1.0E-2, 10.0)
    ax.semilogy()
    ax.set_xticklabels(labels)
    ax.set_ylabel(r'Standard Deviation')
    plt.tight_layout()
    
    fig.savefig('onezone_std.png')
    plt.close()
    return

def plot_masses(sim_data, onez_data):

    fig, ax = plt.subplots()
    fig.set_size_inches(12, 6)

    labels = sort_labels( list(onez_data.keys()) )
    x      = np.arange(1, len(labels) + 1)

    yz = np.array([ onez_data[l]['mean'] for l in labels])
    ax.scatter(x, yz,  c = 'orange', s = 20,  label = 'Onezone Average')

    ylow = np.array([onez_data[l]['min'] for l in labels])
    yhigh = np.array([onez_data[l]['max'] for l in labels])

    ax.set_yscale('log')

    ax.errorbar(x, yz, yerr = [ np.abs(yz-ylow), np.abs(yz-yhigh)], fmt='o', c = 'orange', lw = 3)

    ysim = [ sim_data[l] for l in labels]
    ax.scatter(x, ysim,  c = 'black', s = 20,  label = 'Simulation Data')

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel(r'Mass (M$_{\odot}$)')
    plt.tight_layout()
    fig.savefig('total_mass.png')

    ax.set_ylim(1.0E-3, 1.0E4)
    fig.savefig('total_mass_compact.png')
    plt.close()


    fig, ax = plt.subplots()
    fig.set_size_inches(12,6)

    error = np.abs(ysim-yz)/yz
    ax.scatter(x, error, s = 20, color = 'black')
    print(error)
    ax.set_yscale('log')
    ax.set_ylabel("Fractional Error (sim - onez)/onez")

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    plt.tight_layout()
    fig.savefig('fractional_error.png')
    plt.close()

    return


if __name__ == "__main__":

    plot_average = False

    sim, onez = compute_all_data()

    plot_variance(sim, onez)
    plot_masses(sim, onez)
