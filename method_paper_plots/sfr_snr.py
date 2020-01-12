from galaxy_analysis.plot.plot_styles import *
from galaxy_analysis.analysis.compute_time_average import compute_time_average
from galaxy_analysis.utilities import utilities
import deepdish as dd
import numpy as np
import matplotlib.pyplot as plt
import glob
import sys


def plot(work_dir, separate_axis = False, TMAX = 500.0, outdir = './'):

    data_list = np.sort(glob.glob(work_dir + 'DD*_galaxy*.h5'))
    data = dd.io.load( data_list[-1] )

    all_data = {}
    times = data['time_data']['time'] / 1.0E6
    sfr   = (data['time_data']['SFR'])
    snr   = (data['time_data']['SNII_snr'])
    agbr  = data['time_data']['AGB_rate']

    times = times - times[0]

    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    plot_histogram(ax, times, sfr, lw = line_width, color = 'black', label = 'SFR')

    if separate_axis:
        ax2 = ax.twinx()
        norm = 1.0
    else:
        ax2 = ax
        norm = 100.0
    plot_histogram(ax2, times, snr*norm, lw = line_width, color = 'C3', label = 'SNRx100', ls = '--')

    if separate_axis:
        ax2.semilogy()
        ax2.set_ylabel(r'SNR (yr$^{-1}$)')

    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(r'SFR (M$_{\odot}$ yr$^{-1}$)')
    ax.semilogy()

    ax.set_xlim(np.min(times),np.min([TMAX,np.max(times)]))
    #ax.legend(loc='best')
    plt.tight_layout()
    plt.minorticks_on()

    if separate_axis:
        ax.set_ylim(2E-7,5.0E-4)
        ax2.set_ylim(ax.get_ylim())
    else:
        ax.set_ylim(1.0E-5, 1.0E-3)

    if separate_axis:
        outname = 'sfr_snr_2ax.png'
    else:
        ax.legend(loc = 'upper right')
        outname = 'sfr_snrx100.png'

    fig.savefig(outdir + outname)


    plot_histogram(ax2, times, agbr*10, lw = line_width, color = 'C2', label = 'AGB Rate x 10', ls = '-')
    ax.legend(loc = 'upper right')

    fig.savefig(outdir + 'sfr_snrx100_agb.png')

    print(np.min(agbr), np.max(agbr), np.average(agbr))
    return

if __name__ == "__main__":
    if len(sys.argv) > 1:
        work_dir = sys.argv[1]
    else:
        work_dir      = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11_30km/final_sndriving/'

    plot(work_dir, TMAX = 1.0E3)
