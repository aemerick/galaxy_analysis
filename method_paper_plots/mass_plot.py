from galaxy_analysis.plot.plot_styles import *
from galaxy_analysis.analysis.compute_time_average import compute_time_average
from galaxy_analysis.utilities import utilities
import deepdish as dd
import numpy as np
import matplotlib.pyplot as plt
import glob

work_dir      = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11_30km/final_sndriving/'
#work_dir      = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11/corrected_sndriving/'
data_list, times = utilities.select_data_by_time(dir = work_dir,
                                                 tmin=0.0,tmax=1000)

all_data = {}
M_HI    = np.ones(np.size(data_list))
M_star  = np.ones(np.size(data_list))
M_total = np.ones(np.size(data_list))
M_H2    = np.ones(np.size(data_list))
for i,k in enumerate(data_list):
    M_HI[i] = dd.io.load(k, '/meta_data/M_HI')
    M_star[i] = dd.io.load(k, '/meta_data/M_star')
    M_total[i] = dd.io.load(k, '/meta_data/M_H_total') + dd.io.load(k,'/meta_data/M_He_total')
    M_H2[i] = dd.io.load(k, '/meta_data/M_H2I')



def plot_mass_evolution( t_f = None, image_num = 0):

    selection = (times == times) # all vals
    plot_times = times
    if not (t_f is None):
        selection  = times <= t_f
        plot_times = times[selection]

    plot_times = plot_times - plot_times[0]

    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    plot_histogram(ax, plot_times, M_total[selection], lw = line_width, color = 'black', label = r'M$_{\rm total}$')
    plot_histogram(ax, plot_times, M_HI[selection], lw = line_width, color = 'C0', label = r'M$_{\rm HI}$')
    plot_histogram(ax, plot_times, M_H2[selection], lw = line_width, color = 'C1', label = r'M$_{\rm H_2 }$')
    plot_histogram(ax, plot_times, M_star[selection], lw = line_width, color = 'C3', label = r'M$_*$')

    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(r'Mass in Disk (M$_{\odot}$)')
    ax.semilogy()

    ax.set_xlim(np.min(times-times[0]), np.max(times - times[0]))
    ax.legend(loc='lower right')
    plt.tight_layout()
    plt.minorticks_on()

    if t_f is None:
        fig.savefig('mass_evolution.png')
    else:
        fig.savefig('./mass_evolution_movie/mass_evolution_%0004i.png'%(image_num))

    plt.close()

    return


if __name__ == "__main__":
    plot_mass_evolution()

    if True:
        for i in np.arange(np.size(times)):
            print i
            plot_mass_evolution(t_f = times[i], image_num = i)
