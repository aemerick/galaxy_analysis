from galaxy_analysis.plot.plot_styles import *
from galaxy_analysis.analysis.compute_time_average import compute_time_average
from galaxy_analysis.utilities import utilities
import sys
import numpy as np
import matplotlib.pyplot as plt

#filepath = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11_30km/final_sndriving'


def plot(workdir = './', t_min = 250.0, t_max = 350.0,
         dv = 10, outdir = './'):

    phase_colors = {'cold' : 'C0', 'warm' : 'C1', 'hot' : 'C3',
                    'WNM'  : 'C0', 'WIM' : 'C1', 'HIM'  : 'C3'}

#   override with global
    for k in phase_colors:
        if k in color_dict.keys():
            phase_colors[k] = color_dict[k]

    labels = {'cold' : 'Cold' , 'warm' : 'Warm', 'hot' : 'Hot',
              'WNM' : "WNM", "WIM" : "WIM", "HIM" : "HIM"}

    fig, ax  = plt.subplots()
    fig.set_size_inches(8,8)

    sum = None
    for phase in ['WNM','WIM','HIM']: #['cold','warm','hot']:

        x,avg,min,max,std = compute_time_average(['gas_profiles','velocity','halo',phase], tmin = t_min, tmax = t_max,
                                                 dir = workdir, x_field = 'vbins')
        print np.min(x), np.max(x)
        x, avg = utilities.simple_rebin(x, avg, dv) # re-bin in 10 km/s
        print np.min(x), np.max(x)
        plot_histogram(ax, x, avg, color = phase_colors[phase], lw = line_width,
                                  ls = '-', label = labels[phase])
        if sum is None:
            sum = 1.0 * avg
        else:
            sum += avg

    plot_histogram(ax, x, sum, color = 'black', lw = line_width, ls = '-', label = 'Total')

    ax.set_xlabel(r'Radial Velocity (km s$^{-1})$')
    ax.set_ylabel(r'Mass (M$_{\odot}$)')
    ax.semilogy()
    ymin = 0.001
    ax.set_xlim(np.min(x[:-1][sum>=ymin]) ,np.max(  x[1:][sum>=ymin] ))
    ax.set_ylim(ymin, 4.0E5)
    ax.plot([0.0,0.0], ax.get_ylim(), lw = 2.0, color = 'black', ls = '--')
    plt.minorticks_on()
    plt.tight_layout()
    ax.legend(loc='best')

    fig.savefig(outdir + 'velocity_distribution_time_average.png')
    plt.close()

    f = open(outdir + 'velocity_percentiles.dat','w')
    cum_sum = np.cumsum(sum)
    percent = cum_sum / (cum_sum[-1]) * 100

    f.write("#percentile bin val\n")
    for q in np.arange(0,100,5):
        bin = np.max( [len(percent[percent <= q]) - 1, 0])
        f.write("%3i percentile: %3.3E %3.3E\n"%(q, bin, x[bin]))

    f.write("#outflowing gas ONLY\n")
    xcent = 0.5 * (x[1:] + x[:-1])
    x     = x[x>0]
    sum = sum[ xcent > 0.0]
    cum_sum = np.cumsum(sum)
    percent = cum_sum / (cum_sum[-1]) * 100.0
    for q in np.arange(0,100,5):
        bin = np.max( [ len(percent[percent <= q]) - 1, 0])
        f.write("%3i percentile: %3.3E %3.3E\n"%(q, bin, x[bin]))
    f.close()


    return

if __name__ == "__main__":
    work_dir = './'
    if len(sys.argv) > 1:
        work_dir = sys.argv[1]
    out_dir = './'
    if len(sys.argv) > 2:
        out_dir = sys.argv[2]

    tmin, tmax = 250.0, 350.0
    if len(sys.argv) > 3:
        tmin = float(sys.argv[3])
    if len(sys.argv) > 4:
        tmax = float(sys.argv[4])

    plot(workdir = work_dir, outdir = out_dir, t_min = tmin, t_max = tmax)
