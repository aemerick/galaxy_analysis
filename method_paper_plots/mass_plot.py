from galaxy_analysis.plot.plot_styles import *
from galaxy_analysis.analysis.compute_time_average import compute_time_average
from galaxy_analysis.utilities import utilities
import deepdish as dd
import numpy as np
import matplotlib.pyplot as plt
import glob
import sys

rc('font',size=22)
#data_list, times = utilities.select_data_by_time(dir = work_dir,
#                                                 tmin=0.0,tmax=650.0)
#M_HI    = np.ones(np.size(data_list))
#M_star  = np.ones(np.size(data_list))
#M_total = np.ones(np.size(data_list))
#M_H2    = np.ones(np.size(data_list))
#for i,k in enumerate(data_list):
#    M_HI[i] = dd.io.load(k, '/meta_data/M_HI')
#    M_star[i] = dd.io.load(k, '/meta_data/M_star')
#    M_total[i] = dd.io.load(k, '/meta_data/M_H_total') + dd.io.load(k,'/meta_data/M_He_total')
#    M_H2[i] = dd.io.load(k, '/meta_data/M_H2I')
#
#
def plot_mass_resolution(work_dir = './', output_dir = None, comparison = None, new_color = False, colors = None):

    if output_dir is None:
        output_dir = work_dir

    if comparison is None:
        labels = {'3pcH2' : '3.6 pc' , '6pcH2' : '7.2 pc', 'Fiducial' : 'Fiducial'}
        lstyle = {'3pcH2' : '--', '6pcH2' : '-.', 'Fiducial' : '-'}
        dirs   = {'3pcH2' : '../3pc_H2/' , '6pcH2' : '../6pc_H2/', 'Fiducial' : work_dir}
    
        comparison = {}
        for k in labels.keys():
            comparison[k] = (dirs[k],labels[k],lstyle[k])

    else:
        dirs = {}
        labels = {}
        lstyle = {}
        for k in comparison.keys():
            dirs[k]   = work_dir + comparison[k][0]
            labels[k] = comparison[k][1]
            lstyle[k] = comparison[k][2]

#    labels = {'3pc_hsn' : '3.6 pc - SNx2', '3pc' : '3.6 pc', 'final_sndriving' : 'Fiducial', '6pc_hsn' : '7.2 pc'}
#    lstyle     = {'3pc_hsn' : '--', '3pc' : ':', 'final_sndriving' : '-', '6pc_hsn' : '-.'}

    

    all_data = {}
    for k in comparison.keys():
        all_data[k] = {}

#        if k == 'final_sndriving':
#            all_data[k]['times'] = times
#            all_data[k]['M_HI'] = M_HI
#            all_data[k]['M_star'] = M_star
#            all_data[k]['M_total'] = M_total
#            all_data[k]['M_H2I']    = M_H2
#
        if True:
            dl, t = utilities.select_data_by_time(dir =  dirs[k], tmin = 0.0, tmax=1000.0)

            all_data[k]['times'] = t
            all_data[k]['M_HI'] = np.zeros(np.size(dl))
            all_data[k]['M_star'] = np.zeros(np.size(dl))
            all_data[k]['M_total'] = np.zeros(np.size(dl))
            all_data[k]['M_H2I'] = np.zeros(np.size(dl))


            for i,d in enumerate(dl):
                all_data[k]['M_HI'][i] = dd.io.load(d, '/meta_data/M_HI')
                all_data[k]['M_star'][i] = dd.io.load(d, '/meta_data/M_star')
                all_data[k]['M_total'][i] = dd.io.load(d, '/meta_data/M_H_total') + dd.io.load(d,'/meta_data/M_He_total')
                all_data[k]['M_H2I'][i] = dd.io.load(d, '/meta_data/M_H2I')

    # done collecting data rom all data sets

    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    if colors is None:
        colors = {}
        i = 0
        for k in comparison.keys():
            if not (k in colors.keys()):
                colors[k] = 'C%1i'%(i)
                i = i + 1
    
    
    lstyle['M_total'] = '-'
    lstyle['M_HI'] = '--'
    lstyle['M_H2I'] = '-.'
    lstyle['M_star'] = ':'
    
    for k in comparison.keys():


        for field,color in [('M_total','black'), ('M_HI','C0'), ('M_H2I','C1'), ('M_star','C3')]:
#        for field, ls in [('M_total','-'), ('M_HI', '--'), ('M_H2I', ':')]:  # , ('M_star' , '-.')]:

            if field == 'M_total':
                label = labels[k]
            else:
                label = None
            # print k, field, np.size(all_data[k]['times']), np.size(all_data[k][field])
            if new_color:
                ax.plot(all_data[k]['times'] - all_data[k]['times'][0], all_data[k][field], 
                        ls = lstyle[field],
                        lw = line_width, color = colors[k])
            else:
                        
                plot_histogram(ax, all_data[k]['times'] - all_data[k]['times'][0], all_data[k][field], ls = lstyle[k], lw = line_width, color = color)

    if new_color:
        ax.plot((-1,-1), (-2,-2), ls = '-', lw = 3, color = 'black', label = r'M$_{\rm total}$')
        ax.plot((-1,-1), (-2,-2), ls = '--', lw = 3, color = 'black', label = r'M$_{\rm HI}$')
        ax.plot((-1,-1), (-2,-2), ls = '-.', lw = 3, color = 'black', label = r'M$_{\rm H_2}$')
        ax.plot((-1,-1), (-2,-2), ls = ':', lw = 3, color = 'black', label = r'M$_{\rm *}$')        
    else:
        ax.plot((-1,-1), (-2,-2), ls = '-', lw = 3, color = 'black', label = r'M$_{\rm total}$')
        ax.plot((-1,-1), (-2,-2), ls = '-', lw = 3, color = 'C0', label = r'M$_{\rm HI}$')
        ax.plot((-1,-1), (-2,-2), ls = '-', lw = 3, color = 'C1', label = r'M$_{\rm H_2}$')
        ax.plot((-1,-1), (-2,-2), ls = '-', lw = 3, color = 'C3', label = r'M$_{\rm *}$')


    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(r'Mass in Disk (M$_{\odot}$)')
    ax.semilogy()

    ax.set_xlim(0.0, 500.0)
    ax.legend(loc='lower right')
#    ax.set_ylim(6E4, 3E6)
    plt.tight_layout()
    plt.minorticks_on()

    fig.savefig(work_dir + output_dir + 'mass_evolution_resolution.png')
    plt.close()
    return

def plot_mass_evolution(work_dir, t_f = None, image_num = 0, outdir = './',
                        TMAX = 500.0):

    data_list, times = utilities.select_data_by_time(dir = work_dir,
                                                     tmin=0.0,tmax=650.0)
    M_HI    = np.ones(np.size(data_list))
    M_star  = np.ones(np.size(data_list))
    M_total = np.ones(np.size(data_list))
    M_H2    = np.ones(np.size(data_list))
    for i,k in enumerate(data_list):
        M_HI[i] = dd.io.load(k, '/meta_data/M_HI')
        M_star[i] = dd.io.load(k, '/meta_data/M_star')
        M_total[i] = dd.io.load(k, '/meta_data/M_H_total') + dd.io.load(k,'/meta_data/M_He_total')
        M_H2[i] = dd.io.load(k, '/meta_data/M_H2I')

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
    ax.set_xlim(np.min(times-times[0]), np.min([TMAX, np.max(times - times[0])]) )
    ax.legend(loc='lower right')
    plt.tight_layout()
    plt.minorticks_on()

    if t_f is None:
        fig.savefig(outdir + 'mass_evolution.png')
    else:
        fig.savefig('./mass_evolution_movie/mass_evolution_%0004i.png'%(image_num))

    plt.close()

    return


if __name__ == "__main__":

    if len(sys.argv) > 1:

        wdir = ''
        work_dir = sys.argv[1]

    else:
        wdir = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11_30km/'
        work_dir = wdir + 'final_sndriving/'

    plot_mass_evolution(work_dir, TMAX = 1.0E3)

    plot_mass_resolution(work_dir)


    if False:
        for i in np.arange(np.size(times)):
            print(i)
            plot_mass_evolution(t_f = times[i], image_num = i)
