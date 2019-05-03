import numpy as np
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
import glob
import deepdish as dd
import sys
import os
import yt

# some general plot styles for consistency
from galaxy_analysis.analysis import Galaxy
from galaxy_analysis.plot import plot_styles as ps
from galaxy_analysis.utilities import utilities
from galaxy_analysis.static_data import anum_to_asym


phases = ['Disk','CNM','WNM','WIM','HIM','FullBox']


def get_data(element, directory = './', t0 = None):

    all_output = np.sort(glob.glob(directory + '/DD*.h5'))

    data = {}
    t         = np.zeros(len(all_output))

    for s in phases:
        data[s] = np.zeros(len(all_output))

    for i in np.arange(len(all_output)):
        t[i] = dd.io.load(all_output[i], '/meta_data/Time')
        x    = dd.io.load(all_output[i], '/gas_meta_data/masses')

        for s in phases:
            data[s][i] = x[s][element]


    if t0 is None:
        t0 = t[0]

    t = t - t0

    return t, data

def compute_local_environment(filename = './mixing_events.in',
                              outfile = "event_data_table.dat"):

    mix_data   = load_mixing_file(filename = './mixing_events.in')

    #
    # first, need to know datafile to check for each event
    #
    all_files = np.sort(glob.glob('DD*.h5'))
    t = np.zeros(np.size(all_files))
    for i in np.arange(len(all_files)):
        t[i] = dd.io.load(all_files[i], '/meta_data/Time')

    index = np.zeros(np.size(mix_data['time']))
    for i in np.arange(np.size(mix_data['time'])):
        diff = mix_data['time'][i] - t

        if np.size( diff[diff>0]) == 0:
            index[i] = int(-9999) # hack - means need to look at "original" file
        else:
            diff[diff<0] = 1.0E10
            index[i] = int(np.argmin(diff))

    for field in ['avg_n','vol_avg_n','max_n','min_n','local_mass','local_volume','avg_T','vol_avg_T','N_massive_stars']:
        mix_data[field] = np.zeros(np.size(mix_data['time']))
#    mix_data['avg_n'] = np.zeros(np.size(mix_data['time']))
#    mix_data['vol_avg_n'] = np.zeros(np.size(mix_data['time']))
#    mix_data['max_n'] = np.zeros(np.size(mix_data['time']))
#    mix_data['min_n'] = np.zeros(np.size(mix_data['time']))
#    mix_data['local_mass']         = np.zeros(np.size(mix_data['time']))
#    mix_data['local_volume']       = np.zeros(np.size(mix_data['time']))
#    mix_data['avg_T']  = np.zeros(np.size(mix_data['time']))
#    mix_data['vol_avg_T'] = np.zeros(np.size(mix_data['time']))
    for i in np.arange(np.size(mix_data['time'])):
        # now loop through and compute
        if int(index[i]) == -9999:
            dsname = glob.glob('./original/DD????/DD????')[0]
        else:
            dsname = all_files[int(index[i])].strip('_galaxy_data.h5')
            dsname = dsname + '/' + dsname

        print dsname
        ds   = yt.load(dsname)
        gal  = Galaxy(dsname.split('/')[0])
        data = ds.all_data()

        center  = np.array([mix_data['x'][i],
                            mix_data['y'][i],
                            mix_data['z'][i]])

        center  = center / ds.length_unit.to('pc').value + ds.domain_center.value

        r       = np.min(data['dx']) * 4

        sp = ds.sphere(center, r)
        large_sp = gal.ds.sphere(center, 50.0*yt.units.pc)


        v  = sp['cell_volume']
        m  = sp['cell_mass']
        T  = sp['Temperature']
        n  = sp['number_density']

        mix_data['avg_n'][i] = np.sum(m * n) / np.sum(m)
        mix_data['vol_avg_n'][i] = np.sum(v*n)/np.sum(v)
        mix_data['max_n'][i] = np.max(n)
        mix_data['min_n'][i] = np.min(n)
        mix_data['avg_T'][i]    = np.sum(m * T) / np.sum(m)
        mix_data['vol_avg_T'][i] = np.sum(v * T) / np.sum(v)
        mix_data['local_mass'][i]         = np.sum(m.to('Msun'))
        mix_data['local_volume'][i]       = np.sum(v.to('pc**3'))

        t_o   = large_sp['creation_time'].to('Myr')
        death = t_o + large_sp[('io','particle_model_lifetime')].to('Myr')
        bm     = large_sp['birth_mass']
        tnow = ds.current_time.to('Myr')
        N_SN_stars = np.size( bm[ (bm > 8.0)*(bm<25.0)*((death - tnow) < 10.0*yt.units.Myr)*((death-tnow) > -5.0*yt.units.Myr)])

        mix_data['N_massive_stars'][i]    = N_SN_stars


    dd.io.save(outfile + '.h5', mix_data)

    outf = open(outfile,'w')

    outf.write("r_cyl z r_sph E51 M_ej M_Metal n_m n_v n_min n_max T_m T_v local_mass local_volume N_SN_stars time Anum Element\n")
    for i in np.arange(np.size(mix_data['time'])):
        outf.write("%.3f %.3f %.3f %.3E %.3E %.3E %.4E %.4E %.4E %.4E %.4E %.4E %.4E %.4E %5i %.2f %2i %2s\n"%(
                    mix_data['r_cyl'][i], mix_data['z'][i], mix_data['r_sph'][i],
                    mix_data['E_ej'][i], mix_data['M_ej'][i], mix_data['M_metal'][i],
                    mix_data['avg_n'][i], mix_data['vol_avg_n'][i], mix_data['min_n'][i], mix_data['max_n'][i],
                    mix_data['avg_T'][i], mix_data['vol_avg_T'][i], mix_data['local_mass'][i], mix_data['local_volume'][i],
                    mix_data['N_massive_stars'][i],
                    mix_data['time'][i], mix_data['Anum'][i], mix_data['element'][i]))
    outf.close()


    return

def get_all_data(mixing_filename = './mixing_events.in', directory = './'):

    mix_data = load_mixing_file(mixing_filename)

    all_data = {}
    for i, e in enumerate(mix_data['element']):

        t, data     = get_data(e, directory = directory, t0 = mix_data['time'][i])
        all_data[e] = data
        all_data[e]['time'] = t

    return all_data

def plot_average(cut_field = None, field_cuts = [],
                 annotation = None, show_legend = False):


    mix_data  = load_mixing_file('./mixing_events.in')
    #
    # plots the average of the evolution
    #

    all_data = get_all_data()

    if cut_field is None:

        #
        # Loop over elements, average fractions relative to disk
        #
        average = {}
        fullbox_average = {}
        individual_result = {}

        for e in all_data.keys():
            individual_result[e] = {}

        for field in phases:
            average[field] = np.zeros(np.size(all_data[ all_data.keys()[0] ][field]))
            fullbox_average[field] = np.zeros(np.size(all_data[ all_data.keys()[0] ][field]))

        average['CGM'] = np.zeros(np.size(average['Disk']))
        fullbox_average['CGM'] = np.zeros(np.size(average['Disk']))

        for f in phases + ['CGM']:

            if f == 'CGM':
                for e in all_data.keys():
                    average[f] += (all_data[e]['FullBox'] - all_data[e]['Disk']) / all_data[e]['FullBox']
                    fullbox_average[f] += (all_data[e]['FullBox'] - all_data[e]['Disk']) / all_data[e]['FullBox']
                    individual_result[e][f] = (all_data[e]['FullBox'] - all_data[e]['Disk']) / all_data[e]['FullBox']
            else:
                for e in all_data.keys():
                    average[f] += all_data[e][f] / all_data[e]['Disk']
                    fullbox_average[f] += (all_data[e][f] / all_data[e]['FullBox'])
                    individual_result[e][f] = all_data[e][f] / all_data[e]['FullBox']

            average[f] = average[f] / (1.0 * np.size(all_data.keys()))
            fullbox_average[f] = fullbox_average[f] / (1.0 * np.size(all_data.keys()))

        t = all_data[e]['time']

        fig, ax = plt.subplots()
        fig.set_size_inches(8,8)


        for s in phases:
            if s == 'FullBox' or s == 'Disk':
                continue
            ax.plot(t, average[s], color = ps.color_dict[s], ls = '-', lw = 3, label = s)
#        ax.plot(t, average['CGM'], color = 'black', ls = '--', lw = 3, label = 'CGM')
#        ax.plot(t, average['FullBox'] - plot_data['Disk'], color = 'black', ls = '--', label = 'CGM', lw = 3)

        ax.set_ylim(0.0,1.1)
        ax.set_xlim(0.0, np.max( [t[-1], 150.0]))
        ax.plot( ax.get_xlim(), [1.0,1.0], lw = 1, color = 'black', ls = '-')

        plt.minorticks_on()
        ax.legend(loc='upper right', ncol=2)

        ax.set_xlabel('Time (Myr)')
        ax.set_ylabel('Fraction of Enrichment Metal')


        val = np.average( average['CGM'][-4:-1] )

        if annotation is None:
            annotation = r"log(E$_{\rm ej}$ / [10$^{51}$ erg]) = %.1f"%(np.log10(np.average(mix_data['E_ej'])))

        xytext = (0.02,0.94)
        ax.annotate(annotation, xy=xytext, xytext=xytext,
                    xycoords = 'axes fraction',
                    textcoords = 'axes fraction')
        annotation2 = r"f$_{\rm outflow}$ = %0.3f"%(val)
        xytext = (0.02,0.86)
        ax.annotate(annotation2, xy=xytext, xytext=xytext,
                    xycoords = 'axes fraction',
                    textcoords = 'axes fraction')

        plt.tight_layout()

        outname = "enrichment_evolution_average"

        fig.savefig(outname + '.png')
        plt.close()
#
#
#       
#
#
        fig, ax = plt.subplots()
        fig.set_size_inches(8,8)


        for s in phases:
            if s == 'FullBox' or s == 'Disk':
                continue
            ax.plot(t, fullbox_average[s], color = ps.color_dict[s], ls = '-', lw = 3, label = s)

        ax.plot(t, fullbox_average['CGM'],   color = 'black', ls = '-', lw = 3, label = 'CGM')
        ax.plot(t, fullbox_average['Disk'] , color = 'black', ls = '--', label = 'Disk', lw = 3)

        ax.set_ylim(0.0,   1.0)
        ax.set_xlim(0.0, 150.0)
#        ax.plot( ax.get_xlim(), [1.0,1.0], lw = 1, color = 'black', ls = '-')

        plt.minorticks_on()

        if show_legend:
            ax.legend(loc='upper right', ncol=2)

        ax.set_xlabel('Time (Myr)')
        ax.set_ylabel('Fraction of Enrichment Metal')


        val = np.average( average['CGM'][-4:-1] )

#        annotation = r"log(E$_{\rm ej}$ / [10$^{51}$ erg]) = %.1f"%(np.log10(np.average(mix_data['E_ej'])))
        xytext = (0.02,0.94)
        ax.annotate(annotation, xy=xytext, xytext=xytext,
                    xycoords = 'axes fraction',
                    textcoords = 'axes fraction')
        #annotation = r"f$_{\rm outflow}$ = %0.3f"%(val)
        #xytext = (0.02,0.86)
        #ax.annotate(annotation, xy=xytext, xytext=xytext,
        #            xycoords = 'axes fraction',
        #            textcoords = 'axes fraction')

        plt.tight_layout()

        outname = "enrichment_evolution_average_CGM"

        fig.savefig(outname + '.png')
        plt.close()
#
#       --------------------------------------------------------
#
        header = "#time CNM WNM WIM HIM Disk CGM\n"
        f = open("average_enrichment_evolution.dat","w")
        f.write(header)

        for i in np.arange(np.size(t)):

            f.write("%5.5E"%(t[i]))

            for phase in ["CNM","WNM","WIM","HIM","Disk","CGM"]:
                f.write(" %5.5E"%( fullbox_average[phase][i]))
            f.write("\n")

        f.close()


        for e in all_data.keys():

            f = open(e + "_enrichment_evolution.dat","w")
            f.write(header)
            for i in np.arange(np.size(t)):

                f.write("%5.5E"%(t[i]))

                for phase in ["CNM","WNM","WIM","HIM","Disk","CGM"]:
                    f.write(" %5.5E"%( individual_result[e][phase][i]))
                f.write("\n")

            f.close()

    return



def enrichment_evolution(element,
                         directory = "./",
                         normalization = 1.0,
                         annotation = None, xytext = None,
                         t0 = None):


    t, plot_data = get_data(element, directory = directory, t0 = t0)

    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    if normalization == 'Disk' or normalization == 'disk':
        norm = plot_data['Disk']
    else:
        norm = normalization

    for s in phases:
        if s == 'FullBox':
            continue

        if ((normalization == "Disk" or normalization == 'disk') and  (s == 'Disk')):
            continue

        ax.plot(t, plot_data[s] / norm, color = ps.color_dict[s], ls = '-', lw = 3, label = s)

    if (not (normalization == "Disk" or normalization == 'disk')):
        ax.plot(t, plot_data['FullBox'] - plot_data['Disk'], color = 'black', ls = '--', label = 'CGM', lw = 3)

    ax.set_ylim(0.0,1.1)
    ax.set_xlim(t[0], np.max( [t[-1],40.0]))
    ax.plot( ax.get_xlim(), [1.0,1.0], lw = 1, color = 'black', ls = '-')

    plt.minorticks_on()
    ax.legend(loc='best', ncol=2)

    ax.set_xlabel('Time (Myr)')
    ax.set_ylabel('Fraction of Enrichment Metal')

    if (not (annotation is None)):
       if (xytext is None):
           xytext = (0.02,0.925)

       ax.annotate(annotation, xy=xytext, xytext=xytext,
                   xycoords = 'axes fraction',
                   textcoords = 'axes fraction')

    plt.tight_layout()

    outname = "enrichment_evolution_" + element
    if normalization == 'Disk' or normalization == 'disk':
        outname += '_disk'

    fig.savefig(outname + '.png')

    plt.close()

    return

def load_mixing_file(filename = './mixing_events.in'):

    # load mixing file data
    _data = np.genfromtxt(filename)

    LU = 1.82E23 # BAD TO HARD CODE - length unit of sim in cm

    data = {}
    single_event   = np.size(np.shape(_data)) == 1
    for i,k in enumerate(['time','x','y','z','M_ej','E_ej','Anum','M_metal']):

        if single_event:
            data[k] = _data[i]
        else:
            data[k] = _data[:,i]

    # convert to pc and normalize positions
    for k in ['x','y','z']:
        data[k] = (data[k] - 0.5) * LU * 3.24078E-19 # convert cm to pc

    data['r_cyl'] = np.sqrt(data['x']**2 + data['y']**2)
    data['r_sph'] = np.sqrt(  data['r_cyl']**2 + data['z']**2)

    if single_event:
        data['element'] = anum_to_asym[data['Anum']]
    else:
        data['element'] = [ anum_to_asym[x] for x in data['Anum'] ]

    if os.path.isfile("./event_data_table.dat"):
        d2 = np.genfromtxt("./event_data_table.dat", names = True)
        for name in d2.dtype.names:
            if not (name in data.keys()):
                data[name] = d2[name]

    return data

def plot_from_mixing_file(filename = './mixing_file',
                          norm     = 'disk',
                          annotate = True):
    """
    norm -> disk   or total (reads from file)
    """

    data = load_mixing_file(filename)

    annotation = None
    xytext     = None
    for i in np.arange(np.size(data['x'])):

        if annotate:
            if 'n_m' in data.keys():
                annotation = "<n> = %2.2E  -  <T> = %2.2E  -   E = %.1E  -  r_cyl = %.1f"%(data['n_m'][i],
                                                                   data['T_m'][i],
                                                                   data['E_ej'][i],
                                                                   data['r_cyl'][i])
            else:
                annotation = "r = %0.1f, z = %0.1f, E = %.1E"%(data['r_cyl'][i],
                                                               data['z'][i],
                                                               data['E_ej'][i])
            xytext     = (0.02,0.95)

        if norm == 'disk' or norm == 'Disk':
            x = 'disk'
        elif norm == 'total' or norm == 'Total':
            x = data['M_ej'][i]
        else:
            x = 1.0
            print "WARNING: Requested normalization not recognized", norm

        enrichment_evolution(data['element'][i], t0 = data['time'][i] + 0.1,
                             normalization = x,
                             annotation = annotation, xytext = xytext)

    return

if __name__ == "__main__":

    if len(sys.argv) > 1:

        if sys.argv[1] == 'all':
            mixing_file = './mixing_events.in'

            plot_from_mixing_file(mixing_file)

        elif sys.argv[1] == 'average':

            annotation   = None
            show_legend  = False

            if len(sys.argv) >= 3:
                annotation = sys.argv[2]
            elif len(sys.argv) >= 4:
                show_legend = bool( sys.argv[3] )

            plot_average(annotation = annotation, show_legend = show_legend)


        else:
            element = sys.argv[1]
            norm    = float(sys.argv[2])

            enrichment_evolution(element, norm = norm)

    else:
        enrichment_evolution('Na',norm = 0.24)
