from galaxy_analysis.plot.plot_styles import *

import numpy as np
import matplotlib.pyplot as plt
import deepdish as dd
import h5py, glob, sys

from galaxy_analysis.utilities import utilities
from galaxy_analysis.analysis import Galaxy

from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import binned_statistic_2d

# temporary
import time as cpu_time


def load_abundance_data(data, data_list,
                        fields,
                        property_list,
                        phases = ['CNM','WNM','WIM','HIM'],
                        field_types = None):

    time_data = {}
    for i, field in enumerate(fields):
        time_data[field] = {}

        if field_types is None:
            if '_Fraction' in field or '_over_' in field:
                ft = 'mass_fraction'
         #   elif '_over_' in field:
         #       ft = 'abundance'
            else:
                print("Cannot properly determing field type for " + field)
                raise RunTimeError

        elif hasattr(field_types, 'keys') :
            ft = field_types[field]
        else:
            ft = field_types[i]

        for phase in phases:
            time_data[field][phase] = {}

            for property in property_list:

                if property == 'IQR' or property == 'inner_quartile_range':
                    for p in ['Q1','Q3']:
                        field_path = [phase,ft,field,p]
                        # need to load the data here into a time array
                        time_data[field][phase][p] = utilities.get_property(field_path,
                                                                               file_list = data.filename,
                                                                               data_list = data_list,
                                                                               self_contained = True)
                        time_data[field][phase][p] = np.array(time_data[field][phase][p], dtype = np.float)

                    if '_Fraction' in field: # do not do for abundance ratios
                        time_data[field][phase][property] = np.log10(time_data[field][phase]['Q3']) -\
                                                            np.log10(time_data[field][phase]['Q1'])
                    else:
                        time_data[field][phase][property] = time_data[field][phase]['Q3'] -\
                                                            time_data[field][phase]['Q1']

                elif property == 'mean_median_distance':
                    for p in ['median','mean']:
                        if not p in list(time_data[field][phase].keys()):
                            field_path = [phase,ft,field,p]
                            time_data[field][phase][p] = utilities.get_property(field_path,
                                                                                file_list = data.filename,
                                                                                data_list = data_list,
                                                                                self_contained = True)
                            time_data[field][phase][p] = np.array(time_data[field][phase][p], dtype = np.float)
                    if '_Fraction' in field: # do not do for abundance ratios
                        time_data[field][phase]['mean_median_distance'] =\
                                                   np.log10(time_data[field][phase]['mean']) -\
                                                   time_data[field][phase]['median']
                    else:
                        time_data[field][phase]['mean_median_distance'] = time_data[field][phase]['mean']-\
                                                                          time_data[field][phase]['median']

                elif property == 'inner_decile_range' or property == 'd9_d1_range':
                    for p in ['decile_1','decile_9']:
                        field_path = [phase,ft,field,p]
                        # need to load the data here into a time array
                        time_data[field][phase][p] = utilities.get_property(field_path,
                                                                               file_list = data.filename,
                                                                               data_list = data_list,
                                                                               self_contained = True)
                        time_data[field][phase][p] = np.array(time_data[field][phase][p], dtype = np.float)

                    if '_Fraction' in field: # do not do for abundance ratios
                        time_data[field][phase][property] = np.log10(time_data[field][phase]['decile_9']) -\
                                                            np.log10(time_data[field][phase]['decile_1'])
                    else:
                        time_data[field][phase][property] = time_data[field][phase]['decile_9']-\
                                                            time_data[field][phase]['decile_1']
                else:
                    field_path = [phase,ft,field,property]
                    time_data[field][phase][property] = utilities.get_property(field_path,
                                                                               file_list = data.filename,
                                                                               data_list = data_list,
                                                                               self_contained = True)
                    time_data[field][phase][property] = np.array(time_data[field][phase][property], dtype = np.float)

                    if '_Fraction' in field: # do not do for abundance ratios
                        if property in ['median','mean']:
                            time_data[field][phase][property] = np.log10(time_data[field][phase][property])

    return time_data

def plot_stellar_2d_hist(galaxy, field, time_bins = np.arange(0.0,20.0,0.2),
                         ybins = np.arange(-20,-6,0.1)):

    if '_Fraction' in field:
        yval = np.log10(galaxy.df[('io','particle_' + field.strip('_Fraction') + '_fraction')].value)
    else:
        yval = galaxy.df[('io','particle_' + field)].value

    creation_time = galaxy.df['creation_time'].to('Myr').value

    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)
    statistic_data = np.ones(np.size(creation_time))
    N, x_edge, y_edge, binnum = binned_statistic_2d(creation_time - np.min(creation_time), yval, statistic_data,
                                                    statistic = 'count', bins =  (time_bins,ybins))

    fraction = N / (1.0 * np.size(creation_time))
    fraction[fraction <= 0] = -99
    fraction[fraction >  0] = np.log10(fraction[fraction > 0])
    fraction = np.log10(N / (1.0 * np.size(creation_time)))
    plot_val = fraction

    xmesh, ymesh = np.meshgrid(x_edge, y_edge)
    img1 = ax.pcolormesh(xmesh, ymesh, plot_val.T,
                              cmap = 'magma', vmin = -4, vmax = -1)

    divider = make_axes_locatable(ax)
    cax1 = divider.append_axes('right', size = '5%', pad = 0.05)
    fig.colorbar(img1, cax=cax1, label = "Fraction")
    ax.set_xlabel(r'Time (Myr)')

    plt.minorticks_on()
    plt.tight_layout()
    fig.savefig('stellar_O_2d_hist.png')

    return

def plot_stellar_separation(time, data, galaxy,
                            field    = 'O_Fraction',
                            property = 'median',
                            phases   = ['CNM'],
                            figdim=None,
                            field_types = None,
                            labels = None, line_styles = {},
                            xlim = None, ylim = None,
                            annotate_text = None, ylabels = None):
    if labels is None:
        labels = {}

    for k in phases + [field] + [property]:
        if not (k in list(labels.keys())):
            labels[k] = k

    if figdim is None:
        if len(phases) == 1:
            nrow, ncol = 1, 1
        else:
            ncol = 3
            nrow = 2
    else:
        nrow, ncol = figdim

    data_list = np.array(np.sort([x for x in list(data.keys()) if 'DD' in x]))
    data_list = data_list[:len(time)]

    time_data = load_abundance_data(data, data_list, [field], [property],
                                    phases = phases, field_types = field_types)

    # now need to load the data output to get the stellar values

    #
    if nrow*ncol > 1:
        fig, all_axes = plt.subplots(nrow,ncol,sharex=True,sharey=True)
        fig.set_size_inches(ncol*5, nrow*5)
        fig.subplots_adjust(hspace=0.0,wspace=0.0)
    else:
        fig, ax = plt.subplots()
        fig.set_size_inches(6,6)

    axi = axj = 0

    creation_time = galaxy.df['creation_time'].convert_to_units('Myr').value

    for i, phase in enumerate(phases):
        if nrow*ncol > 1:
            ax = all_axes[(axi,axj)]

        property_value = time_data[field][phase][property]

        if '_Fraction' in field:
            star_val = np.log10(galaxy.df['particle_' + field.strip('_Fraction') + '_fraction'])
        else:
            star_val = galaxy.df['particle_' + field]

        select = creation_time > 0.0 #do all for now (np.min(creation_time) + 100.0)

        distance = star_val[select] - np.interp(creation_time[select],time,property_value)

        ax.scatter(creation_time[select] - np.min(creation_time),
                   distance, alpha = 0.75, color = 'black', s = 20)
        ax.set_xlim(0.0, 500.0)
        ax.set_ylim(-2.5,2.5)
        ax.plot( ax.get_xlim(), [0.0,0.0], color = 'black', lw = line_width, ls = '--')
        plt.minorticks_on()

        xy = (200.0, ax.get_ylim()[1] - 0.35)
        #q3dist = q3 - median_distance
        #q1dist = median_distance - q1
#        ax.annotate(labels[phase] + " : %.2f + %.2f - %.2f"%(median_distance,q3dist,q1dist), xy=xy,xytext=xy)
        #d1, d9 = np.percentile(distance, [0.1,0.9])
        select = creation_time > (np.min(creation_time) + 120.0)

        ax.annotate("Median Dist. = %.2f"%(np.median(np.abs(distance[select]))), xy=xy,xytext=xy)
        xy = (200.0, ax.get_ylim()[1] - 0.7)
        q1, q3 = np.percentile(np.abs(distance[select]), [25,75])
        print(q3, q1, q3-q1, np.size(distance[select][distance[select]<0])/(1.0*np.size(distance[select])), np.size(distance[select][distance[select]>0]))
        ax.annotate("IQR = %.2f"%(q3-q1), xy=xy,xytext=xy)



        axj = axj + 1
        if axj >= ncol:
            axj = 0
            axi = axi + 1

    #plt.tight_layout()
    if nrow*ncol>1:
        for i in np.arange(ncol):
            all_axes[(nrow-1,i)].set_xlabel(r'Time (Myr)')
        for i in np.arange(nrow):
            all_axes[(i,0)].set_ylabel(r'Distance From ' + labels[property] + '[dex]')
    else:
        
        ax.set_xlabel(r'Time (Myr)')
        if ylabels is None:
            ax.set_ylabel(r'Distance to ' + labels[phase] + ' ' + labels[property] + ' [dex]')
        else:
            ax.set_ylabel(ylabels[0])
    if nrow*ncol == 1:
        plt.tight_layout()
    fig.savefig("stellar_distance_to_median.png")

    return

def plot_abundace_resolution_study():

    fields = ['O_Fraction','Ba_Fraction']
    property_list = ['d9_d1_range']
    phases = ['CNM','WNM','WIM','HIM']

    simulations = {'3pcH2' : './3pc_H2/abundances/gas_abundances.h5',
                   '6pcH2' : './6pc_H2/abundances/gas_abundances.h5'}

    fig, ax = plt.subplots(2,2,sharex=True,sharey=True)
    fig.set_size_inches(12,12)
    fig.subplots_adjust(hspace=0.0,wspace=0.0)

    all_time_data = {}
    for sim in simulations:
        data = h5py.File(simulations[sim])

        data_list = np.sort([x for x in list(data.keys()) if 'DD' in x])
        times     = np.array([float(x.strip('DD')) for x in data_list])

        all_time_data[sim] = load_abundance_data(data, data_list,
                                                 fields, property_list)
        all_time_data[sim]['time'] = times - times[0]

    ls = {'CNM':'-','WNM':'-','WIM':'-','HIM':':'}

    for phase in phases:

        ax[(0,0)].plot(all_time_data['3pcH2']['time'],
                       all_time_data['3pcH2']['O_Fraction'][phase]['d9_d1_range'],
                       color = color_dict[phase], ls = ls[phase], lw = line_width)
        ax[(0,1)].plot(all_time_data['6pcH2']['time'],
                       all_time_data['6pcH2']['O_Fraction'][phase]['d9_d1_range'],
                       color = color_dict[phase], ls = ls[phase], lw = line_width)
        ax[(1,0)].plot(all_time_data['3pcH2']['time'],
                       all_time_data['3pcH2']['Ba_Fraction'][phase]['d9_d1_range'],
                       color = color_dict[phase], ls = ls[phase], lw = line_width)
        ax[(1,1)].plot(all_time_data['6pcH2']['time'],
                       all_time_data['6pcH2']['Ba_Fraction'][phase]['d9_d1_range'],
                       color = color_dict[phase], ls = ls[phase], lw = line_width)

    for a1 in ax:
        for a2 in a1:
            a2.set_xlim(0, 500)
            a2.set_ylim(0, 3)
            a2.minorticks_on()

    for i in [0,1]:
        ax[(i,0)].set_ylabel(r'Inner Decile Range [dex]')
        ax[(1,i)].set_xlabel(r'Time (Myr)')

    x = 300
    y = 2.7
    size = 20
    ax[(0,0)].text(x, y, r'O  - 3.6 pc', color = 'black', size = size)
    ax[(1,0)].text(x, y, r'Ba - 3.6 pc', color = 'black', size = size)
    ax[(0,1)].text(x, y, r'O  - 7.2 pc', color = 'black', size = size)
    ax[(1,1)].text(x, y, r'Ba - 7.2 pc', color = 'black', size = size)

    plt.minorticks_on()
    fig.savefig('O_Ba_resolution_study.png')

    return

def plot_abundance_evolution(time, data,
                             fields = ['O_Fraction','Ba_Fraction'],
                             property_list = ['median','IQR'],
                             phases = ['CNM','WNM','WIM','HIM'],
                             figdim=None,
                             field_types = None,
                             labels = None, line_styles = None,
                             xlim = None, ylim = None,
                             annotate_text = None, fsize=5):

    if figdim is None:
        ncol = len(property_list)
        nrow = len(fields)
    else:
        nrow, ncol = figdim

    data_list = np.array(np.sort([x for x in list(data.keys()) if 'DD' in x]))
    data_list = data_list[:len(time)]



    start = cpu_time.time()

    time_data = load_abundance_data(data, data_list, fields, property_list,
                                    phases = phases, field_types = field_types)

    print("DATA LOADING TOOK %2.2E"%(cpu_time.time()- start))

    fig, ax = plt.subplots(nrow,ncol)
    fig.set_size_inches(fsize*ncol, fsize*nrow)

    axi = 0
    axj = 0
    xval = np.array(time)
    xval = xval - xval[0]
    for field in fields:

        for property in property_list:
            if len(fields) == 1:
                axindex = axj
            else:
                axindex = (axi,axj)

            for phase in phases:
                yval = time_data[field][phase][property] # - time_data[field]['CNM'][property]


                if phase in list(line_styles.keys()):
                    ls = line_styles[phase]
                else:
                    ls = '-'

                ax[axindex].plot(xval, yval, lw = line_width,
                                 ls = ls, color = color_dict[phase],
                                 label = phase)

            ax[axindex].set_xlabel(r'Time (Myr)')

            if hasattr(labels[property],"keys"):
                ax[axindex].set_ylabel(labels[property][field])
            else:
                ax[axindex].set_ylabel(labels[property])

            if xlim is None:
                ax[axindex].set_xlim(xval[0],xval[-1])
            else:
                ax[axindex].set_xlim(xlim)

            if not (ylim is None):
                if np.size(fields) == 1:
                    ymin,ymax = ylim[axindex]
                else:
                    ymin,ymax = ylim[axindex[0]][axindex[1]]
                ax[axindex].set_ylim(ymin,ymax)


            axj = axj + 1
            if axj >= ncol:
                axj = 0
                axi = axi + 1

    if len(fields) == 1:
        indexzero=0
    else:
        indexzero=(0,0)

    ax[indexzero].legend(loc='lower right', ncol=2)

    if not annotate_text is None:
        for axi in len(fields):
            if len(fields) == 1:
                axindex0 = axi
                axindex1 = axi
            else:
                axindex0 = (axi,1)
                axindex1 = (axi,0)

            xy = annotate_text[axindex0]
            ax[axindex1].annotate( annotate_text[axindex0], xy,xy)

    plt.minorticks_on()
    #plt.tight_layout()
    fig.savefig( '_'.join(fields + property_list) + '_evolution.png')

    return fig, ax

if __name__ == "__main__":


#    plot_abundace_resolution_study()

    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = './gas_abundances.h5'

    data = h5py.File(filename)

    #time = np.arange(0.0, 10.0*(len( [x for x in data.keys() if 'DD' in x])-1) + 0.1, 10.0)

    time = [float(x.strip('DD')) for x in np.sort(list(data.keys())) if 'DD' in x]

    galaxy = Galaxy('DD0619')

    plot_stellar_2d_hist(galaxy, 'O_Fraction')


    plot_stellar_separation(time, data, galaxy,
                             field    = 'O_Fraction',
                           property = 'median',
                            phases   = ['CNM'],
                            labels = {'median' : 'Median'},
                            ylabels = [r"[O/H] - [O/H]$_{\rm CNM}$"])


    plot_abundance_evolution(time, data,
                             fields = ['O_Fraction','Ba_Fraction'],
                             phases = ['CNM','WNM','WIM','HIM'],
                             property_list = ['median','mean_median_distance','IQR','d9_d1_range'],
                             labels = {'median' : r'log(Median)',
                                       'IQR'    : r'Inner Quartile Range [dex]',
                                       'mean_median_distance' : r'log(Mean) - log(Median) [dex]',
                                       'd9_d1_range' : r'Inner Decile Range [dex]'},
                             line_styles = {'HIM' : ':'},
                              ylim = [ [(-8,-2),(0.,1.5),(0.0,1.5),(0.0,3.0)], [(-16,-10),(0.,1.5),(0.0,1.5),(0.0,3.0)]],
                            # ylim = [ [(-1,8),(-3,3),(-3,3.0)], [(-1,8),(-3,3),(-3,3.0)]],
                             annotate_text = [ ('Oxygen', (20.0,-2.5)),   ('Barium',(20,-10.5))],
                             xlim = (0.0, 500.0))

    plot_abundance_evolution(time, data,
                             fields = ['O_Fraction','N_Fraction'],
                             phases = ['CNM','WNM','WIM','HIM'],
                             property_list = ['median','IQR','d9_d1_range'],
                             labels = {'median' : r'log(Median)',
                                       'IQR'    : r'Inner Quartile Range [dex]',
                                       'd9_d1_range' : r'Inner Decile Range [dex]'},
                             line_styles = {'HIM' : ':'},
                              ylim = [ [(-8,-2),(0.0,1.5),(0.0,3.0)], [(-9,-3),(0.0,1.5),(0.0,3.0)]],
                           # ylim = [ [(-1,8),(-3,3),(-3,3.0)], [(-1,8),(-3,3),(-3,3.0)]],
                             annotate_text = [ ('Oxygen', (20.0,-2.5)),   ('Nitrogen',(-8.5,-10.5))])

    
