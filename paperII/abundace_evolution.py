from galaxy_analysis.plot.plot_styles import *

import numpy as np
import matplotlib.pyplot as plt
import deepdish as dd
import h5py, glob, sys

from galaxy_analysis.utilities import utilities

# temporary
import time as cpu_time

def plot_abundance_evolution(time, data,
                             fields = ['O_Fraction','Ba_Fraction'],
                             property_list = ['median','IQR'],
                             phases = ['CNM','WNM','WIM','HIM'],
                             figdim=None,
                             field_types = None,
                             labels = None, line_styles = None,
                             xlim = None, ylim = None,
                             annotate_text = None):

    if figdim is None:
        ncol = len(property_list)
        nrow = len(fields)
    else:
        nrow, ncol = figdim

    data_list = np.array(np.sort([x for x in data.keys() if 'DD' in x]))
    data_list = data_list[:len(time)]



    start = cpu_time.time()

    time_data = {}
    for i, field in enumerate(fields):
        time_data[field] = {}

        if field_types is None:
            if '_Fraction' in field:
                ft = 'mass_fraction'
            elif '_over_' in field:
                ft = 'abundance'
            else:
                print "Cannot properly determing field type for " + field
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
                    time_data[field][phase][property] = np.log10(time_data[field][phase]['Q3']) -\
                                                        np.log10(time_data[field][phase]['Q1'])
                elif property == 'mean_median_distance':
                    for p in ['median','mean']:
                        if not p in time_data[field][phase].keys():
                            field_path = [phase,ft,field,p]
                            time_data[field][phase][p] = utilities.get_property(field_path,
                                                                                file_list = data.filename,
                                                                                data_list = data_list,
                                                                                self_contained = True)
                            time_data[field][phase][p] = np.array(time_data[field][phase][p], dtype = np.float)
                    time_data[field][phase]['mean_median_distance'] =\
                                               np.log10(time_data[field][phase]['mean']) -\
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

                    time_data[field][phase][property] = np.log10(time_data[field][phase]['decile_9']) -\
                                                        np.log10(time_data[field][phase]['decile_1'])
                else:
                    field_path = [phase,ft,field,property]
                    time_data[field][phase][property] = utilities.get_property(field_path,
                                                                               file_list = data.filename,
                                                                               data_list = data_list,
                                                                               self_contained = True)
                    time_data[field][phase][property] = np.array(time_data[field][phase][property], dtype = np.float)

                    if property in ['median','mean']:
                        time_data[field][phase][property] = np.log10(time_data[field][phase][property])

    print "DATA LOADING TOOK %2.2E"%(cpu_time.time()- start)

    fig, ax = plt.subplots(nrow,ncol)
    fig.set_size_inches(5*ncol, 5*nrow)

    axi = 0
    axj = 0
    xval = np.array(time)
    xval = xval - xval[0]
    for field in fields:

        for property in property_list:
            axindex = (axi,axj)

            for phase in phases:
                yval = time_data[field][phase][property] # - time_data[field]['CNM'][property]


                if phase in line_styles.keys():
                    ls = line_styles[phase]
                else:
                    ls = '-'

                ax[axindex].plot(xval, yval, lw = line_width,
                                 ls = ls, color = color_dict[phase],
                                 label = phase)

            ax[axindex].set_xlabel(r'Time (Myr)')
            ax[axindex].set_ylabel(labels[property])

            if xlim is None:
                ax[axindex].set_xlim(xval[0],xval[-1])
            else:
                ax[axindex].set_xlim(xlim)

            if not (ylim is None):
                ymin,ymax = ylim[axi][axj]
                ax[axindex].set_ylim(ymin,ymax)


            axj = axj + 1
            if axj >= ncol:
                axj = 0
                axi = axi + 1

    ax[(0,0)].legend(loc='lower right', ncol=2)

    if not annotate_text is None:
        for axi in [0,1]:
            xy = annotate_text[axi][1]
            ax[(axi,0)].annotate( annotate_text[axi][0], xy,xy)

    plt.minorticks_on()
    plt.tight_layout()
    fig.savefig( '_'.join(fields + property_list) + '_evolution.png')

    return

if __name__ == "__main__":

    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = './gas_abundances.h5'

    data = h5py.File(filename)

    #time = np.arange(0.0, 10.0*(len( [x for x in data.keys() if 'DD' in x])-1) + 0.1, 10.0)

    time = [float(x.strip('DD')) for x in np.sort(data.keys()) if 'DD' in x]

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
