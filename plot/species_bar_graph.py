import deepdish as dd
import numpy as np

from galaxy_analysis.plot.plot_styles import *
import matplotlib.pyplot as plt
import glob

import sys

from joblib import Parallel, delayed
import multiprocessing


###
###
from galaxy_analysis.utilities import utilities
###
###

def species_bar_graph(name, data, fraction = True, disk_only = False, outname = None,
                      display_total = None, show_individual_amounts = False, **kwargs):
   """
   Uses an analysis output to give a bar graph
   """

   if not 'alpha' in kwargs:
       kwargs['alpha'] = 0.8

   # set default output name
   if outname is None:
       s1, s2 = '', ''
       if fraction:
           s1 = '_fraction'
       else:
           s1 = '_total'

       if disk_only:
           s2 = '_disk_only'

       outname = name + '_species_bar' + s1 + s2 + '.png'

   # set defaults for displaying total amount
   if display_total is None:
       display_total = fraction

   # construct the plot
   exclude = ['Total','HI','HeI','HeII','HII','H2','Metals','Total Tracked Metals']

   masses = data['gas_meta_data']['masses']
   all_fields  = masses['FullBox'].keys()
   species     = [x for x in all_fields if (not (x in exclude))]
   ordered_species = utilities.sort_by_anum(species)
   N = len(ordered_species)

   total  = {}

   # total normalizes the data - set to 1 if not plotting fraction
   if fraction:
       if disk_only: # use disk mass + stars
           for species in ordered_species:
               total[species] = masses['Disk'][species] + masses['stars'][species]
       else: # full box, stars, and mass loss in box
           for species in ordered_species:
               total[species] = 0.0
               for k in ['FullBox', 'OutsideBox', 'stars']:
                   total[species] += masses[k][species]
   else: # otherwise set normalization to 1.0
       for species in ordered_species:
           total[species] = 1.0

   fig, ax = plt.subplots()
   fig.set_size_inches(16,8) # may need to make this longer

   index  = np.arange(N)
   width  = 0.75

   if not disk_only:
       fields = ['Disk','stars','Halo']

       colors = {'Disk' : 'purple', 'stars' : 'gold',
                 'Halo' : 'black'}
       colors = {'Disk': plasma(0.0), 'stars': plasma(1.0),
                 'Halo': plasma(0.5)}
       barplot = {}
       bottom  = np.zeros(N)
       sum     = np.zeros(N)
       # plot!
       for f in fields:
           bar_values = np.array([masses[f][k]/total[k] for k in ordered_species])

           barplot[f] = ax.bar(index, bar_values, width,
                               color = colors[f], bottom = bottom, label = f, **kwargs)
           bottom += bar_values * 1.0
           sum    = sum + bottom

   else:
       fields = ['Molecular', 'CNM', 'WNM', 'WIM', 'HIM', 'stars']
       colors = {'stars': 'gold', 'CNM' : 'blue', 'WNM' : 'green',
                 'WIM' : 'orange', 'HIM'  : 'red', 'Molecular' : 'black'}

       colors = {'Molecular' : plasma(0.0), 'CNM' : plasma(1.0/5.0),
                 'WNM' : plasma (2.0/5.0), 'WIM' : plasma(3.0/5.0),
                 'HIM' : plasma(4.0/5.0), 'stars' : plasma(1.0)}

       barplot = {}
       bottom  = np.zeros(N)
       # plot!
       for f in fields:
           bar_values = np.array([masses[f][k]/total[k] for k in ordered_species])

           barplot[f] = ax.bar(index, bar_values, width,
                               color = colors[f], bottom = bottom,
                               label = f, **kwargs)
           bottom += bar_values

   nfields = len(fields)

   if display_total:
       rects = barplot[np.array(fields)[-1]]
       totals = np.array([ total[k] for k in ordered_species])
       for i, rect in enumerate(rects):
           ax.text(rect.get_x() + rect.get_width()/2.0, 0.85*bottom[i],
                   '%.1E' % float(totals[i]), ha='center', va='bottom', color = 'white',
                   rotation='vertical')

   if show_individual_amounts:
       # put code here to add amounts to each field
       print "cannot include individual amounts yet"
   #
   #
   #
   #
   #
   plt.tight_layout()
   ax.set_xticks(index)
   ax.set_xticklabels(ordered_species)

   # make legend, reverse label ordering
   handles, labels = ax.get_legend_handles_labels()
   ax.legend(handles[::-1], labels[::-1], loc='best')

   if fraction:
       ax.set_ylim(0.0,1.0)
   else:
       ax.semilogy()

   # turn on minorticks, but keep x axis ticks off
   plt.minorticks_on()
   ax.tick_params(axis='x',which='minor',bottom='off')

   fig.savefig(outname)
   plt.close()

   return

if __name__ == "__main__":
    # example useage

    def _plot_both(name, data):
        for fraction in [True, False]:
            species_bar_graph(name, data,
                              fraction      = fraction,
                              disk_only     = False,
                              display_total = True,
                              show_individual_amounts = False)

            species_bar_graph(name, data,
                              fraction      = fraction,
                              disk_only     = True,
                              display_total = True,
                              show_individual_amounts = False)
        return

    if False:
        name = 'DD0126'
        data = dd.io.load(name + '_galaxy_data.h5')
        _plot_both(name, data)

    else:
        files = np.sort(glob.glob('DD????_galaxy_data.h5'))
        names = [x[:6] for x in files]

        if len(sys.argv) == 1:
            # assume not parallel
            _plot_both(names[-1], dd.io.load(names[-1] + '_galaxy_data.h5'))
        elif len(sys.argv) == 3:

            i,j = int(sys.argv[1]), int(sys.argv[2])
            n_jobs = multiprocessing.cpu_count()

            Parallel(n_jobs = n_jobs)(\
                    delayed(_plot_both)(x, dd.io.load(x+'_galaxy_data.h5')) for x in names[i:j])
#            for f in names:
#                data = dd.io.load(f + '_galaxy_data.h5')
#            _plot_both(f,data)
#            del(data)
