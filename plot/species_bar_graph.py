import deepdish as dd
import numpy as np

from matplotlib import rc

fsize = 17
rc('text', usetex=False)
rc('font', size=fsize)#, ftype=42)
line_width = 3
point_size = 30

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
import glob


###
###
from galaxy_analysis.utilities import utilities
###
###

def species_bar_graph(name, data, fraction = True, disk_only = False, outname = None,
                      display_total = None, show_individual_amounts = False):
   """
   Uses an analysis output to give a bar graph
   """

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
   exclude = ['Total','HI','HeI','HeII','HII','H2','Metals']

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
   width  = 0.4

   if not disk_only:
       colors = {'Disk' : 'purple', 'stars' : 'gold',
                 'Halo' : 'black'}
       barplot = {}
       bottom  = np.zeros(N)
       sum     = np.zeros(N)
       # plot!
       for f in ['Disk','stars','Halo']:
           bar_values = np.array([masses[f][k]/total[k] for k in ordered_species])

           barplot[f] = ax.bar(index, bar_values, width,
                               color = colors[f], bottom = bottom)
           bottom += bar_values * 1.0
           sum    = sum + bottom

   else:
       colors = {'stars': 'gold', 'CNM' : 'blue', 'WNM' : 'green',
                 'WIM' : 'orange', 'HIM'  : 'red', 'Molecular' : 'black'}

       barplot = {}
       bottom  = np.zeros(N)
       # plot!
       for f in ['Molecular','CNM','WNM','WIM','HIM','stars']:
           bar_values = np.array([masses[f][k]/total[k] for k in ordered_species])

           barplot[f] = ax.bar(index, bar_values, width,
                               color = colors[f], bottom = bottom)
           bottom += bar_values
   #
   #
   #
   #
   #
   plt.tight_layout()
   ax.set_xticks(index)
   ax.set_xticklabels(ordered_species)
   if fraction:
       ax.set_ylim(0.0,1.0)
   else:
       ax.set_ylim(0.0, 50.0)
#   plt.minorticks_on()
   fig.savefig(outname)
   plt.close()

   return

if __name__ == "__main__":
    # example useage

    def _plot_both(name, data):
        species_bar_graph(name, data,
                              fraction      = True,
                              disk_only     = False,
                              display_total = True,
                              show_individual_amounts = False)

        species_bar_graph(name, data,
                          fraction      = True,
                          disk_only     = True,
                          display_total = True,
                          show_individual_amounts = False)
        return

    if False:
        name = 'DD0114'
        data = dd.io.load(name + '_galaxy_data.h5')
        _plot_both(name, data)

    else:
        files = np.sort(glob.glob('DD????_galaxy_data.h5'))
        names = [x[:6] for x in files]

        for f in names:
            data = dd.io.load(f + '_galaxy_data.h5')
            _plot_both(f,data)
            del(data)
