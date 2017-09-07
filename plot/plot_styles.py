"""
   Defining some constants for plots. Colors and
   general styles.
"""
import matplotlib as mpl
mpl.use('Agg')

from matplotlib import rc, cm
viridis = cm.get_cmap('viridis')
magma   = cm.get_cmap('magma')
plasma  = cm.get_cmap('plasma')

fsize = 17
rc('text', usetex=False)
rc('font', size=fsize)#, ftype=42)
line_width = 3.5
point_size = 30

purple  = '#7d03a8'
magenta = '#cb4679'
blue    = '#0c0887'
orange  = '#fdc328'
black   = 'black'
green   = 'green'
lw      = 4.5
