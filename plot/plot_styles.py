"""
   Defining some constants for plots. Colors and
   general styles.
"""
import matplotlib as mpl
import numpy as np
mpl.use('Agg')
mpl.rcParams['hatch.linewidth'] = 8
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


#
# Define consistent colors for simulation names or
# other things that can be easily loaded
# here to make life easy....
#
color_dict = {
              # the four dwarf galaxy runs - to use when comp to each other
              '17' : orange, '25' : magenta,
              '30' : purple, '40' : black
              # ----
             }




def plot_histogram(ax, x, y, *args, **kwargs):
    """
    Wrapper on ax.plot from matplotlib to handle making
    a histogram plot with step locations that make sense.
    This works by taking in the actual bins as the x value,
    meaning that len(x) = len(y) + 1 must be true.
    """

    if len(x) == len(y):
        xnew = np.zeros(len(x) + 1)
        dx   = np.zeros(len(x) + 1)
        dx[1:-1]   = x[1:] - x[:-1]
        dx[0] = dx[1]
        dx[-1] = dx[-2]
        xnew[:-1] = x - 0.5*dx[:-1]
        xnew[-1] = x[-1] + 0.5*dx[-1]
        x = xnew

    if len(x) == len(y) + 1:
        ynew = np.zeros(len(x))
        ynew[:len(y)] = y
        ynew[-1] = y[-1]

    if 'drawstyle' in kwargs.keys():
        print "'drawstyle' included as kwarg as " +\
              kwargs['drawstyle'] +\
              "; are you sure? Changing this to 'steps-post'"

    kwargs['drawstyle'] = 'steps-post'

    p = ax.plot(x, ynew, *args, **kwargs)

    return p
