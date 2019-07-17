import yt
import numpy as np

from galaxy_analysis import Galaxy
from galaxy_analysis.utilities import utilities

from astroML.time_series import ACF_EK
from astroML.time_series import ACF_scargle

from matplotlib import rc

fsize = 17
rc('text', usetex=False)
rc('font', size=fsize)#, ftype=42)
line_width = 3
point_size = 30

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt



def simple_plot(x,y,name):

    fig, ax = plt.subplots()
    ax.plot(x,y, lw  = 3 , color = 'black')
    plt.minorticks_on()
    plt.tight_layout()
    fig.set_size_inches(8,8)
    plt.savefig(name)
    plt.close()

    return
#
# testing script for now
#
gal  = Galaxy('DD0126')

#    select a random SF region
n    = gal.df['number_density']
T    = gal.df['temperature']
select = (n >= 200.0) * (T < 200.0)
x,y,z = gal.df['x'][select], gal.df['y'][select], gal.df['z'][select]
pos   = np.array([x[0],y[0],z[0]])

rmax   = 50.0
dr     =   5.0
sphere = gal.ds.sphere(pos, (rmax, 'pc')) # do small for now


x      = sphere['spherical_radius'].convert_to_units('pc').value
y      = sphere['Fe_Fraction'].value

p = yt.ProfilePlot(sphere, "radius", ["Fe_Fraction",'Fe_over_H','O_over_Fe'],
                       weight_field = 'cell_volume', accumulation = False)
p.set_unit('radius','pc')
p.save()

bins   = np.arange(0.0, rmax+dr, dr)
#acf, bins = ACF_scargle(x, y, dy = 0.0000001, n_omega = 2**12, omega_max = np.pi/5.0) #, bins = bins)
acf, err, bins = ACF_EK(x,y,dy=1.0E-8, bins = bins)

print(acf)
print(bins)
simple_plot(0.5*(bins[1:]+bins[:-1]), acf, 'Fe_Fraction_acf.png')

print('----------------------------------------------------------')
print('----------------------------------------------------------')


x      = sphere['spherical_radius'].convert_to_units('pc').value
y      = sphere['Fe_over_H'].value

bins   = np.arange(0.0, rmax+dr, dr)

#acf, bins = ACF_scargle(x, y, dy = 0.0001, n_omega = 2**12, omega_max = np.pi/5.0) #, bins = bins)
acf, err, bins = ACF_EK(x, y, dy = 0.00001, bins = bins)

simple_plot(0.5*(bins[1:]+bins[:-1]), acf, 'Fe_over_H_acf.png')

print(acf)
print(bins)

print("-------------------------------------------------------------------")
print("-------------------------------------------------------------------")

