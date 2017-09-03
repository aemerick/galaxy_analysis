import yt
import numpy as np

from galaxy_analysis import Galaxy
from galaxy_analysis.utilities import utilities

from astroML.time_series import ACF_EK
from astroML.time_series import ACF_scargle



#
# testing script for now
#

gal  = Galaxy('DD0114')

pos    = gal.ds.domain_center

rmax   = 50.0
sphere = gal.ds.sphere(pos, (rmax, 'pc')) # do small for now

x      = sphere['spherical_radius'].convert_to_units('pc').value
y      = sphere['Fe_over_H'].value

bins   = np.linspace(0.0, rmax, 25.0)

#acf, bins = ACF_scargle(x, y, dy = 0.0001, n_omega = 2**12, omega_max = np.pi/5.0) #, bins = bins)
acf, err, bins = ACF_EK(x, y, dy = 0.00001, bins = bins)

print acf
print bins

print "-------------------------------------------------------------------"
print "-------------------------------------------------------------------"

x      = sphere['spherical_radius'].convert_to_units('pc').value
y      = sphere['Fe_Fraction'].value

bins   = np.linspace(0.0, rmax, 25.0)
#acf, bins = ACF_scargle(x, y, dy = 0.0000001, n_omega = 2**12, omega_max = np.pi/5.0) #, bins = bins)
acf, err, bins = ACF_EK(x,y,dy=1.0E-8, bins = bins)

print acf
print bins

