import numpy as np

#
# load data from plot digitizer
#
vr0   = np.genfromtxt('prantzos_f4_v0.dat')
vr150 = np.genfromtxt('prantzos_f4_v150.dat')
vr300 = np.genfromtxt('prantzos_f4_v300.dat')


#
# grab all and subsample with linear interpolation
#
npoints = 100
Fe_H = np.linspace(-4.0, 0.5, npoints)
vr0_resample   = np.interp(Fe_H, vr0[:,0], vr0[:,1])
vr150_resample = np.interp(Fe_H, vr150[:,0], vr150[:,1])
vr300_resample = np.interp(Fe_H, vr300[:,0], vr300[:,1])

#
# ensure this sums to 1
#
vr_sum           = vr0_resample + vr150_resample + vr300_resample
vr0_resample     = vr0_resample / vr_sum
vr150_resample   = vr150_resample / vr_sum
vr300_resample   = vr300_resample / vr_sum

#
# save to file
#
f = open("prantzos_fig4_fractions.dat",'w')

f.write("#[Fe/H] v0 v150 v300\n")

for i in np.arange(npoints):

    f.write("%5.3f %6.4f %6.4f %6.4f\n"%(Fe_H[i], vr0_resample[i],
                                         vr150_resample[i], vr300_resample[i]))
f.close()

