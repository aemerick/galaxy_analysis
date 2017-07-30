from matplotlib import rc

fsize = 17
rc('text', usetex=False)
rc('font', size=fsize)#, ftype=42)
line_width = 3
point_size = 30

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
import yt
import numpy as np
import glob
from onezone import star

#
# Step 1: load "final" and "initial" simulations
# Step 2: Load all massive star particle remnants, final mass, initial mass, etc.
# Step 3: Run each through the star particle class in onezone
#         (make sure sn ejecta for m > 25 is zero)
# Step 4: Compute for each SN ejecta, wind ejecta, and total ejecta
# Step 5: Compute:
#          % error = ((Birth - Current) - (sn_ej+wind_ej)) / (sn_ej+wind_ej)
#
# Step 6: Plot cumulative distribution of SNII remnants error


def generate_model_stars(m, z):
    """
    Makes a list of star objects from one zone model
    """
    all_star = [None]*np.size(m)
    for i in np.arange(np.size(m)):
        all_star[i] = star.Star(M=m[i],Z=z[i])
    return all_star


def compute_SNII_error(ds, data, uselog = True):

    pm = data['particle_mass'].convert_to_units('Msun').value
    bm = data['birth_mass'].value
    pt = data['particle_type']

    # select all particles that could have gone supernova
    select = (pt == 13) * (bm > 8.0) * (bm < 25.0)

    pm = pm[select]
    bm = bm[select]
    z  = data['metallicity_fraction'][select]

    all_stars = generate_model_stars(bm, z)

    total_ejecta = np.zeros(np.size(bm))
    error        = np.zeros(np.size(bm))
    wind_error   = np.zeros(np.size(bm))
    sn_error     = np.zeros(np.size(bm))
    for i,s in enumerate(all_stars):
        s.set_SNII_properties()
        wind = s.wind_ejecta_masses()
        sn   = s.sn_ejecta_masses

        total_ejecta[i] = wind['m_tot'] + sn['m_tot']
        error[i] = ( -1.0*(bm[i]-pm[i]) + total_ejecta[i]) / (total_ejecta[i])

        wind_error[i] = ( wind['m_tot'] / total_ejecta[i] )
        sn_error[i]   = ( sn['m_tot']   / total_ejecta[i] )

    snavg  , snstd   = np.average(sn_error), np.std(sn_error)
    windavg, windstd = np.average(wind_error), np.std(wind_error)

    # now plot cumulative distribution of positive error (error > 0 = missing mass)
    pos_error = error[error>0]

    fig, ax = plt.subplots()
    if uselog:
        xdata = np.log10(pos_error)
        bins  = np.arange(-2, 0.05, 0.025)
    else:
        xdata = pos_error
        bins  = np.linspace(0.0, 1.0, 200)

    hist,bins = np.histogram(xdata, bins = bins)
    cent = (bins[1:] + bins[:-1])*0.5

    ax.plot(cent, np.cumsum(hist)/(1.0*np.sum(hist)), lw = 3, color = 'black')
    ylim = [0.0, 1.05]
    ax.set_ylim(ylim)

    def _plot_line(x, color, ls, log, label):
        if log:
            if x <= 0:
                return
            x = np.log10(x)

        ax.plot([x,x],ylim, color = color, ls = ls, label = label, lw = 2)
        return

#    _plot_line(snavg, 'blue', '-', uselog, 'SN fraction')
#    _plot_line(snavg-snstd, 'blue', '-', uselog, None)
#    _plot_line(snavg+snstd, 'blue', '-', uselog, None)
#    _plot_line(windavg, 'purple', '-', uselog, 'Wind fraction')
#    _plot_line(windavg - windstd, 'purple', '--', uselog, None)
#    _plot_line(windavg + windstd, 'purple', '--', uselog, None)

    ax.set_xlabel('Error in Ejected Mass')
    ax.set_ylabel('Fraction of SN')
    fig.set_size_inches(8,8)
    plt.tight_layout()
    plt.minorticks_on()

    fig.savefig('sn_cum_mass_error.png')
    plt.close()
#
#
#   histogram
#
#
    fig, ax = plt.subplots()
    if uselog:
	xdata = np.log10(pos_error)
        bins  = np.arange(-2, 0.05, 0.025)
    else:
	xdata = pos_error
        bins  = np.linspace(0.0, 1.0, 200)

    hist,bins = np.histogram(xdata, bins = bins)
    cent = (bins[1:] + bins[:-1])*0.5

    ax.plot(cent, hist, lw = 3, color = 'black')

    energy_error  = ( np.sum(pos_error)) / (np.size(pos_error)*1.0)

    ax.plot([np.average(pos_error),np.average(pos_error)], [0,np.max(hist)], color = 'black' ,ls = '--', lw = 3)
    ax.annotate("Energy Error = %0.2f percent"%(100*energy_error), xy=(0.5,0.9*np.max(hist)),
                                           xytext=(0.5,0.9*np.max(hist)))
    ax.set_ylim([0,np.max(hist)])
    ax.set_xlabel('Error in Ejected Mass')
    ax.set_ylabel('Counts')
    fig.set_size_inches(8,8)
    plt.tight_layout()
    plt.minorticks_on()

    fig.savefig('sn_mass_error.png')

    return error, fig, ax


if __name__=="__main__":
    name = np.sort(glob.glob('DD????/DD????'))[-1]
    ds = yt.load(name)
    data = ds.all_data()

    error, fig, ax = compute_SNII_error(ds,data, uselog=False)
