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

from galaxy_analysis.utilities import utilities as util

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


def generate_model_stars(m, z, abund = ['m_tot','m_metal']):
    """
    Makes a list of star objects from one zone model
    """
    all_star = [None]*np.size(m)

    if not 'm_tot' in abund:
        abund.extend(['m_tot'])
    if not 'm_metal' in abund:
        abund.extend(['m_metal'])

    ele = {}
    for k in abund:
        if k is 'm_tot':
            val = 1.0
        elif k is 'm_metal':
            val = z[0]
        else:
            val = 0.0

        ele[k] = val

    for i in np.arange(np.size(m)):
        s = star.Star(M=m[i],Z=z[i], abundances=ele)
        s.set_SNII_properties()
        all_star[i] = s
    return all_star

def check_all_masses(ds, data, d0 = None, time_cut = -1.0):

    bm = data['birth_mass'].value
    pm = data['particle_mass'].convert_to_units('Msun').value
    z  = data['metallicity_fraction'].value
    pt = data['particle_type']

    elements = util.species_from_fields(ds.field_list)

    all_stars = generate_model_stars(bm,z, abund = elements)

    lifetime = data['dynamical_time'].convert_to_units('Myr')
    birth    = data['creation_time'].convert_to_units('Myr')
    age      = ds.current_time.convert_to_units('Myr') - birth

    model_wind_ejecta = {} # total_wind_ejecta
    for k in all_stars[0].wind_ejecta_masses().keys():
        model_wind_ejecta[k] = np.array([x.wind_ejecta_masses()[k] for x in all_stars])
    model_sn_ejecta = {}
    for k in all_stars[0].sn_ejecta_masses.keys():
        model_sn_ejecta[k] = np.array([x.sn_ejecta_masses[k] for x in all_stars])

    # correct for AGB stars that haven't died
    AGB = (bm < 8.0) * (pt == 11)
    select = (bm > 8.0) * (pt == 11)
    factor = age / lifetime
    factor[factor>1.0] = 1.0

    time_select = birth > time_cut

    for k in model_wind_ejecta.keys():
        model_wind_ejecta[k][AGB]        = 0.0
        model_sn_ejecta[k][ (pt == 11) ] = 0.0

        model_wind_ejecta[k][select] = model_wind_ejecta[k][select]*factor[select]

    total_model_ejecta = {}
    for k in model_wind_ejecta.keys():
        total_model_ejecta[k] = np.sum(model_sn_ejecta[k][time_select]) + np.sum(model_wind_ejecta[k][time_select])

    # now do this for the individual abundances on grid:
    grid_masses = {}
    for k in model_wind_ejecta.keys():
        if k is 'm_tot' or k is 'm_metal':
            continue

        grid_masses[k] = np.sum(data[k + '_Density'] * ds.mass_unit / ds.length_unit**3 *\
                                   data['cell_volume']).convert_to_units('Msun').value

        if not (d0 is None):
            grid_masses[k] = grid_masses[k] - np.sum(d0[k + '_Density'] * ds.mass_unit / ds.length_unit**3 *\
                                       d0['cell_volume']).convert_to_units('Msun').value

    print total_model_ejecta
    print grid_masses

    print grid_masses.keys()
    print "Element Total_on_Grid Total_model_mass Percent_error"
    for k in grid_masses.keys():
        error =100 *  (grid_masses[k] - total_model_ejecta[k] ) / total_model_ejecta[k]
        print "%2s     %8.8E %8.8E %3.3f"%(k,grid_masses[k], total_model_ejecta[k], error)

    return all_stars, model_sn_ejecta, model_wind_ejecta, total_model_ejecta




def check_wind_ejecta(ds, data):

    bm = data['birth_mass'].value
    pm = data['particle_mass'].convert_to_units('Msun').value
    z  = data['metallicity_fraction'].value
    pt = data['particle_type']

    elements = util.species_from_fields(ds.field_list)

    all_stars = generate_model_stars(bm,z, abund = elements)

    lifetime = data['dynamical_time'].convert_to_units('Myr')
    birth    = data['creation_time'].convert_to_units('Myr')
    age      = ds.current_time.convert_to_units('Myr') - birth

    # total wind ejecta over entire lifetime
    total_wind_ejecta = np.array([x.wind_ejecta_masses()['m_tot'] for x in all_stars])

    # correct for AGB stars that haven't died
    AGB = (bm < 8.0)
    model_wind_ejecta = total_wind_ejecta * 1.0
    model_wind_ejecta[ AGB * (pt == 11)] = 0.0

    # adjust wind to correct fraction given lifetime
    select = (bm > 8.0) * (pt == 11)
    factor = age / lifetime
    factor[factor>1.0] = 1.0
    model_wind_ejecta[select] = model_wind_ejecta[select] * factor[select]

    # load actual injection from simulation
    actual_wind_ejecta = data['wind_mass_ejected'].value

    # compute percent error
    model_wind_ejecta = model_wind_ejecta[age>1]
    actual_wind_ejecta = actual_wind_ejecta[age>1]
    error = (model_wind_ejecta - actual_wind_ejecta)
    error[model_wind_ejecta>0] = error[model_wind_ejecta>0]/model_wind_ejecta[model_wind_ejecta>0]

    error_mass = error[model_wind_ejecta>0]
    all = 1.0 * np.size(error_mass)
    print np.size( error_mass[ (np.abs(error_mass) < 0.05) ])/all
    print np.size( error_mass[ (np.abs(error_mass) < 0.10) ])/all
    print np.size( error_mass[ (np.abs(error_mass) < 0.15) ])/all
    print np.size( error_mass[ (np.abs(error_mass) < 0.20) ])/all
    print np.size( error_mass[ (np.abs(error_mass) < 0.25) ])/all

    #error_mass = error_mass[birth[model_wind_ejecta>0] > 110]
    #error_mass = error_mass[error_mass>0]


    print np.min(error_mass), np.max(error_mass), np.average(error_mass), np.median(error_mass)
    print error_mass
    select = (age>1)
    bm = bm[select]
    pm = pm[select]
    age = age[select]
    lifetime = lifetime[select]
    total_wind_ejecta = total_wind_ejecta[select]
    select = (model_wind_ejecta>0)
    bm = bm[select]
    pm = pm[select]
    age = age[select]
    lifetime = lifetime[select]
    model_wind_ejecta = model_wind_ejecta[select]
    actual_wind_ejecta = actual_wind_ejecta[select]
    total_wind_ejecta = total_wind_ejecta[select]
    print "BM   PM   Percent_error    Model_wind    actual_wind    lifetime_wind"
    for i in np.arange(np.size(error_mass)):
        print "%5.5f %3.3f %5.5f %5.5E %5.5E %5.5E"%(bm[i],pm[i],error_mass[i]*100,model_wind_ejecta[i], actual_wind_ejecta[i], total_wind_ejecta[i])
    print np.min(error_mass), np.max(error_mass), np.average(error_mass), np.median(error_mass)

   

#    print bm[error > 0.9], pm[error>0.9], pt[error>0.9]
#    print age[error>0.9]
#    print actual_wind_ejecta[error>0.9]
#    print model_wind_ejecta[error>0.9]
    #print actual_wind_ejecta[birth > 110]
    #print model_wind_ejecta[birth > 110]

    return

def compute_SNII_error(ds, data, uselog = True):

    pm = data['particle_mass'].convert_to_units('Msun').value
    bm = data['birth_mass'].value
    pt = data['particle_type']

    # select all particles that could have gone supernova
    select = (pt == 13) * (bm > 8.0) * (bm < 25.0)

    pm = pm[select]
    bm = bm[select]
    z  = data['metallicity_fraction'][select]

    elements = util.species_from_fields(ds.field_list)

    all_stars = generate_model_stars(bm, z, abund = elements)

    total_ejecta = np.zeros(np.size(bm))
    error        = np.zeros(np.size(bm))
    wind_error   = np.zeros(np.size(bm))
    sn_error     = np.zeros(np.size(bm))
    ej_frac = np.zeros(np.size(bm))
    for i,s in enumerate(all_stars):
        s.set_SNII_properties()
        wind = s.wind_ejecta_masses()
        sn   = s.sn_ejecta_masses

        total_ejecta[i] = wind['m_tot'] + sn['m_tot']
        error[i] = ( -1.0*(bm[i]-pm[i]) + total_ejecta[i]) / (total_ejecta[i])
        ej_frac[i] = (bm[i]-pm[i]) / total_ejecta[i]

        wind_error[i] = ( wind['m_tot'] / total_ejecta[i] )
        sn_error[i]   = ( sn['m_tot']   / total_ejecta[i] )

    snavg  , snstd   = np.average(sn_error), np.std(sn_error)
    windavg, windstd = np.average(wind_error), np.std(wind_error)

    # now plot cumulative distribution of positive error (error > 0 = missing mass)
    pos_error = error[error>0]

    fig, ax = plt.subplots()
    if uselog:
        xdata = np.log10(pos_error)
        bins  = np.arange(-4, 1.0, 0.025)
    else:
        xdata = pos_error
        bins  = np.linspace(0.0, 1.0, 200)

    hist,bins = np.histogram(np.log10(ej_frac), bins = bins)
    cent = (bins[1:] + bins[:-1])*0.5

    ax.plot(cent, np.cumsum(hist) / (1.0*np.sum(hist)), lw = 3, color = 'black')
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

    ax.set_xlabel('Fraction of Mass Actually Injected')
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
    print energy_error
    ax.set_ylim([0,np.max(hist)])
    ax.set_xlabel('Error in Ejected Mass')
    ax.set_ylabel('Counts')
    fig.set_size_inches(8,8)
    plt.tight_layout()
    plt.minorticks_on()

    fig.savefig('sn_mass_error.png')

    return error, fig, ax


if __name__=="__main__":
    name_list = np.sort(glob.glob('DD????/DD????'))
    try:
        ds = yt.load(name_list[-1])
    except:
        print "Could not load ", name_list[-1], " trying the next one"
        ds = yt.load(name_list[-2])
    data = ds.all_data()

    if ('enzo','wind_mass_ejected') in ds.field_list or\
       ('io','wind_mass_ejected') in ds.field_list:
        try:
            check_wind_ejecta(ds,data)
        except:
            print "failing in wind ejecta"
        try:
            error, fig, ax = compute_SNII_error(ds,data, uselog=True)
        except:
            print "failing in SNII check"
#    ds0 = yt.load('./../lowres/DD0035/DD0035')
#    d0  = ds0.all_data()


    check_all_masses(ds,data) #, d0 = d0)
