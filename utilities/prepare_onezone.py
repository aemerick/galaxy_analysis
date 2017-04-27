__author__ = "Andrew Emerick"
"""
   Prepare a onezone chemical evolution simulation
   run script. Does the following:

       1) Loads parameter file from most recent data
          output (via yt)

       2) Computes sfr from this output, using 
          pre-defined time-step sizes (10 Myr)

       3) Checks if galaxy analysis file exists,
          if not make some assumptions about time

       4) Outputs text file with:
           # t sim_time SFR
          where "t" stars at zero, and is time 
          of first star formation onward, sim_time
          is translation between current time and 
          sim time. SFR is the star formation rate. 
          Time in Myr, SFR in Msun/yr
"""

import galaxy_analysis as ga
import glob
import os


import yt
import numpy as np
import deepdish as dd


_none_parameters = ['config.zone.initial_gas_mass',
                    'config.zone.initial_metallicity',
                    'config.zone.initial_stellar_mass']

def check_analysis_output(dsname, make_analysis = True):
    """
    Check if latest analysis output exists. If not, make
    it.
    """

    if os.path.isfile(dsname + '_galaxy_data.h5'):
        return True

    #
    # make the output if it doesn't exist
    # 

    if make_analysis:
        gal = ga.Galaxy(dsname)
        gal.compute_everything()
        gal.save()
        del(gal)

    return make_analysis

def get_parameters_from_analysis(galname, onez = None, gal_files = None):
    """
    Get parameters from analysis output
    """

    if onez is None:
        onez = {}

    if os.path.isfile('DD0001_galaxy_data.h5'):
        data = dd.io.load('DD0001_galaxy_data.h5')
        m_o  = data['particle_meta_data']['total_birth_mass']

        onez['config.zone.initial_stellar_mass'] = m_o

    else:
        onez['config.zone.initial_stellar_mass'] = 0.0


    if gal_files is None:
        return onez


    data = dd.io.load(galname)
    t_first_star = data['particle_meta_data']['t_first_star']

    all_t = np.zeros(len(gal_files))
    for i,f in enumerate(gal_files):
        all_t[i] = dd.io.load(f, '/meta_data/Time')

    imin = np.argmin( np.abs(all_t - t_first_star))

    onez['config.zone.initial_gas_mass'] = dd.io.load(gal_files[imin], '/meta_data/M_gas')
    onez['config.zone.initial_metallicity'] = dd.io.load(gal_files[imin], '/meta_data/Z_avg')

    onez['config.zone.t_final'] = np.max(all_t)

    return onez, t_first_star

def enzo_to_onezone(parameters, onez = None):

    if onez is None:
        onez = {}

    onez['config.units.time'] = parameters['TimeUnits']


    onez['config.zone.stochastic_sample_mass']  = parameters['IndividualStarSFGasMassThreshold']/\
                                                  parameters['IndividualStarMassFraction']

    onez['config.zone.M_min'] = parameters['IndividualStarIMFLowerMassCutoff']
    onez['config.zone.M_max'] = parameters['IndividualStarIMFUpperMassCutoff']

    onez['config.stars.SNII_mass_threshold']        = parameters['IndividualStarSNIIMassCutoff']
    onez['config.stars.SNIa_candidate_mass_bounds'] = [parameters['IndividualStarSNIaMinimumMass'],
                                                       parameters['IndividualStarWDMaximumMass']]
    onez['config.stars.DTD_slope']                  = parameters['IndividualStarDTDSlope']
    onez['config.stars.NSNIa']                      = parameters['IndividualStarSNIaFraction']

    onez['config.stars.AGB_wind_mass_threshold']    = parameters['IndividualStarAGBThreshold']
    onez['config.stars.AGB_wind_velocity']          = parameters['IndividualStarAGBWindVelocity'] 

    onez['config.stars.direct_collapse_mass_threshold'] = parameters['IndividualStarDirectCollapseThreshold']

#    onez['config.zone.initial_gas_mass'] = parameters['
#
#

    return onez


def set_default_parameters(onez = None):
    """
    Set defaults that will always be true and may not have
    Enzo analogs. This will be somewhat redundant with
    onez defaults.
    """


    if onez is None:
        onez = {}

    onez['config.zone.star_formation_method'] = 4
#    onez['config.zone.SFR_filename']          = "SFR.in"
    onez['config.zone.constant_SFR']          = 0.0

    onez['config.zone.cosmological_evolution'] = False
    onez['config.zone.use_stochastic_mass_sampling'] = True    


    # change these if getting fancy - aka read in mass loading from file as well
    onez['config.zone.inflow_factor'] = 0.0
    onez['config.zone.mass_loading_factor'] = 0.0
    onez['config.zone.mass_loading_index']  = 0.0
    onez['config.zone.mass_outflow_method'] = 2   # read from table
    onez['config.zone.SFR_efficiency']      = 0.0 # zero since using SFH
    onez['config.zone.SFR_dyn_efficiency']  = 0.0 # zero since using SFH

    onez['config.zone.t_o'] = 0.0
    onez['config.zone.dt']              = 1.0
    onez['config.zone.adaptive_timestep'] = True
    onez['config.zone.timestep_safety_factor'] = 16

    onez['config.stars.use_snII'] = True
    onez['config.stars.use_snIa'] = True
    onez['config.stars.use_stellar_winds'] = True
    onez['config.stars.use_AGB_wind_phase'] = True

    onez['config.io.dt_summary']    = 10.0
    onez['config.io.cycle_summary'] = 0
#    onez['config.io.summary_output_filename'] = 'summary_output.txt'
    onez['config.io.dt_dump']       = 100.0
    onez['config.io.radiation_binned_output'] = 1

    onez['config.zone.species_to_track'] = ['m_tot', 'm_metal', 'H', 'He', 'C', 'N', 'O', 'Mg',
                                            'Si', 'S', 'K', 'Ca', 'Ti', 'V', 'Mn', 'Fe',
                                            'Co', 'Ni', 'Cu', 'Sr', 'Y', 'Ba', 'Eu']

    #
    # set these to none to ensure they are set somehow
    #
    for k in _none_parameters:
        onez[k] = None

    return onez

def gather_mass_flow(all_files, t_o = 0.0, r = 0.25, mode = 'outflow'):
    """
    Goes through all available analysis outputs and gathers
    the mass loading factor for all species at the defined radius bin.
    """

    mass_loading = False
    if mode == 'mass_loading':
        mass_loading = True
        mode = 'outflow'

    Nfiles = len(all_files)
    t   = np.zeros(Nfiles)
    sfr = np.zeros(Nfiles)
    m   = np.zeros(Nfiles)

    # open up the first file, check which species we are dealing with
    gal    = dd.io.load(all_files[0])
    raw_fields = gal['gas_profiles'][mode]['sphere'].keys()
    fields = [x for x in raw_fields if ( not ('center' in x)) and (not ('dL' in x)) and ( not ('bin' in x))]
    fields = [x[1] for x in fields if (not ('_p0_' in x[1])) and (not ('_p1_' in x[1]))]
    ele    = [x.replace('_Mass','') for x in fields]
    ele    = [x.replace('_total','') for x in fields]
    ele    = [x.replace('cell_mass','total_mass') for x in fields]

    centers_rvir = gal['gas_profiles'][mode]['sphere']['centers_rvir']
    centers      = gal['gas_profiles'][mode]['sphere']['centers']

    flowbin = np.argmin( np.abs( r - centers_rvir) )

    pos  = centers[flowbin]

    # clean up field names a bit to just individual species, and not ionization states
    # e.g. just total, metals, H, He, C, N, O.....

    flow = {}
    mass = {}
    norm = {}
    for k in ele:
        flow[k] = np.zeros(Nfiles)
        mass[k] = np.zeros(Nfiles)
        norm[k] = np.zeros(Nfiles)

    for i, galfile in enumerate(all_files):
        t[i]   = dd.io.load(galfile, '/meta_data/Time')
        sfr[i] = dd.io.load(galfile, '/meta_data/SFR')   # 4/26 need to put this as a meta data point

        flow_data = dd.io.load(galfile, '/gas_profiles/' + mode + '/sphere')
        mass_data = dd.io.load(galfile, '/gas_profiles/accumulation/sphere')

        massbin   = np.argmin( np.abs( pos - mass_data['xbins']) )


        for j, name in enumerate(ele):
            flow[name][i] = flow_data[('gas',fields[i])][flowbin]              # outflow rate in Msun / yr at desired r
            mass[name][i] = np.sum( mass_data[('gas',fields[i])][0:massbin] )  # total mass of species interior to r

        if sfr[i] > 0:
            for name in enumerate(ele):
                norm[name][i] = flow[name][i] / mass[name][i] / sfr[i]             # normalized result - intput to onezone
        else:
            for name in enumerate(ele):
                norm[name][i] = 0.0


    f = open('./onez_model/mass_outflow.in', 'w')
    f.write("# This output contains time (Myr) and fractional mass outflow rate (f_x) through a shell centered on %.2f R_vir\n"%(r))
    f.write("# these columns are computed for species x as f_x = (dM_x/dt) / SFR / M_x, where dM_x/dt is the mass outflow\n")
    f.write("# rate in units of Msun/yr, SFR is the star formation rate in units of Msun/yr, \n")
    f.write("# and M_x is the total mass of species x in the galaxy interior to R\n")
    f.write("# To extend to other simulations / models, compute the outflow rate as f_x * SFR * M_x\n")
    f.write("#t ") # time in first column
    for e in ele:
        f.write(e + " ")
    f.write("\n")

    # fudge the last time a tiny bit to make sure no interpolation issues later on
    t[-1] = t[-1]*1.0001

    for i in np.arange(Nfiles):
        f.write("%3.3E"%(t[i] ))

        for e in ele:
            f.write(" %5.5E"%(norm[e][i]))
        f.write("\n")

    f.close()

    return

def generate_sfr(ds, t_o = 0.0):

    data = ds.all_data()

    selection = data['creation_time'].convert_to_units('Myr').value > t_o

    t, SFR = ga.particle_analysis.sfrFromParticles(ds, data = ds.all_data(), selection = selection, times = 5.0, t_o = 0.0)

    t = t / 1.0E6 # convert to Myr

    f = open('./onez_model/SFR.in', 'w')

    f.write("#t SFR\n")
    f.write("%8.8E %8.8E\n"%(t[0], SFR[0]))

    for i in np.arange(np.size(SFR)):
        f.write("%8.8E %8.8E\n"%( 0.5 * (t[i] + t[i+1]), SFR[i]))

    f.write("%8.8E %8.8E\n"%(t[-1], SFR[-1]))

    f.close()

    return

def save_script(onez, outname = 'onez_param_file.py'):

    f = open(outname,'w')

    f.write("# This python script is the parameter file and\n" +\
            "# initialization routine for a onez chemical evolution\n" +\
            "# model matching the initial star formation conditions\n" +\
            "# from an isolated dwarf galaxy feedback + chemical\n" +\
            "# evolution simulation run with Enzo. This can be used to\n" +\
            "# generate multiple (simplistic) realizations of the chemical\n" +\
            "# evolution based on the galaxy IC's and SFR from the\n" +\
            "# more expensive radiative hydrodynamics simulations.\n" +\
            "# Currently, the onez model is a closed-box model, but can\n" +\
            "# be extended to account for the mass loading factors produced\n" +\
            "# in the simulation (future work).\n\n" +\
            "# The onezone model can be run simply as:\n" +\
            "#      python onez_parameter_file.py\n")

    f.write("from onezone import zone\n")
    f.write("from onezone import imf\n")
    f.write("from onezone import config\n\n\n\n")

    f.write("# -------------------- Set parameters ---------------------\n")
    # organize the output parameters slightly
    for x in ['config.units','config.zone','config.stars','config.io']:

        for k in onez.keys():
            if x in k:
                f.write(k + "          = " + str(onez[k]) + "\n")

        f.write("\n\n")

    f.write("# -------------------- Set up and run the simulation -------------------\n")
    f.write("istart = len(glob.glob('run????_summary_output.txt'))\n")
    f.write("for i in np.arange(100):\n")
    f.write("    np.random.seed(i)\n")
    f.write("    config.io.summary_output_filename = 'run%0004i_summary_output.txt'%(istart + i)\n")
    f.write("    config.io.dump_output_basename    = 'run%0004i_dump'%(istart + i)\n")
    f.write("    sim = zone.Zone()\n")
    f.write("    sim.set_initial_abundances(config.zone.species_to_track)\n")
    f.write("    try:\n")
    f.write("        sim.evolve()\n")
    f.write("    except:\n")
    f.write("        print 'failing in run ' + str(i) + ' skipping and going to the next one'\n")
    f.write("    del(sim)\n")

    f.close()

    return


if __name__=="__main__":

    # glob for both data files AND galaxy hdf5 files

    make_analysis = False # change to command line parameter

    ds_files = glob.glob('DD????')
    ds_files = np.sort(ds_files)

    gal_files = glob.glob('DD????_galaxy_data.h5')
    gal_files = np.sort(gal_files)

    if len(ds_files) > 1:
        dsname = ds_files[-1]

        if check_analysis_output(ds_files[-1], make_analysis = make_analysis):
            analysis_to_use = ds_files[-1] + '_galaxy_data.h5'
        else:
            if len(gal_files) < 1:
                print "Need to allow for generating the analysis file, or have files already generated"
                raise ValueError

            # if we aren't generating file, use the latest one made
            print "Warning: Not using most up to date output dump"
            analysis_to_use = gal_files[-1]

    else:
        print "Currently, requires data dump to function - None found"
        raise ValueError 

        if len(gal_files) < 1:
            print "No data dumps or analysis files found - what are you doing?"
            raise ValueError

        analysis_to_use = gal_files[-1]

    ds = yt.load(dsname + '/' + dsname)

    onez = set_default_parameters()
    onez = enzo_to_onezone(ds.parameters, onez = onez)
    onez, t_o = get_parameters_from_analysis(analysis_to_use, onez = onez, gal_files = gal_files)

    generate_sfr(ds, t_o = t_o)
    gather_mass_flow(gal_files, r = 0.25)

    save_script(onez, outname = "./onez_model/onez_parameter_file.py")
