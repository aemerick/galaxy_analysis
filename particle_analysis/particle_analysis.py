import yt
import numpy as np




def compute_binned_sfr(ds, data,
                       particle_types = ['all_stars','all_popIII_stars','all_popII_stars'],
                       tmin = 0.0, tmax = None, dt = 10.0):
    """
    Compute the SFR for each particle type and returns
    a dictionary containing this information.

    This assumes that stellar particle types have been defined
    using galaxy_analysis.field_generators.define_particle_types

    Parameters
    -----------
    ds   :   yt ds object
    data :   yt data region object to compute on
    particle_types : yt-defined / registered particle types to compute SFR for . Optional
    tmin : start time for bins in Myr . Optional, 0.0
    tmax : end time for bins in Myr . Optional, None. Taken to be curent time rouded up to be even in dt
    dt   : bin size in Myr. Optional. 10.0
    """

    if tmax is None:
        # ceil to nearest dt
        tmax = np.ceil(ds.data.current_time.to('Myr').value / dt) * dt

    tbins = np.arange(tmin, tmax + dt*0.5, dt) * yt.units.Myr

    # particle_types =

    sfr_data = {}

    sfr_data['bins']  = tbins
    sfr_data['cbins'] = 0.5 * (tbins[1:] + tbins[:-1])

    for pt in particle_types:

        sfr_data[pt + '_SFR'] = np.zeros(np.size(tbins)-1)

        t0   = data[(pt,'creation_time')].to('Myr')

        if np.size(t0) == 0:
            continue

        bm = data[(pt,'birth_mass')].value

        sfr_data[pt + '_SFR'], edges = np.hist(t0, bins=tbins, weights = bm)

    return sfr_data



def compute_end_of_life(ds, data,
                        particle_types = None,
                        dt = 10):
    """
    Computes the end of life behavior of each particle type and returns
    a dictionary containing this information. Options are:

        SNR      : Pop II core collapse SNe rate
        SNR_III  : Pop III core collapse SNe rate
        PISNe    : Pop III PISNe rate
        AGB      : Pop II AGB star rate


    This assumes that stellar particle types have been defined
    using field_generators.define_particle_types

    Parameters
    -----------
    ds   :   yt ds object
    data :   yt data region object to compute on
    particle_types : yt-defined / registered particle types to compute SFR for . Optional
    tmin : start time for bins in Myr . Optional, 0.0
    tmax : end time for bins in Myr . Optional, None. Taken to be curent time rouded up to be even in dt
    dt   : bin size in Myr. Optional. 10.0
    """

    if tmax is None:
        # ceil to nearest dt
        tmax = np.ceil(ds.data.current_time.to('Myr').value / dt) * dt

    tbins = np.arange(tmin, tmax + dt*0.5, dt) * yt.units.Myr

    # particle_types =

    eol_data = {}

    eol_data['bins']  = tbins
    eol_data['cbins'] = 0.5 * (tbins[1:] + tbins[:-1])

    # translate particle type labels to something sensible for the purpose
    # of this code
    labels = {'ccsne_remnant' : 'SNR_II', 'popIII_ccsne_remnant' : 'SNR_III',
              'snia_progenitor' : 'SNR_Ia', 'snia_dds_progenitor' : 'SNR_Ia_DDS',
              'snia_hers_progenitor' : 'SNR_Ia_HERS',
              'snia_sch_progenitor' : 'SNR_Ia_SCH', 'snia_sds_progenitor' : 'SNR_Ia_SDS',
              'popIII_pisne_remnant' : 'SNR_PISNe',
              "white_dwarf" : "AGB_rate",
              "direct_collapse_remnant" : "direct_collapse_rate",
              'popIII_direct_collapse_remnant' : "popIII_direct_collapse_rate"}

    if particle_type is None:
        particle_type = list(labes.keys())

    ds.add_particle_filter("snia_dds_progenitor")
    ds.add_particle_filter("snia_hers_progenitor")
    ds.add_particle_filter("snia_sds_progenitor")
    ds.add_particle_filter("snia_sch_progenitor")

    for pt in particle_types:
        eol_data[pt + '_'] = np.zeros(np.size(tbins)-1)

        t0   = data[(pt,'creation_time')].to('Myr')


    return eol_data
