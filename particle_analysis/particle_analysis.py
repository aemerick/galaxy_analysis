import yt
import numpy as np




def compute_binned_sfr(data,
                       particle_types = ['all_stars','all_popIII_stars','all_popII_stars'],
                       tmin = 0.0, tmax = None, dt = 10.0):
    """
    Compute the SFR for each particle type and returns
    a dictionary containing this information.

    This assumes that stellar particle types have been defined
    using galaxy_analysis.field_generators.define_particle_types

    SFR rate is returned in Msun / Myr

    Parameters
    -----------
    data :   yt data region object to compute on
    particle_types : yt-defined / registered particle types to compute SFR for . Optional
    tmin : start time for bins in Myr . Optional, 0.0
    tmax : end time for bins in Myr . Optional, None. Taken to be curent time rouded up to be even in dt
    dt   : bin size in Myr. Optional. 10.0
    """

    if tmax is None:
        # ceil to nearest dt
        tmax = np.ceil(data.ds.current_time.to('Myr').value / dt) * dt

    tbins = np.arange(tmin, tmax + dt*0.5, dt) * yt.units.Myr

    # particle_types =

    sfr_data = {}

    sfr_data['bins']  = tbins
    sfr_data['cbins'] = 0.5 * (tbins[1:] + tbins[:-1])

    for pt in particle_types:

        sfr_data[pt + '_SFR'] = np.zeros(np.size(tbins)-1)

        t0   = data[(pt,'creation_time')].to('Myr').value

        if np.size(t0) == 0:
            continue

        bm = data[(pt,'birth_mass')].value

        print(np.shape(bm), np.shape(t0))
        print(pt)
        print(np.shape(tbins))

        sfr_data[pt + '_SFR'] = np.histogram(t0, bins=tbins, weights = bm)[0] / (tbins[1:]-tbins[:-1])

    return sfr_data



def compute_end_of_life(data,
                        particle_types = None,
                        tmin = 0.0, tmax = None, dt = 10):
    """
    Computes the end of life behavior of each particle type and returns
    a dictionary containing this information. Options are:

    By default this computes for all of the particle types below (left) and
    outputs as the names on right

          'ccsne_remnant' :                  'SNR_II'
          'popIII_ccsne_remnant' :           'SNR_popiii'
          'snia_progenitor' :                'SNR_Ia'
          'snia_dds_progenitor' :            'SNR_Ia_DDS'
          'snia_hers_progenitor' :           'SNR_Ia_HERS'
          'snia_sch_progenitor' :            'SNR_Ia_SCH'
          'snia_sds_progenitor' :            'SNR_Ia_SDS'
          'popIII_pisne_remnant' :           'SNR_PISNe'
          "white_dwarf" :                    "AGB_rate"
          "direct_collapse_remnant" :        "direct_collapse_rate"
          'popIII_direct_collapse_remnant' : "popIII_direct_collapse_rate



    This assumes that stellar particle types have been defined
    using field_generators.define_particle_types.

    Rate is returned in events / Myr

    Parameters
    -----------
    data :   yt data region object to compute on
    particle_types : yt-defined / registered particle types to compute SFR for . Optional
    tmin : start time for bins in Myr . Optional, 0.0
    tmax : end time for bins in Myr . Optional, None. Taken to be curent time rouded up to be even in dt
    dt   : bin size in Myr. Optional. 10.0
    """

    if tmax is None:
        # ceil to nearest dt
        tmax = np.ceil(data.ds.current_time.to('Myr').value / dt) * dt

    tbins = np.arange(tmin, tmax + dt*0.5, dt) * yt.units.Myr

    # particle_types =

    eol_data = {}

    eol_data['bins']  = tbins
    eol_data['cbins'] = 0.5 * (tbins[1:] + tbins[:-1])

    # translate particle type labels to something sensible for the purpose
    # of this code
    labels = {'ccsne_remnant' :                  'SNR_II',
              'popIII_ccsne_remnant' :           'SNR_popiii',
              'snia_progenitor' :                'SNR_Ia',
              'snia_dds_progenitor' :            'SNR_Ia_DDS',
              'snia_hers_progenitor' :           'SNR_Ia_HERS',
              'snia_sch_progenitor' :            'SNR_Ia_SCH',
              'snia_sds_progenitor' :            'SNR_Ia_SDS',
              'popIII_pisne_remnant' :           'SNR_PISNe',
              "white_dwarf" :                    "AGB_rate",
              "direct_collapse_remnant" :        "direct_collapse_rate",
              'popIII_direct_collapse_remnant' : "popIII_direct_collapse_rate"}

    if particle_types is None:
        particle_types = list(labels.keys())

    for pt in particle_types:
        name = labels[pt]

        eol_data[name] = np.zeros(np.size(tbins)-1)

        t0         = data[(pt,'creation_time')].to('Myr')
        lifetime   = data[(pt,'particle_model_lifetime')].to('Myr')

        if 'snia' in pt:
            # need to ADD the actual lifetime of the particle
            # which represents the time AFTER WD formation for explosion
            lifetime = lifetime + data[(pt,'dynamical_time')].to('Myr')

        death_time = (t0 + lifetime).value

        eol_data[name] = np.histogram(death_time, bins = tbins)[0] / (tbins[1:]-tbins[:-1])


    return eol_data
