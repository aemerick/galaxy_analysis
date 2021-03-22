import yt
import numpy as np


from galaxy_analysis.utilities import utilities





def _define_particle_filter_functions():
    """
    Hidden function to define filtering functions for particles to be used
    with yt data. This is a way to do yt's particle filtering (sort of)
    a bit more manually which is nice because it DOES NOT recompute data
    when filtering.
    """

    def all_stars(dobj, filtered_type = 'all'):
        filter = dobj[(filtered_type, "particle_type")] >= 11
        return filter

    def all_popIII_stars(dobj, filtered_type = 'all'):
        filter = dobj[(filtered_type,"particle_is_popiii")].astype(np.bool)
        return filter

    def all_popII_stars(dobj, filtered_type = 'all'):
        filter = np.logical_not(dobj[(filtered_type,"particle_is_popiii")])

        return filter

    def main_sequence_stars(dobj, filtered_type = 'all'):
        filter = (dobj[(filtered_type, "particle_type")] == 11) +\
                 (dobj[(filtered_type, "particle_type")] == 15)

        return filter

    def main_sequence_popIII_stars(dobj, filtered_type = 'all'):
        filter = dobj[(filtered_type, "particle_type")] == 14
        return filter

    def remnant_stars(dobj, filtered_type = 'all'):
        filter = dobj[(filtered_type, "particle_type")] == 13
        return filter

    def low_mass_stars(dobj, filtered_type = 'all'):
        filter = dobj[(filtered_type, "particle_type")] == 11
        filter = filter * (dobj[(filtered_type,"birth_mass")] > 2.0) * (dobj[(filtered_type,"birth_mass")] < 8.0)
        return filter

    def low_mass_unresolved_stars(dobj, filtered_type = 'all'):
        filter = dobj[(filtered_type, "particle_type")] == 15
        return filter

    def white_dwarf(dobj, filtered_type = 'all'):
        filter = dobj[(filtered_type, "particle_type")] == 12
        return filter


    #
    #
    # End of life filteres for non-snia
    #
    #

    def all_remnants(dobj, filtered_type = 'all'):
        filter = dobj[(filtered_type,"particle_type")] == 13
        return filter

    def popIII_remnant(dobj, filtered_type = 'all'):

        if ('IndividualStarPopIIIFormation' in data.ds.parameters) and\
           ('PopIIIMetalCriticalFraction' in dobj.ds.parameters):
            if dobj.ds.parameters['IndividualStarPopIIIFormation'] > 0:

                if dobj.ds.parameters['PopIIIMetalCriticalFraction'] > 0:
                    filter = dobj[(filtered_type,'metallicity_fraction')] < dobj.ds.parameters['PopIIIMetalCriticalFraction']
                else: # use the Chiaki threshold:

                    filter = np.logical_not( dobj[(filtered_type,'particle_above_chiaki_threshold')] )

        else:
            filter = np.logical_not(dobj[(filtered_type, "birth_mass")] == dobj[(filtered_type, "birth_mass")])


        return filter

    def popIII_ccsne_remnant(dobj, filtered_type = 'all'):

        if ('IndividualStarPopIIIFormation' in dobj.ds.parameters) and\
           ('PopIIIMetalCriticalFraction' in dobj.ds.parameters):
            if dobj.ds.parameters['IndividualStarPopIIIFormation'] > 0:

                if dobj.ds.parameters['PopIIIMetalCriticalFraction'] > 0:
                    filter = dobj[(filtered_type,'metallicity_fraction')] < dobj.ds.parameters['PopIIIMetalCriticalFraction']
                else: # use the Chiaki threshold:

                    filter = np.logical_not( dobj[(filtered_type,'particle_above_chiaki_threshold')] )

                filter = filter * ((dobj[(filtered_type,'birth_mass')] >= dobj.ds.parameters['TypeIILowerMass']) *\
                               (dobj[(filtered_type,'birth_mass')] <= dobj.ds.parameters['TypeIIUpperMass']))
        else:
            filter = np.logical_not(dobj[(filtered_type, "birth_mass")] == dobj[(filtered_type, "birth_mass")])


        return filter


    def popIII_pisne_remnant(dobj, filtered_type = 'all'):

        if ('IndividualStarPopIIIFormation' in dobj.ds.parameters) and\
           ('PopIIIMetalCriticalFraction' in dobj.ds.parameters):
            if dobj.ds.parameters['IndividualStarPopIIIFormation'] > 0:

                if dobj.ds.parameters['PopIIIMetalCriticalFraction'] > 0:
                    filter = dobj[(filtered_type,'metallicity_fraction')] < dobj.ds.parameters['PopIIIMetalCriticalFraction']
                else: # use the Chiaki threshold:

                    filter = np.logical_not( dobj[(filtered_type,'particle_above_chiaki_threshold')] )

                filter = filter * ((dobj[(filtered_type,'birth_mass')] >= dobj.ds.parameters['PISNLowerMass']) *\
                               (dobj[(filtered_type,'birth_mass')] <= dobj.ds.parameters['PISNUpperMass']))
        else:
            filter = np.logical_not(dobj[(filtered_type, "birth_mass")] == dobj[(filtered_type, "birth_mass")])


        return filter

    def popIII_direct_collapse_remnant(dobj, filtered_type = 'all'):

        if ('IndividualStarPopIIIFormation' in dobj.ds.parameters) and\
           ('PopIIIMetalCriticalFraction' in dobj.ds.parameters):
            if dobj.ds.parameters['IndividualStarPopIIIFormation'] > 0:

                if dobj.ds.parameters['PopIIIMetalCriticalFraction'] > 0:
                    filter = dobj[(filtered_type,'metallicity_fraction')] < dobj.ds.parameters['PopIIIMetalCriticalFraction']
                else: # use the Chiaki threshold:

                    filter = np.logical_not( dobj[(filtered_type,'particle_above_chiaki_threshold')] )

                filter = filter *  (np.logical_not((dobj[(filtered_type,'birth_mass')] >= dobj.ds.parameters['PISNLowerMass']) *\
                                                  (dobj[(filtered_type,'birth_mass')] <= dobj.ds.parameters['PISNUpperMass'])) *\
                                    np.logical_not((dobj[(filtered_type,'birth_mass')] >= dobj.ds.parameters['TypeIILowerMass']) *\
                                                  (dobj[(filtered_type,'birth_mass')] <= dobj.ds.parameters['TypeIIUpperMass'])))
        else:
            filter = np.logical_not(dobj[(filtered_type, "birth_mass")] == dobj[(filtered_type, "birth_mass")])


        return filter

    def ccsne_remnant(dobj, filtered_type = 'all'):

        filter = ((dobj[(filtered_type, "birth_mass")] <= dobj.ds.parameters['IndividualStarDirectCollapseThreshold']) *\
                  (dobj[(filtered_type, "birth_mass")] >= dobj.ds.parameters['IndividualStarAGBThreshold']))

        if ('IndividualStarPopIIIFormation' in dobj.ds.parameters) and\
           ('PopIIIMetalCriticalFraction' in dobj.ds.parameters):
            if dobj.ds.parameters['IndividualStarPopIIIFormation'] > 0:

                if dobj.ds.parameters['PopIIIMetalCriticalFraction'] > 0:
                    filter = filter * dobj[(filtered_type,'metallicity_fraction')] < dobj.ds.parameters['PopIIIMetalCriticalFraction']
                else: # use the Chiaki threshold:
                    filter = filter * dobj[(filtered_type,'particle_above_chiaki_threshold')].astype(np.bool)

        return filter

    def direct_collapse_remnant(dobj, filtered_type = 'all'):

        filter = dobj[(filtered_type, "birth_mass")] > dobj.ds.parameters['IndividualStarDirectCollapseThreshold']

        if ('IndividualStarPopIIIFormation' in dobj.ds.parameters) and\
           ('PopIIIMetalCriticalFraction' in dobj.ds.parameters):
            if dobj.ds.parameters['IndividualStarPopIIIFormation'] > 0:

                if dobj.ds.parameters['PopIIIMetalCriticalFraction'] > 0:
                    filter = filter * dobj[(filtered_type,'metallicity_fraction')] < dobj.ds.parameters['PopIIIMetalCriticalFraction']
                else: # use the Chiaki threshold:
                    filter = filter * dobj[(filtered_type,'particle_above_chiaki_threshold')].astype(np.bool)

        return filter

    def snia_progenitor(dobj, filtered_type = 'all'):
        filter = ((dobj[(filtered_type, "birth_mass")] >= dobj.ds.parameters['IndividualStarSNIaMinimumMass']) *\
                  (dobj[(filtered_type, "birth_mass")] <= dobj.ds.parameters['IndividualStarSNIaMaximumMass']))

        filter = filter * ( (dobj[(filtered_type, 'snia_sch_metal_fraction')] < 0) +\
                            (dobj[(filtered_type, 'snia_sds_metal_fraction')] < 0) +\
                            (dobj[(filtered_type, 'snia_hers_metal_fraction')] < 0) +\
                            (dobj[(filtered_type, 'snia_metal_fraction')] < 0) )
        return filter

    def snia_dds_progenitor(dobj, filtered_type = 'all'):
        filter = ((dobj[(filtered_type, "birth_mass")] >= dobj.ds.parameters['IndividualStarSNIaMinimumMass']) *\
                   (dobj[(filtered_type, "birth_mass")] <= dobj.ds.parameters['IndividualStarSNIaMaximumMass']))

        filter = filter * (dobj[(filtered_type, 'snia_metal_fraction')] < 0)

        return filter

    def snia_sch_progenitor(dobj, filtered_type = 'all'):
        filter = ((dobj[(filtered_type, "birth_mass")] >= dobj.ds.parameters['IndividualStarSNIaMinimumMass']) *\
                  (dobj[(filtered_type, "birth_mass")] <= dobj.ds.parameters['IndividualStarSNIaMaximumMass']))

        filter = filter * (dobj[(filtered_type, 'snia_sch_metal_fraction')] < 0)

        return filter

    def snia_hers_progenitor(dobj, filtered_type = 'all'):
        filter = ((dobj[(filtered_type, "birth_mass")] >= dobj.ds.parameters['IndividualStarSNIaMinimumMass']) *\
                  (dobj[(filtered_type, "birth_mass")] <= dobj.ds.parameters['IndividualStarSNIaMaximumMass']))

        filter = filter * (dobj[(filtered_type, 'snia_hers_metal_fraction')] < 0)

        return filter

    def snia_sds_progenitor(dobj, filtered_type = 'all'):
        filter = ((dobj[(filtered_type, "birth_mass")] >= dobj.ds.parameters['IndividualStarSNIaMinimumMass']) *\
                  (dobj[(filtered_type, "birth_mass")] <= dobj.ds.parameters['IndividualStarSNIaMaximumMass']))

        filter = filter * (dobj[(filtered_type, 'snia_sds_metal_fraction')] < 0)

        return filter


    function_dict = { 'all_stars' : all_stars,
                      'all_popIII_stars' : all_popIII_stars,
                      'all_popII_stars'  : all_popII_stars,
                      'main_sequence_stars' : main_sequence_stars,
                      'main_sequence_popIII_stars' : main_sequence_popIII_stars,
                      'remnant_stars' : remnant_stars,
                      'low_mass_stars' : low_mass_stars,
                      'low_mass_unresolved_stars': low_mass_unresolved_stars,
                      'white_dwarf' : white_dwarf,
                      # ----------
                      'all_remnants': all_remnants,
                      'popIII_remnant': popIII_remnant,
                      'popIII_ccsne_remnant' : popIII_ccsne_remnant,
                      'popIII_pisne_remnant' : popIII_pisne_remnant,
                      'popIII_direct_collapse_remnant' : popIII_direct_collapse_remnant,
                      'ccsne_remnant' : ccsne_remnant,
                      'direct_collapse_remnant' : direct_collapse_remnant,
                      'snia_progenitor' : snia_progenitor,
                      'snia_dds_progenitor' : snia_dds_progenitor,
                      'snia_sch_progenitor' : snia_sch_progenitor,
                      'snia_hers_progenitor' : snia_hers_progenitor,
                      'snia_sds_progenitor' : snia_sds_progenitor}

    return function_dict

_particle_filter_function_dict = _define_particle_filter_functions()



def particle_filter(particle_type,
                    data = None,
                    inherited_type = 'all',
                    join_method = np.logical_or,
                    return_function = False):
    """
    Do particle filtering using defined particle filters. This is not the ideal
    way to be doing this, but is faster than routing through yt's particle
    filter interface (at least for now).

    Makes no safety checks for fields being defined or not.

    Parameters
    -----------
    particle_type :  name of particle type. If "help" is supplied, returns
                     a list of particle types this function can filter.

                     Can pass a list and return the join_method of them

    data          : yt data object. Almost always should be passed, but this
                    is the data object to do the filtering on.
                    Optional (default None)

    inherited_type : Type of particle to provide a filter for, where this must
                     be a known particle type in yt. Not quite sure why this
                     would be useful, but godo to have. Default : 'all'

    return_function : Bool. If True, will return the function that does the
                    filtering on the data object instead of the filter itself.
                    data must be None for this to work.
                    Default : False

    Returns
    --------
    filter   : array that filters the dataset for the desired particles
    """

    if particle_type == 'help':
        return list(_particle_filter_function_dict.keys())

    if return_function:
        return lambda data_x : _particle_filter_function_dict[particle_type](data_x, filtered_type=inherited_type)
    else:

        if isinstance(particle_type,list):
            if join_method == np.logical_or :
                join_method = np.sum
            elif join_method == np.logical_and:
                join_method = np.multiply
            else:
                print("join method not supported")
                raise ValueError

            filter = join_method([particle_filter_function_dict[pt](data,filtered_type=inherited_type) for pt in particle_types], axis=0).astype(np.bool)
        else:

            filter = _particle_filter_function_dict[particle_type](data, filtered_type = inherited_type)

        return filter


def compute_stellar_MDFs(data,
                         particle_types = ['main_sequence_stars','main_sequence_popIII_stars'],
                         species = 'all',
                         bins=None, amin = -20, amax = -1,
                         talive = None):
    """
    Compute the stellar MDFs for given particle types and elements
    in a given yt data object.

    MDFs are all computed as [X/H], which can be used to derive aother ratios,
    (e.g. [X/Y] = [X/H] - [Y/H]). Solar normalization follows abundances from
    Asplund+2009.

    In addition, we compute the mass fraction MDFs for each fraction type
    followed (for convenience this is only done if elements = 'all' at the
    moment).
    """

    popiii_metals = None


    #
    # organized
    #
    MDFs = {}


    if species == 'all':

        metals = utilities.species_from_fields(data.ds.field_list)
        metal_fields = ['particle_' + x + '_fraction' for x in metals]

        popiii_metal_fields = [x[1] for x in data.ds.field_list if ('popIII_particle_' in x[1]) and ('_fraction' in x[1])]
        source_fractions = ['intermediate_wind','agb','massive_wind','popIII','rprocess','snii','snia',
                             'snia_sch','snia_hers','snia_sds']

        source_fraction_fields = []
        for name in source_fractions:
            if ('io',name) in data.ds.field_list:
                source_fraction_fields.append(name + '_fraction')

        #
        # now gather these together
        #
        bracket_species  = metals + popiii_metals
        fraction_species = ['metallicity_fraction'] + source_fraction_fields

        print("Not complete 'all' does not yet work")
        raise RuntimeError

    elif species == 'metals_only':

        metals = utilities.species_from_fields(data.ds.field_list)

        # we want to make [X/H] for all of these
        for pt in particle_types:
            MDF[pt] = {}
            bm = data[(pt,'birth_mass')].value
            for e in metals:
                abund = data[(pt,'particle_' + e + '_over_H')].value

                MDF[pt][e + '_over_H'] = np.histogram(abund, bins=bins, weights=weights)[0] / db






        #for i in np.arange(np.size(elements)):
        #    if ('all','particle_' + e + '_popiii_fraction') in data.ds.field_list:
        #        elements =


    return

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

        filter = particle_filter(pt,data)

        if not np.any(filter):
            continue

        t0 = data[('all','creation_time')][filter].value
        bm = data[('all','birth_mass')][filter].value

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

        filter = particle_filter(pt,data)

        if not np.any(filter):
            continue

        t0         = data[('all','creation_time')][filter].to('Myr')
        lifetime   = data[('all','particle_model_lifetime')][filter].to('Myr')

        if 'snia' in pt:
            # need to ADD the actual lifetime of the particle
            # which represents the time AFTER WD formation for explosion
            lifetime = lifetime + data[('all','dynamical_time')][filter].to('Myr')

        death_time = (t0 + lifetime).value

        eol_data[name] = np.histogram(death_time, bins = tbins)[0] / (tbins[1:]-tbins[:-1])


    return eol_data
