def select_formed_stars(ds, data, t_min, t_max):
    """
    Select stars formed in a given time range - useful for determining
    'instantaneous' properties of formed stars...
    """
    lifetime = data[('io','particle_model_lifetime')].convert_to_units('Myr')
    t_o      = data[('io','creation_time')].convert_to_units('Myr')

    pcut     = (t_o >= t_min) * (t_o <= t_max)

    return pcut

def select_alive_stars(ds, data, t_min, t_max,
                       include_partial = True):
    """
    Given a list of all stars, returns boolean list that
    can be used to select the ones that were 'alive'
    during the desired time range, t_min and t_max. This is very
    useful for post-processing without having to acess many datasets

    include_partial : include stars that were alive for only part of the time
                      range, rather than the entire time range. Defualt True
    """

    lifetime = data[('io','particle_model_lifetime')].convert_to_units('Myr')
    t_o      = data[('io','creation_time')].convert_to_units('Myr')

    if include_partial: # include stars that were alive at some point during range
        pcut = (t_o + lifetime >= t_min) * ( t_o <= t_max )
    else: # only stars that were alive during the entire time range
        pcut = (t_o <= t_min) * (t_o + lifetime >= t_max)

    return pcut



def white_dwarfs(ds, data):
    """
    Returns slicing array to determine white dwarf particles
    """
    pcut = data['particle_type'] == 12
    pcut = pcut * (data['particle_mass'] > 0.0)

    return pcut

def snIa(ds, data):

    pcut = data['particle_type'] == 12
    pcut = pcut * (data['particle_mass'] == 0.0)

    return pcut

def main_sequence(ds, data):
    return data['particle_type'] == 11

def core_collapse(ds, data):

    pcut = data['particle_type'] == 13

    pcut = pcut * (data['birth_mass'] < ds.parameters['IndividualStarDirectCollapseThreshold'])*\
                  (data['birth_mass'] > ds.parameters['IndividualStarAGBThreshold'])

    return pcut

def direct_collapse(ds,data):

    pcut = data['particle_type'] == 13
    pcut = pcut * (data['birth_mass'] > ds.parameters['IndividualStarDirectCollapseThreshold'])

    return pcut

