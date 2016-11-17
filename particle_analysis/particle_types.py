def white_dwarfs(ds, data):
    """
    Returns slicing array to determine white dwarf particles
    """
    pcut = data['particle_type'] == 12
    pcut = pcut * data['particle_mass'] > 0.0

    return pcut

def snIa(ds, data):

    pcut = data['particle_type'] == 12
    pcut = pcut * data['particle_mass'] == 0.0

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

