import yt.mods as yt
import numpy as np
from collections import Iterable, OrderedDict

import glob

# -- internal --
from galaxy_analysis.utilities import convert_abundances

def compute_aratio(ds, data, ratios, particle_type = 11):
    """
    Generate abundances for all particles of given particle type.
    If particle_type is None, does this for all particles
    """

    birth_mass = data['birth_mass'].value * yt.units.Msun
    ptype      = data['particle_type']

    if not isinstance(ratios, Iterable):
        ratios = [ratios]
    #
    # split string
    #
    aratios = OrderedDict()

    for ratio in ratios:
        if '/' in ratio:
            ele1, ele2 = ratio.rsplit('/')
        else:
            print "Must provide abundance ratio string of style: Fe/H"
            return

        enzo_name1 = ('io','particle_' + ele1 + '_fraction')
        enzo_name2 = ('io','particle_' + ele2 + '_fraction')

        mass1 = data[enzo_name1] * birth_mass
        mass2 = data[enzo_name2] * birth_mass

        aratios[ratio] = convert_abundances.abundance_ratio( (ele1, mass1), 
                                                             (ele2, mass2), 
                                                             'mass' )
    return aratios


if __name__=='__main__':
    
    #
    # do this for all
    #
    outdir  = './abundances/'
    outfile = 'abundances.h5'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if not os.path.isfile(outdir + outfile):
        hf = h5py.File(outdir + outfile, 'w')
    else:
        hf = h5py.File(outdir + outfile, 'r')

    ds_list = np.sort( glob.glob('./DD????/DD????') )
    times   = np.zeros(np.size(ds_list))
 
    # get elements present:
    ds              = yt.load(dsname[-1])
    fields          = ds.field_list
    element_fields  = [x[1] for x in fields if ('particle_' in x[1] and\
                                                '_fraction' in x[1]) and\
                                               ('all' in x[0])]
    elements = np.sort([x.rsplit('_')[1] for x in element_fields])
    metals   = [x for x in elements if x != 'H']
    ratios   = [ x +'/H' for x in metals]

    for i, dsname in enumerate(ds_list):
        ds   = yt.load(dsname)
        data = ds.all_data()

        groupname = dsname.rsplot('/')[0]
        g = hf.create_group(groupname)
        g.create_dataset('Time'  , data = ds.current_time.convert_to_units('Myr').value)
    
        if ('io', 'particle_type') in ds.field_list:
            
            aratios = compute_aratio(ds, data, ratios)

            g.create_dataset('Nstars', data = np.size(data['particle_mass'][ data['particle_type'] == 11]))
            g.create_dataset('Mstars', data = np.sum( data['particle_mass'][ data['particle_type'] == 11].convert_to_units('Msun').value))

            sg = hf.create_group(groupname + '/abundances')
            for abundance in aratios.keys():
                sg.create_dataset( abundance, data = aratios[abundance])
       else:
           g.create_dataset('Nstars', data = 0.0)
           g.create_dataset('Mstars', data = 0.0)
           sg = hf.create_group(groupname + '/abundances')
           for abundance in aratios.keys():
               sg.create_dataset( abundance, data = 0.0)

    hf.close()
