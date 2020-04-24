"""
   Convert the yields generated from Benoit (in dictionary form) into
   an HDF5 table of the style needed for Enzo simulations.

   These tables have the following groups:

        SN     : CCSNe yields (stars above 8 Msun)
        Wind   : Wind yields (at the moment, must have same grid points as above)
        AGB    : AGB yields (stars below 8 Msun)
        PopIII : PopIII yields

    As noted, SN and wind table must have the same dimensions at the moment,
    even though this is not best to rely on. But this makes the code in Enzo
    easier for the moment. The intent is that SN and wind yields come from the
    same model, but are just split in type.
"""

import numpy as np
import h5py

from galaxy_analysis import static_data as const  # for elements and things

available_yields = {'FRUITY' : 'fruity_elem_stable.npy',
                    'LC18_R_0' : 'LC18_R_0_elem_stable.npy',
                    'LC18_R_150': 'LC18_R_150_elem_stable.npy',
                    'LC18_R_300': 'LC18_R_300_elem_stable.npy',
                    'LC18_winds_R_0' : 'LC18_winds_R_0_elem_stable.npy',
                    'LC18_winds_R_150': 'LC18_winds_R_150_elem_stable.npy',
                    'LC18_winds_R_300': 'LC18_winds_R_300_elem_stable.npy',
                    'Monash' : 'monash_elem_stable.npy',
                    'Nomoto' : 'nomoto_elem_stable.npy',
                    'Nomoto-HNe' : 'nomoto_HNe_elem_stable.npy',
                    'PopIII' : 'popIII_yields.in'}

info = {}
# temporary for now. better description later
for k in available_yields.keys()
    info[k] = available_yields[k]


def massage_dataset(dataset):
    """
    Take kwargs and things from Benoits dict to what I need
    """

    outdict = {}

    outdict['M'] = dataset['M_list']
    outdict['Z'] = dataset['Z_list']

    outdict['num_M'] = np.size(outdict['M'])
    outdict['num_Z'] = np.size(outdict['Z'])

    outdict['colnames'] = np.array(['Mtot','Metal_Mtot'] + dataset['elem_list'])
    outdict['anums']    = np.array([-1,0] + [const.asym_to_anum[x] for x in dataset['elem_list']])
    outdict['num_Y']    = np.size(outdict['anums'])

    # probably a few different ways to do this, but most fool-proof is to loop:
    yields = np.zeros( (outdict['num_M'], outdict['num_Z'], outdict['num_Y']))
    for i,Mkey in enumerate(list(dataset['yields'][ list(dataset.keys())[0] ].keys())):
        for j,zkey in enumerate(list(dataset['yields'].keys())):
            for k,e in enumerate(dataset['elem_list']):
                yields[i][j][k+2] = dataset['yields'][Zkey][Mkey][e]

                yields[i][j][0] += yields[i][j][k+2]
                if (not (e=='H')) and (not (e=='He')):
                    yields[i][j][1] += yields[i][j][k+2]

    #
    outdict['yields'] = yields

    return outdict

def dict_to_dataset(grp, outdict, info_str = ''):

    grp.create_dataset("M", (outdict['num_M'],), dtype='f', data=outdict['M'])
    grp.create_dataset("Z", (outdict['num_Z'],), dtype='f', data=outdict['Z'])

    grp.create_dataset("yield_names", (outdict['ncol'],), dtype='S10',
                            data = outdict['colnames'].astype('S10'))
    grp.create_dataset("atomic_numbers", (outdict['ncol'],), dtype='i',
                            data = outdict['anums'])

    grp.create_dataset("yields", (outdict['num_M'],outdict['num_Z'],outdict['num_Y'],),
                           dtype='f', data=outdict['yields'])

    grp.attrs['Info'] = info_str

    return

def construct_table(SN_model = 'LC18_R_0', wind_model = 'LC18_winds_R_0',
                    AGB_model = 'FRUITY', PopIII_model = 'PopIII',
                    wdir = './', outname = None):
    """
    Construct table
    """

    if outname is None:
        outname = SN_model + '-' + wind_model + '-' + AGB_model + '-' + PopIII_model + '.h5'

    hf = h5py.File(wdir + outname, 'w')

    #
    # Load each set one-by-one
    #
    hf.create_group('SN')
    SN_yields = np.load( wdir + available_yields[SN_model], allow_pickle=True).item()
    outdict   = massage_dataset(SN_yields)
    dict_to_dataset('SN',outdict,info_str=info[SN_model])



    return
