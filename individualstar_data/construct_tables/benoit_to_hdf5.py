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
import copy

from galaxy_analysis import static_data as const  # for elements and things

from galaxy_analysis.individualstar_data.construct_tables import compile_popIII

available_yields = {'FRUITY' : 'fruity_elem_stable.npy',
                    'LC18_R_0' : 'LC18_R_0_elem_stable.npy',
                    'LC18_R_150': 'LC18_R_150_elem_stable.npy',
                    'LC18_R_300': 'LC18_R_300_elem_stable.npy',
                    'LC18_mixture': 'LC18_mixture_elem_stable.npy',
                    'LC18_winds_R_0' : 'LC18_winds_R_0_elem_stable.npy',
                    'LC18_winds_R_150': 'LC18_winds_R_150_elem_stable.npy',
                    'LC18_winds_R_300': 'LC18_winds_R_300_elem_stable.npy',
                    'LC18_winds_mixture': 'LC18_winds_mixture_elem_stable.npy',
                    'Monash' : 'monash_elem_stable.npy',
                    'Nomoto' : 'nomoto_elem_stable.npy',
                    'Nomoto-HNe' : 'nomoto_HNe_elem_stable.npy',
                    'PopIII' : 'popIII_yields.in'}

info = {}
# temporary for now. better description later
for k in available_yields.keys():
    info[k] = available_yields[k]
verbose = True


def add_benoit_yields(setA, setB, factor = 1.0):
    """
    returns yields = setA + factor*setB.

    By default, factor is 1.0 (addition),
    but this can be easily switched to subtraction (factor = -1.0)
    """

    yields = copy.deepcopy(setA)

    for zkey in setA['yields'].keys():
        for Mkey in setA['yields'][zkey].keys():
            for e in setA['yields'][zkey][Mkey].keys():
                yields['yields'][zkey][Mkey][e] = setA['yields'][zkey][Mkey][e] +\
                                                  factor*setB['yields'][zkey][Mkey][e]

                if yields['yields'][zkey][Mkey][e] < 0.0:
                    print("Getting negative yields in subtract operation")
                    print(zkey,Mkey,e)
                    print(setA['yields'][zkey][Mkey][e])
                    print(setB['yields'][zkey][Mkey][e])
                    print(yields['yields'][zkey][Mkey][e])
                    raise RuntimeError


    return yields

def subtract_benoit_yields(setA, setB):
    """
    returns yields = setA - setB
    """

    return add_benoit_yields(setA,setB,factor = -1.0)


def massage_benoit_dataset(dataset):
    """
    Take kwargs and things from Benoits dict to what I need
    """

    outdict = {}

    outdict['M'] = dataset['M_list']
    outdict['Z'] = dataset['Z_list']

    flip_z = False
    if outdict['Z'][0] > outdict['Z'][-1]:
        # need to flip to be in ascending order
        outdict['Z'] = outdict['Z'][::-1]
        flip_z = True

    outdict['num_M'] = np.size(outdict['M'])
    outdict['num_Z'] = np.size(outdict['Z'])

    outdict['colnames'] = np.array(['Mtot','Metal_Mtot'] + dataset['elem_list'])
    outdict['ncol']     = np.size(outdict['colnames'])
    outdict['anums']    = np.array([-1,0] + [const.asym_to_anum[x] for x in dataset['elem_list']])
    outdict['num_Y']    = np.size(outdict['anums'])

    # probably a few different ways to do this, but most fool-proof is to loop:
    yields = np.zeros( (outdict['num_M'], outdict['num_Z'], outdict['num_Y']))

    # get metallicities:
    z_list = np.array(list(dataset['yields'].keys()))
    if flip_z:
        z_list = z_list[::-1]

    for i,Mkey in enumerate(list(dataset['yields'][ list(dataset['yields'].keys())[0] ].keys())):

        for j,zkey in enumerate(z_list):

            for k,e in enumerate(dataset['elem_list']):

                yields[i][j][k+2] = dataset['yields'][zkey][Mkey][e]

                yields[i][j][0] += yields[i][j][k+2]

                if (not (e=='H')) and (not (e=='He')):
                    yields[i][j][1] += yields[i][j][k+2]

    #
    outdict['yields'] = yields

    return outdict

def massage_ascii_dataset(data):
    """
    Given data loaded as ASCII table in old format,
    returns (in dictionary) the things needed to feed
    into the H5 file
    """

    outdict = {}

    outdict['M'] = np.unique( data[:,0])
    outdict['Z'] = np.unique( data[:,1])

    outdict['num_M'] = np.size(outdict['M'])
    outdict['num_Z'] = np.size(outdict['Z'])

    yields = data[:,2:]

    outdict['num_Y'] = np.shape(yields)[1]

    final_shape = (outdict['num_M'],outdict['num_Z'],outdict['num_Y'])
    outdict['yields'] = yields.reshape(final_shape)

    # all element names to Bi
    elements = list(const.anum_to_asym.values())
    outdict['colnames'] = np.array(['Mtot','Metal_Mtot']+elements)
    outdict['ncol'] = np.size(outdict['colnames'])

    outdict['anums'] = [-1,0] + list(const.anum_to_asym.keys())

    return outdict

def massage_popIII_dataset(low_mass_model='heger_woosley_2010',
                           pisne_model = 'heger_woosley_2010'):
    """
    Given data from popIII functions, convert ot format needed in a dict
    """

    #yields, M = compile_popIII.generate_table(low_mass_model, pisne_model)

    data   = np.genfromtxt('popIII_yields.in')
    M      = np.unique(data[:,0])
    yields = data[:,2:]

    outdict = {}

    outdict['M'] = M
    outdict['Z'] = np.array([0])

    outdict['num_M'] = np.size(M)
    outdict['num_Z'] = 1

#    outdict['num_Y'] = np.shape(yields)[0]
    outdict['num_Y'] = np.shape(yields)[1]

    final_shape = (outdict['num_M'], outdict['num_Z'], outdict['num_Y'])
    outdict['yields'] = yields.reshape(final_shape)

    # all element names to Bi
    elements = list(const.anum_to_asym.values())
    outdict['colnames'] = np.array(['Mtot','Metal_Mtot']+elements)
    outdict['ncol'] = np.size(outdict['colnames'])
    outdict['anums'] = [-1,0] + list(const.anum_to_asym.keys())

    return outdict

def dict_to_dataset(grp, outdict, info_str = ''):
    """
    Take the converted dictionary and write it to the HDF5 group
    """

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
    model = {'SN' : SN_model, 'AGB' : AGB_model, 'Wind' : wind_model,
             'PopIII' : PopIII_model}

    # special cases
    group_list = ['Wind','AGB']
    if not ('LC18' in SN_model): # this needs special treatment

        group_list = group_list + ['SN']

    # loop over all the ones that don't need special treatment
    for gname in group_list:
        grp       = hf.create_group(gname)
        yields    = np.load( wdir + available_yields[model[gname]], allow_pickle=True).item()
        outdict   = massage_benoit_dataset(yields)

        print("Name = %8s - Model = %20s: Nm = %i Nz = %i Ny = %i"%(gname, model[gname],
                                                                    outdict['num_M'],
                                                                    outdict['num_Z'],
                                                                    outdict['num_Y']))
        if verbose:
            print("M: ", outdict['M'])
            print("Z: ", outdict['Z'])
            print("--------------------")

        dict_to_dataset(grp,outdict,info_str=info[ model[gname]])

    #
    # If LC18 models are used for SN, wind yields need to be subtracted from total
    #
    if ('LC18' in SN_model):
        print("Generating LC18 SN table from total - winds")
        gname = 'SN'
        grp = hf.create_group(gname)
        total_yields = np.load(wdir + available_yields[model[gname]], allow_pickle=True).item()
        wind_yields  = np.load(wdir + available_yields[model["Wind"]], allow_pickle=True).item()

        yields  = subtract_benoit_yields(total_yields, wind_yields)
        outdict = massage_benoit_dataset(yields)

        print("Name = %8s - Model = %20s: Nm = %i Nz = %i Ny = %i"%(gname, model[gname],
                                                                    outdict['num_M'],
                                                                    outdict['num_Z'],
                                                                    outdict['num_Y']))
        if verbose:
            print("M: ", outdict['M'])
            print("Z: ", outdict['Z'])
            print("--------------------")

        dict_to_dataset(grp, outdict,info_str=info[model[gname]])

    #
    # PopIII table
    #
    gname   = 'PopIII'
    grp     = hf.create_group(gname)
    outdict = massage_popIII_dataset('heger_woosley_2010','heger_woosley_2002')
    print("Name = %8s - Model = %20s: Nm = %i Nz = %i Ny = %i"%(gname, model[gname],
                                                                outdict['num_M'],
                                                                outdict['num_Z'],
                                                                outdict['num_Y']))
    if verbose:
        print("M: ", outdict['M'])
        print("Z: ", outdict['Z'])
        print("--------------------")

    popIII_info = "Yields from Heger Woosley 2002 for PISNe and Heger Woosley 2010 for SNe"
    dict_to_dataset(grp,outdict,info_str = popIII_info)

    hf.close()

    return


if __name__ == "__main__":

    construct_table()
    construct_table(SN_model = 'LC18_R_150', wind_model = 'LC18_winds_R_150')
    construct_table(SN_model = 'LC18_R_300', wind_model = 'LC18_winds_R_300')
    construct_table(SN_model = 'LC18_mixture', wind_model = 'LC18_winds_mixture')
