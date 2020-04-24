"""

   Convert the old-style python tables into the new
   HDF5 format. This is done for backwards compatability
   purposes with previous runs.

"""
import numpy as np
import h5py
from galaxy_analysis import static_data as const


def load_ascii_data(fname):
    """
    Wrapper on genfromtxt in case we need to do fancy
    things later
    """
    return np.genfromtxt(fname)

def massage_dataset(data):
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


#
# A bad experiment... won't work without a lot of work 
#
#def generate_new_form_table(outname="IndividualStarYields_Merged"):
#    """
#    Produce a new-style yield table using the old yields. The
#    new style merges together all yields of individual types and removes
#    the need for contructing uniform M,Z across all yield types.
#
#    SN : core collapse SNe
#    winds : stellar winds, both AGB and massive stars
#    PopIII : popIII
#    """

#    outname = outname + ".h5"

#    hf = h5py.File(outname, 'w')

#    group_names = ["SN","winds","popIII"]

#    info = {'SN' : "NuGrid yields for CCSNe. Used in Emerick+2018-2020",
#            'AGB_winds' :
#            'winds' : "NuGrid yields for massive stars to 25 Msun"+\
#                      " combined with yields from Slemer+ (unpublished). Used in Emerick2018-2020",
#            'popIII' : "Heger+Woosley 2002 and 2010 yields for PopIII SNe and PISNe."}

#    fnames = {'SN' : 'stellar_yields_sn.in', 'winds' : 'stellar_yields_wind.in',
#              'massive_star' : 'stellar_yields_massive_star.in', 'popIII' : 'popIII_yields.in'}


    # do SNE

#    SN_data = load_ascii_data(fnames['SN'])

    # remove stars < 12
#    SN_data = SN_data[ SN_data[:,0] > 10 ]

#    grp = hf.create_group('SN')
#    outdict = massage_dataset(SN_data)
#    dict_to_dataset('SN',outdict, info_str=info['SN'])

    #
    # do winds
    #
#    winds_NuGrid = load_ascii_data(fnames['winds'])
#    winds_Slemer = load_ascii_data(fnames['massive_star'])
#
#    # stitch tables together. NuGrid 1 - 25, slemer above
#    winds_Slemer
#
#    return

def generate_basic_table(outname="IndividualStarYields"):


    outname = outname + ".h5"

    hf = h5py.File(outname, 'w')

    # need to load and create a few groups
    group_names = ["SN", "winds", "massive_star", "popIII"]

    info = {'SN' : "NuGrid Yields for core collapse SNe. Used in Emerick+2018-2020",
            'winds' : 'NuGrid wind yields for massive stars. Used in Emerick+2018-2020',
            'massive_star' : "Slemer+ (unpublished) yields for star winds above NuGrid tables. Used in Emerick+2018-2020",
            'popIII' : "Heger+Woosley 2002 and 2010 yields for PopIII SNe and PISNe."}

    fnames = {'SN' : 'stellar_yields_sn.in',
              'winds' : 'stellar_yields_wind.in',
              'massive_star' : 'stellar_yields_massive_star.in',
              'popIII' : 'popIII_yields.in'}

    # do basic things here
    for gn in group_names:

        grp = hf.create_group(gn)
        print(gn)
        # basic load data
        ascii_data = load_ascii_data(fnames[gn]) #, filetype = fname, cut_masses = merge_winds)

        # basic convert to dict
        outdict = massage_dataset(ascii_data)

        # write to H5 group
        dict_to_dataset(grp, outdict, info_str = info[gn])
    hf.close()

    print("Generated new yield table as %30s"%(outname))

    return



if __name__ == "__main__":

    generate_basic_table()

    generate_new_form_table()
