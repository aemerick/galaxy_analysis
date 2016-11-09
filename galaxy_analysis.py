from __future__ import division


import numpy as np
import glob


def create_h5_file(directory, data, dataset_name):
    output_name = os.path.sep.join([directory, dataset_name+'.h5'])
    if os.path.exists(output_name):
        return output_name
    with h5py.File(output_name) as f:
        f.create_dataset(dataset_name, data = data, compression='lzf')
    return output_name

_all_fields = ['density', 'temperature', 'cell_mass']

class Galaxy(object):
    """
    Analyzes time evolution properties of a galaxy
    """

    def __init__(self, dir = './'):

        self.ds_list = np.sort( glob.glob(dir + 'DD????/DD????') )

        ds      = yt.load(ds_list[-1])
        species = [x[1] for x in fields if ('_Density') in x[1]]
        self.species = np.sort( [x.rsplot('_')[0] for x in species])

        return


    def create_radial_profile(self, fields = None, dsi = None):

        if dsi is None:
            data_files = ds_list
        else:
            data_files = ["DD%0004i/DD%0004i"%(i)]

        #
        # loop
        #
        for dsname in ds_list:
            ds   = yt.load(dsname)
            data = ds.all_data()

            profile_1d = yt.create_profile(data, ['density', 'temperature',
                                                  '
