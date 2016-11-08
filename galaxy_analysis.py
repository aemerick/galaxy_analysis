from __future__ import division


import numpy as np



def create_h5_file(directory, data, dataset_name):
    output_name = os.path.sep.join([directory, dataset_name+'.h5'])
    if os.path.exists(output_name):
        return output_name
    with h5py.File(output_name) as f:
        f.create_dataset(dataset_name, data = data, compression='lzf')
    return output_name


class GalaxyEvolution(object):
    """
    Analyzes time evolution properties of a galaxy
    """

    def __init__(self):
        return
