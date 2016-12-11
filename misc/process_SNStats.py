import numpy as np

_col_names = ['grid_id', 'pid', 'time', 'mass', 'birth_mass', 'metallicity',
             'm_eject', 'cells', 'volume', 'fraction', 'injected_mass',
             'energy_injected', 'ISM_mass', 'max_rho', 'avg_rho',
             'total_metal_mass', 'z_avg']

_dtypes = [int, int, float, float, float, float,
          float, int, float, float, float,
          float, float, float, float,
          float, float]

_ndtype = []
for i in np.arange(len(_col_names)):
    _ndtype += (col_names[i], dtypes[i])


def read_data(filename, filter = True):
    """
    Read in a single file and filter to individual SN events if
    filter is true.
    """


    return

def read_all_data(directory):
    """
    Go through all files in directory and combine
    """

    return

def plot_density(data, plot_type = 'max'):

    return


if __name__ == "__main__":

    directory = 

    read_all_data(directory)

