import numpy as np
import h5py
from scipy import interpolate

# Grackle data file is 3D: n, T, cooling-rate
#    This is split between a primordial cooling
#    rate and a metal cooling rate
#
#    Table is at solar metallicity and assumed
#    to scale linearly with Z.
#
# File is located in grackle repo:
#   $GRACKLE_INSTALL/input/CloudyData_noUVB.h5
#

class CoolingTable:

    def __init__(file = './CloudyData_noUVB.h5'):
        self.file_path = file
        self.load_data()

        # hard-code n and T sample values in table
        self.n = np.linspace(-10,4,29)
        self.T = np.logspace(1,9,161)

        return

    def load_data(file = './CloudyData_noUVB.h5'):

        data = h5.File(self.file_path)
        self.metal_cooling = data['CoolingRates']['Metals']['Cooling']
        self.prim_cooling  = data['CoolingRates']['Metals']['Cooling']

        return



