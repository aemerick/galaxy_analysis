"""
   Take the LC18 yields and average them together using an optional upsampling
   in metallicity. By default, use Prantzos+18 weighting scheme to do the
   averaging.
"""

import numpy as np
import h5py
import copy


from galaxy_analysis.individualstar_data.construct_tables.benoit_to_hdf5 import available_yields
