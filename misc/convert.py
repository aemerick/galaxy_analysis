"""
Fixes issues in how I stored strings / complex keys in
dictionaries in HDF5 files in python2 that can be weird
in python3.
"""


def convert_keys_to_string(dictionary):
    """Recursively converts dictionary keys to strings."""
    if not isinstance(dictionary, dict):
        return dictionary
    return dict((str(k), convert_keys_to_string(v)) 
        for k, v in dictionary.items())

import glob, os
import numpy as np
import deepdish as dd
import copy

files = np.sort(glob.glob('DD????_galaxy_data.h5'))

if not (os.path.exists('py3temp/')):
    os.mkdir('py3temp/')

for f in files:
    if os.path.exists('py3temp/'+f):
        continue
    try:
        data = dd.io.load(f)
    except:
        print("--------------Loading failed for file", f)
    data2 = copy.deepcopy(data)
    data3 = convert_keys_to_string(data2)
    dd.io.save('py3temp/' + f,data3)
    del(data)
    del(data2)
    del(data3)
    print(f)

