"""
HDF5 in memory object

Britton Smith <brittonsmith@gmail.com>
"""

import h5py
import numpy as np
import sys

class H5InMemory(object):
    def __init__(self, fh):
        self.attrs = {}
        if fh is None:
            self.data = {}
            return
        
        if isinstance(fh, str):
            fh = h5py.File(fh, "r")

        if hasattr(fh, "attrs"):
            self.attrs = dict([(attr, fh.attrs[attr]) \
                               for attr in fh.attrs])

        if isinstance(fh, h5py.Dataset):
            self.value = fh.value
        elif isinstance(fh, np.ndarray):
            self.value = fh
        else:
            self.data = dict([(field, H5InMemory(fh[field])) \
                               for field in fh.keys()])

        if isinstance(fh, h5py.File):
            fh.close()

    def create_group(self, key):
        self.data[key] = H5InMemory(None)
        return self.data[key]

    def create_dataset(self, key, data=None):
        self.data[key] = H5InMemory(data)
        return self.data[key]
            
    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key, value):
        self.data[key] = value

    def __delitem__(self, key):
        del self.data[key]

    def __iter__(self):
        for key in self.data:
            yield key

    def iterkeys(self):
        for key in self.data:
            yield key
        
    def __contains__(self, key):
        return key in self.data

    def __repr__(self):
        if hasattr(self, "value"):
            return "<H5InMemory Data object %s>" % \
              str(self.value.shape)
        else:
            return "<H5InMemory Group object (%d items)>" % \
              len(self.keys())

    def __str__(self):
        return self.__repr__()

    def __iter__(self):
        for field in self.keys():
            yield field

    def keys(self):
        if hasattr(self, "data"):
            return self.data.keys()
        return None

    def save(self, fh):
        top = False
        if isinstance(fh, str):
            top = True
            fh = h5py.File(fh, "w")

        for attr in self.attrs:
            fh.attrs[attr] = self.attrs[attr]

        if hasattr(self, "data"):
            for field in self:
                if hasattr(self.data[field], "data"):
                    self.data[field].save(fh.create_group(field))
                else:
                    dfh = fh.create_dataset(field,
                                            data=self.data[field].value)
                    self.data[field].save(dfh)

        if top:
            fh.close()
