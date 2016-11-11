from __future__ import division


import numpy as np
import glob

_hdf5_compression = 'lzf'


def create_h5_file(directory, data, dataset_name):
    output_name = os.path.sep.join([directory, dataset_name+'.h5'])
    if os.path.exists(output_name):
        return output_name
    with h5py.File(output_name) as f:
        f.create_dataset(dataset_name, data = data, compression='lzf')
    return output_name

_all_fields = ['density', 'temperature', 'cell_mass']

def _create_data_structure(f, ds_list):
    """
    Create the general data structure
    """

    f.create_group("TimeEvolution")
    grp = f.create_group("Spatial")

    f.create_dataset('Time',

    for dsname in ds_list:
        dsname = dsname.rsplit('/')[0]
        subgrp = grp.create_group(dsname)

        _create_ds_subgroups(subgrp)


    return

def _create_ds_subgroups(grp):
    """
    Create the necessary sub groups
    """

    for subgrp_name in ['Gas', 'Stars']:
        sbgrp = grp.create_group(subgrp_name)
        sbgrp.create_group('RadialProfiles')
        sbgrp.create_group('Profiles')

    return


def _update_data_structure(f, ds_list):
    """
    Updates data structure in HDF5 file if there are 
    new data sets present that are not currently included
    """

    ds_in_hdf5 = f['Spatial'].keys()

    if len(ds_in_hdf5) == len(ds_list):
        # nothing to do here
        return False

    else:
        # find the missing keys
        if len(ds_list) < len(ds_in_hdf5):
            print "Available data sets less than those in HDF5 - is this OK"
            return False

        missing_ds = [x for x in ds_list if x not in ds_in_hdf5]

        for dsname in missing_ds:
            _create_ds_subgroups( f.create_group(dsname) )

    return True


class Galaxy(object):
    """
    Analyzes time evolution properties of a galaxy
    """

    def __init__(self, dir = './', h5file = None):

        self.ds_list = np.sort( glob.glob(dir + 'DD????/DD????') )

        self.ds0        = yt.load(ds_list[0])
        self.dsf        = yt.load(ds_list[-1])
        self.parameters = ds0.parameters

        if h5file is None:
            hf5file = dir + 'galaxy_data.h5'

        #
        # check if h5file exists -- if not, create it
        # and construct data structure groups
        #
        new_data = False
        if not os.path.isfile( h5file ):

            f = h5py.File( hf5file, 'r+')
            _create_data_structure(f, ds_list)

            new_data = True
        else:

            f = h5py.File( h5file, 'r+')
            new_data = _update_data_structure(f, ds_list)

        self.data = f

       # maybe store time values from all data sets here            

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

    def particle_radial_profile(self, fields = None, dsi= None):

    def sfr_from_last():

        if 'SFR' in self.data['TimeEvolution'].keys():
            if self.ds1.current_time > self.data['TimeEvolution']['SFR'][lastval]:
                self.data.__delitem__('TimeEvolution/SFR')
                del self.data['TimeEvolution/SFR']

            else:
                # up to date
                return
        
        # now compute the SFR
        times, sfr = particle_analysis.sfrFromParticles(self.dsf, self.dsf.all_data())


        self.data['TimeEvolution'].create_dataset('SFR', data=[times,sfr],
                                                  compression=_hdf5_compression)

        return

    def snr_from_last():

        for sn_type in ['SNII, SNIa']:
            rate = sn_type + '_rate'

            if rate in self.data['TimeEvolution'].keys():
                if self.dsf.current_time > self.data['TimeEvolution'][rate][lastval]:
                    self.data.__delitem__('TimeEvolution/' + rate)
                    del self.data['TimeEvolution/' + rate]

            else:
                return

            times, snr = particle_analysis.snrFromParticles(
    
