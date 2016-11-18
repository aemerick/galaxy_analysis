from __future__ import division

import yt
import numpy as np
import glob
import os
import h5py

from collections import Iterable


# --------- internal imports --------------
from utilities import utilities as util
from static_data import LABELS,\
                        FIELD_UNITS,\
                        IMAGE_COLORBAR_LIMITS,\
                        PLOT_LIMITS

from particle_analysis import particle_types as pt
from particle_analysis import IMF

from yt_fields import field_generators as fg

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

    f.create_dataset('Time')

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


class GalaxyAnalyzer(object):
    """
    Analyzes time evolution properties of a galaxy
    """

    def __init__(self, dir = './', h5file = None):

        self.ds_list = np.sort( glob.glob(dir + 'DD????/DD????') )

        self.ds0        = yt.load(ds_list[0])
        self.dsf        = yt.load(ds_list[-1])
        self.parameters = ds0.parameters

        _set_data_region_properties()


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
        raise NotImplementedError

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

        return
            #profile_1d = yt.create_profile(data, ['density', 'temperature',
            #                                      '

    def particle_radial_profile(self, fields = None, dsi= None):
        raise NotImplementedError
        return

    def sfr_from_last(self):
        """ Calculate SFR from final data dump """

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

    def snr_from_last(self):
        """ Calculate the SNR from final data dump """

        raise NotImplementedError

        for sn_type in ['SNII, SNIa']:
            rate = sn_type + '_rate'

            if rate in self.data['TimeEvolution'].keys():
                if self.dsf.current_time > self.data['TimeEvolution'][rate][lastval]:
                    self.data.__delitem__('TimeEvolution/' + rate)
                    del self.data['TimeEvolution/' + rate]

            else:
                return

#            times, snr = particle_analysis.snrFromParticles(
        return

class Galaxy(object):

    def __init__(self, dsname, dir = './'):
        """
        Object for individual data dump analysis
        """

        self.dir    = dir
        self.dsname = dsname

        # load, generate fields, reload
        self.ds     = yt.load(self.dir + '/' + self.dsname + '/' + self.dsname)
        fg.generate_derived_fields(self.ds)
        self.ds     = yt.load(self.dir + '/' + self.dsname + '/' + self.dsname)

        self.df     = self.ds.all_data()

        hdf5_file   = self.dir + '/' + self.dsname + '_galaxy_data.h5'

        self._set_data_region_properties()
        self.species_list = util.species_from_fields(self.ds.field_list)

        self._set_accumulation_fields()

        self.particle_meta_data = {}
        self.gas_meta_data      = {}
        self.meta_data          = {}

        self.construct_regions()

        self.total_quantities = {}
        self._total_quantities_calculated = False


        if not os.path.isfile( hdf5_file ):
            f = h5py.File(hdf5_file, 'w')

            self._create_data_structure()

        else:

            f = h5py.File(hdf5_file, 'a')
            self._update_data_structure()

        self.data = f

        return

    def _create_data_structure(self):
        print 'does nothing for now'
        return

    def _update_data_structure(self):
        print "does nothing for now"
        return

    def calculate_projected_disk_field(self, field, **kwargs):
        """
        Calculate the projected surface density of a given field.
        Examples:
           For total gas density:
           > galaxy.calculate_surface_density(('gas','Density'))
           For H2 density:
           > galaxy.calculate_surface_density(('gas','H2I_Density'))
        """

        proj = ds.proj( field, 'z', data_source = self.disk, **kwargs)

        # return raw set to handle elsewhere
        return proj

    def calculate_total_quantities(self, fields = None, *args, **kwargs):
        """
        Computes and saves total quantities
        """

        if fields is None:
            # calculate all
            fields = self._accumulation_fields

        for field in fields:
            self.total_quantities[field] = np.sum(self.df[field].convert_to_units(FIELD_UNITS[field].units))

        self._total_quantities_calculated = True

        return

    def calculate_mass_fraction_profile(self, fields = None, *args, **kwargs):

        if fields is None:
            fields = self._accumulation_fields

        rbins, centers, profiles = self.calculate_mass_profile(fields = fields)

        compute_total_fields = [x for x in fields if x not in self.total_quantities]
        self.calculate_total_quantities(fields = compute_total_fields)

        for field in fields:
            profiles[field] = profiles[field] / self.total_quantities[field]

        return rbins, centers, profiles

    def calculate_mass_profile(self, fields = None, *args, **kwargs):
        """
        Compute mass fraction of given species contained within spherical
        bins out to the halo radius. Used to estimate metal retention fractions
        of given species.
        """

        if fields is None:
            fields = self._accumulation_fields

        rbins = self.rbins_halo_sphere

        # data
        r = self.halo_sphere['spherical_r'].convert_to_units(rbins.units)

        profiles = {}

        for field in fields:
            profiles[field] = np.zeros(np.size(rbins)-1)

        for i in np.arange(np.size(rbins)-1):
            radial_filter = (r >= rbins[i]) * (r < rbins[i+1])

            for field in fields:
                profiles[field][i] = np.sum(\
                      self.halo_sphere[field][radial_filter].convert_to_units(FIELD_UNITS[field].units))

        centers = 0.5 * (rbins[1:] + rbins[:-1])

        return rbins, centers, profiles


    def calculate_surface_density_profile(self, fields, *args, **kwargs):
        """
        Computes a 1D surface density profile in a cylindrical region of the galaxy
        using already set disk selection region.
        """

        if not isinstance(fields, Iterable):
            fields = [fields]

        # project the fields
        proj = self.ds.proj(fields, 'z', data_source = self.disk, **kwargs)

        r    = np.sqrt(( (proj['px']-self.disk.center[0])**2 + (self.disk.center[1] -proj['py'])**2)).convert_to_units('pc')
        A    = (proj['pdx']*proj['pdy']).convert_to_units('pc**2')

        # set up bins
        rbins = self.rbins_disk

        profiles = {}

        # construct the profiles
        for field in fields:
            profiles[field] = np.zeros(np.size(rbins) - 1)

        for i in np.arange(np.size(rbins)-1):
            radial_filter = (r >= rbins[i]) * (r < rbins[i+1])

            annulus_area  = 2.0 * np.pi * (rbins[i+1]**2 - rbins[i]**2)

            if np.size(proj['dx'][radial_filter]) > 0:

                for field in fields:

                    projection = np.array(proj[field][radial_filter].convert_to_units('Msun/pc**2'))

                    profiles[field][i] =\
                       (np.sum(projection * np.array(A[radial_filter]))\
                                 / annulus_area.value) * yt.units.Msun / yt.units.pc**2

        centers = 0.5 * (rbins[:-1] + rbins[1:])

        return rbins, centers, profiles

    def compute_all_meta_data(self):
        """
        Wrapper to compute all meta data
        """

        self.compute_meta_data()
        self.compute_particle_meta_data()
        self.compute_gas_meta_data()

        return

    def compute_particle_meta_data(self):
        """
        Computes and saves all particle meta data for this time step.
        This is:
            1) Total mass, and total mass of each particle type
            2) Total number, and number of each particle type
            3) Number of each: core collapse, SNIa, AGB, direct collapse
            4)

        """

        MS_stars = pt.main_sequence(self.ds, self.df)
        WD_stars = pt.main_sequence(self.ds, self.df)
        SNII     = pt.core_collapse(self.ds, self.df)
        SNIa     = pt.snIa(self.ds, self.df)

        particle_mass = self.df['birth_mass'].convert_to_units(UNITS['Mass'].units)

        self.particle_meta_data['total_mass']    = np.sum(particle_mass)
        self.particle_meta_data['total_mass_MS'] = np.sum(particle_mass[MS_stars])
        self.particle_meta_data['total_number']  = np.size(particle_mass)
        self.particle_meta_data['total_number_MS'] = np.size(particle_mass[MS_stars])

        self.particle_meta_data['total_mass_WD'] = np.sum(particle_mass[WD_stars])
        self.particle_meta_data['total_number_WD'] = np.size(particle_mass[WD_stars])

        self.particle_meta_data['total_number_SNII'] = np.size(particle_mass[SNII])
        self.particle_meta_data['total_number_SNIa'] = np.size(particle_mass[SNIa])

        self.particle_meta_data['N_OTRAD']           = np.size(particle_mass[(particle_mass>ds.parameters['IndividualStarOTRadiationMass'])*\
                                                                             (MS_stars)])
        self.particle_meta_data['N_ionizing']        = np.size(particle_mass[(particle_mass>ds.parameters['IndividualStarRadiationMinimumMass'])*\
                                                                             (MS_stars)])

        self.particle_meta_data['metallicity_stas'] = util.compute_stats(self.df['metallicity_fraction'])

        # compute theoretical total IMF and 'observed' IMF of just MS stars
        self.particle_meta_data['IMF_obs']        = IMF.compute_IMF(ds,data, mode='mass',       bins=25)
        self.particle_meta_data['IMF_birth_mass'] = IMF.compute_IMF(ds,data, mode='birth_mass', bins=25)

        # in principle store abundance information here, but might just leave this as separate

        return

    def compute_gas_meta_data(self):


        return

    def compute_meta_data(self):
        """
        Computes general meta data information and stores it as a
        dictionary. Really should just be looking into parameter file
        """

        self.meta_data['Time']  = self.ds.current_time.convert_to_units(UNITS['Time'].units)
        self.meta_data['dx']    = np.min(self.df['dx'].convert_to_units(UNITS['Length'].units)

        return

    def compute_time_evolution(self):
        """
        Computes current SFR, SNR, and SFH from particles
        """

        self.time_data['time'], self.time_data['SFR'] = sfrFromParticles(ds,data)
        self.time_data['time'], self.time_data['SFH'] = sfhFromParticles(ds, data, times=self.time_data['times'])
        self.time_data['time'], self.time_data['SNII_snr'] = snr(ds, data,times=self.time_data['times'], sn_type ='II')
        self.time_data['time'], self.time_data['SNIa_snr'] = snr(ds, data,times=self.time_data['times'], sn_type ='Ia')

        return

    def compute_everything(self):
        """
        Compute everything we need to add to the HDF5 files:

          - Meta quantities:
            1) Total masses of all species
            2) Total stellar mass
            3) 

            1) Mass profiles of all species
            2)
        """


        return

    def construct_regions(self, disk_kwargs = None, sphere_kwargs = None,
                                halo_sphere_kwargs = None):

        if disk_kwargs is None:
            disk_kwargs = {}

            for n in ['normal','radius','height','center']:
                disk_kwargs[n] = self.disk_region[n]

        if sphere_kwargs is None:
            sphere_kwargs = {}

            for n in ['radius','center']:
                sphere_kwargs[n] = self.spherical_region[n]

        if halo_sphere_kwargs is None:
            halo_sphere_kwargs = {}

            for n in ['radius','center']:
                halo_sphere_kwargs[n] = self.halo_spherical_region[n]


        self.disk   = self.ds.disk(**disk_kwargs)

        self.sphere = self.ds.sphere(**sphere_kwargs)

        self.halo_sphere = self.ds.sphere(**halo_sphere_kwargs)

        return

    def _set_data_region_properties(self):
        """
        Sets parameters used to define disk and spherecal data regions
        Separate function currently useless, but allows this to be
        expanded so one can guess sizes from parameter file settings,
        rather than hard coding as is done below.
        """

        # dr and dz set spacing for cylindrical or spherical regions
        # to use to construct shells

        self.disk_region = {'normal' : np.array([0.0, 0.0, 1.0]),
                            'radius' : 2.0 * yt.units.kpc,
                            'height' : 500.0 * yt.units.pc,
                            'center' : self.ds.domain_center,
                            'dr'     : 25.0 * yt.units.pc,
                            'dz'     : 100.0 * yt.units.pc }

        self.spherical_region = {'center' : self.ds.domain_center,
                                 'radius' : 2.0 * yt.units.kpc,
                                 'dr'     : 25.0 * yt.units.pc   }

        self.halo_spherical_region = {'center' : self.ds.domain_center,
                                      'radius' : 5.0 * yt.units.kpc,
                                      'dr'     : 50.0 * yt.units.pc}

        return

    def _set_accumulation_fields(self):

        self._accumulation_fields = [('gas','H_total_mass'), ('gas','He_total_mass'),
                                     ('gas','metal_mass')]

        for e in self.species_list:
            self._accumulation_fields += [('gas', e +'_Mass')]

        return


    @property
    def rbins_sphere(self):

        rmin = 0.0
        rmax = self.sphere.radius
        dr   = self.spherical_region['dr']
        rmax = rmax.convert_to_units(dr.units)

        return np.arange(rmin, rmax + dr, dr) * dr.unit_quantity

    @property
    def rbins_halo_sphere(self):
        rmin = 0.0
        rmax = self.halo_sphere.radius
        dr   = self.halo_spherical_region['dr']
        rmax = rmax.convert_to_units(dr.units)

        return np.arange(rmin, rmax + dr, dr) * dr.unit_quantity

    @property
    def zbins_disk(self):
        zmin = 0.0
        zmax = self.disk.height * 0.5
        dz   = self.disk_region['dz']
        zmax = zmax.convert_to_units(dz.units)

        return np.arange(zmin, zmax + dz, dz) * dz.unit_quantity


    @property
    def rbins_disk(self):
        rmin = 0.0
        rmax = self.disk.radius
        dr   = self.disk_region['dr']
        rmax = rmax.convert_to_units(dr.units)

        return np.arange(rmin, rmax + dr, dr) * dr.unit_quantity
