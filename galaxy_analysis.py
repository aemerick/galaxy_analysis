from __future__ import division

import yt
import numpy as np
from scipy import optimize
import glob
import os
import h5py
import copy

from collections import Iterable

# this is used to save / load generated
# hierarchical ditionaries of data very easily
# in / out of HDF5 files
import deepdish as dd


# --------- internal imports --------------
from utilities import utilities
from static_data import LABELS,\
                        FIELD_UNITS,\
                        IMAGE_COLORBAR_LIMITS,\
                        PLOT_LIMITS,\
                        UNITS,\
                        ISM

from particle_analysis import particle_types as pt
from particle_analysis import IMF

# need to have better API
from particle_analysis.sfrFromParticles import sfrFromParticles
from particle_analysis.sfhFromParticles import sfhFromParticles
from particle_analysis.sn_rate          import snr

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
        times, sfr = sfrFromParticles(self.dsf, self.dsf.all_data())


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
        with utilities.nooutput(): # silence yt for now
            self.ds     = yt.load(self.dir + '/' + self.dsname + '/' + self.dsname)
        fg.generate_derived_fields(self.ds)
        self.ds     = yt.load(self.dir + '/' + self.dsname + '/' + self.dsname)

        self.df     = self.ds.all_data()

        self.hdf5_filename   = self.dir + '/' + self.dsname + '_galaxy_data.h5'

        self._set_data_region_properties()
        self.species_list = utilities.species_from_fields(self.ds.field_list)

        self._set_accumulation_fields()

        self.particle_meta_data = {}
        self.gas_meta_data      = {}
        self.meta_data          = {}
        self.gas_profiles       = {}
        self.time_data          = {}

        self.construct_regions()

        self.total_quantities = {}
        self._total_quantities_calculated = False

        if os.path.isfile( self.hdf5_filename ):
            self.load()

#        if not os.path.isfile( hdf5_file ):
#            f = h5py.File(hdf5_file, 'w')
#
#            self._create_data_structure()
#
#        else:
#
#            f = h5py.File(hdf5_file, 'a')
#            self._update_data_structure()
#
#        self.data = f

        return

    def load(self, filename = None, nocopy = False):
        if filename is None:
            filename = self.hdf5_filename
        else:
            self.hdf5_filename = filename

        if os.path.isfile(filename):
            if nocopy:
                return dd.io.load(self.hdf5_filename)
            else:
                self._output_data_dict = dd.io.load(self.hdf5_filename)
                self._map_output_to_class()

        else:
            print "No hdf5 output file exists at " + filename

        return

    def save(self, filename = None):
        if filename is None:
            filename = self.hdf5_filename
        else:
            self.hdf5_filename = filename

        self._map_class_to_output()
 
        dd.io.save(filename, self._output_data_dict)

        return

    def _map_class_to_output(self):
        """
        Map the class structure to the output file
        """

        self._output_data_dict['gas_meta_data']      = self.gas_meta_data
        self._output_data_dict['particle_meta_data'] = self.particle_meta_data
        self._output_data_dict['time_data']          = self.time_data

        return

    def _map_output_to_class(self):

        def _verify_and_add(x, name):
            if name in self._output_data_dict.keys():
                return self._output_data_dict[name]
            else:
                return x

        self.gas_meta_data      = _verify_and_add(self.gas_meta_data, 'gas_meta_data')
        self.particle_meta_data = _verify_and_add(self.particle_meta_data, 'particle_meta_data')
        self.time_data          = _verify_and_add(self.time_data, 'time_data')

        return

    def _update_data_structure(self):
        print "does nothing for now"
        return

    def calculate_projected_disk_field(self, field, axis = 'z', **kwargs):
        """
        Calculate the projected surface density of a given field.
        Examples:
           For total gas density:
           > galaxy.calculate_surface_density(('gas','Density'))
           For H2 density:
           > galaxy.calculate_surface_density(('gas','H2I_Density'))
        """

        proj = ds.proj( field, axis, data_source = self.disk, **kwargs)

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


    def calculate_dMdt_profile(self, fields = None, mode = 'sphere', n_cell = 4,
                               outflow = True, *args, **kwargs):
        """
        Returns the mass inflow or outflow rate as a function of radius. This can be used
        to compute mass loading factor by dividing by the current SFR.

        outflow == True (default) computes outflow rate (requiring that v > 0.0). False
        computes inflow rate (requiring that v < 0.0).

        if no fields are provided, computes total mass flow rate and rate for
        all species.
        """

        if fields is None:
            fields = self._accumulation_fields

        if mode == 'sphere':
            xbins  =  self.rbins_halo_sphere
            xdata  =  self.halo_sphere['spherical_radius'].convert_to_units(UNITS["Length"].units)
            vel    =  self.halo_sphere['velocity_spherical_radius'].convert_to_units(UNITS['Velocity'].units)
            data   =  self.halo_sphere

        elif mode == 'disk':
            xbins  =  self.zbins_disk
            xdata  =  (self.disk['z'] - self.disk.center[2]).convert_to_units(UNITS["Length"].units)
            vel    =  self.disk['velocity_cylindrical_z'].convert_to_units(UNITS['Velocity'].units)
            data   =  self.disk

        else:
            raise ValueError("Must choose either disk or sphere for mass outflow profile mode")

        profile = {}

        for field in fields:
            profile[field] = np.zeros(np.size(xbins)-1)

        center = 0.5 * (xbins[1:] + xbins[:-1])

        dx = n_cell * np.min(data['dx'].convert_to_units(xdata.units))

        if outflow: # compute outflow
            v_filter = vel > 0.0
        else:       # compute inflow
            v_filter = vel < 0.0

        for i in np.arange(np.size(xbins)-1):
            x_filter = ( xdata >= (center[i] - 0.5*dx)) * ( xdata < (center[i] + 0.5*dx))

            filter = x_filter * v_filter
            for field in fields:
                M    = data[field].convert_to_units(UNITS['Mass'].units)
                Mdot     = np.sum( M[filter] * vel[filter] ) / dx
                Mdot     = Mdot.convert_to_units('Msun/yr')

                profile[field][i] = Mdot

        return xbins, center, profile

    def calculate_mass_fraction_profile(self, fields = None, *args, **kwargs):
        """
        Computes fractional radial mass profiles for all species. Can be easily
        used to make cumulative profiles
        """

        #
        # check if profiles exist already and don't recalculate unless ordered to
        # 
        if fields is None:
            fields = self._accumulation_fields

        rbins, centers, profiles = self.calculate_mass_profile(fields = fields)

        compute_total_fields = [x for x in fields if x not in self.total_quantities]
        self.calculate_total_quantities(fields = compute_total_fields)

        for field in fields:
            profiles[field] = profiles[field] / self.total_quantities[field]

        # save fractional profiles here?

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

        #
        # save and store profiles here?
        #

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
        self.compute_time_evolution()

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

        particle_mass = (self.df['birth_mass'].value *
                         yt.units.Msun).convert_to_units(UNITS['Mass'].units)

        self.particle_meta_data['total_mass']    = np.sum(particle_mass)
        self.particle_meta_data['total_mass_MS'] = np.sum(particle_mass[MS_stars])
        self.particle_meta_data['total_number']  = np.size(particle_mass)
        self.particle_meta_data['total_number_MS'] = np.size(particle_mass[MS_stars])

        self.particle_meta_data['total_mass_WD'] = np.sum(particle_mass[WD_stars])
        self.particle_meta_data['total_number_WD'] = np.size(particle_mass[WD_stars])

        self.particle_meta_data['total_number_SNII'] = np.size(particle_mass[SNII])
        self.particle_meta_data['total_number_SNIa'] = np.size(particle_mass[SNIa])

        self.particle_meta_data['N_OTRAD']           = np.size(particle_mass[(particle_mass>self.ds.parameters['IndividualStarOTRadiationMass'])*\
                                                                             (MS_stars)])
        self.particle_meta_data['N_ionizing']        = np.size(particle_mass[(particle_mass>self.ds.parameters['IndividualStarIonizingRadiationMinimumMass'])*\
                                                                             (MS_stars)])

        self.particle_meta_data['metallicity_stars'] = utilities.compute_stats(self.df['metallicity_fraction'])

        # compute theoretical total IMF and 'observed' IMF of just MS stars
        self.particle_meta_data['IMF_obs']        = IMF.compute_IMF(self.ds, self.df, mode='mass',       bins=25)
        self.particle_meta_data['IMF_birth_mass'] = IMF.compute_IMF(self.ds, self.df, mode='birth_mass', bins=25)

        self.compute_half_light_radius()

        # in principle store abundance information here, but might just leave this as separate

        return

    def compute_gas_meta_data(self):

        self.gas_meta_data['masses'] = self.compute_gas_sequestering()

        return

    def compute_gas_profiles(self):

        return

    def compute_gas_sequestering(self):
        """
        Computes the sequestering of gas into various categ ories (geometric and phase) for 
        all species as well as the total gas mass. Returns a dictionary of dictionaries
        containing all of this information
        """

        mdict = {}
        
        cut_region_names = ['Molecular', 'CNM', 'WNM', 'WIM', 'HIM']        
        fields = {'H':'H_total_mass','He':'He_total_mass','Total':'cell_mass'}

        # do this for the disk ISM regions
        for crtype in cut_region_names:
            mdict[crtype] = {}
            for s in fields:
                mdict[crtype][s] = np.sum(self.disk.cut_region( ISM[crtype])[fields[s]]).convert_to_units('Msun')

            for s in self.species_list:
                mdict[crtype][s] = np.sum(self.disk.cut_region( ISM[crtype])[('gas',s + '_Mass')]).convert_to_units('Msun')

        # now do this for the whole disk
        mdict['Disk'] = {}
        for s in fields:
            mdict['Disk'][s] = np.sum(self.disk[fields[s]]).convert_to_units('Msun')
        for s in self.species_list:
            mdict['Disk'][s] = np.sum(self.disk[('gas',s + '_Mass')]).convert_to_units('Msun')
   
        # now do this for the halo
        mdict['Halo'] = {}
        for s in fields:
            mdict['Halo'][s] = np.sum(self.halo_sphere[fields[s]]).convert_to_units('Msun')
        for s in self.species_list:
            mdict['Halo'][s] = np.sum(self.halo_sphere[('gas',s + '_Mass')]).convert_to_units('Msun')

        # now do this for full box
        mdict['FullBox'] = {}
        for s in fields:
            mdict['FullBox'][s] = np.sum(self.df[fields[s]]).convert_to_units('Msun')
        for s in self.species_list:
            mdict['FullBox'][s] = np.sum(self.df[('gas', s + '_Mass')]).convert_to_units('Msun')

        # now we need to do some subtraction of the fields
        mdict['OutsideHalo'] = {}
        for s in fields.keys() + self.species_list:
            mdict['OutsideHalo'][s] = mdict['FullBox'][s] - mdict['Halo'][s]
            mdict['Halo'][s]        = mdict['Halo'][s]    - mdict['Disk'][s]

        # now we compute the gravitationally bound gas IF potential is present
        if 'PotentialField' in self.ds.field_list or 'GravPotential' in self.ds.field_list:
            mdict['GravBound'] = {}
            for s in fields:
                mdict['GravBound'][s] = np.sum( self.ds.cut_region(self.df, "obj[('gas','gravitationally_bound')] > 0" )[fields[s]]).convert_to_units('Msun')
            for s in self.species_list:
                mdict['GravBound'][s] = np.sum(self.ds.cut_region(self.df, "obj[('gas','gravitationally_bound')] > 0")[('gas', s + '_Mass')]).convert_to_units('Msun')

        self.gas_sequestering = mdict
        return mdict

    @property
    def fractional_gas_sequestering(self):
        if not hasattr(self, 'gas_sequestering'):
            discard = self.compute_gas_sequestering()
        
        fields = self.gas_sequestering['Disk'].keys()
        x      = copy.deepcopy(self.gas_sequestering)
        for region in self.gas_sequestering.keys():
            for s in fields:
                x[region][s] /= self.gas_sequestering['FullBox'][s]

        return x        

    def compute_meta_data(self):
        """
        Computes general meta data information and stores it as a
        dictionary. Really should just be looking into parameter file
        """

        self.meta_data['Time']  = self.ds.current_time.convert_to_units(UNITS['Time'].units)
        self.meta_data['dx']    = np.min(self.df['dx'].convert_to_units(UNITS['Length'].units))

        return

    def compute_time_evolution(self):
        """
        Computes current SFR, SNR, and SFH from particles
        """

        if not hasattr(self, 'time_data'):
            self.time_data = {}

        self.time_data['time'], self.time_data['SFR'] = sfrFromParticles(self.ds, self.df)
        self.time_data['time'], self.time_data['SFH'] = sfhFromParticles(self.ds, self.df, times=self.time_data['time'])
        self.time_data['time'], self.time_data['SNII_snr'] = snr(self.ds, self.df ,times=self.time_data['time'], sn_type ='II')
        self.time_data['time'], self.time_data['SNIa_snr'] = snr(self.ds, self.df ,times=self.time_data['time'], sn_type ='Ia')

        self.time_data['time'] = 0.5 * (self.time_data['time'][1:] + self.time_data['time'][:-1]) # bin centers

        return

    def instantaneous_SFR(self):
        """
        Returns instantaneous SFR as computed from the global SFR from
        particle formation times. Uses linear interpolation on sampled
        SFR to get exact instantaneous SFR
        """

        if hasattr(self, 'time_data'):
            if not 'SFR' in self.time_data:
                self.compute_time_evolution()
        else:
            self.compute_time_evolution()

        return np.interp(self.ds.current_time.convert_to_units(UNITS['Time'].units),
                         0.5*(self.time_data['time'][:-1]+self.time_data['time'][1]), self.time_data['SFR']) * yt.units.Msun / self.time_data['time'].unit_quantity

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
        raise NotImplementedError

        return


    def compute_half_light_radius(self):
        """
        Compute the radial Luminosity profile as determined
        from stellar evolution model used to described stars, then
        calculate the half light radius.
        """

        if hasattr(self, 'particle_profiles'):
            if not 'luminosity' or ('io','particle_model_luminosity') in self.particle_profiles:
                out = self.particle_profile([('io','particle_model_luminosity')], pt = 11)
                centers = out[1]; prof = out[2]
                self.particle_profiles[('io','particle_model_luminosity')] = prof[('io','particle_model_luminosity')]

        else:
            self.particle_profiles = {}
            self.particle_profiles['r'], centers, prof = self.particle_profile([('io','particle_model_luminosity')], pt = 11)
            self.particle_profiles[('io','particle_model_luminosity')] = prof[('io','particle_model_luminosity')]

        cum_luminosity = np.cumsum(self.particle_profiles[('io','particle_model_luminosity')])

        frac_luminosity = cum_luminosity / cum_luminosity[-1]

        if hasattr(frac_luminosity, 'value'):
            frac_luminosity = frac_luminosity.value

        x      = centers.value

        func   = lambda xval : np.interp(xval, x, frac_luminosity) - 0.5

        r_half = optimize.brentq(func, x[0], x[-1])

        self.particle_meta_data['half_light_radius'] = r_half * centers.unit_quantity

        return r_half * centers.unit_quantity


    def get_star_model_properties(self):
        """
        Go to stellar model and compute stellar model properties of the stars
        as well as radiation properties of the stars. Used to compute more fun things
        like total luminosity and half light radius
        """

        self.star_model_properties = {}


        return

    def particle_profile(self, fields, mode = 'sphere', xtype = 'radial',
                         accumulate=True, weight_field = None, pt=None):
        """
        Constructs a radial profile of the corresponding field. xtype = 'radial' for
        mode 'sphere' ONLY. For mode = 'disk', xtype = 'z' or xtype = 'radial'.
        """

        if (not weight_field is None) and accumulate:
            raise ValueError("Cannot have weight field and accumulation True")

        if not isinstance(fields, Iterable):
            fields = [fields]

        if mode == 'sphere':
            xbins = self.rbins_sphere
            x     = self.sphere['particle_position_spherical_radius'].convert_to_units(xbins.units)
            data  = self.sphere

        elif mode == 'disk':
            if xtype == 'radial':
                xbins = self.rbins_disk
                x     = self.disk['particle_position_cylindrical_radius'].convert_to_units(xbins.units)
            else:
                xbins = self.zbins_disk
                x     = np.abs(self.disk['particle_position_cylindrical_height']).convert_to_units(xbins.units)

            data = self.disk

        if pt == None:
            particle_filter = [True] * np.size(x)
        else:
            particle_filter = data['particle_type'] == pt

        profiles = {}
        for field in fields:
            profiles[field] = np.zeros(np.size(xbins) - 1)

        for field in fields:
            for i in np.arange(np.size(xbins)-1):
                x_filter   = (x < xbins[i]) * (x >= xbins[i-1])
                filter     = x_filter * particle_filter
                field_data = data[field][filter]

                if accumulate:
                    profiles[field][i] = np.sum( field_data )
                elif weight_field is None:
                    profiles[field][i] = np.average( field_data )
                else:
                    weights = data[weight_field][filter]
                    profiles[field][i] = np.average( field_data, weights = weights)


        centers = 0.5 * (xbins[1:] + xbins[:-1])
        return xbins, centers, profiles

    def construct_regions(self, disk_kwargs = None, sphere_kwargs = None,
                                halo_sphere_kwargs = None):
        """
        Defines the pre-defined (or user modified) geometric regions to 
        perform analysis on. These are the galaxy disk, a sphere around the
        galaxy, and a halo sphere (out to virial radius).
        """

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
                            'height' : 400.0 * yt.units.pc,
                            'center' : self.ds.domain_center,
                            'dr'     : 25.0 * yt.units.pc,
                            'dz'     : 50.0 * yt.units.pc }

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
