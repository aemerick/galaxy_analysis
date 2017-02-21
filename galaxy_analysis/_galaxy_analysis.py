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
from ..utilities import utilities
from ..static_data import LABELS,\
                        FIELD_UNITS,\
                        IMAGE_COLORBAR_LIMITS,\
                        PLOT_LIMITS,\
                        UNITS,\
                        ISM

from galaxy_analysis import particle_analysis as pa

from ..particle_analysis import particle_types as pt
#from ..particle_analysis import IMF

# need to have better API
# from ..particle_analysis.sfrFromParticles import sfrFromParticles
# from ..particle_analysis.sfhFromParticles import sfhFromParticles
# from ..particle_analysis.sn_rate          import snr

from ..yt_fields import field_generators as fg

from ..misc import process_boundary_flux

_hdf5_compression = 'lzf'

_all_fields = ['density', 'temperature', 'cell_mass']

__all__ = ['Galaxy']

class Galaxy(object):

    def __init__(self, dsname, wdir = './'):
        """
        Object for individual data dump analysis
        """

        self.wdir    = wdir
        self.dsname = dsname

        # load, generate fields, reload
        with utilities.nooutput(): # silence yt for now
            self.ds     = yt.load(self.wdir + '/' + self.dsname + '/' + self.dsname)

        if not fg.FIELDS_DEFINED:
            dfiles = glob.glob(self.wdir + '/' + 'DD????/DD????')
            dfiles = np.sort(dfiles)
            fg.generate_derived_fields(yt.load(dfiles[-1]))
            fg.FIELDS_DEFINED = True

        self.ds     = yt.load(self.wdir + '/' + self.dsname + '/' + self.dsname)
        self.current_time = self.ds.current_time.convert_to_units(UNITS['Time'].units).value

        self.df     = self.ds.all_data()

        self._has_particles = False
        if ('io','particle_position_x') in self.ds.field_list:
            self._has_particles = True

        self.hdf5_filename   = self.wdir + '/' + self.dsname + '_galaxy_data.h5'

        self._set_data_region_properties()
        self.species_list = utilities.species_from_fields(self.ds.field_list)

        # define some fields that will be used for automatic 
        # computation of profiles when no user supplied field is given
        self._set_accumulation_fields()
        self._set_projection_fields()
        self._set_radiation_fields()

        self.particle_meta_data = {}
        self.gas_meta_data      = {}
        self.meta_data          = {}
        self.gas_profiles       = {}
        self.particle_profiles  = {}
        self.time_data          = {}

        self.construct_regions()

        self.total_quantities = {}
        self._total_quantities_calculated = False

        self._has_boundary_mass_file, self.boundary_mass_flux =\
                    process_boundary_flux(data = None, wdir = self.wdir)

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
        Map the class parameter structure to the output dictionary
        """
        if not hasattr(self, '_output_data_dict'):
            self._output_data_dict = {}

        self._output_data_dict['meta_data']          = self.meta_data
        self._output_data_dict['gas_meta_data']      = self.gas_meta_data
        self._output_data_dict['particle_meta_data'] = self.particle_meta_data
        self._output_data_dict['time_data']          = self.time_data
        self._output_data_dict['gas_profiles']       = self.gas_profiles
        self._output_data_dict['particle_profiles']  = self.particle_profiles

        return

    def _map_output_to_class(self):
        """
        Analyzed data on disk is read in and stored to self._output_data_dict.
        Map this dictionary to the associated class parameters
        """

        def _verify_and_add(x, name):
            if name in self._output_data_dict.keys():
                return self._output_data_dict[name]
            else:
                return x

        self.meta_data          = _verify_and_add(self.meta_data, 'meta_data')
        self.gas_meta_data      = _verify_and_add(self.gas_meta_data, 'gas_meta_data')
        self.particle_meta_data = _verify_and_add(self.particle_meta_data, 'particle_meta_data')
        self.time_data          = _verify_and_add(self.time_data, 'time_data')
        self.gas_profiles       = _verify_and_add(self.gas_profiles, 'gas_profiles')
        self.particle_profiles  = _verify_and_add(self.particle_profiles, 'particle_profiles')

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


    def calculate_dMdt_profile(self, fields = None, mode = 'large_disk', n_cell = 4,
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

        # set up bin values and data region
        xbins, xdata, data = self._get_bins_and_data(mode = mode)

        # get velocity corresponding to outflow / inflow
        if mode == 'sphere':
            velname = 'velocity_spherical_radius'
        else:
            velname = 'velocity_cylindrical_z'
        vel = data[velname].convert_to_units(UNITS['Velocity'].units)

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

    def _get_bins_and_data(self, mode = None, axis='z'):

        if mode == 'sphere':
            xbins  =  self.rbins_halo_sphere
            xdata  =  self.halo_sphere['spherical_radius'].convert_to_units(UNITS["Length"].units)
            data   =  self.halo_sphere

        elif mode == 'disk':
            if axis == 'z':
                xbins  =  self.zbins_disk
                xdata  =  (self.disk['z'] - self.disk.center[2]).convert_to_units(UNITS["Length"].units)


            elif axis == 'r':
                xbins = self.rbins_disk
                xdata = self.disk['cylindrical_r'].convert_to_units(UNITS["Length"].units)

            data  = self.disk

        elif mode == 'large_disk':
            if axis == 'z':
                xbins = self.zbins_large_disk
                xdata = (self.large_disk['z'] - self.large_disk.center[2]).convert_to_units(UNITS["Length"].units)
            elif axis == 'r':
                xbins = self.rbins_large_disk
                xdata = self.large_disk['cylindrical_r'].convert_to_units(UNITS["Length"].units)

            data  = self.large_disk

        else:
            raise ValueError("Must choose disk, large_disk, or sphere for region")

        return xbins, xdata, data

    def calculate_mass_fraction_profile(self, fields = None, mode = 'sphere', axis = 'r', *args, **kwargs):
        """
        Computes fractional radial mass profiles for all species. Can be easily
        used to make cumulative profiles
        """

        #
        # check if profiles exist already and don't recalculate unless ordered to
        # 
        if fields is None:
            fields = self._accumulation_fields

        rbins, centers, profiles = self.calculate_mass_profile(fields = fields, mode = mode, axis = axis)

        compute_total_fields = [x for x in fields if x not in self.total_quantities]
        self.calculate_total_quantities(fields = compute_total_fields)

        for field in fields:
            profiles[field] = profiles[field] / self.total_quantities[field]

        # save fractional profiles here?

        return rbins, centers, profiles

    def calculate_mass_profile(self, fields = None, mode = 'sphere', axis = 'r', *args, **kwargs):
        """
        Compute mass fraction of given species contained within spherical
        bins out to the halo radius. Used to estimate metal retention fractions
        of given species.
        """

        if fields is None:
            fields = self._accumulation_fields

        xbins, xdata, data = self._get_bins_and_data(mode, axis)

        profiles = {}

        for field in fields:
            profiles[field] = np.zeros(np.size(xbins)-1)

        for i in np.arange(np.size(xbins)-1):
            x_filter = (xdata >= xbins[i]) * (xdata < xbins[i+1])

            for field in fields:
                profiles[field][i] = np.sum(\
                      data[field][x_filter].convert_to_units(FIELD_UNITS[field].units))

        centers = 0.5 * (xbins[1:] + xbins[:-1])

        #
        # save and store profiles here?
        #
        prof_type = 'accumulation'
        if not prof_type in self.gas_profiles.keys():
            self.gas_profiles[prof_type] = {}

        if not mode in self.gas_profiles[prof_type].keys():
            self.gas_profiles[prof_type][mode] = {}
        
        self.gas_profiles[prof_type][mode].update( profiles )
        self.gas_profiles[prof_type][mode]['xbins'] = xbins

        return xbins, centers, profiles


    def calculate_surface_density_profile(self, fields = None, data_source = None, rbins = None,
                                          *args, **kwargs):
        """
        Computes a 1D surface density profile in a cylindrical region of the galaxy
        using already set disk selection region.
        """

        if fields is None:
            fields = self._projection_fields

        if not isinstance(fields, Iterable):
            fields = [fields]

        if data_source is None:
            data_source = self.disk

        if rbins is None:
            rbins       = self.rbins_disk

        # project the fields
        proj = self.ds.proj(fields, 'z', data_source = data_source, **kwargs)

        r    = np.sqrt(( (proj['px']-data_source.center[0])**2 + (data_source.center[1] -proj['py'])**2)).convert_to_units('pc')
        A    = (proj['pdx']*proj['pdy']).convert_to_units('pc**2')

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

        #
        # save profiles
        #
        prof_type = 'surface_density'
        mode = 'disk'
        if not prof_type in self.gas_profiles.keys():
            self.gas_profiles[prof_type] = {}

        if not mode in self.gas_profiles[prof_type].keys():
            self.gas_profiles[prof_type][mode] = {}
        
        self.gas_profiles[prof_type][mode].update( profiles )
        self.gas_profiles[prof_type][mode]['xbins'] = rbins


        return rbins, centers, profiles

    def compute_radiation_profiles(self, fields = None, mode = 'disk', axis = 'r'):
        
        if fields is None:
            fields = self._radiation_fields()

        xbins, xdata, data = self._get__bins_and_data(mode, axis)

        prof = {}

        return xbins, centers, prof

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
        m = (self.df['particle_mass'].convert_to_units(UNITS['Mass'].units))

        self.particle_meta_data['total_mass'] = np.sum(m)
        self.particle_meta_data['total_mass_MS'] = np.sum(m[MS_stars])
        self.particle_meta_data['total_birth_mass'] = np.sum(particle_mass)
        self.particle_meta_data['total_birth_mass_MS'] = np.sum(particle_mass[MS_stars])

        self.meta_data['M_star'] = (self.particle_meta_data['total_mass_MS'] * 1.0)

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
        self.particle_meta_data['IMF_obs']        = pa.compute_IMF(self.ds, self.df, mode='mass',       bins=25)
        self.particle_meta_data['IMF_birth_mass'] = pa.compute_IMF(self.ds, self.df, mode='birth_mass', bins=25)

        self.compute_half_light_radius()

        # in principle store abundance information here, but might just leave this as separate

        return

    def compute_gas_meta_data(self):

        self.gas_meta_data['masses'] = self.compute_gas_sequestering()

        self.meta_data['M_HI']  = np.sum(self.disk['HI_Density'] *\
                                  self.disk['cell_volume']).convert_to_units(UNITS['Mass'].units)
        self.meta_data['M_gas'] = np.sum((self.disk['cell_mass']).convert_to_units(UNITS['Mass'].units))
        self.meta_data['Z_avg'] = np.sum( (self.disk[('gas','Metal_Density')]*\
                                           self.disk['cell_volume']).convert_to_units(UNITS['Mass'].units))/\
                                                  self.meta_data['M_gas']

        return



    def compute_gas_profiles(self):

        junk = self.calculate_dMdt_profile()
        junk = self.calculate_surface_density_profile()
        junk = self.calculate_mass_profile(mode = 'disk')
        junk = self.calculate_mass_profile(mode = 'sphere')

        return

    def compute_gas_sequestering(self):
        """
        Computes the sequestering of gas into various categ ories (geometric and phase) for 
        all species as well as the total gas mass. Returns a dictionary of dictionaries
        containing all of this information
        """

        mdict = {}
        
        cut_region_names = ['Molecular', 'CNM', 'WNM', 'WIM', 'HIM']        
        fields = {'H':'H_total_mass','He':'He_total_mass','Total':'cell_mass','Metals':'metal_mass'}

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
        if 'PotentialField' in self.ds.field_list or ('enzo','GravPotential') in self.ds.field_list:
            mdict['GravBound'] = {}
            for s in fields:
                mdict['GravBound'][s] = np.sum( self.ds.cut_region(self.df, "obj[('gas','gravitationally_bound')] > 0" )[fields[s]]).convert_to_units('Msun')
            for s in self.species_list:
                mdict['GravBound'][s] = np.sum(self.ds.cut_region(self.df, "obj[('gas','gravitationally_bound')] > 0")[('gas', s + '_Mass')]).convert_to_units('Msun')

        # and finally add up the mass in stars
        mdict['stars'] = {}
        for s in ['H','He'] + self.species_list:
            if self._has_particles:
                mdict['stars'][s] = np.sum( (self.df['birth_mass'].value *
                                             self.df['particle_' + s + '_fraction'])[self.df['particle_type'] == 11])
            else:
                mdict['stars'][s] = 0.0

        mdict['OutsideBox'] = {}
        if self._has_boundary_mass_file:
            index = np.argmin( np.abs(self.current_time - self.boundary_mass_flux['Time']) )
            diff = np.abs(self.boundary_mass_flux['Time'][index] - self.current_time)
            if diff > 1.0: 
                print "WARNING: Nearest boundary mass flux data point is > 1 Myr from current simulation time"
                print "T_now = %5.5E T_file = %5.5E diff = %5.5E"%(self.current_time, self.boundary_mass_flux['Time'][index], diff)


            if np.size(index) > 1:
                index = index[0]

            for s in self.species_list:
                mdict['OutsideBox'][s] = self.boundary_mass_flux[s + '_Density'][index]

            _fields = ['HI_Density','HII_Density','H2I_Density','H2II_Density','HM_Density']
            mdict['OutsideBox']['H'] = np.sum([ self.boundary_mass_flux[field][index] for field in _fields])
            _fields = ['HeI_Density','HeII_Density','HeIII_Density']
            mdict['OutsideBox']['He'] = np.sum([ self.boundary_mass_flux[field][index] for field in _fields])
            mdict['OutsideBox']['Metals'] = self.boundary_mass_flux['Metal_Density'][index]
        else:
            for s in ['H','He','Metals'] + self.species_list:
                mdict['OutsideBox'][s] = 0.0
                                       

        if self._has_particles:
            mdict['stars']['metals'] = np.sum( (self.df['birth_mass'].value * self.df['metallicity_fraction'])[self.df['particle_type'] == 11])
        else:
            mdict['stars']['metals'] = 0.0

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

        self.time_data['time'], self.time_data['SFR'] = pa.sfrFromParticles(self.ds, self.df)
        self.time_data['time'], self.time_data['SFH'] = pa.sfhFromParticles(self.ds, self.df, times=self.time_data['time'])
        self.time_data['time'], self.time_data['SNII_snr'] = pa.snr(self.ds, self.df ,times=self.time_data['time'], sn_type ='II')
        self.time_data['time'], self.time_data['SNIa_snr'] = pa.snr(self.ds, self.df ,times=self.time_data['time'], sn_type ='Ia')

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
        self.compute_all_meta_data()
        self.compute_gas_profiles()
        self.compute_gas_sequestering()

        return


    def compute_half_light_radius(self):
        """
        Compute the radial Luminosity profile as determined
        from stellar evolution model used to described stars, then
        calculate the half light radius.
        """

        if hasattr(self, 'particle_profiles'):
            if (not 'luminosity' in self.particle_profiles) or (not ('io','particle_model_luminosity') in self.particle_profiles):
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

        if pt is None:
            particle_filter = [True] * np.size(data['particle_type'])
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

    def construct_regions(self, disk_kwargs = None, large_disk_kwargs = None,
                                sphere_kwargs = None, halo_sphere_kwargs = None):
        """
        Defines the pre-defined (or user modified) geometric regions to 
        perform analysis on. These are the galaxy disk, a sphere around the
        galaxy, and a halo sphere (out to virial radius).
        """

        if disk_kwargs is None:
            disk_kwargs = {}

            for n in ['normal','radius','height','center']:
                disk_kwargs[n] = self.disk_region[n]

        if large_disk_kwargs is None:
            large_disk_kwargs = {}

            for n in ['normal','radius','height','center']:
                large_disk_kwargs[n] = self.large_disk_region[n]

        if sphere_kwargs is None:
            sphere_kwargs = {}

            for n in ['radius','center']:
                sphere_kwargs[n] = self.spherical_region[n]

        if halo_sphere_kwargs is None:
            halo_sphere_kwargs = {}

            for n in ['radius','center']:
                halo_sphere_kwargs[n] = self.halo_spherical_region[n]


        self.disk   = self.ds.disk(**disk_kwargs)

        self.large_disk = self.ds.disk(**large_disk_kwargs)

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

        self.large_disk_region = {'normal' : np.array([0,0,1]),
                                  'radius' : 2.0 * yt.units.kpc,
                                  'height' : 2.0 * yt.units.kpc,
                                  'center' : self.ds.domain_center,
                                  'dr'     : 25.0*yt.units.pc,
                                  'dz'     : 50.0*yt.units.pc}

        self.spherical_region = {'center' : self.ds.domain_center,
                                 'radius' : 2.0 * yt.units.kpc,
                                 'dr'     : 25.0 * yt.units.pc   }

        self.halo_spherical_region = {'center' :    self.ds.domain_center,
                                      'radius' : 14.0 * yt.units.kpc,
                                      'dr'     : 50.0 * yt.units.pc}

        return


    #
    # TO DO: Move the "set_xxx_fields" functions to a utilities
    #        or static data function, rather than here... but keep hidden.
    #        this may be a little cleaner - Feb 2017
    #
    #
    def _set_accumulation_fields(self):

        self._accumulation_fields = [('gas','H_total_mass'), ('gas','He_total_mass'),
                                     ('gas','metal_mass'), ('gas','cell_mass'),
                                     ('gas','H_p0_mass')]

        for e in self.species_list:
            self._accumulation_fields += [('gas', e +'_Mass')]

        return

    def _set_radiation_fields(self):

        self._radiation_fields = [('enzo','HI_kph'), ('enzo','HeI_kph'),
                                  ('enzo','OTLW_kdissH2I'), ('gas','Pe_heating_rate_masked'),
                                  ('gas','G_o')]

        return

    def _set_projection_fields(self):

        self._projection_fields = [('enzo','HI_Density'), ('enzo','H2I_Density'),
                                   ('enzo','Density'), ('enzo','HII_Density'),
                                   ('gas','Metal_Density')]

        for e in self.species_list:
            self._projection_fields += [('gas',e + '_Density')]

        return

    def _load_boundary_mass_flux(self):

        self.boundary_mass_flux = {}

        if not os.path.isfile(self.wdir + 'boundary_mass_flux.dat'):
            return

        f = np.genfromtxt(self.wdir + 'boundary_mass_flux.dat')

        # do data filtering --- code this somewhere else and do it with function call

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
        zmax = self.disk.height
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

    @property
    def zbins_large_disk(self):
        zmin = 0.0
        zmax = self.large_disk.height
        dz   = self.large_disk_region['dz']
        zmax = zmax.convert_to_units(dz.units)

        return np.arange(zmin, zmax + dz, dz) * dz.unit_quantity


    @property
    def rbins_large_disk(self):
        rmin = 0.0
        rmax = self.large_disk.radius
        dr   = self.large_disk_region['dr']
        rmax = rmax.convert_to_units(dr.units)

        return np.arange(rmin, rmax + dr, dr) * dr.unit_quantity

