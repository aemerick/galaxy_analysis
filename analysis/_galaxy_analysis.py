from __future__ import division

import yt
yt.funcs.mylog.setLevel(40)

#
# Nov 2017: Units bug in Enzo domain mass flux that needs to be corrected
#           by root grid dx^2. However, to avoid confusion, for now,
#           just make this a global post-process correction. Leaving as bug
#           in Enzo.
#
APPLY_CORRECTION_TO_BOUNDARY_MASS_FLUX_VALUES = True

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
from galaxy_analysis.utilities import utilities as util
#from ..utilities import utilities
from ..static_data import LABELS,\
                        FIELD_UNITS,\
                        IMAGE_COLORBAR_LIMITS,\
                        PLOT_LIMITS,\
                        UNITS,\
                        ISM, CUT_REGION, ISM_FILTER

from galaxy_analysis import particle_analysis as pa

from ..particle_analysis import particle_types as pt
#from ..particle_analysis import IMF

# need to have better API
# from ..particle_analysis.sfrFromParticles import sfrFromParticles
# from ..particle_analysis.sfhFromParticles import sfhFromParticles
# from ..particle_analysis.sn_rate          import snr

from ..yt_fields import field_generators as fg

from ..misc import process_boundary_flux
from ..misc import dm_halo as dmprof

_hdf5_compression = 'lzf'

_all_fields = ['density', 'temperature', 'cell_mass']

__all__ = ['Galaxy']


class Galaxy(object):

    def __init__(self, dsname, wdir = './'):
        """
        Galaxy object to run full analysis on individual data dumps
        in yt. Defines uniform set of galaxy disk an halo regions, species
        and abundance fields given field list, and functions for running 
        a variety of analysis.
        """

        self.wdir    = wdir
        self.dsname = dsname

        # load, generate fields, reload
        with util.nooutput(): # silence yt for now
            self.ds     = yt.load(self.wdir + '/' + self.dsname + '/' + self.dsname)

        # define fields if they have not yet been defined
        if not fg.FIELDS_DEFINED:
            dfiles = glob.glob(self.wdir + '/' + 'DD????/DD????')
            dfiles = np.sort(dfiles)

            if len(dfiles) == 1:
                dstemp = yt.load(dfiles[0])
            else:
                for i in np.arange(-1, -np.size(dfiles), -1):
                    try:
                        dstemp = yt.load(dfiles[i])
                    except:
                        print i
                        continue

                    break

            fg.generate_derived_fields(dstemp)
            fg.FIELDS_DEFINED = True

        # load data set and data
        self.ds     = fg.load_and_define(self.wdir + '/' + self.dsname + '/' + self.dsname)
        self.current_time = self.ds.current_time.convert_to_units(UNITS['Time'].units).value
        self.df     = self.ds.all_data()

        # check for particles
        self._has_particles = False
        if ('io','particle_position_x') in self.ds.field_list:
            self._has_particles = True

        self.hdf5_filename   = self.wdir + '/' + self.dsname + '_galaxy_data.h5'

        self._compute_virial_parameters()

        self._set_data_region_properties()
        self.species_list = util.species_from_fields(self.ds.field_list)

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
        self.observables        = {}

        self.construct_regions()

        self.total_quantities = {}
        self._total_quantities_calculated = False

#        self._has_boundary_mass_file, self.boundary_mass_flux =\
#                    process_boundary_flux(data = None, wdir = self.wdir)
#
#        self._has_boundary_mass_file = os.path.isfile(self.wdir + "filtered_boundary_mass_flux.dat")
#
#        if self._has_boundary_mass_file:
#            self.boundary_mass_flux = np.genfromtxt(self.wdir + "filtered_boundary_mass_flux.dat", names = True)

        self._load_boundary_mass_flux() # load directly from parameter file

        if os.path.isfile( self.hdf5_filename ):
            self.load()

        return

    def _compute_virial_parameters(self):
        """
        Assuming standard cosmology and z = 0, compute virial mass
        and radius from input parameters. Burkert profile only
        """
        r_s = self.ds.parameters['DiskGravityDarkMatterR'] * yt.units.Mpc
        r_s = r_s.convert_to_units('kpc')

        rho_o = self.ds.parameters['DiskGravityDarkMatterDensity'] * yt.units.g / yt.units.cm**3

        self.M_vir, self.R_vir = \
                   dmprof.burkert_solve_virial(r_s, rho_o)

        return

    def load(self, filename = None, nocopy = False):
        """
        Load the specified hdf5 file.
        """
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
        """
        Save all of the generated data to file using deepdish.
        This constructs a nested dictionary which is then outputted to file.
        """
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
        self._output_data_dict['observables']        = self.observables

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
        self.observables        = _verify_and_add(self.observables, 'observables')

        if 'R_vir' in self.meta_data.keys():
            self.R_vir = self.meta_data['R_vir']
        if 'M_vir' in self.meta_data.keys():
            self.M_vir = self.meta_data['M_vir']

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

    def calculate_velocity_distributions(self):
        """
        Compute dM / dv vs. v_r distributions (mass moving at certain velocity)
        outside the disk of the galaxy.
        """

        # use the halo sphere. First compute this for ALL gas within the halo,
        # then do it again in the same bins as the outflow profiles
        self.gas_profiles['velocity'] = {}
        self.gas_profiles['velocity']['halo']   = {} # for the entire halo
        self.gas_profiles['velocity']['binned'] = {} # binned by radius - empty for now

        vbins = np.arange(-100.0, 1000.5, 5.0) * yt.units.km / yt.units.s
        self.gas_profiles['velocity']['vbins'] = vbins

        local_cr = self.cut_region
        not_disk = self.halo_sphere.cut_region(local_cr['not_disk'])

        # do this for the entire halo, everything except the disk
        self.gas_profiles['velocity']['halo'] = {}
        for k in CUT_REGION.keys():  # loop through all cut regions
            data  = self.halo_sphere.cut_region(CUT_REGION[k] + "&" + local_cr['not_disk'] ) # phase + geometric cut
            test  = self.halo_sphere.cut_region(CUT_REGION[k])
            test2 = self.halo_sphere.cut_region(local_cr['not_disk'])

            # loop over velocity bins
            hist = np.zeros(np.size(vbins) - 1)
            for i in np.arange(np.size(vbins) -1):
                v       = data['velocity_spherical_radius'].convert_to_units('km/s')
                hist[i] = np.sum( data['cell_mass'][(v < vbins[i+1]) * (v >= vbins[i])].convert_to_units('Msun') )

            self.gas_profiles['velocity']['halo'][k] = hist

        # compute mass weighted statistics about both the inflow and outflow
        v = not_disk['velocity_spherical_radius'].convert_to_units('km/s')
        m = not_disk['cell_mass'].convert_to_units('Msun').value
        v_out = v[v>0] ; m_out = m[v>0]
        v_in  = v[v<0] ; m_in  = m[v<0]
        if np.size(v_out) <= 3:
            v_out = np.zeros(10); m_in = np.ones(10)
        if np.size(v_in) <= 3:
            v_in = np.zeros(10); m_in = np.ones(10)
        self.gas_profiles['velocity']['halo']['outflow_stats'] =\
                    util.compute_weighted_stats(v_out, m_out, return_dict = True)
        self.gas_profiles['velocity']['halo']['inflow_stats'] =\
                    util.compute_weighted_stats(v_in, m_in, return_dict = True)

        return

    def calculate_radiation_profiles(self, fields = None, mode = 'disk'):

        if not ('radiation' in self.gas_profiles.keys()):
            self.gas_profiles['radiation'] = {}
        if not (mode in self.gas_profiles['radiation'].keys()):
            self.gas_profiles['radiation'][mode] = {}

        # want to make radial profiles of G_o, Q_0, Q_1, F_LW, F_PE, Gamma_PE
        fields = [('gas','G_o'), ('gas','FUV_flux'), ('gas','LW_flux'),
                    ('gas','Pe_heating_rate_masked'), ('gas','Q0_flux'), ('gas','Q1_flux')]

        midplane = self.ds.disk(self.disk.center, self.disk.field_parameters['normal'], self.disk.radius,
                                 2.5 * np.min(self.disk['dx'].convert_to_units('pc')))

        n_bins = np.ceil( self.disk.radius.convert_to_units('pc').value / 5.0) # 5 pc bins
        profiles = yt.create_profile(midplane, 'radius', fields, n_bins = n_bins,
                                        weight_field = 'cell_volume', logs = {'radius':False})
        self.gas_profiles['radiation'][mode]['xbins'] = profiles.x_bins.convert_to_units('pc').value

        for f in fields:
            self.gas_profiles['radiation'][mode][f[1]] = profiles[f]

        return self.gas_profiles['radiation'][mode]['xbins'], self.gas_profiles['radiation'][mode]

    def calculate_dMdt_profile(self, fields = None, mode = 'sphere', n_cell = 4,
                               outflow = True, phase = None, *args, **kwargs):
        """
        Returns the mass inflow or outflow rate as a function of radius. This can be used
        to compute mass loading factor by dividing by the current SFR.

        outflow == True (default) computes outflow rate (requiring that v > 0.0). False
        computes inflow rate (requiring that v < 0.0). If using a disk, v is v_z, but 
        if using a sphere, v is v_r

        if no fields are provided, computes total mass flow rate and rate for
        all species.
        """

        if fields is None:
            fields = self._accumulation_fields

        # set up bin values and data region
        xbins, xdata, data = self._get_bins_and_data(mode = mode)

        # 

        # get velocity corresponding to outflow / inflow
        if mode == 'sphere' or mode == 'halo_sphere':
            velname = 'velocity_spherical_radius'
        else:
            velname = 'velocity_cylindrical_z'
        vel = data[velname].convert_to_units(UNITS['Velocity'].units)

        profile = {}
        profile['mass_profile'] = {} # Bin up the total mass in outflowing material

        #
        # Following typical definitions, construct bins to be centered at 
        # 0.25, 0.5, 0.75, 1.0, and 1.25 R_vir, with a width of 0.1 R_vir
        #

        center = np.arange(0.25, 1.30, 0.25) # in units of R_vir
        center = np.array([0.1, 0.2, 0.25, 0.5, 0.75, 1.0, 1.25]) # in units of R_Vir
        dL     = 0.1                           # in units of R_vir

        for field in fields:
            profile[field] = np.zeros(np.size(center))
            profile['mass_profile'][field] = np.zeros(np.size(center))

        # convert from r_vir to kpc
        center = (center * self.R_vir).convert_to_units('kpc')
        dL     = (dL     * self.R_vir).convert_to_units('kpc')

#        dx = n_cell * np.min(data['dx'].convert_to_units(xdata.units))

        if outflow: # compute outflow
            v_filter = vel > 0.0
        else:       # compute inflow
            v_filter = vel < 0.0

        if phase is None:
            phase_filter = (vel == vel)
        else:
            phase_filter = ISM_FILTER[phase](data)

        for i in np.arange(np.size(center)):
            # define the shell
            x_filter = ( xdata >= (center[i] - 0.5*dL)) * ( xdata < (center[i] + 0.5*dL))

            # filter out the data
            filter = x_filter * v_filter * phase_filter
            for field in fields:
                M        = data[field].convert_to_units(UNITS['Mass'].units)
                Mdot     = np.sum( M[filter] * vel[filter] ) / dL
                Mdot     = Mdot.convert_to_units('Msun/yr')

                profile[field][i] = Mdot
                profile['mass_profile'][field][i] = np.sum(M[filter])

        #
        # save profiles
        #
	prof_type = 'outflow'
        if not outflow:
            prof_type = 'inflow'

        if not prof_type in self.gas_profiles.keys():
            self.gas_profiles[prof_type] = {}

        if not mode in self.gas_profiles[prof_type].keys():
            self.gas_profiles[prof_type][mode] = {}

        if not (phase is None):
            if not (phase in self.gas_profiles[prof_type][mode].keys()):
                self.gas_profiles[prof_type][mode][phase] = {}

            self.gas_profiles[prof_type][mode][phase].update(profile)
        else:

            self.gas_profiles[prof_type][mode].update( profile )
            self.gas_profiles[prof_type][mode]['centers'] = center
            self.gas_profiles[prof_type][mode]['centers_rvir'] = (center.convert_to_units('kpc') / self.R_vir.convert_to_units('kpc')).value
            self.gas_profiles[prof_type][mode]['dL']      = dL
            self.gas_profiles[prof_type][mode]['dL_rvir'] = (dL.convert_to_units('kpc')/self.R_vir.convert_to_units('kpc')).value

        return xbins, center, profile

    def _get_bins_and_data(self, mode = None, axis='z'):

        if mode == 'sphere' or mode == 'halo_sphere':
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

    def compute_observables(self, young_star = 10.0 * yt.units.Myr):
        """
        Does some computing of observational diagnostics that would be
        important to check. When possible, does this using various definitions of
        how one might compute the observable.
        """
        #
        # First, compute the surface densities one would want to compare
        #        against schmidt law.
        #
        # Following Roychowdhury et. al. 2017, 2014, etc., define a "SF region"
        #   rather than the whole galaxy. For our sake, require this to be at least 100 pc
        age  = (self.ds.current_time - self.df['creation_time'])
        r_sf =  self.df['particle_position_cylindrical_radius'][ age <= young_star]
        if np.size(r_sf) <= 5:
            r_sf = self.df['particle_position_cylindrical_radius'][ age <= young_star*2 ]

        if np.size(r_sf) == 0: # no stars formed !!!
            keys = ['r_sf','A_sf','A','SD_HI_sf','SD_gas_sf','SD_gas_sf_obs','SD_HI',
             'SD_gas','SD_gas_obs','SD_H2_sf','SD_H2','SD_SFR','SD_SFR_sf','SD_stellar','SD_stellar_sf']
            for k in keys:
                self.observables[k] = 0.0
        else:
            r_sf = np.max(r_sf.convert_to_units('pc').value) * yt.units.pc
            if r_sf < 100.0 * yt.units.pc:
                r_sf = 100*yt.units.pc

            sf_disk = self.ds.disk([0.5,0.5,0.5],[0,0,1], r_sf, self.disk_region['height'])
            A       = (np.pi * r_sf * r_sf).convert_to_units("pc**2")
            A_disk  = (np.pi * self.disk.radius * self.disk.radius).convert_to_units('pc**2')

            self.observables['r_sf']   = r_sf * 1.0
            self.observables['A_sf']   = A * 1.0
            self.observables['A']      = A_disk * 1.0
            self.observables['SD_HI_sf' ] = np.sum( sf_disk['H_p0_mass'].convert_to_units('Msun') ) / A
            self.observables['SD_gas_sf'] = np.sum( sf_disk['cell_mass'].convert_to_units('Msun') ) / A
            self.observables['SD_gas_sf_obs'] = self.observables['SD_HI_sf'] * 1.34
            self.observables['SD_HI'] = np.sum(self.disk['H_p0_mass'].convert_to_units('Msun')) / A_disk
            self.observables['SD_gas'] = np.sum(self.disk['cell_mass'].convert_to_units('Msun')) / A_disk
            self.observables['SD_gas_obs'] = self.observables['SD_HI'] * 1.34
            self.observables['SD_H2_sf'] = np.sum( (sf_disk['H2_p0_mass'] + sf_disk['H2_p1_mass']).convert_to_units('Msun'))/A
            self.observables['SD_H2']    = np.sum( (self.disk['H2_p0_mass'] + self.disk['H2_p1_mass']).convert_to_units('Msun'))/A_disk

            self.observables['SD_SFR']    = self.meta_data['SFR'] / A_disk.convert_to_units("kpc**2")
            self.observables['SD_SFR_sf'] = self.meta_data['SFR'] / A.convert_to_units("kpc**2")

            self.observables['SD_stellar'] = self.meta_data['M_star'] / A_disk.convert_to_units('kpc**2')
            self.observables['SD_stellar_sf'] = self.meta_data['M_star'] / A.convert_to_units('kpc**2')


        return self.observables

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

#    def compute_radiation_profiles(self, fields = None, mode = 'disk', axis = 'r'):
#        
#        if fields is None:
#            fields = self._radiation_fields()
#
#        xbins, xdata, data = self._get__bins_and_data(mode, axis)
#
#        prof = {}
#
#        return xbins, centers, prof

    def compute_all_meta_data(self):
        """
        Wrapper to compute all meta data
        """

        self.compute_meta_data()
        self.compute_all_particle_profiles()
        self.compute_particle_meta_data()
        self.compute_gas_meta_data()
        self.compute_time_evolution()
        self.compute_observables()

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
        WD_stars = pt.white_dwarfs(self.ds, self.df)
        SNII     = pt.core_collapse(self.ds, self.df)
        SNIa     = pt.snIa(self.ds, self.df)

        particle_mass = (self.df['birth_mass'].value *
                         yt.units.Msun).convert_to_units(UNITS['Mass'].units)
        m = (self.df['particle_mass'].convert_to_units(UNITS['Mass'].units))

        self.particle_meta_data['t_first_star'] = np.min( self.df['creation_time'].convert_to_units(UNITS['Time'].units))

        if 'GalaxySimulationInitialStellarDist' in self.ds.parameters:
            if self.ds.parameters['GalaxySimulationInitialStellarDist'] > 0:
                ct = self.df['creation_time'].convert_to_units('Myr')
                if np.size(ct[ct>1.0]) > 0:
                    self.particle_meta_data['t_first_star'] = np.min( ct[ct > 1.0].convert_to_units(UNITS['Time'].units))


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

        self.particle_meta_data['metallicity_stars'] = util.compute_stats(self.df['metallicity_fraction'])
        # compute theoretical total IMF and 'observed' IMF of just MS stars
        self.particle_meta_data['IMF_obs']        = pa.compute_IMF(self.ds, self.df, mode='mass',       bins=25)
        self.particle_meta_data['IMF_birth_mass'] = pa.compute_IMF(self.ds, self.df, mode='birth_mass', bins=25)

        self.compute_half_light_radius()

        # in principle store abundance information here, but might just leave this as separate

        return

    def compute_gas_meta_data(self):

        self.gas_meta_data['masses'] = self.compute_gas_sequestering()

        # save individual specise masses from non-equillibrium chemistry
        for e in ['HI','HII','HeI','HeII','HeIII','H2I','H2II','HM']:
            self.meta_data['M_' + e]  = np.sum(self.disk[e + '_Density'] *\
                                        self.disk['cell_volume']).convert_to_units(UNITS['Mass'].units)

        # total masses for species where ionization statest are tracked
        self.meta_data['M_H_total'] = self.meta_data['M_HI'] + self.meta_data['M_HII'] + self.meta_data['M_H2I'] +\
                                      self.meta_data['M_H2II'] + self.meta_data['M_HM']
        self.meta_data['M_He_total'] = self.meta_data['M_HeI'] + self.meta_data['M_HeII'] + self.meta_data['M_HeIII']
        self.meta_data['M_H2_total'] = self.meta_data['M_H2I'] + self.meta_data['M_H2II']

        self.meta_data['M_gas'] = np.sum((self.disk['cell_mass']).convert_to_units(UNITS['Mass'].units))
        self.meta_data['Z_avg'] = np.sum( (self.disk[('gas','Metal_Density')]*\
                                           self.disk['cell_volume']).convert_to_units(UNITS['Mass'].units))/\
                                                  self.meta_data['M_gas']

        #
        # compute total mass in disk for all species
        #
        for e in self.species_list:
            fname = e + '_Mass'
            self.meta_data['M_' + e] = np.sum(self.disk[fname]).convert_to_units(UNITS["Mass"].units)


        #
        # compute the volume occupied by each phase in the disk
        #
        cut_region_names = ['Molecular', 'CNM', 'WNM', 'WIM', 'HIM']
        self.gas_meta_data['volume_fractions'] = {}
        self.gas_meta_data['mass_fractions']   = {}

        total_volume = np.sum(self.disk['cell_volume'].convert_to_units('cm**(3)'))
        for crtype in cut_region_names:
            v = self.disk.cut_region(ISM[crtype])['cell_volume'].convert_to_units('cm**(3)')
            self.gas_meta_data['volume_fractions'][crtype] = np.sum(v) / total_volume

            self.gas_meta_data['mass_fractions'][crtype] = self.gas_meta_data['masses'][crtype]['Total']/\
                                                               self.gas_meta_data['masses']['Disk']['Total']

        self.gas_meta_data['volume_fractions']['Total'] = total_volume

        return


    def compute_gas_profiles(self):

        junk = self.calculate_radiation_profiles()
        junk = self.calculate_dMdt_profile()               # mass outflow rate
        junk = self.calculate_dMdt_profile(outflow=False)  # mass inflow rate
        for phase in ['CNM','WNM','WIM','HIM']:
            junk = self.calculate_dMdt_profile(phase = phase)
            junk = self.calculate_dMdt_profile(phase = phase, outflow = False)
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
        fields = {'H':'H_total_mass','He':'He_total_mass','Total':'cell_mass','Metals':'metal_mass', 'H2' : 'H2_mass',
                  'HI':'H_p0_mass', 'HII': 'H_p1_mass'}

        def _sum_tracked_metals(d): # sum tracked metals species only
            return np.sum([d[k] for k in d.keys() if (not any([k in ['Metals','Total','H','H2','He','HI','HeI','HeII','HeIII','H2I','H2II','HII']]))])

        # do this for the disk ISM regions
        for crtype in cut_region_names:
            mdict[crtype] = {}
            for s in fields:
                mdict[crtype][s] = np.sum(self.disk.cut_region( ISM[crtype])[fields[s]]).convert_to_units('Msun')

            for s in self.species_list:
                mdict[crtype][s] = np.sum(self.disk.cut_region( ISM[crtype])[('gas',s + '_Mass')]).convert_to_units('Msun')
            mdict[crtype]['Total Tracked Metals'] = _sum_tracked_metals(mdict[crtype])

        # now do this for the whole disk
        mdict['Disk'] = {}
        for s in fields:
            mdict['Disk'][s] = np.sum(self.disk[fields[s]]).convert_to_units('Msun')
        for s in self.species_list:
            mdict['Disk'][s] = np.sum(self.disk[('gas',s + '_Mass')]).convert_to_units('Msun')
        mdict['Disk']['Total Tracked Metals'] = _sum_tracked_metals(mdict['Disk'])

        # now do this for the halo
        mdict['Halo'] = {}
        for s in fields:
            #print s
            #print fields[s]
            #print self.halo_sphere[fields[s]]
            mdict['Halo'][s] = np.sum(self.halo_sphere[fields[s]]).convert_to_units('Msun')
        for s in self.species_list:
            mdict['Halo'][s] = np.sum(self.halo_sphere[('gas',s + '_Mass')]).convert_to_units('Msun')
        mdict['Halo']['Total Tracked Metals'] = _sum_tracked_metals(mdict['Halo'])

        # now do this for full box
        mdict['FullBox'] = {}
        for s in fields:
            mdict['FullBox'][s] = np.sum(self.df[fields[s]]).convert_to_units('Msun')
        for s in self.species_list:
            mdict['FullBox'][s] = np.sum(self.df[('gas', s + '_Mass')]).convert_to_units('Msun')
        mdict['FullBox']['Total Tracked Metals'] = _sum_tracked_metals(mdict['FullBox'])

        # now we need to do some subtraction of the fields
        mdict['OutsideHalo'] = {}
        for s in fields.keys() + self.species_list + ['Total Tracked Metals']:
            mdict['OutsideHalo'][s] = mdict['FullBox'][s] - mdict['Halo'][s]
            mdict['Halo'][s]        = mdict['Halo'][s]    - mdict['Disk'][s]

        # now we compute the gravitationally bound gas IF potential is present
        if 'PotentialField' in self.ds.field_list or ('enzo','GravPotential') in self.ds.field_list:
            mdict['GravBound'] = {}
            for s in fields:
                mdict['GravBound'][s] = np.sum( self.ds.cut_region(self.df, "obj[('gas','gravitationally_bound')] > 0" )[fields[s]]).convert_to_units('Msun')
            for s in self.species_list:
                mdict['GravBound'][s] = np.sum(self.ds.cut_region(self.df, "obj[('gas','gravitationally_bound')] > 0")[('gas', s + '_Mass')]).convert_to_units('Msun')
            mdict['GravBound']['Total Tracked Metals'] = _sum_tracked_metals(mdict['GravBound'])

        # and finally add up the mass in stars
        mdict['stars'] = {}
        for s in ['H','He'] + self.species_list:
            if self._has_particles:
                mdict['stars'][s] = np.sum( (self.df['birth_mass'].value *
                                             self.df['particle_' + s + '_fraction'])[self.df['particle_type'] == 11])
            else:
                mdict['stars'][s] = 0.0
        mdict['stars']['Total Tracked Metals'] = _sum_tracked_metals(mdict['stars'])

        mdict['OutsideBox'] = {}
#        if self._has_boundary_mass_file:
#            index = np.argmin( np.abs(self.current_time - self.boundary_mass_flux['Time']) )
#            diff = np.abs(self.boundary_mass_flux['Time'][index] - self.current_time)
#            if diff > 1.0: 
#                print "WARNING: Nearest boundary mass flux data point is > 1 Myr from current simulation time"
#                print "T_now = %5.5E T_file = %5.5E diff = %5.5E"%(self.current_time, self.boundary_mass_flux['Time'][index], diff)
#
#
#            if np.size(index) > 1:
#                index = index[0]
#
        for s in self.species_list:
            mdict['OutsideBox'][s] = self.boundary_mass_flux[s + '_Density'] # [index]

        _fields = ['HI_Density','HII_Density','H2I_Density','H2II_Density','HM_Density']
#        mdict['OutsideBox']['H'] = np.sum([ self.boundary_mass_flux[field][index] for field in _fields])
        mdict['OutsideBox']['H'] = np.sum([ self.boundary_mass_flux[field] for field in _fields])
        _fields = ['HeI_Density','HeII_Density','HeIII_Density']
#        mdict['OutsideBox']['He'] = np.sum([ self.boundary_mass_flux[field][index] for field in _fields])
        mdict['OutsideBox']['He']     = np.sum([self.boundary_mass_flux[field] for field in _fields])
        mdict['OutsideBox']['Metals'] = self.boundary_mass_flux['Metal_Density'] # [index]
#        else:
#            for s in ['H','He','Metals'] + self.species_list:
#                mdict['OutsideBox'][s] = 0.0
        mdict['OutsideBox']['Total Tracked Metals'] = _sum_tracked_metals(mdict['OutsideBox'])

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

        #
        # Make a set of SFR and SNR evolutions with default bin spacing (10 Myr)
        #
        # x is temp dummy variable
        self.time_data['time'], self.time_data['SFR'] = pa.sfrFromParticles(self.ds, self.df)
        x, self.time_data['SFH'] = pa.sfhFromParticles(self.ds, self.df, times=self.time_data['time'])
        x, self.time_data['SNII_snr'] = pa.snr(self.ds, self.df, times=x, sn_type ='II')
        x, self.time_data['SNIa_snr'] = pa.snr(self.ds, self.df, times=x, sn_type ='Ia')
        x, self.time_data['AGB_rate'] = pa.snr(self.ds, self.df, times=x, sn_type = 'AGB')

        self.time_data['time'] = 0.5 * (x[1:] + x[:-1]) # bin centers

        self.meta_data['SFR']  = self.time_data['SFR'][-1]

        #
        # Make a set of SFR and SNR evolutions with longer bin spacing (100 Myr)
        #
        self.time_data['time_100'], self.time_data['SFR_100'] = pa.sfrFromParticles(self.ds, self.df, times = 100.0 * yt.units.Myr)
        x, self.time_data['SFH_100'] = pa.sfhFromParticles(self.ds, self.df, times=self.time_data['time_100'])
        x, self.time_data['SNII_snr_100'] = pa.snr(self.ds, self.df ,times=x, sn_type ='II')
        x, self.time_data['SNIa_snr_100'] = pa.snr(self.ds, self.df ,times=x, sn_type ='Ia')

        self.time_data['time_100'] = 0.5 * (x[1:] + x[:-1]) # bin centers

        self.meta_data['SFR_100']  = self.time_data['SFR_100'][-1]

        #
        # Make a set of SFR and SNR evolutions with very short bin spacing (1 Myr)
        #
        self.time_data['time_1'], self.time_data['SFR_1'] = pa.sfrFromParticles(self.ds, self.df, times = 1.0 * yt.units.Myr)
        x, self.time_data['SFH_1'] = pa.sfhFromParticles(self.ds, self.df, times=self.time_data['time_1'])
	x, self.time_data['SNII_snr_1'] = pa.snr(self.ds, self.df ,times=x, sn_type ='II')
        x, self.time_data['SNIa_snr_1'] = pa.snr(self.ds, self.df ,times=x, sn_type ='Ia')

        self.time_data['time_1'] = 0.5 * (x[1:] + x[:-1]) # bin centers

        self.meta_data['SFR_1']  = self.time_data['SFR_1'][-1]



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
                         0.5*(self.time_data['time'][:-1]+self.time_data['time'][1:]), self.time_data['SFR']) * yt.units.Msun / self.time_data['time'].unit_quantity

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
        self.calculate_velocity_distributions()


        return


    def compute_half_light_radius(self):
        """
        Compute the radial Luminosity profile as determined
        from stellar evolution model used to described stars, then
        calculate the half light radius.
        """

        make_profile = False
        if hasattr(self, 'particle_profiles'):
            if not util.nested_haskey(self.particle_profiles, ['disk','radial','sum',('io','particle_model_luminosity')]):
                make_profile = True
        else:
            self.particle_profiles = {}
            make_profile = True

        if make_profile:
            junk = self.compute_particle_profile( [ ('io','particle_model_luminosity'), ], xtype = 'radial',
                                                  accumulate = True, mode = 'disk', pt = 11)

        xbins          = self.particle_profiles['disk']['radial']['xbins']
        centers        = 0.5 * (xbins[1:] + xbins[:-1])
        cum_luminosity = np.cumsum( self.particle_profiles['disk']['radial']['sum'][('io','particle_model_luminosity')] )

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

    def compute_all_particle_profiles(self):
        """
        Computes all particle profiles we want, including abundance profiles
        for the particles. 
        """

        junk = self.compute_particle_profile( [('io','particle_mass'),], mode = 'disk', pt = 11)
        junk = self.compute_particle_profile( [('io','particle_mass'),], mode = 'disk', xtype = 'z', pt = 11)

        junk = self.compute_particle_profile( [('io','particle_age'),
                                              ('io','metallicity_fraction')],
                                              mode = 'disk', accumulate = False, pt = 11)

        abundance_fields = ['Fe_over_H', 'C_over_Fe', 'O_over_Fe', 'Mg_over_Fe',
                            'O_over_Fe']

        for i, a in enumerate( abundance_fields):
            abundance_fields[i] = ('io','particle_' + a)

        junk = self.compute_particle_profile(abundance_fields, mode = 'disk', pt = 11,
                                             xtype = 'radial', accumulate = False)

        junk = self.compute_particle_profile(abundance_fields, mode = 'disk', pt = 11,
                                             xtype = 'z', accumulate = False)

        return

    def compute_particle_profile(self, fields, mode = 'disk', xtype = 'radial',
                         accumulate=True, weight_field = None, pt=None):
        """
        Constructs a radial profile of the corresponding field. xtype = 'radial' for
        mode 'sphere' ONLY. For mode = 'disk', xtype = 'z' or xtype = 'radial'. Here, disk
        is assumed to be stellar disk always
        """

        if (not weight_field is None) and accumulate:
            raise ValueError("Cannot have weight field and accumulation True")

        if not isinstance(fields, Iterable):
            fields = [fields]

        if mode == 'sphere':
            xbins = self.rbins_sphere
            x     = self.sphere[('io','particle_position_spherical_radius')].convert_to_units(xbins.units)
            data  = self.sphere

        elif mode == 'disk':
            if xtype == 'radial':
                xbins = self.rbins_stellar_disk
                x     = self.stellar_disk[('io','particle_position_cylindrical_radius')].convert_to_units(xbins.units)
            else:
                xbins = self.zbins_stellar_disk
                x     = np.abs(self.stellar_disk[('io','particle_position_cylindrical_z')]).convert_to_units(xbins.units)

            data = self.stellar_disk

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

                if field in UNITS:
                    field_data = field_data.convert_to_units(FIELD_UNITS[field].units)

                if accumulate:
                    profiles[field][i] = np.sum( field_data )
                elif weight_field is None:
                    profiles[field][i] = np.average( field_data )
                else:
                    weights = data[weight_field][filter]
                    profiles[field][i] = np.average( field_data, weights = weights)

        #
        # save profiles
        #
        prof_type = mode
        if not prof_type in self.particle_profiles.keys():
            self.particle_profiles[prof_type] = {}

        if not xtype in self.particle_profiles[prof_type].keys():
            self.particle_profiles[prof_type][xtype] = {}

        weight = weight_field
        if weight is None and accumulate:
            weight = "sum"
        elif weight is None and not accumulate:
            weight = "average"

        if not weight in self.particle_profiles[prof_type][xtype].keys():
            self.particle_profiles[prof_type][xtype][weight] = {}

        self.particle_profiles[prof_type][xtype][weight].update( profiles )
        self.particle_profiles[prof_type][xtype]['xbins'] = xbins

        centers = 0.5 * (xbins[1:] + xbins[:-1])
        return xbins, centers, profiles

    def construct_regions(self, disk_kwargs = None, large_disk_kwargs = None, stellar_disk_kwargs = None,
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

        if stellar_disk_kwargs is None:
           stellar_disk_kwargs = {}

           for n in ['normal', 'radius','height','center']:
               stellar_disk_kwargs[n] = self.stellar_disk_region[n]

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

        self.stellar_disk = self.ds.disk(**stellar_disk_kwargs)

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

        self.stellar_disk_region = {'normal' : np.array([0.0, 0.0, 1.0]),
                                    'radius' : 1.0 * yt.units.kpc,
                                    'height' : 400.0 * yt.units.pc,
                                    'center' : self.ds.domain_center,
                                    'dr'     : 10.0 * yt.units.pc,
                                    'dz'     : 10.0 * yt.units.pc}


        self.disk_region = {'normal' : np.array([0.0, 0.0, 1.0]),
                            'radius' : 600.0 * yt.units.pc,
                            'height' : 200.0 * 2 * yt.units.pc, # 200 pc above and below
                            'center' : self.ds.domain_center,
                            'dr'     : 25.0 * yt.units.pc,
                            'dz'     : 50.0 * yt.units.pc }

        # HACK HACK HACK:
        if self.ds.parameters['DiskGravityStellarDiskMass'] > 1.0E7:
            self.disk_region['height'] = 2.0 * yt.units.kpc
            self.disk_region['radius'] = 1.5 * yt.units.kpc


        self.large_disk_region = {'normal' : np.array([0,0,1]),
                                  'radius' : 2.0 * yt.units.kpc,
                                  'height' : 2.0 * yt.units.kpc,
                                  'center' : self.ds.domain_center,
                                  'dr'     : 25.0*yt.units.pc,
                                  'dz'     : 50.0*yt.units.pc}

        self.spherical_region = {'center' : self.ds.domain_center,
                                 'radius' : 2.0 * yt.units.kpc,
                                 'dr'     : 50.0 * yt.units.pc   }

        # 
        # need to not hard code virial radius
        #

        self.halo_spherical_region = {'center' :    self.ds.domain_center,
                                      'radius' :    self.R_vir}
        self.halo_spherical_region['dr'] = self.halo_spherical_region['radius'] * 0.05

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
                                     ('gas','H_p0_mass'), ('gas','H2_mass'),
                                     ('gas','H_p1_mass')]

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
        """
        Previously, this was done by reading an external file containing
        cumulative mass loss from grid at each root grid timestep. Instead,
        this is stored directly in the paramater file. Its gross, but
        just do this.
        """

        # BoundarMassFluxFieldNumbers stores field number which is connected
        # to field name via data label
        num_flux = len( [x for x in self.ds.parameters.keys() if 'BoundaryMassFluxFieldNumbers' in x])

        if hasattr(self, 'boundary_mass_flux'):
            if len(self.boundary_mass_flux.keys()) == num_flux:
                return # don't need to re-make this

        self.boundary_mass_flux = {}

        conv = 1.0
        if APPLY_CORRECTION_TO_BOUNDARY_MASS_FLUX_VALUES: # obnoxius on purpose, see note at top
            conv = np.max( self.df['dx'].convert_to_units('code_length').value )
            conv = 1.0 / (conv**2) # correct by dividing by dx^2

        for i in np.arange(num_flux):
            field = self.ds.parameters["DataLabel[%i]"%(self.ds.parameters['BoundaryMassFluxFieldNumbers[%i]'%(i)])]
            self.boundary_mass_flux[field] = self.ds.parameters['BoundaryMassFluxContainer[%i]'%(i)] * conv


        # easy !
        return

    @property
    def cut_region(self):
        """
        Dictionary of cut region strings to be used in yt's cut region functionality.
        At the moment these are geometric (or kinematic) exclusion regions
       (e.x. 'not_disk') to be used in conjuction with the already defined regions
       (disk, sphere, etc.) in order to easily select part of the volume in between
        two regions (there is no support for this directly in some sort of object).
        Otherwise, the disk, sphere and such objects can be used for continious and
        un-interrupted region selection.

        Other cut regions are defined in static_data that are immutable, like
        the ISM phases; these depend on properties of the galaxy itself that need to be
        defined on the fly.

        Though not in place yet, this could be extended to add in a LOT more convenience
        selection functions.
        """

        all_cr = {}

        #
        # set up cut regions to select bands above and below disk corresponding
        # to regions where one might want to examine select outflow / inflow properties
        #
        center = np.array([0.1, 0.2, 0.25, 0.5, 0.75, 1.0, 1.25]) # in units of R_Vir
        dL     = 0.1                                              # in units of R_vir

        center = (center * self.R_vir).convert_to_units('kpc')
        dL     = (dL * self.R_vir).convert_to_units('kpc')

        rvir_ranges = [None] * len(center)
        for i in np.arange(np.size(center)):
            lower_lim = (center[i] - 0.5 * dL).value
            upper_lim = (center[i] + 0.5 * dL).value

            rvir_ranges[i] = ["(obj['magnitude_cylindrical_z'].in_units('kpc') > %.4E)"%(lower_lim), "&",
                              "(obj['magnitude_cylindrical_z'].in_units('kpc') > %.4E)"%(upper_lim)]

            all_cr['z_rvir_%i'] = rvir_ranges[i]

        all_cr['z_rvir_info'] = {'centers' : center, 'dL' : dL}

        # cut region string for outside the disk region
        disk_r =       self.disk.radius.convert_to_units('pc').value
        disk_z = 0.5 * self.disk.height.convert_to_units('pc').value

        # add this cut region as an attribute of the galaxy
        not_disk = ["(obj['cylindrical_r'].in_units('pc') > %.4E)"%(disk_r), "&", # outside certain radius
                    "(obj['magnitude_cylindrical_z'].in_units('pc') > %.4E)"%(disk_z)]
        all_cr['not_disk'] = ' '.join(not_disk)

        return all_cr

    @property
    def rbins_stellar_disk(self):
        rmin = 0.0
        rmax = self.stellar_disk.radius
        dr   = self.stellar_disk_region['dr']
        rmax = rmax.convert_to_units(dr.units)

        return np.arange(rmin, rmax + dr, dr) * dr.unit_quantity

    @property
    def zbins_stellar_disk(self):
        zmin = 0.0
        zmax = self.stellar_disk.height
        dz   = self.stellar_disk_region['dz']
        zmax = zmax.convert_to_units(dz.units)

        return np.arange(zmin, zmax + dz, dz ) * dz.unit_quantity


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


