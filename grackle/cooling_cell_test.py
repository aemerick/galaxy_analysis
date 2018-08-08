import matplotlib.pyplot as plt
import os
import numpy as np
import yt

from pygrackle import \
    FluidContainer, \
    chemistry_data, \
    evolve_constant_density

from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    sec_per_Myr, \
    cm_per_mpc

import sys

tiny_number = 1e-20

class NoStdStreams(object):
    def __init__(self,stdout = None, stderr = None):
        self.devnull = open(os.devnull,'w')
        self._stdout = stdout or self.devnull or sys.stdout
        self._stderr = stderr or self.devnull or sys.stderr

    def __enter__(self):
        self.old_stdout, self.old_stderr = sys.stdout, sys.stderr
        self.old_stdout.flush(); self.old_stderr.flush()
        sys.stdout, sys.stderr = self._stdout, self._stderr

    def __exit__(self, exc_type, exc_value, traceback):
        self._stdout.flush(); self._stderr.flush()
        sys.stdout = self.old_stdout
        sys.stderr = self.old_stderr
        self.devnull.close()


def cooling_cell(density = 12.2,
                 initial_temperature = 2.0E4,
                 final_time = 200.0,
                 metal_fraction = 4.0E-4,
                 make_plot = False,
                 save_output = False,
                 outname = None, save_H2_fraction = False,
                 return_result = False,
                 verbose = False, H2_converge = None,
                 *args, **kwargs):


    current_redshift = 0.

    # Set solver parameters
    my_chemistry = chemistry_data()
    my_chemistry.use_grackle = 1
    my_chemistry.with_radiative_cooling = 1
    my_chemistry.primordial_chemistry = 2
    my_chemistry.metal_cooling = 1
    my_chemistry.UVbackground = 1
    my_chemistry.self_shielding_method = 3
    my_chemistry.H2_self_shielding = 2
    my_chemistry.h2_on_dust = 1
    my_chemistry.three_body_rate = 4
    grackle_dir = "/home/aemerick/code/grackle-emerick/"
    my_chemistry.grackle_data_file = os.sep.join( #['/home/aemerick/code/grackle-emerick/input/CloudyData_UVB=HM2012.h5'])
        [grackle_dir, "input","CloudyData_UVB=HM2012_shielded.h5"])


    # set the factors
    my_chemistry.LW_factor = kwargs.get("LW_factor", 1.0)
    my_chemistry.k27_factor = kwargs.get("k27_factor", 1.0)
    #if 'LW_factor' in kwargs.keys():
    #    my_chemistry.LW_factor = kwargs['LW_factor']
    #else:
    #    my_chemistry.LW_factor = 1.0

    #if 'k27_factor' in kwargs.keys():
    #    my_chemistry.k27_factor = kwargs['k27_factor']
    #else:
    #    my_chemistry.k27_factor = 1.0

    # Set units
    my_chemistry.comoving_coordinates = 0 # proper units
    my_chemistry.a_units = 1.0
    my_chemistry.a_value = 1. / (1. + current_redshift) / \
        my_chemistry.a_units
    my_chemistry.density_units = mass_hydrogen_cgs # rho = 1.0 is 1.67e-24 g
    my_chemistry.length_units = cm_per_mpc         # 1 Mpc in cm
    my_chemistry.time_units = sec_per_Myr          # 1 Myr in s
    my_chemistry.velocity_units = my_chemistry.a_units * \
        (my_chemistry.length_units / my_chemistry.a_value) / \
        my_chemistry.time_units


    rval = my_chemistry.initialize()

    fc = FluidContainer(my_chemistry, 1)
    fc["density"][:] = density
    if my_chemistry.primordial_chemistry > 0:
        fc["HI"][:] = 0.76 * fc["density"]
        fc["HII"][:] = tiny_number * fc["density"]
        fc["HeI"][:] = (1.0 - 0.76) * fc["density"]
        fc["HeII"][:] = tiny_number * fc["density"]
        fc["HeIII"][:] = tiny_number * fc["density"]
    if my_chemistry.primordial_chemistry > 1:
        fc["H2I"][:] = tiny_number * fc["density"]
        fc["H2II"][:] = tiny_number * fc["density"]
        fc["HM"][:] = tiny_number * fc["density"]
        fc["de"][:] = tiny_number * fc["density"]
    if my_chemistry.primordial_chemistry > 2:
        fc["DI"][:] = 2.0 * 3.4e-5 * fc["density"]
        fc["DII"][:] = tiny_number * fc["density"]
        fc["HDI"][:] = tiny_number * fc["density"]
    if my_chemistry.metal_cooling == 1:
        fc["metal"][:] = metal_fraction * fc["density"] * \
          my_chemistry.SolarMetalFractionByMass

    fc["x-velocity"][:] = 0.0
    fc["y-velocity"][:] = 0.0
    fc["z-velocity"][:] = 0.0
    fc['H2_self_shielding_length'][:] = 1.8E-6

    fc["energy"][:] = initial_temperature / \
        fc.chemistry_data.temperature_units
    fc.calculate_temperature()
    fc["energy"][:] *= initial_temperature / fc["temperature"]

    # timestepping safety factor
    safety_factor = 0.001

    # let gas cool at constant density

    #if verbose:
    print "Beginning Run"
    data = evolve_constant_density(
            fc, final_time=final_time, H2_converge = H2_converge,
            safety_factor=safety_factor, verbose = verbose)
    #else:
    #    print "Beginning Run"

    #    with NoStdStreams():
    #        data = evolve_constant_density(
    #            fc, final_time=final_time, H2_converge = 1.0E-6,
    #            safety_factor=safety_factor)
    #    print "Ending Run"


    if make_plot:
        p1, = plt.loglog(data["time"].to("Myr"), data["temperature"],
                            color="black", label="T")
        plt.xlabel("Time [Myr]")
        plt.ylabel("T [K]")

        data["mu"] = data["temperature"] / \
            (data["energy"] * (my_chemistry.Gamma - 1.) *
             fc.chemistry_data.temperature_units)
        plt.twinx()
        p2, = plt.semilogx(data["time"].to("Myr"), data["mu"],
                              color="red", label="$\\mu$")
        plt.ylabel("$\\mu$")
        plt.legend([p1,p2],["T","$\\mu$"], fancybox=True,
                      loc="center left")
        plt.savefig("cooling_cell.png")


    # save data arrays as a yt dataset
    if outname is None:
        outname = 'cooling_cell_%.2f_%.2f'%(my_chemistry.k27_factor,
                                           my_chemistry.LW_factor)

    if save_output:

        yt.save_as_dataset({}, outname + '.h5', data)

    H2_fraction = (data['H2I'] + data['H2II']) / data['density']

    if save_H2_fraction:
        #np.savetxt(outname + ".dat", [data['time'], H2_fraction])

        f = open("all_runs_d_%.2f.dat"%(density),"a")
#        f.write("# k27 LW f_H2 T time\n")
        f.write("%8.8E %8.8E %8.8E %8.8E %8.8E \n"%(my_chemistry.k27_factor,
                                              my_chemistry.LW_factor,
                                              H2_fraction[-1], data['temperature'][-1],
                                              data['time'][-1] ))
        f.close()

    if return_result:
        return data
    else:
        return


def cooling_cell_grid(k27_factors = None, LW_factors = None,
                      fmin = 1.0, fmax = 100.0, npoints = 100):

    if k27_factors is None:
        k27_factors = np.logspace(np.log10(fmin),
                                  np.log10(fmax), npoints)
    if LW_factors is None:
        LW_factors  = 1.0 * k27_factors

    call_cell = lambda x, y : cooling_cell(k27_factor = x,
                                           LW_factor  = y, save_H2_fraction = True)


    for i,k27 in enumerate(k27_factors):

        print (i)*np.size(LW_factors)

        temp_cell = lambda y : call_cell(k27,y)
        map(temp_cell, LW_factors)



    return


if __name__ == "__main__":

    # test this
    #cooling_cell(k27_factor = 0.99, LW_factor = 0.99,
    #             save_output = False, save_H2_fraction=True)

    cooling_cell_grid(npoints = 32)
