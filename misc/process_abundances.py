"""
    process_abundances

    Author: A. Emerick

    Notes: script and functions to process stellar abundane output from a
           simulation where chemical tags of stars are written to stdout. This
           loops through all files that could contain this info in the given
           directory (or files that match chosen string) and collates the
           abundance information into a single file. Stars are identified
           their particle IDs, which can be readily tagged with data output
           info from Enzo.

           In order to work, this MUST have a data output associated with the
           simulation, in order to read in species names for what is followed
           in the run.
"""
import numpy as np
import glob
import yt
import subprocess
from galaxy_analysis.utilities.utilities import species_from_fields

_base_col_names = ["grid_id","pid","ptype","x","y","z","particle_mass","birth_mass","metallicity"]

_base_dtypes = [int, int, float, float, float, float]

def filter_data(data):
    """
    Filter data, returning only the non-repeating values
    """


    return

def _read_sf_data(directory = '.'):

    bash_commands = ["grep --no-filename -e '^P(' ./enzo_job*.out  > SA_temp.dat",
                     'sed -e "s/P(//g" -i SA_temp.dat',
                     'sed -e "s/individual_star_maker//g" -i SA_temp.dat',
                     'sed -e "s/ \[add\]\://g" -i SA_temp.dat',
                     'sed -i "/CFRF/d" SA_temp.dat',
                     'sed -e "s/new star particles//g" -i SA_temp.dat',
                     "awk " + "'{$1=" + '""; print $0}' + "' SA_temp.dat > sf.dat",
                     "rm SA_temp.dat"]

    for bc in bash_commands:
        subprocess.call(bc, shell=True)

    data = np.genfromtxt('sf.dat')

    print(np.sum(data[:]))

    return

def read_all_data(directory = '.'):

    # First lets do some hacky BS
    bash_commands = ["grep --no-filename -e '^StellarAbundances P' ./enzo_job*.out  > SA_temp.dat",
                     'sed -e "s/StellarAbundances P(//g" -i SA_temp.dat',
                     "awk " + "'{$1=" + '""; print $0}' + "' SA_temp.dat > StellarAbundances.dat",
                     "rm SA_temp.dat"]

    for bc in bash_commands:
        subprocess.call(bc, shell=True)

    try:
        # need simulation parameter file to get column names
        ds_names = glob.glob(directory + "/DD????/DD????")
        ds = yt.load(ds_names[-1])
        species = ['H','He'] + species_from_fields(ds.field_list)
    except:
        species = ['H','He','C','N','O','K','Fe','Zn','Sr','Ba','AGB','PopIII','SNIa','SNII']

    _col_names = _base_col_names + species
    _dtypes = _base_dtypes + [float] * len(species)
    _ndtypes = [(x,y) for x,y in zip(_col_names,_dtypes)]

    data = np.genfromtxt(directory +'/StellarAbundances.dat', dtype = _ndtypes, invalid_raise=False)



    return data


if __name__ == "__main__":

    _read_sf_data()

    read_all_data()
