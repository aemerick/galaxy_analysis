


from galaxy_analysis.cosmology import find_halos
from galaxy_analysis.cosmology import filter_tree
import glob, os, sys

#
# Load a bunch of paths to various modules
#
from galaxy_analysis  import __path__ as ga_path
from rockstar_galaxies import __path__ as rockstar_path # requires rockstar_galaxies folder to be soft linked to rockstar-galaxies
ga_path       = [x for x in ga_path][0]
rockstar_path = [x for x in rockstar_path][0]


#
# An attempt to better automate all the halo finding junk
#

def run_rockstar(mpicall = "ibrun", nproc = 3, overwrite=False):
    """
    Run rockstar halo finder using my wrapper. This searches for the enzo
    parameter file in a directory (*.enzo) to run a yt-based call of rockstar
    halo finder using a simulation series in yt. Will default to restarting
    the rockstar run if rockstar_galaxies/ folder is present. MUST be run
    in simulation directory
    """

    #
    # search for enzo parameter file
    #
    parfiles = np.sort(glob.glob("*.enzo"))
    if len(parfiles) == 0:
        print("No parameter files found")
        raise RuntimeError
    elif len(partifles) > 1:
        print("More than one parameter file found, using last in sorted list")
        print(parfiles)
        parfile = parfiles[-1]
    else:
        parfile = parfiles[0]
        print("Using parameter file: ", parfile)

    # assume True
    if (not overwrite ) and (os.path.isdir("rockstar_halos/")):
        restart = True
    else:
        restart = False

    if mpicall == "mpirun":
        procstring = "-np"
    elif mpicall == "ibrun":
        procstring = '-n'

    bash_command = mpicall + " " + procstring + " " + ga_path[0] +\
                   "/cosmology/find_halos.py " + parfile + " " + str(restart)


    print("Running rockstar with command : ", bash_command)
    call(bash_command, shell=True)
    print("Finished with run_rockstar")
    return

def run_consistent_trees():
    """
    Wrapper to run calls on consistent trees.
    """

    rockstar_file = './rockstar_halos/rockstar.cfg'
    if not os.path.isfile(rockstar_file):
        print(rockstar_file, " cannot be found. Aborting run_consistent_trees")
        raise RuntimeError

    bash_command = "perl " + rockstar_path + '/scripts/gen_merger_cfg.pl' +\
                   rockstar_file

    print("Running consistent trees set up with ", bash_command)
    call(bash_command, shell=True)

    full_workdir_path = os.getcwd()

    bash_command = "cd " + consistent_trees_path
    call(bash_command, shell=True)

    bash_command = "make"
    call(bash_command, shell = True)

    bash_command = "perl do_merger_tree.pl " + full_workdir_path +\
                   "/rockstar_halos/outputs/merger_tree.cfg"

    print("Generating merger tree with " + bash_command)
    call(bash_command, shell=True)

    bash_command = "cd " + full_workdir_path
    call(bash_command, shell=True)

    print("Finishing run_consistent_trees")

    return


def run_filter_halos():
    """
    Run my particular halo filter
    """

    tree_file = "./rockstar_halos/trees/tree_0_0_0.dat"

    if (not (os.path.isfile(tree_file))):
        print("Cannot find tree file at ", tree_file)
        raise RuntimeError

    filter_tree.halo_filter(tree_file)


    return

if __name__ == "__main__":

    run_rockstar()
    run_consistent_trees()
    run_filter_halos()

    #
    # and after this, select the halo you want, make the changes
    # to the _mrp file and run
    #
    # mpirun -n <number-of-cores> python enzo-mrp.music.py <mrp_file> <level>
    #
    # then run enzo again and repeat if needed 
