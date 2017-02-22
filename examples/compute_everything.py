import galaxy_analysis as ga
import glob
import os
import numpy as np
import sys

from galaxy_analysis.misc import process_boundary_flux as pbf

from joblib import Parallel, delayed
import multiprocessing

def _parallel_loop(i):
    if not os.path.isfile('DD%0004i/DD%0004i'%(i,i)):
        print "File path not found : DD%0004i/DD%0004i"%(i,i)
        return

    galaxy = ga.Galaxy('DD%0004i'%(i))
    galaxy.compute_everything()
    galaxy.save()
    del(galaxy)

    return




def run_analysis(imin, imax, di, n_jobs = None):

    if n_jobs is None:
        n_jobs = multiprocessing.cpu_count()

    pbf(filename = 'boundary_mass_flux.dat')

    print "Beginning analysis on %i files on %i processors"%((imax+di-imin)/(1.0*di),n_jobs)

    Parallel(n_jobs=n_jobs)(delayed(_parallel_loop)(i) for i in np.arange(imin,imax+di,di))

    return


if __name__ == "__main__":

    n_jobs = None

    if len(sys.argv) == 5:
        imin, imax, di, n_jobs = [ int(x) for x in sys.argv[1:] ]

    elif len(sys.argv) == 4:
        imin, imax, di = [ int(x) for x in sys.argv[1:] ]

    elif len(sys.argv) == 3:
        imin, imax = int(sys.argv[1]), int(sys.argv[2])
        di = 1

    elif len(sys.argv) == 1:
        files = glob.glob('DD????')
        files = np.sort(files)

        ivals = np.array( [ int(re.search(r'\d+', x)) for x in files ] )
        imin  = np.min(ivals)
        imax  = np.max(ivals)
        di    = np.min(ivals[1:] - ivals[:-1])

        
    run_analysis(imin, imax, di, n_jobs = n_jobs)
