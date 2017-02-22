import galaxy_analysis as ga
import glob
import os
import numpy as np

if True:

    imin = 100
    imax = 101
    di   = 1

    for i in np.arange(imin,imax+di,di):

        if not os.path.isfile('DD%0004i/DD%0004i'%(i,i)):
            print "File path not found : DD%0004i/DD%0004i"%(i,i)
            continue

        galaxy = ga.Galaxy('DD%0004i'%(i))
        galaxy.compute_everything()
        galaxy.save()
        del(galaxy)


