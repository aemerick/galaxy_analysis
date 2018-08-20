import yt
import numpy as np
import glob
import os
import h5py
import deepdish as dd

from galaxy_analysis.analysis import Galaxy

# parallel
from multiprocessing import Pool
from contextlib import closing
import itertools
import sys

nproc = 4
nmax  = 200

def _parallel_loop(dsname):

    groupname = dsname.rsplit('/')[1]
    print "starting computation on ", groupname
    gal = Galaxy(groupname)

    M    = gal.df['birth_mass']
    L    = gal.df[('io','particle_model_L_1_3eV')].value
    L_LW = gal.df[('io','particle_model_L_LW')].value
    pt   = gal.df['particle_type']

    L_low  = np.sum(  L[ (M < 8.0) * (pt == 11)])
    L_high = np.sum(  L[ (M > 8.0) * (pt == 11)])

    L_LW_low  = np.sum( L_LW[ (M < 8.0) * (pt == 11)])
    L_LW_high = np.sum( L_LW[ (M > 8.0) * (pt == 11)])

    dict = {}
    dict[groupname] = [L_low, L_high, L_LW_low, L_LW_high]

    del(gal)

    return dict

# only do at most 100
dsimin = int(sys.argv[1])
dsimax = int(sys.argv[2])

ds_list = []
dsi     = []
for i in np.arange(dsimin,dsimax+1,1):
    dsname = 'DD%0004i/DD%0004i'%(i,i)
    if os.path.exists(dsname):
        ds_list = ds_list + [dsname]
        dsi = dsi + [i]

#ds_list = np.sort(glob.glob('./DD????/DD????'))
ds_list = ds_list[: np.min([nmax,np.size(ds_list)])]
times   = np.array(dsi)

fulldict = {}
for sub_list in itertools.izip_longest(*(iter(ds_list),) * nproc):
    sub_list = list(sub_list)
    sub_list = [s for s in sub_list if s is not None]
    reduced_nproc = np.min( [len(sub_list), nproc] )

    pool = Pool(reduced_nproc)
    results = pool.map_async(_parallel_loop, sub_list)

    pool.close()
    pool.join()

    for r in results.get():
        fulldict[ r.keys()[0] ] = r[r.keys()[0]]
    del(results)

f = open("low_eV_luminosity_%0004i_%0004i.dat"%(dsimin,dsimax),'w')
f.write("#Time L_low L_high L_low_LW L_high_LW\n")
for i, k in enumerate(np.sort(fulldict.keys())):

    f.write("%2.2E %8.8E %8.8E %8.8E %8.8E\n"%(times[i], fulldict[k][0], fulldict[k][1], fulldict[k][2], fulldict[k][3]))

f.close()
