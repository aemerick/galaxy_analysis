import yt
import numpy as np
import glob
import os
import h5py
import deepdish as dd

# parallel
from multiprocessing import Pool
from contextlib import closing
import itertools

nproc = 14
nmax  = 100

def _parallel_loop(dsname):

    groupname = dsname.rsplit('/')[1]
    print "starting computation on ", groupname
    #gal = Galaxy(groupname)
    ds = yt.load(dsname)
    data  = ds.all_data()

    M = data['cell_mass'].to('Msun').value
    n = data['number_density'].to('cm**(-3)')

    dict = {}
    dict[groupname] = np.sum( M[n > 200.0] )

    del(data)
    del(ds)

    return dict

# only do at most 100
dsimin = 119
dsimax = 219

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

f = open("number_density_cut_mass.dat",'w')
f.write("#Time Mass\n")
for i, k in enumerate(np.sort(fulldict.keys())):

    f.write("%2.2E %8.8E\n"%(times[i], fulldict[k]))

f.close()
