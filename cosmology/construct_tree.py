import yt
#yt.enable_parallelism()
from treefarm import TreeFarm
import numpy as np


def build_merger_tree():

    ds = yt.load("rockstar_halos/halos_RD0023/halos_0.0.bin")
    i_max = np.argmax(ds.r["halos", "particle_mass"])
    my_id = ds.r['particle_identifier'][i_max]

    ts      = yt.DatasetSeries("RD????/RD????")
    my_tree = TreeFarm(ts)
    my_tree.trace_ancestors("halos", my_id, filename = "test_halo/")

    return



if __name__=="__main__":

    build_merger_tree()
