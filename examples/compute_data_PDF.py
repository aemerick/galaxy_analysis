import numpy as np
from galaxy_analysis.analysis import gas_abundance as GA
from galaxy_analysis.analysis import Galaxy

from galaxy_analysis.static_data import ISM # some pre-defined masks

#
# Test code to load a data set and then compute the gas
# phase mass fraction / abundance distribution for an arbitrarily
# defined mask
#
gal = Galaxy('DD0401')        # load galaxy


# mask is just a boolean mask over the cells. This hould
# be easy to make arbitrarily complex things
data_source = gal.disk
#mask        = np.ones(np.shape(gal.disk['x']))

abundances = {}
for n in [50, 100, 150]:
    mask        = (data_source['number_density'] > n)
    abundances[n] = GA.compute_abundance_stats(gal.ds, data_source, mask = mask,
                               fraction_fields = ['O_Fraction','N_Fraction'])


