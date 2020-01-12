import yt
import numpy as np
from galaxy_analysis.plot.plot_styles import *
import matplotlib.pyplot as plt

from scipy.stats import pearsonr

def correlation(field, ds, dmax = 600.0, N = 100, Ndist = 60,
                           Niter = 1):
    """
    Computes the Pearson correlation statistic for a given field
    as a function of distance. This is a spatial autocorrelation.
    Currently samples points in the midplane of a disk, then computes
    correlation as a function of spherical radius around those points.

    Updates needed for this function:
        1) Ability to request random initial points over arbitrary geometry
        2) Ability to use other distances (cylindrical, etc.)
        3) Use of Niter parameter to sample Niter times around initial
           points and average, rather than just once
           (makes for Niter x N x Ndist total samples.. so.. a lot..)
        4) Parallelization? Primarily if iterations used. Read from
           disk may be issue here though (probably not... just preload?)
        5) If above two, sampling is inefficient constructing the
           pt_x and pt_x_prime arrays of point objects. Could do a faster
           job with a FRB (maybe?) with data loaded in memory already
    """

    domain_width = ds.domain_width
    center = ds.domain_center.to('pc') / domain_width

    # get initial random points
    theta = np.random.rand(N) * np.pi * 2.0
    r_cyl = np.random.rand(N) * dmax**2
    xvec  = (np.sqrt(r_cyl) * np.cos(theta)) * yt.units.pc / domain_width[0] + center[0]
    yvec  = (np.sqrt(r_cyl) * np.cos(theta)) * yt.units.pc / domain_width[1] + center[1]
    zvec  = (np.zeros(N))                    * yt.units.pc / domain_width[2] + center[2]

    # list of yt point objects
    pt_x  = [ ds.r[ [xvec[i].value, yvec[i].value, zvec[i].value] ] for i in np.arange(N)]

    # field values
    S_x   = np.array([p[metal_field][0] for p in pt_x])

    rsample   = np.linspace(0.0, dmax, Ndist)
    corrcoeff = np.zeros(np.size(rsample))
    for i in np.arange(np.size(rsample)):

        # get new random points around each
        theta = np.random.rand(N) * np.pi * 2.0
        r_cyl = np.random.rand(N) * (rsample[i])**2.0

        xprime  = (np.sqrt(r_cyl) * np.cos(theta)) * yt.units.pc / domain_width[0] + xvec
        yprime  = (np.sqrt(r_cyl) * np.cos(theta)) * yt.units.pc / domain_width[1] + yvec
        zprime  = (np.zeros(N))                    * yt.units.pc / domain_width[2] + zvec

        # new list of point objects and new values
        pt_x_prime = [ ds.r[ [xprime[j].value,yprime[j].value,zprime[j].value]] for j in np.arange(N)]
        S_x_prime  = np.array([p[metal_field][0] for p in pt_x_prime])

        # compute coefficient
        corrcoeff[i] = pearsonr(S_x, S_x_prime)[0]

    return corrcoeff
