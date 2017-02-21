"""
    process_boundary_flux

    Author: A. Emerick

    Notes: script and functions to process output 
           from domain boundary mass flux computation
"""

import numpy as np
import os

__all__ = ['process_boundary_flux']

def process_boundary_flux(data = None, filename = None, wdir = ''):
    """
    Given a set of boundary mass flux data, loop through
    and stitch this together so that there is no double
    counting of timesteps.
    """

    if data is None:

        if filename is None:
            filename = wdir + 'boundary_mass_flux.dat'
        
        if not os.path.isfile(filename):
            print 'boundary mass flux file not found at ' + filename
            return False, 0

        data = np.genfromtxt(filename)
        with open(filename, 'r') as f:
            header = f.readline()

    data = data[data[:,1].argsort()]

    unique_time = np.unique(data[:,1])

    filtered_data = [None]*np.size(unique_time)
    
    for i, t in enumerate(unique_time):

        selection = (data[:,1] == t)
        filtered_data[i] = np.mean(data[selection], axis= 0)
    
    filtered_data = np.array(filtered_data)

    for i in np.arange(2, np.size(filtered_data[0])):
        filtered_data[:,i] = np.cumsum(filtered_data[:,i])

    outfile = 'filtered_boundary_mass_flux.dat'

    np.savetxt(outfile, filtered_data, fmt = ('%0.6E'), header = header)

    data = np.genfromtxt(outfile, names = True)

    return True, data

    
    
