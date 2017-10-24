from galaxy_analysis.plot.plot_styles import *

import numpy as np
import deepdish as dd
import matplotlib.pyplot as plt
import h5py

def get_index_number(n, z = None, z_index = None):
    """
    Get correct index number for grackle data files given 
    number density and redshift
    """

    all_n = np.arange(-10, 4.0 , 0.5)

    all_z = np.array([0.0000e+00, 1.2202e-01, 2.5893e-01, 4.1254e-01, 5.8489e-01, 7.7828e-01, 9.9526e-01, 1.2387e+00,
             1.5119e+00, 1.8184e+00, 2.1623e+00, 2.5481e+00, 2.9811e+00, 3.4668e+00, 4.0119e+00, 4.6234e+00,
             5.3096e+00, 6.0795e+00, 6.9433e+00, 7.9125e+00, 9.0000e+00, 1.0220e+01, 1.1589e+01, 1.3125e+01,
             1.4849e+01])

    
    index1 = np.argmin( np.abs(all_n - n))
    
    if z_index is None and z is None:
        index2 = 0 # redshift zero
    elif z_index is None:
        index2 = np.argmin( np.abs(all_z - z))
    elif z is None:
        if z_index < 0:
            z_index = np.size(all_z) + z_index
            
        index2 = z_index
    else:
        print "Cannot provide both a redshift value and redshift bin value"
        raise ValueError
        
    
    run_num = (index1)*len(all_z) + (index2 + 1)
    
    return index1, index2, run_num
#
# Load data
#
thin   = h5py.File('./../grackle/CloudyData_UVB=HM2012.h5')
shield = h5py.File('./../grackle/CloudyData_UVB=HM2012_shielded.h5')
    
thin_metal_c = thin['CoolingRates']['Metals']['Cooling']
thin_metal_h = thin['CoolingRates']['Metals']['Heating']

shield_metal_c = shield['CoolingRates']['Metals']['Cooling']
shield_primordial_c = shield['CoolingRates']['Primordial']['Cooling']

shield_metal_h = shield['CoolingRates']['Metals']['Heating']
shield_primordial_h = shield['CoolingRates']['Primordial']['Heating']


def plot_cooling_model_comparison():
    T = np.logspace(1,9,161)
    T = np.log10(T)

    plot_n = [-2, -1, 0, 1, 2, 3] # densities to plot
    Z = 0.1                       # metallicity (in solar)
    fig, ax = plt.subplots(2,3, sharex=True, sharey=True)

    axis_tuple = [(0,0),(0,1),(0,2),(1,0),(1,1),(1,2)]
    z = 0.0 # redshift
    z_index = None # index of redshift (found for you)

    def _plot_data(ax, x, y, color, label):
        if np.size(y[y<0]) > 0:

            ax.plot(x[y<0], np.log10(np.abs(y[y<0])), lw = 3, color = color, ls = '--')
            ax.plot(x[y>0], np.log10(np.abs(y[y>0])), lw = 3, color = color, ls = '-', label = label)
    
            # fill in the gaps between the two regions if there are any - don't attempt if there are multiple
            test = np.where(y<0)[0]
            if all( test[1:] == (test[:-1] + 1)):
                heat_max_index = np.argmax(y[y<0])
                cool_min_index = np.argmin(y[y>0]) + np.size(y[y<0]) + 1
                ax.plot(x[heat_max_index:cool_min_index], np.log10(np.abs(y[heat_max_index:cool_min_index])), 
                    lw = 3, color = color, ls = '-')  
            else:
                print 'cannot plot the separate lines'
        else:
            ax.plot(x, np.log10(np.abs(y)), lw = 3, color = color, ls = '-', label = label)
        return
    
    
    for i in [0,1,2,3,4,5]:
    
        n = plot_n[i]
    
        index1, index2, run_num = get_index_number(n, z = z, z_index = z_index)
    
        #
        # Plot Forbes et al. type cooling
        # this should be shielded primordial + optically thin metal
        #
        total_cool = thin_metal_c[index1][index2] * Z + shield_primordial_c[index1][index2]
        total_heat = thin_metal_h[index1][index2] * Z + shield_primordial_h[index1][index2]
        net = total_cool - total_heat
    
        i = axis_tuple[i]
        _plot_data(ax[i], T, net, 'orange', 'Incorrect')
    
        # plot consistent model using shielding data
        total_cool = shield_metal_c[index1][index2] * Z + shield_primordial_c[index1][index2]
        total_heat = shield_metal_h[index1][index2] * Z + shield_primordial_h[index1][index2]
        net = total_cool - total_heat    
    
        _plot_data(ax[i], T, net, 'black', 'Our Model')
        
        x_ann = 1.2
        y_ann = -21.5
        ax[i].annotate(r'n = %.3f'%(10**n), xy=(x_ann,y_ann),xytext=(x_ann,y_ann))

        ax[i].xaxis.set_ticks(np.arange(1,9.1,1))
    
        ax[i].set_xlim(0.75,9.25)
        ax[i].set_ylim(-30.25, -20.75)

    xlabel = r'log(T [K])'
    for index in [(1,0),(1,1),(1,2)]:
        ax[index].set_xlabel(xlabel)
        ax[index].set_xticklabels([1,2,3,4,5,6,7,8,9])

    ylabel = r'log($\Lambda$ [(erg cm$^{3}$ s$^{-1}$])'
    ax[(0,0)].set_ylabel(ylabel)
    ax[(1,0)].set_ylabel(ylabel)
    ax[(0,0)].legend(loc = 'best')
    plt.minorticks_on()
    fig.set_size_inches(12,8)
    plt.tight_layout()
    plt.subplots_adjust(wspace=0,hspace=0)
    fig.savefig('cooling_model_comparison')


if __name__ == '__main__':
    plot_cooling_model_comparison()
