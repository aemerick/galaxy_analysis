from galaxy_analysis.plot.plot_styles import *
import matplotlib.pyplot as plt
import os
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

def bins_from_centers(x):
    xnew = np.zeros(len(x) + 1)
    dx   = np.zeros(len(x) + 1)
    dx[1:-1]   = x[1:] - x[:-1]
    dx[0] = dx[1]
    dx[-1] = dx[-2]
    xnew[:-1] = x - 0.5*dx[:-1]
    xnew[-1] = x[-1] + 0.5*dx[-1]

    return xnew

def plot_2d_histogram(datafile = 'all_runs_d_12.20.dat'):

    data = np.genfromtxt(datafile) # names = True)


    k27      = data[:,0]
    LW       = data[:,1]
    k27_vals = np.linspace(np.log10(np.min(k27)), np.log10(np.max(k27)),
                           int(np.sqrt(np.size(k27) )))
    k27_vals = bins_from_centers(k27_vals)
    LW_vals  = np.linspace(np.log10(np.min(LW)), np.log10(np.max(LW)),
                           int(np.sqrt(np.size(LW))))
    LW_vals = bins_from_centers(LW_vals)

    k27_mesh, LW_mesh = np.meshgrid(LW_vals, k27_vals)

    #f_H2[data['k27'] == 1.58489319] = 100.0 # flag to figure out orientation
    f_H2 = data[:,2]
    z_mesh = f_H2.reshape( int(np.sqrt(np.size(k27))), int(np.sqrt(np.size(LW))))
    #z_mesh = z[:-1,:-1]

    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    img1 = ax.pcolormesh(10.0**(LW_mesh),
                         10.0**(k27_mesh),
                         np.log10(z_mesh.T), cmap = 'magma',
                         vmin = np.min(np.log10(z_mesh)),
                         vmax = np.max(np.log10(z_mesh)))
    ax.semilogx()
    ax.semilogy()
    ax.set_xlabel('log(LW Factor)')
    ax.set_ylabel("log(k27 Factor)")

    divider = make_axes_locatable(ax)
    cax1 = divider.append_axes('right', size = '5%', pad = 0.05)
    fig.colorbar(img1, cax=cax1, label = r'log(f$_{\rm H_2}$)')

    plt.minorticks_on()

    plt.tight_layout(h_pad = 0, w_pad = 0.05)
    fig.savefig("test_fH2.png")
    plt.close()

    f_H2 = data[:,3]
    z_mesh= f_H2.reshape( int(np.sqrt(np.size(k27))), int(np.sqrt(np.size(LW))))
    #z_mesh = z[:-1,:-1]

    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    img1 = ax.pcolormesh(10.0**(LW_mesh),
                         10.0**(k27_mesh),
                         np.log10(z_mesh.T), cmap = 'RdYlBu_r',
                         vmin = np.min(np.log10(z_mesh)),
                         vmax = np.max(np.log10(z_mesh)))
    ax.semilogx()
    ax.semilogy()
    ax.set_xlabel('log(LW Factor)')
    ax.set_ylabel("log(k27 Factor)")

    divider = make_axes_locatable(ax)
    cax1 = divider.append_axes('right', size = '5%', pad = 0.05)
    fig.colorbar(img1, cax=cax1, label = r'log(Temperature [K])')

    plt.minorticks_on()

    plt.tight_layout(h_pad = 0, w_pad = 0.05)
    fig.savefig("test_T.png")
    plt.close()

    return

if __name__ == "__main__":
    plot_2d_histogram()
