from galaxy_analysis.plot.plot_styles import *
import matplotlib.pyplot as plt
import numpy as np
import sys


def plot_lowrad(file = 'low_eV_luminosity.dat'):

    data = np.genfromtxt(file, names = True)

    t = data['Time'] - data['Time'][0]

    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    Ltot = data['L_low'] + data['L_high']
    Ltot = 1.0

    ax.plot(t,  data['L_low'] / Ltot, lw = 3, ls = '-', color = 'C0', label = 'M < 8')
    ax.plot(t,  data['L_high'] / Ltot, lw = 3, ls = '-', color = 'C1', label = 'M > 8')

    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(r'H- Photodetachment Luminosity')
    ax.set_xlim(0.0, 500.0)
    ax.semilogy()

    ax.legend(loc='best')
    plt.tight_layout()
    fig.savefig('low_eV_plot.png')
    plt.close()
############
    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    Ltot = data['L_low_LW'] + data['L_high_LW']
    Ltot = 1.0

    ax.plot(t,  data['L_low_LW'] / Ltot, lw = 3, ls = '-', color = 'C0', label = 'M < 8')
    ax.plot(t,  data['L_high_LW'] / Ltot, lw = 3, ls = '-', color = 'C1', label = 'M > 8')

    ax.set_xlabel(r'Time (Myr)')
    ax.set_ylabel(r'LW Luminosity')
    ax.set_xlim(0.0, 500.0)
    ax.semilogy()

    ax.legend(loc='best')
    plt.tight_layout()
    fig.savefig('LW_L_plot.png')
    plt.close()
    return

if __name__ == "__main__":
    plot_lowrad() 


