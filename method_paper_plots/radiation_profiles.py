import yt
from galaxy_analysis.plot.plot_styles import *
import deepdish as dd

from galaxy_analysis.utilities import utilities
import numpy as np
import matplotlib.pyplot as plt

# time average
from galaxy_analysis.analysis import compute_time_average as ta
import sys

def plot(work_dir = './', t_min = 300.0, t_max = 350.0, xmax = 600.0,
         fields = ['G_o','LW_flux','FUV_flux','Q0_flux','Q1_flux','Pe_heating_rate_masked'],
         outdir = './'):
    """

    """


    # gather time averaged profiles for all
    all_profiles = {}
    for f in fields:
        radius, avg, min, max, std, q1,q2,q3= ta.compute_time_average(['gas_profiles','radiation','disk',f], dir = work_dir,
                                                  tmin = t_min, tmax = t_max, x_field = 'xbins', return_quartiles =True)
        all_profiles[f] = {'avg':avg,'min':min,'max':max,'std':std, 'radius' : radius,
                           'centers' : 0.5*(radius[1:] + radius[:-1]), 'q1' : q1, 'q2' : q2, 'q3' : q3, 'median' : q2}

    # Plot for G_o #
    if 'G_o' in fields:
        fields = [x for x in fields if not (x=='G_o')]
        fig, ax = plt.subplots()
        fig.set_size_inches(8,8)
        color = 'black'
        # now plot
        plot_histogram(ax, all_profiles['G_o']['radius'], all_profiles['G_o']['avg'],
                        color = color, label = r'G$_{\rm o}$', lw = line_width)
        fillmin = all_profiles['G_o']['min'] #- all_profiles['G_o']['std']
        fillmax = all_profiles['G_o']['avg'] + np.abs( all_profiles['G_o']['std'])
        fillmax = all_profiles['G_o']['max']

        fillmin[fillmin <= 0.0] = 0.00324

        ax.fill_between(all_profiles['G_o']['centers'], fillmin, fillmax,
                           facecolor = color, interpolate=False, edgecolor = color, alpha = 0.5)
        ax.set_xlabel(r'Radius (pc)')
        ax.set_xlim(0.0, xmax)
        ax.set_ylim( 3.0E-3, 5.0)
        ax.set_ylabel(r'G$_{\rm o}$')
        ax.semilogy()
        plt.minorticks_on()
        plt.tight_layout()
        fig.savefig(outdir + 'G_o_profile.png')
        plt.close()

    if 'Pe_heating_rate_masked' in fields:
        fields = [x for x in fields if not (x=='Pe_heating_rate_masked')]
        fig, ax = plt.subplots()
        fig.set_size_inches(8,8)
        color = 'C7'
        # now plot
        plot_histogram(ax, all_profiles['Pe_heating_rate_masked']['radius'], all_profiles['Pe_heating_rate_masked']['avg'],
                            color = color, label = r'G$_{\rm o}$', lw = line_width)
        fillmin = all_profiles['Pe_heating_rate_masked']['min'] #- all_profiles['Pe_heating_rate_masked']['std']
        fillmax = all_profiles['Pe_heating_rate_masked']['avg'] + np.abs(all_profiles['Pe_heating_rate_masked']['std'])

        ax.fill_between(all_profiles['Pe_heating_rate_masked']['centers'], fillmin, fillmax,
                               facecolor = color, interpolate=False, edgecolor = color, alpha = 0.5)
        ax.set_xlabel(r'Radius (pc)')
        ax.set_xlim(0.0, xmax)
    #    ax.set_ylim( 3.0E-3, 1.0)
        ax.set_ylabel(r'$\Gamma_{\rm PE}$')
        ax.semilogy()
        plt.minorticks_on()
        plt.tight_layout()
        fig.savefig(outdir + 'pe_heating_rate_profile.png')
        plt.close()


    if 'Q0_flux' in fields:
        fields = [x for x in fields if not (x=='Q0_flux')]
        fig, ax = plt.subplots()
        fig.set_size_inches(8,8)
        color1 = 'black'
        color2 = 'C1'

        # now plot
        q0_conv = 1.0 # (13.6 * yt.units.eV).to('erg').value
        q1_conv = 1.0 # (24.6 * yt.units.eV).to('erg').value
        lower_lim = 1.0E-8

        plot_histogram(ax, all_profiles['Q0_flux']['radius'], all_profiles['Q0_flux']['avg'] * q0_conv,
                            color = color1, label = r'Q$_{\rm 0}$', lw = line_width)

        fillmin = all_profiles['Q0_flux']['q1']
        fillmax = all_profiles['Q0_flux']['max']

        fillmin[ fillmin <= 0.0] = np.min(fillmin[fillmin>0.0]) # lower_lim # * np.min(all_profiles['Q0_flux']['avg'])
        fillmax[ fillmax <= 0.0] = np.min(fillmax[fillmax>0.0]) #lower_lim

        # shade in std
        ax.fill_between(all_profiles['Q0_flux']['centers'], fillmin*q0_conv, fillmax*q0_conv,
                               facecolor = color1, interpolate=False, edgecolor = color1, alpha = 0.5)

    #    plot_histogram(ax, all_profiles['Q1_flux']['radius'], all_profiles['Q1_flux']['avg']*q1_conv,
    #                        color = color2, label = r'Q$_{\rm 1}$', lw = line_width)

        fillmin = all_profiles['Q1_flux']['q1']
        fillmax = all_profiles['Q1_flux']['max']

        fillmin[ fillmin <= 0.0] = np.min(fillmin[fillmin>0.0]) # lower_lim # * np.min(all_profiles['Q0_flux']['avg'])
        fillmax[ fillmax <= 0.0] = np.min(fillmax[fillmax>0.0]) #lower_lim


        # shade in std
    #    ax.fill_between(all_profiles['Q1_flux']['centers'], fillmin*q1_conv, fillmax*q1_conv,
    #                           facecolor = color2, interpolate=False, edgecolor = color2, alpha = 0.5)

        ax.set_xlabel(r'Radius (pc)')
        ax.set_xlim(0.0, xmax)
        ax.set_ylabel(r'Ionizing Radiation Flux (erg s$^{-1}$ cm$^{-2}$)')
        ax.semilogy()
        ax.set_ylim(1.0E-9, 1.0E-3)
        plt.minorticks_on()
        plt.tight_layout()
        fig.savefig(outdir + 'ionizing_photon_profile.png')
        plt.close()


    if 'LW_flux' in fields:
        fields = [x for x in fields if not (x=='LW_flux')]
        fig, ax = plt.subplots()
        fig.set_size_inches(8,8)
        color1 = 'C1'
        color2 = 'C2'

        # now plot
        plot_histogram(ax, all_profiles['LW_flux']['radius'], all_profiles['LW_flux']['avg'],
                            color = color1, label = r'Lyman-Werner', lw = line_width)
        fillmin = all_profiles['LW_flux']['avg'] #- all_profiles['LW_flux']['std']
        fillmax = all_profiles['LW_flux']['avg'] + all_profiles['LW_flux']['std']
        #    fillmin[ fillmin < 0 ] = 1.0E-10 * np.min(all_profiles['LW_flux']['avg'])
        # shade in std
        ax.fill_between(all_profiles['LW_flux']['centers'], fillmin, fillmax,
                               facecolor = color1, interpolate=False, edgecolor = color1, alpha = 0.5)

        plot_histogram(ax, all_profiles['FUV_flux']['radius'], all_profiles['FUV_flux']['avg'],
                            color = color2, label = r'FUV', lw = line_width)
        fillmin = all_profiles['FUV_flux']['avg'] # - all_profiles['FUV_flux']['std']
        fillmax = all_profiles['FUV_flux']['avg'] + all_profiles['FUV_flux']['std']
        #    fillmin[ fillmin < 0 ] = 1.0E-6

        # shade in std
        ax.fill_between(all_profiles['FUV_flux']['centers'], fillmin, fillmax,
                           facecolor = color2, interpolate=False, edgecolor = color2, alpha = 0.5)

        ax.set_xlabel(r'Radius (pc)')
        ax.set_xlim(0.0, xmax)
        ax.set_ylabel(r'Optically Thin Radiation Flux (erg s$^{-1}$ cm$^{-2}$')
        ax.semilogy()
        ax.set_ylim(1.0E-5, 1.0E-3)
        plt.minorticks_on()
        ax.legend(loc='best')
        plt.tight_layout()
        fig.savefig(outdir + 'optically_thin_radiation_profile.png')
        plt.close()

    if len(fields) > 0:
        print "Failure in radiation plots. Currently no defined methods to plot:"
        print fields

    return


if __name__ == "__main__":

    work_dir = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11_30km/final_sndriving/'
    plot(work_dir = workd_dir)
