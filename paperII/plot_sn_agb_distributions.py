from galaxy_analysis.plot.plot_styles import *
import deepdish as dd
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u


colors = {'AGB': 'C0', 'SNII' : 'C1', 'SNIa' : 'C2'}
ls     = {'AGB': '-', 'SNII' : '-', 'SNIa' : '-'}

def plot_data(data, field, log_field = True, field_label = None, field_bins = None, outname = None, abs_field = False):

    fig, ax = plt.subplots()
    fig.set_size_inches(8,8)

    #
    # plot this separately for AGB stars, SN, SNIa
    #
    select = {}
    select['AGB']   = (data['M_o'] > 3.0) * (data['M_o'] < 8.0) #* (data['ptype'] == 11)
    select['SNIa']  = np.array([False]*np.size(data['M_o'])) # (data['M_o'] > 3.0) * (data['M_o'] < 8.0) * (data['ptype'] != 11)
    select['SNII']  = (data['M_o'] > 8.0) * (data['M_o'] < 25.0)

    maxval = -1
    for k in ['SNII','AGB']:

        if field == 'M':
            volume = (4.0 * np.pi / 3.0 * (   ((25.0 * u.pc).to(u.cm))**3)).value
            xdata = ((data['M_tot'][select[k]]*u.Msun).to(u.g)).value / (volume)
        else:
            xdata = data[field][select[k]]

        if abs_field:
            xdata = np.abs(xdata)

        if log_field:
            xdata = np.log10(xdata)

        hist, bins = np.histogram( xdata, bins = field_bins)

        if np.sum(hist) > 0:
            hist = hist / (1.0 * np.sum(hist))
        else:
            hist = np.zeros(np.size(hist))

        plot_histogram(ax, bins, hist, lw = line_width, color = colors[k], ls = ls[k], label = k)

        maxval = np.max([maxval, np.max(hist)])


    ax.set_ylabel('Fraction of Events')

    if field_label is None:
        field_label = field

        if log_field:
            field_label = r'log(' + field_label + ')'

    ax.legend(loc='best')
    ax.set_xlabel(field_label)
    plt.minorticks_on()
    ax.set_xlim(field_bins[0], field_bins[-1])
    ax.set_ylim(0, maxval + 0.02)
    plt.tight_layout()
    fig.savefig(outname)

    plt.close()
    return


if __name__ == '__main__':
    data = np.genfromtxt('stellar_environment.dat',names=True)

    plot_data(data, 'n_v_avg', log_field = True, field_bins = np.arange(-6,3,0.2), field_label = r'log(n [cm$^{-3}$])', outname = 'SE_density.png')
    plot_data(data, 'n_max', log_field = True, field_bins = np.arange(-6,3,0.2), field_label = r'log(n max [cm$^{-3}$])', outname = 'SE_n_max.png')
    plot_data(data, 'n_min', log_field = True, field_bins = np.arange(-6,3,0.2), field_label = r'log(n min [cm$^{-3}$])', outname = 'SE_n_min.png')


    plot_data(data, 'T_m_avg', log_field = True, field_bins = np.arange(1,3,0.2), field_label = r'log(T [K])', outname = 'SE_T_m.png')
    plot_data(data, 'T_v_avg', log_field = True, field_bins = np.arange(0,8,0.2), field_label = r'log(T [K])', outname = 'SE_T_v.png')

    plot_data(data, 'r',   log_field = False, field_bins = np.arange(0, 601, 10), field_label = r'R (pc)', outname = 'SE_r.png')
    plot_data(data, 'z',   log_field = False, abs_field = True, field_bins = np.arange(0, 101, 5),  field_label = r'|z| (pc)', outname = 'SE_z.png')
    plot_data(data, 'M', log_field = True, field_bins = np.arange(-30, -20, 0.2), field_label = r'log(rho)', outname = 'SE_M_density.png')

