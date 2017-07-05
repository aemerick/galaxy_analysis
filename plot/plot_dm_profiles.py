import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

#
from dwarfs.analysis.initial_conditions import ic_list as icl

purple  = '#7d03a8'
magenta = '#cb4679'
blue    = '#0c0887'
orange  = '#fdc328'
black   = 'black'
green   = 'green'

colors = [orange, magenta, purple, black]
ls     = ['-','-','-','-']

#
# style defines
#
fsize = 17
rc('text', usetex=False)
rc('font', size=fsize)#, ftype=42)
line_width = 4
point_size = 30


#
# Get the 4 Leo P galaxy IC's
#
run11_largebox = icl.ic_object_dict['Leo_P']
run11_2rs      = icl.ic_object_dict['Leo_P_2rs']
run11_30km     = icl.ic_object_dict['Leo_P_30km']
run11_40km     = icl.ic_object_dict['Leo_P_40km']

galaxies = [run11_largebox, run11_2rs, run11_30km, run11_40km]
labels    = ['17', '25', '30', '40']
for i in np.arange(len(labels)):
    labels[i] = r'v$_{\rm c,max}$ = ' + labels[i] + r' km s$^{-1}$'


print labels

def overplot_scale_radius(ax, gal, color = 'black', ls = '--',
                              virial = False):
    if virial:
        norm = gal.ic['R200']
    else:
        norm = icl.cgs.kpc

    ylim = ax.get_ylim()
    x = [ gal.ic['b'] / norm,  gal.ic['b'] / norm ]
    ax.plot(x, ylim, lw = line_width*0.75, color = color, ls = ls)
    ax.set_ylim(ylim)

    return

def plot_dm_profiles(virial = False):
    #
    # Plot dark matter density as a function of radius
    #
    fig, ax = plt.subplots()

    if virial:
        r = np.logspace(-4, 0) # in virial radius
    else:
        r = np.logspace(-3,2)

    ax.loglog()
    for i,g in enumerate(galaxies):
    
        if virial:
            x = r * g.ic['R200']
        else:
            x = r * icl.cgs.kpc

        ax.plot( r, g.DM_density(x),
                    lw = line_width, color = colors[i], ls = ls[i],
                    label = labels[i])

    for i,g in enumerate(galaxies):
        overplot_scale_radius(ax, g, color = colors[i], ls = '--', virial = virial)

    #ax.set_xlim(-4, 0)
    #ax.set_ylim()

    if virial:
        ax.set_xlabel(r'log[R / R$_{\rm vir}$]')
    else:
        ax.set_xlabel(r'R (kpc)')
    ax.set_ylabel(r'Dark Matter Density (g cm$^{-3}$)')

    ax.legend(loc = 'best')
    fig.set_size_inches(8,8)
    plt.tight_layout()
    
    outname = 'dark_matter_density'
    if virial:
        outname += '_virial'

    plt.savefig(outname +'.png')

    return


def plot_mass_profiles(virial):
    fig, ax = plt.subplots()

    if virial:
        r = np.logspace(-4, 0) # in virial radius
    else:
        r = np.logspace(-3,2)

    ax.loglog()

    for i,g in enumerate(galaxies):

        if virial:
            x = r * g.ic['R200']
        else:
            x = r * icl.cgs.kpc

        ax.plot( r, g.M_r(x) / icl.cgs.Msun,
                    lw = line_width, color = colors[i], ls = ls[i],
                    label = labels[i])

    for i, g in enumerate(galaxies):
        overplot_scale_radius(ax, g, color = colors[i], ls = '--', virial = virial)

    #ax.set_xlim(-4, 0)
    #ax.set_ylim()

    if virial:
        ax.set_xlabel(r'log[R / R$_{\rm vir}$]')
    else:
        ax.set_xlabel(r'R (kpc)')
    ax.set_ylabel(r'Cumulative Mass (M$_{\odot}$)')

    ax.legend(loc = 'best')
    fig.set_size_inches(8,8)
    plt.tight_layout()

    outname = 'dark_matter_mass'
    if virial:
        outname += '_virial'

    plt.savefig(outname + '.png')

    return

def plot_circular_velocity_profiles(virial):
    fig, ax = plt.subplots()

    if virial:
        r = np.logspace(-4, 0) # in virial radius
    else:
        r = np.logspace(-3,2)

    ax.loglog()
    for i,g in enumerate(galaxies):

        if virial:
            x = r * g.ic['R200']
        else:
            x = r * icl.cgs.kpc

        ax.plot( r, g.circular_velocity(x) / 1.0E5,
                    lw = line_width, color = colors[i], ls = ls[i],
                    label = labels[i])

    for i,g in enumerate(galaxies):
        overplot_scale_radius(ax, g, color = colors[i], ls = '--',
                              virial = virial)

    #ax.set_xlim(-4, 0)
    #ax.set_ylim()

    if virial:
        ax.set_xlabel(r'log[R / R$_{\rm vir}$]')
    else:
        ax.set_xlabel(r'R (kpc)')
    ax.set_ylabel(r'V$_{\rm c}$ (km s$^{-1}$)')

    ax.legend(loc = 'best')
    fig.set_size_inches(8,8)
    plt.tight_layout()

    outname = 'circular_velocity'
    if virial:
         outname += '_virial'

    plt.savefig(outname + '.png')


    return

if __name__ == '__main__':

    for virial in [True, False]:

        plot_dm_profiles(virial)
        plot_mass_profiles(virial)
        plot_circular_velocity_profiles(virial)
