import numpy as np
import yt
from matplotlib import rc
fsize = 17
rc('text', usetex=False)
rc('font', size=fsize)#, ftype=42)
line_width = 3
point_size = 30
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from galaxy_analysis.particle_analysis import particle_types as pdef

def plot_dtd(ds):
    data = ds.all_data()

    snIa = pdef.snIa(ds, data)
    WD   = pdef.white_dwarfs(ds, data)

    WD_death   =  data['dynamical_time'][WD] # + data['creation_time'][WD] 
    SNIa_death =  data['dynamical_time'][snIa] # + data['creation_time'][snIa]

    WD_death = list(WD_death.convert_to_units("Gyr").value)
    SNIa_death = list(SNIa_death.convert_to_units("Gyr").value)

    fig, ax = plt.subplots()
    all = np.array( WD_death + SNIa_death)

    hist, bins = np.histogram(all, bins = np.arange(0,14.25,0.5))

    x = 0.5 * (bins[1:] + bins[:-1])
    ax.plot(x, hist, lw = 3, color = 'black', ls = '-')
    y = x**(-1.0* ds.parameters['IndividualStarDTDSlope'])
    norm = hist[0] / y[0]
    ax.plot(x, norm*y, lw = 3, color = 'black', ls='--')
    ax.plot(x, hist[0]/((x[0])**(-1.01)) * x**(-1.01),lw =3, color = 'black',ls=':')
    ax.set_xlabel(r'Time (Gyr)')
    ax.set_ylabel(r'Binned SNIa (counts)')
    ax.loglog()

    fig.set_size_inches(8,8)
    plt.tight_layout()
    plt.minorticks_on()
    fig.savefig('dtd.png')
    plt.close()

    return

if __name__ == "__main__":
    ds   = yt.load('DD0205/DD0205')
    data = ds.all_data()

    plot_dtd(ds)
