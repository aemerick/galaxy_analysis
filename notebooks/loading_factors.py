import numpy as np
import yt
from galaxy_analysis import Galaxy
import matplotlib.pyplot as plt
from galaxy_analysis.plot.plot_styles import *
import deepdish as dd
import glob
import os

from galaxy_analysis.utilities import utilities



wdir = "./"

#
# Define region
#

def compute_mass_outflow(ds, data):
    
    dL = 200.0 * yt.units.pc
    above_region = data.cut_region(["(obj['cylindrical_z'].in_units('kpc') > 0.9) & (obj['cylindrical_z'].in_units('kpc') < 1.1) & (obj['cylindrical_radius'].in_units('kpc') < 2)"])
    above_hot_region = above_region.cut_region(["obj['temperature'].in_units('K') > 3.0E5"])
    above_cold_region = above_region.cut_region(["obj['temperature'].in_units('K') < 3.0E5"])
    above_colder_region = above_region.cut_region(["obj['temperature'].in_units('K') < 1.0E4"])

    below_region = data.cut_region(["(obj['cylindrical_z'].in_units('kpc') < -0.9) & np.abs(obj['cylindrical_z'].in_units('kpc') > -1.1) & (obj['cylindrical_radius'].in_units('kpc') < 2)"])
    below_hot_region = below_region.cut_region(["obj['temperature'].in_units('K') > 3.0E5"])
    below_cold_region = below_region.cut_region(["obj['temperature'].in_units('K') < 3.0E5"])
    below_colder_region = below_region.cut_region(["obj['temperature'].in_units('K') < 1.0E4"]) 
    
    i = 0
    
    total_vals = np.zeros(3)
    metal_vals = np.zeros(3)
    
    for above_reg, below_reg in [(above_region, below_region),
                                 (above_hot_region,below_hot_region),
                                 (above_cold_region,below_cold_region)]:
                                # (above_colder_region,below_colder_region)]:
        
        M = above_reg['cell_mass'].to('Msun')
        O = above_reg['O_Density'].value / above_reg['Density'].value * M    
        v = above_reg['z-velocity'].to('km/s')
        select = v > 0
    
        total_outflow_rate = np.sum(M[select]*v[select]/dL).to('Msun/yr')
        metal_outflow_rate = np.sum(O[select]*v[select]/dL).to('Msun/yr')
    
        M = below_reg['cell_mass'].to('Msun')
        O = below_reg['O_Density'].value / below_reg['Density'].value * M
        v = below_reg['z-velocity'].to('km/s')
        select = v < 0
    
        total_outflow_rate += np.sum(-M[select]*v[select]/dL).to('Msun/yr')
        metal_outflow_rate += np.sum(-O[select]*v[select]/dL).to('Msun/yr')
    
        total_vals[i] = total_outflow_rate
        metal_vals[i] = metal_outflow_rate
        
        i = i + 1
    
    # now I need to compute the SFR and metal production rate
    
    
    return total_vals[0], total_vals[1], total_vals[2], metal_vals[0], metal_vals[1], metal_vals[2]
   
def compute_energy_outflow(ds, data):
    
    birth_mass = data['birth_mass']
    birth_time = data['creation_time'].to('Myr')
    pt = data['particle_type']

    t = ds.current_time.to('Myr')


    select = ((birth_mass > 8.0) * (birth_mass < 25.0)) * (pt > 11) * ((t - birth_time ) < 40*yt.units.Myr)

    N_SN = np.size(birth_mass[select])

    dL = 200 * yt.units.pc
    region = None
    above_region = data.cut_region(["(obj['cylindrical_z'].in_units('kpc') > 0.9) & (obj['cylindrical_z'].in_units('kpc') < 1.1) & (obj['cylindrical_radius'].in_units('kpc') < 2)"])
    above_hot_region = above_region.cut_region(["obj['temperature'].in_units('K') > 3.0E5"])
    above_cold_region = above_region.cut_region(["obj['temperature'].in_units('K') < 3.0E5"])
    above_colder_region = above_region.cut_region(["obj['temperature'].in_units('K') < 1.0E4"])

    below_region = data.cut_region(["(obj['cylindrical_z'].in_units('kpc') < -0.9) & np.abs(obj['cylindrical_z'].in_units('kpc') > -1.1) & (obj['cylindrical_radius'].in_units('kpc') < 2)"])
    below_hot_region = below_region.cut_region(["obj['temperature'].in_units('K') > 3.0E5"])
    below_cold_region = below_region.cut_region(["obj['temperature'].in_units('K') < 3.0E5"])
    below_colder_region = below_region.cut_region(["obj['temperature'].in_units('K') < 1.0E4"])


    vals = np.zeros(4)
    
    if N_SN == 0:
        print("No Supernova this time step")
        return vals
    
    
    i = 0
    for above_reg, below_reg in [(above_region, below_region),
                                 (above_hot_region,below_hot_region),
                                 (above_cold_region,below_cold_region),
                                 (above_colder_region,below_colder_region)]:
    
        M  = above_reg['cell_mass'].to('Msun')
        v  = above_reg['z-velocity'].to('km/s')
        cV = above_reg['cell_volume'].to('cm**(3)')
        KE = 0.5*M*above_reg['velocity_magnitude']**2
        
        select = v > 0                  
        above_E_out_therm = (v / dL * above_reg['thermal_energy'].to('erg/g') * M).to('erg/s')[select]
        above_E_out_kin   = (v / dL * KE).to('erg/s')[select]
    
        M  = below_reg['cell_mass'].to('Msun')
        v  = below_reg['z-velocity'].to('km/s')
        cV = below_reg['cell_volume'].to('cm**(3)')
        KE = 0.5*M*below_reg['velocity_magnitude']**2
        
        select = v < 0           
        below_E_out_therm = (-v / dL * below_reg['thermal_energy'].to('erg/g') * M).to('erg/s')[select]
        below_E_out_kin   = (-v / dL * KE).to('erg/s')[select]    
    
    
        E_out_therm = np.sum(above_E_out_therm) + np.sum(below_E_out_therm)
        E_out_kin   = np.sum(above_E_out_kin  ) + np.sum(below_E_out_kin  )
        E_out_tot   = E_out_therm + E_out_kin
    
        #print("Out-rate      :  %10.3E %10.3E %10.3E"%(E_out_therm.value, E_out_kin.value, E_out_tot.value))
        #x = N_SN*yt.units.erg*1.0E51 / ((40.0*yt.units.Myr).to('s'))
        x = 1.0
        #print("Loading Factor: %10.3E %10.3E %10.3E"%((E_out_therm/x).value, (E_out_kin/x).value, (E_out_tot/x).value))
    

        vals[i] = (E_out_tot/x).value * 3.154E7 # convert to erg / yr
        i = i + 1
    
    return vals



#
#
#
#
#
#
def compute_loading_denominators(directory):
    
    files     = np.sort(glob.glob(directory + '/*galaxy_data*.h5'))
    
    all_times  = np.zeros(np.size(files))
    #SFR    = np.zeros(np.size(files))
    #SNR    = np.zeros(np.size(files))
    #O_rate = np.zeros(np.size(files))
    
    for i,f in enumerate(files):
        meta_data = dd.io.load(f,'/meta_data')
        all_times[i] = meta_data['Time']
    #    SFR[i]   = meta_data['SFR_100']
    #    SNR[i]   = dd.io.load(f,'/time_data/SNII_snr_100')

    tdata = dd.io.load(files[-1], '/time_data')

    times = tdata['time_100']
    SNR   = tdata['SNII_snr_100']
    SFR   = tdata['SFR_100']

    Omass = np.zeros(np.size(SFR))

    for i,t in enumerate(times):
        index = np.argmin(np.abs(50.0 + t/1.0E6-all_times - all_times[0]))
        Omass[i] = dd.io.load(files[index],
                             '/gas_meta_data/masses/FullBox')['O'] +\
                   dd.io.load(files[index],
                              '/gas_meta_data/masses/OutsideBox')['O']

    Orate = np.zeros(np.size(SFR))
    Orate[:-1] = (Omass[1:] - Omass[:-1]) / (times[1:] - times[:-1])
    Orate[-1]  = Orate[-2]

    return times, SFR, SNR, Orate


t_100, SFR_100, SNR_100, Omass_100 = compute_loading_denominators(wdir)

t_100   = t_100 / 1.0E6 # in Myr

f_SFR   = lambda x: np.interp(x,t_100, SFR_100)
f_SNR   = lambda x: np.interp(x,t_100, SNR_100)
f_Omass = lambda x: np.interp(x,t_100, Omass_100)


#
#
#
#
#
#


dsnames = np.sort(np.array(glob.glob(wdir + "DD???5/DD???5") +glob.glob(wdir+"DD???0/DD???0")))

if len(dsnames) > 120: # only select first 100 assuming only grabbing every 5 Myr (first 500 Myr)
    dsnames = dsnames[:120]

j = 0
ds_select = dsnames[:]

hot_loading = np.zeros(np.size(ds_select))
cold_loading = np.zeros(np.size(ds_select))
colder_loading = np.zeros(np.size(ds_select))

mass_load = np.zeros(np.size(ds_select))
mass_load_hot= np.zeros(np.size(ds_select))
mass_load_cold = np.zeros(np.size(ds_select))

metal_load = np.zeros(np.size(ds_select))
metal_load_hot = np.zeros(np.size(ds_select))
metal_load_cold = np.zeros(np.size(ds_select))


sfr = np.zeros(np.size(mass_load))
snr = np.zeros(np.size(sfr))
orate = np.zeros(np.size(snr))

loading = np.zeros(np.size(ds_select))

if os.path.isfile(wdir + "loading_factors.dat"):
    previous_data = np.genfromtxt(wdir + "loading_factors.dat", names=True)

    if len(previous_data['time']) == 0: # hack
        previous_data = {}
        previous_data['time'] = np.ones(2)*-1

    file = open(wdir + "loading_factors.dat","a")

else:
    previous_data = {'time' : np.array([-1,-1])}
    file = open(wdir + "loading_factors.dat","w")
    file.write("#time E_out E_hot_out E_cold_out E_colder_out M_out M_out_hot M_out_cold Metal_out Metal_out_hot Metal_out_cold SFR SNR O_rate\n")

times = np.zeros(np.size(ds_select))
for j,dsname in enumerate(ds_select):

    ds = yt.load(dsname)
    
    times[j] = ds.current_time.to('Myr')

    if times[j] <= previous_data['time'][-1]:
        del(ds)
        continue

    if times[j] > (500.0 + 120.0):
        break

    data = ds.all_data()
    
    try:
        x,y,z,w = compute_energy_outflow(ds,data)
    except:
        continue
    
    loading[j]        = x
    hot_loading[j]    = y
    cold_loading[j]   = z
    colder_loading[j] = w
    
    mass_load[j], mass_load_hot[j], mass_load_cold[j], metal_load[j], metal_load_hot[j], metal_load_cold[j] = compute_mass_outflow(ds,data)
    
    sfr[j] = f_SFR(times[j])
    snr[j] = f_SNR(times[j])
    orate[j] = f_Omass(times[j])
    
    
    print("%5.2f Myr:   %10.3E  %10.3E  %10.3E  %10.3E"%(times[j],loading[j],hot_loading[j],cold_loading[j],colder_loading[j]))
    file.write("%5.2f %10.3E  %10.3E  %10.3E  %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E\n"%(times[j],loading[j],hot_loading[j],cold_loading[j],colder_loading[j],
                                                              mass_load[j],mass_load_hot[j], mass_load_cold[j], metal_load[j], metal_load_hot[j], metal_load_cold[j]
                                                                                                   ,sfr[j],snr[j],orate[j]))
    file.flush()

    del(data)
    del(ds)
    
file.close()
