import numpy as np
import glob
import deepdish as dd

#
# This may win awards for quickest-dirtiest scripts I've made 
# recently......
#

# compute averages:

# convert to loading factors - compute averages

#file_list = np.sort(glob.glob('loading_factors*.dat'))

#data_list = [None] * len(file_list)

#for name in data_list:
#    data_list[i] = np.genfromtxt(name,names=True)

data   = np.genfromtxt('loading_factors.dat',names=True)
#data2   = np.genfromtxt('loading_factors2.dat',names=True)

#data = {}
#for k in data1.dtype.names:
#    data[k] = np.array(list(data1[k]) + list(data2[k])) # ew


#
#correct_orate = np.genfromtxt('orate.dat',names=True)
#
#data['O_rate'] = correct_orate['O_rate'] #####


loading_data = {}
for k in ['M_out','M_out_hot','M_out_cold']:
    select = data['M_out'] > 0
    loading_data[k] =   np.average( data[k][select] / data['SFR'][select])   
for k in ['Metal_out','Metal_out_hot','Metal_out_cold']:
    select = data['Metal_out'] > 0
    loading_data[k] =   np.average( data[k][select] / data['O_rate'][select])      
for k in ['E_out','E_hot_out','E_cold_out','E_colder_out']:
    select = data['E_out'] > 0
    loading_data[k] = np.average( data[k][select] / (data['SNR'][select]*1.0E51))

#
# compute SFR for gas and stars
#
files = np.sort(glob.glob("./*galaxy_data*.h5"))
gassdens = np.zeros(np.size(files))

for i,f in enumerate(files):
    try:
        r = dd.io.load(f,"/gas_profiles/surface_density/disk/xbins")
    except:
        gassdens[i] = np.average(gassdens[gassdens>0])
        continue

    gassdens[i] =dd.io.load(f, "/meta_data/M_gas")
A = np.pi * 600.0 * 600.0 # surface area of galaxy disk
gassdens  = np.average(gassdens) / A                                                   # in Msun / pc^2
starsdens = np.average(   dd.io.load( files[-1], '/time_data/SFR_100')) / A  * 1.0E6   # in Msun / yr / kpc^2

# compute other values
loading_data['Eta_h-Eta_c'] = loading_data['E_hot_out'] / loading_data['E_cold_out']
loading_data['e_s'] = np.average(  data['E_out'] / (data['M_out']*1.989E33 )) # erg / g
select = data['M_out_hot'] > 0
loading_data['e_s_hot'] = np.average(  data['E_hot_out'][select] / (data['M_out_hot'][select]*1.989E33 )) # erg / g
loading_data['e_s_cold'] = np.average(  data['E_cold_out'] / (data['M_out_cold']*1.989E33 )) # erg / g
loading_data['e_s_h-e_s_c'] = loading_data['e_s_hot'] / loading_data['e_s_cold']


loading_data['Sigma_gas'] = gassdens
loading_data['Sigma_sfr'] = starsdens

loading_data['E_h-Metal_h'] = loading_data['E_hot_out'] / loading_data['Metal_out_hot']

names = {  'M_out'      : 'Eta_{mass}',
           'M_out_hot'  : "Eta_{mass,hot}",
           'M_out_cold' : "Eta_{mass,cold}",
           'Metal_out'  : "Eta_{metal}",
           'Metal_out_hot' : "Eta_{metal,hot}",
           'Metal_out_cold' : "Eta_{metal,cold}",
           'E_out' : "Eta_{E}",
           'E_hot_out' : "Eta_{E,hot}",
           'E_cold_out' : "Eta_{E,cold}",
           'Eta_h-Eta_c' : "Eta_{E,hot} / Eta_{E,cold}",
           'E_h-Metal_h' : 'Eta_{E,hot} / Eta_{metal,hot}',
           'e_s'        : 'e_s (erg/g)',
           'e_s_hot'    : 'e_s,h (erg/g)',
           'e_s_cold'   : 'e_s,c (erg/g)',
           'e_s_h-e_s_c' : 'e_{s,h} / e_{s,c}',
           'Sigma_gas'   : 'Sigma_gas (Msun / pc^2)',
           'Sigma_sfr'   : 'Sigma_sfr (Msun / yr / kpc^2)'}


for k in loading_data.keys():
    if 'colder' in k:
        continue
    print("%30s  %8.3E"%(names[k],loading_data[k]))


