from galaxy_analysis.plot.plot_styles import *
from galaxy_analysis.analysis.compute_time_average import compute_time_average
from galaxy_analysis.utilities import utilities
import deepdish as dd
import numpy as np
import matplotlib.pyplot as plt
import glob

# filepath = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11_30km/final_sndriving'
work_dir      = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11_30km/final_sndriving/'
#work_dir      = '/mnt/ceph/users/emerick/enzo_runs/pleiades/starIC/run11/corrected_sndriving/'

data_list, times = utilities.select_data_by_time(dir = work_dir,
                                                 tmin=0.0,tmax=1000)
print data_list[0], times[0]
print data_list[-1], times[-1]
print times[-1] - times[0]
all_data = {}

# gather metallicities for:
#   1) ISM
#   2) Halo
#   3) Outflowing gas (at all radii)
#   4) Inflowing gas  (at all radii)
#data_list = data_list[:60]
#times     = times[:60]

for k in ['ISM','Halo','outflow','inflow']:
    all_data[k] = [None]*np.size(data_list)

x = dd.io.load(data_list[-1], '/gas_profiles/outflow/sphere')
all_data['R'] = x['centers_rvir']

for i, dname in enumerate(data_list):
    all_data['ISM'][i]  = dd.io.load( dname, '/meta_data/Z_avg')
    halo                = dd.io.load( dname, '/gas_meta_data/masses/Halo')
    all_data['Halo'][i] = halo['Metals'] / halo['Total']

    temp_outflow = dd.io.load( dname, '/gas_profiles/outflow/sphere')
    temp_inflow  = dd.io.load( dname, '/gas_profiles/inflow/sphere')
    all_data['outflow'][i] = temp_outflow['mass_profile'][('gas','metal_mass')] / temp_outflow['mass_profile'][('gas','cell_mass')]
    all_data['inflow'][i] = temp_inflow['mass_profile'][('gas','metal_mass')]  / temp_inflow['mass_profile'][('gas','cell_mass')]
for k in all_data.keys():
    all_data[k] = np.array(all_data[k])

### now plot ###
fig, ax = plt.subplots()
fig.set_size_inches(8,8)
times = times - times[0]
ax.plot( times, all_data['ISM'], lw = line_width, color = 'black', label = 'ISM')
ax.plot( times, all_data['Halo'], lw = line_width, color = 'C0',  label = 'Halo')
ax.plot( times, all_data['outflow'][:,0], lw = line_width, color = 'C1', label = 'Outflow', ls = '-') #
ax.plot( times, all_data['outflow'][:,2], lw = line_width, color = 'C1',  ls = '--')
ax.plot( times, all_data['outflow'][:,5], lw = line_width, color = 'C1',  ls = ':')

#ax.plot( times, all_data['inflow'][:,0], lw = line_width, color = 'C2', label = r'0.10 $R_{\rm vir}$', ls = '-') #
#ax.plot( times, all_data['inflow'][:,2], lw = line_width, color = 'C2', label = r'0.25 $R_{\rm vir}$', ls = '--')
#ax.plot( times, all_data['inflow'][:,5], lw = line_width, color = 'C2', label = r'1.00 $R_{\rm vir}$', ls = ':')


ax.set_ylabel(r'Metallicity')
ax.set_xlabel(r'Time (Myr)')
ax.set_xlim(times[0], times[-1])
ax.semilogy()
ax.set_ylim(1.0E-4, 5.0E-3)

plt.tight_layout()
plt.minorticks_on()
ax.legend(loc='best')
fig.savefig('metallicity_evolution.png')
