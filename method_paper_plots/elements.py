import yt
import numpy as np

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from galaxy_analysis.yt_fields import field_generators as fg
from galaxy_analysis.utilities import utilities

import galaxy_analysis as ga

yt.funcs.mylog.setLevel(50)

from joblib import Parallel, delayed
import multiprocessing
import deepdish as dd

temp = dd.io.load('DD0100_galaxy_data.h5', '/gas_meta_data/masses/CNM')
elements = utilities.sort_by_anum([x for x in temp.keys() if len(x) <= 2 and (not (x in ['H','He','H2','HI','HII','HeI']))])

print(np.size(elements), '-----------')


def plot_single_frames(dsi, fpath = '.'):

    fields = ['N_Number_Density','O_Number_Density','Fe_Number_Density','Ba_Number_Density'] #

    fields = ['O_over_Fe', 'Fe_over_H', 'N_over_H', 'O_over_H', 'N_over_O', 'Ba_over_H']
    ds = yt.load(fpath + '/DD%0004i/DD%0004i'%(dsi,dsi))
    fg.generate_derived_fields(ds)
    ds = yt.load(fpath + '/DD%0004i/DD%0004i'%(dsi,dsi))
    data = ds.all_data()

    data = ds.disk([0.5,0.5,0.5],[0,0,1],(2.0,'kpc'), (500,'pc'))

    proj = yt.ProjectionPlot(ds, 'z', center =[0.5,0.5,0.5], fields = fields, weight_field = 'Density', width = (1.0,'kpc'))
    proj.set_buff_size(2048)


    if True:
        for f in fields:
            proj.set_log(f, False)

            proj.set_cmap(f, 'plasma')

        proj.set_zlim('O_over_Fe', -1.5, 1.5)
        proj.set_zlim('Fe_over_H', -3, -1.5)
        proj.set_zlim('N_over_O',  -0.5, -1.0)
        proj.set_zlim('O_over_H',  -3, -1.5)
        for f in fields:
            label = '[' + f.split('_over_')[0]
            label += f.split('_over_')[1] + ']'
            proj.set_colorbar_label(f, label)

    else:
        proj.set_cmap('N_Number_Density', 'viridis')
        proj.set_cmap('O_Number_Density', 'magma')
        proj.set_cmap('Fe_Number_Density', 'plasma')
        proj.set_cmap('Ba_Number_Density', 'RdYlBu_r')

        for f in fields:
            fmin = np.min(proj.data_source[('gas',f)])
            fmax = np.max(proj.data_source[('gas',f)])
            plot_min = 1.0E-4 * fmax
            plot_max = fmax

            proj.set_zlim(f, plot_min, plot_max)  #1.0E9,1.0E16)
#    proj.annotate_quiver('velocity_y','velocity_z', 64)
#    proj.annotate_streamlines('velocity_y','velocity_z')


    proj.save('./elements3/')

def plot_slice_frames(dsi, fpath = '.'):

    fields = ['N_Fraction','O_Fraction','Fe_Fraction','Ba_Fraction', 'velocity_spherical_radius'] #
    ds = yt.load(fpath + '/DD%0004i/DD%0004i'%(dsi,dsi))
    fg.generate_derived_fields(ds)
    ds = yt.load(fpath + '/DD%0004i/DD%0004i'%(dsi,dsi))
    data = ds.all_data()

    proj = yt.SlicePlot(ds, 'z', center =[0.5,0.5,0.5], fields = fields, width = 2.0*yt.units.kpc)
    proj.set_buff_size(2048)

#    proj.set_cmap('N_N', 'viridis')
#    proj.set_cmap('O_Number_Density', 'magma')
#    proj.set_cmap('Fe_Number_Density', 'plasma')
#    proj.set_cmap('Ba_Number_Density', 'RdYlBu_r')

    for f in fields:
        if f == 'velocity_spherical_radius':
            proj.set_cmap(f, 'cubehelix')
            proj.set_log(f, False)
            proj.set_unit(f, 'km/s')
            proj.set_zlim(f, -400.0, 400.0)
        else:
            proj.set_cmap(f, 'magma')
#        fmin = np.min(proj.data_source[('gas',f)])
#        fmax = np.max(proj.data_source[('gas',f)])
            fmax = np.max( data[f][ data['spherical_radius'] < 1.0*yt.units.kpc])
            plot_min = 1.0E-4 * fmax
            plot_max = fmax
            proj.set_zlim(f, plot_min, plot_max)  #1.0E9,1.0E16)

#    proj.annotate_quiver('velocity_y','velocity_z', 64)
    proj.annotate_streamlines('velocity_y','velocity_z', density = 4)

    proj.save('./elements3/')


def plot_individual_panels(dsi, fpath = './', fration = False):

    fraction = True

    fields = ['H_p0_number_density','H2_p0_number_density',
              'O_Number_Density','Mg_Number_Density','Fe_Number_Density']
    try:
        ds = yt.load(fpath + '/DD%0004i/DD%0004i'%(dsi,dsi))
    except:
	print('load failed on ', dsi)
        return

    fg.generate_derived_fields(ds)
    ds = yt.load(fpath + '/DD%0004i/DD%0004i'%(dsi,dsi))
    data = ds.all_data()

    fig = plt.figure()
    grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (2, 5),
                axes_pad = 0.7,
                label_mode = "1",
                share_all = False,
                cbar_location="right",
                cbar_mode="each",
                cbar_size="3%",
                cbar_pad="0%")

    sp = yt.ProjectionPlot(ds, 'x', center = [0.5,0.5,0.5],
                           fields = fields, weight_field = None,
                           width = (6.8,'kpc'))

    disk = ds.disk([0.5,0.5,0.5],[0,0,1],(1,'kpc'), (75,'pc'))
    sp2 = yt.ProjectionPlot(ds, 'z',center = [0.5,0.5,0.5],
                            fields = fields, weight_field = None,
                            width = (1.5,'kpc'), data_source = disk)
    disk2 = ds.disk([0.5,0.5,0.5],[0,0,1],(1,'kpc'), (50,'pc'))

    sp3 = yt.ProjectionPlot(ds, 'z', center = [0.5,0.5,0.5],
                            fields=['O_Fraction','Mg_Fraction','Fe_Fraction'],
                            weight_field = 'density',
                            width = (1.5,'kpc'), data_source =disk2)

    sp.set_buff_size(512)
    sp2.set_buff_size(512)
    sp3.set_buff_size(512)

    sp3.set_cmap('O_Fraction','magma')
    sp3.set_cmap('Mg_Fraction','magma')
    sp3.set_cmap('Fe_Fraction','magma')

    for x in ['O','Mg','Fe']:
        sp3.set_zlim(x + '_Fraction', 1.0E-6, 1.0E-3)

    for x in [sp,sp2]:
        x.set_cmap('H_p0_number_density', 'Greys')
        x.set_cmap('H2_p0_number_density', 'arbre')
        x.set_cmap('O_Number_Density','magma')
        x.set_cmap('Mg_Number_Density','magma')
        x.set_cmap('Fe_Number_Density','magma')

        x.set_zlim('H_p0_number_density', 1.0E16, 1.0E22)
        x.set_zlim('H2_p0_number_density',1.0E16, 1.0E20)
        x.set_zlim('O_Number_Density',1.0E11, 1.0E16)
        x.set_zlim('Mg_Number_Density',1.0E11, 1.0E16)
        x.set_zlim('Fe_Number_Density',1.0E11,1.0E16)

        x.set_colorbar_label('H_p0_number_density',r'HI Column (cm$^{-2}$)')
        x.set_colorbar_label('H2_p0_number_density',r'H$_2$ Column (cm$^{-2}$)')
        x.set_colorbar_label('O_Number_Density',r'O Column (cm$^{-2}$)')
       	x.set_colorbar_label('Mg_Number_Density',r'Mg Column (cm$^{-2}$)')
       	x.set_colorbar_label('Fe_Number_Density',r'Fe Column (cm$^{-2}$)')



    for x in [sp,sp2,sp3]:
        x.set_font({'family': 'sans-serif','size': 10})

    for i, field in enumerate(fields):
        plot = sp.plots[field]
        plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]

    if not fraction:
        for i, field in enumerate(fields):
            plot = sp2.plots[field]
            plot.figure = fig
            plot.axes = grid[i + len(fields)].axes
            plot.cax = grid.cbar_axes[i + len(fields)]
    else:
        for i, field in enumerate(['H_p0_number_density','H2_p0_number_density']):
            plot = sp2.plots[field]
            plot.figure = fig
            plot.axes = grid[i + len(fields)].axes
            plot.cax = grid.cbar_axes[i + len(fields)]

        for i, field in enumerate(['O_Fraction','Mg_Fraction','Fe_Fraction']):
            plot = sp3.plots[field]
            plot.figure = fig
            plot.axes = grid[i + len(fields) + 2].axes
            plot.cax = grid.cbar_axes[i + len(fields) + 2]

    sp._setup_plots()
    sp2._setup_plots()
    sp3._setup_plots()
#    sp2.hide_colorbar()
#    sp3.hide_colorbar()

    mstar = 0.0
    if ('io','particle_mass') in ds.field_list:
        mstar = np.sum(data['particle_mass'][ data['particle_type'] == 11].convert_to_units('Msun').value)
    fig.suptitle(r"Time = %1.1f Myr       -       M$_*$ = %2.2E M$_{\odot}$"%(ds.current_time.convert_to_units('Myr').value, mstar), fontsize=10)
    fig.set_size_inches(11,9)

    plt.savefig(fpath + '/elements2/DD%0004i_elements.png'%(dsi))



def plot_element_panels(dsi, fpath = './', axis = 'x', width = (14.35,'kpc'),
                        cmap = 'magma'):
    try:
        ds = yt.load(fpath + '/DD%0004i/DD%0004i'%(dsi,dsi))
    except:
        print('load failed on ', dsi)
        return

    fg.generate_derived_fields(ds)
    ds = yt.load(fpath + '/DD%0004i/DD%0004i'%(dsi,dsi))
    
    #
    # plot 8 elements: HI, Z, C, N, O, Mg, Fe, Ni, Eu
    #
    fig = plt.figure()
    #1 2 6 7 8 12 14 16 20 25 26 28 39 56 63  
#    metals = [ 'C', 'N','O','Mg','Si','S','Ca','Mn','Fe','Ni','Ba','Y','Eu']
    metals = elements
    
    fields = ['Metal_Number_Density']
    for ele in metals:
        fields += [ele + '_Number_Density']
        
    grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (4, 4),
                axes_pad = 1.0,
                label_mode = "1",
                share_all = True,
                cbar_location="right",
                cbar_mode="each",
                cbar_size="3%",
                cbar_pad="0%")
    
    sp = yt.ProjectionPlot(ds, axis, center = [0.5,0.5,0.5],
                           fields = fields, weight_field = None,
                           width = width)
    for f in fields:
        sp.set_cmap(f, 'viridis')
        sp.set_unit(f, 'cm**(-2)')
        sp.set_zlim(f, 1.0E16, 1.0E21)
        
    sp.set_cmap('Metal_Number_Density', 'arbre')
    sp.set_zlim('Metal_Number_Density', 1.0E13, 1.0E18)


    for met in metals:
        f = met + '_Number_Density'
        fmin = np.min(sp.data_source[('gas',f)])
        fmax = np.max(sp.data_source[('gas',f)])
        plot_min = 1.0E-7 * fmax
        plot_max = fmax

        sp.set_cmap(f, cmap)
        sp.set_zlim(f, plot_min, plot_max)  #1.0E9,1.0E16)

#    for met in ['Ba','Y','Sr']:
#        f = met + '_Number_Density'
#        sp.set_cmap(f, 'YlOrRd')
#        sp.set_zlim(f, 1.0E4, 1.0E9)

    sp.set_buff_size(512)

    for i, field in enumerate(fields):
        plot = sp.plots[field]
        plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
        
    sp._setup_plots()
    fig.set_size_inches(32,32)
    
    plt.savefig(fpath + '/elements/DD%0004i_elements.png'%(dsi))
    
    return
    
    


#for i in np.arange(64, 231, 1):
#    try:
#        ds = yt.load('DD%0004i/DD%0004i'%(i,i))
#    except:
#        print 'file ', i, ' does not load'
#        continue
#plot_single_frames(401)
for i in [200, 400, 500, 550]:
    plot_slice_frames(i)
    plot_single_frames(i)

#n_jobs = 28
#Parallel(n_jobs = n_jobs)(\
#    delayed(plot_element_panels)(i,'./') for i in np.arange(481,532,1))
#    delayed(plot_element_panels)(i,'./') for i in np.arange(350, 401, 2))

#    plot_element_panels(i, './')
