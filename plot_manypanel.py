import yt
import numpy as np
import matplotlib.pyplot as plt
#
#
from galaxy_analysis.analysis import Galaxy
import sys

imin = 167
imax = 168

if len(sys.argv) > 1:
    imin = int(sys.argv[1])
    imax = int(sys.argv[2]) 

for i in np.arange(imin, imax, 1):
	dsname = 'DD%0004i'%(i)
	gal  = Galaxy(dsname)
	ds   = gal.ds
	data = ds.all_data()
	disk = ds.disk(ds.domain_center, [0,0,1], (800,'pc'), (300,'pc'))
	x    = np.array([800,800,800]) * yt.units.pc
	reg_vert  = ds.region(ds.domain_center, ds.domain_center - x, ds.domain_center+x)
	x    = np.array([800,800,15]) * yt.units.pc
	reg_plane = ds.region(ds.domain_center, ds.domain_center - x, ds.domain_center+x)
	width = (800,'pc')



	text_args={'color':'white', 'size':40}
	inset_box_args={'boxstyle':'square,pad=0.3',
		                         'facecolor':'black',
		                         'linewidth':1,
		                         'edgecolor':'black', 'alpha':0.3}
	xyannotate=(-350,-350)





	field = 'G_o'
	pp = yt.ProjectionPlot(ds, 'z', field, weight_field ='number_density', 
		               data_source = reg_plane, width = width)
	pp.set_cmap(field,'cubehelix')
	pp.set_zlim(field,0.05, 80)
	pp.annotate_text(xyannotate, r'ISRF (G$_o$)', 
		         text_args=text_args,
		         coord_system='plot',inset_box_args=inset_box_args)
	#pp.show()
	#pp.annotate_particles(0.9)
	#pp.show()
	pp.hide_axes()
	pp.hide_colorbar()
	pp.save("multipanel/" + dsname + '_' + field + '_face_on_square.png',mpl_kwargs={'dpi':100})
	#
	#
	#
	#
	field = 'Q0_flux'
	pp = yt.ProjectionPlot(ds, 'z', field, weight_field ='number_density', 
		               data_source = reg_plane, width = width)
	pp.set_cmap(field,'magma')
	pp.set_zlim(field, 1.0E-6, 1.0E-1)
	pp.annotate_text(xyannotate, 'HI Ionizing Radiation', 
		         text_args=text_args,
		         coord_system='plot',inset_box_args=inset_box_args)
	#pp.annotate_particles(0.9)
	#pp.show()
	pp.hide_axes()
	pp.hide_colorbar()
	pp.save("multipanel/" + dsname + '_' + field + '_face_on_square.png',mpl_kwargs={'dpi':100})
	#
	#
	#
	#
	field = 'Q1_flux'
	pp = yt.ProjectionPlot(ds, 'z', field, weight_field ='number_density', 
		               data_source = reg_plane, width = width)
	pp.set_cmap(field,'magma')
	pp.set_zlim(field, 1.0E-6, 1.0E-1)
	pp.annotate_text(xyannotate, 'HeI Ionizing Radiation', 
		         text_args=text_args,
		         coord_system='plot',inset_box_args=inset_box_args)
	#pp.annotate_particles(0.9)
	#pp.show()
	pp.hide_axes()
	pp.hide_colorbar()
	pp.save("multipanel/" + dsname + '_' + field + '_face_on_square.png',mpl_kwargs={'dpi':100})
	#
	#
	#
	#
	field = 'number_density'
	pp = yt.ProjectionPlot(ds, 'z', field, weight_field ='number_density', 
		               data_source = reg_plane, width = width)
	pp.set_cmap(field,'viridis')
	pp.set_zlim(field, 0.001, 100.0)
	pp.annotate_text(xyannotate, 'Number Density', 
		         text_args=text_args,
		         coord_system='plot',inset_box_args=inset_box_args)
	#pp.annotate_particles(0.9)
	#pp.show()
	pp.hide_axes()
	pp.hide_colorbar()
	pp.save("multipanel/" + dsname + '_' + field + '_face_on_square.png',mpl_kwargs={'dpi':100})
	#
	#
	#

	#
	#
	#
	field = 'Temperature'
	pp = yt.ProjectionPlot(ds, 'z', field, weight_field ='number_density', 
		               data_source = reg_plane, width = width)
	pp.set_cmap(field,'RdYlBu_r')
	pp.set_zlim(field, 1.0, 3.0E7)
	#pp.annotate_particles(0.9)
	pp.annotate_text(xyannotate, field, 
		         text_args=text_args,
		         coord_system='plot',inset_box_args=inset_box_args)
	#pp.annotate_particles(0.9)
	#pp.show()
	pp.hide_axes()
	pp.hide_colorbar()
	pp.save("multipanel/" + dsname + '_' + field + '_face_on_square.png',mpl_kwargs={'dpi':100})
	#
	#
	#


	field = 'number_density'
	pp = yt.ProjectionPlot(ds, 'z', field, weight_field ='number_density', 
		               data_source = reg_plane, width = width)
	pp.set_cmap(field,'viridis')
	pp.set_zlim(field, 0.001, 100.0)
	#pp.show()

	xpos = (reg_plane['particle_position_x'] - ds.domain_center[0]).to('pc').value
	ypos = (reg_plane['particle_position_y'] - ds.domain_center[1]).to('pc').value
	m    = (reg_plane['birth_mass'].value)

	pp.set_buff_size(800)

	pp.annotate_text(xyannotate, 'Number Density', 
		         text_args=text_args,
		         coord_system='plot',inset_box_args=inset_box_args)
	pp[field].image.axes.scatter( xpos[m>8.0], ypos[m>8.0], s = 80, 
		                      color = 'C1', marker = '*')

	pp.hide_axes()
	pp.hide_colorbar()
	pp.save("multipanel/" + dsname + '_' + field + '_particles_face_on_square.png',mpl_kwargs={'dpi':100})
	#pp[field].image.write_png(field + '_particles_face_on_square.png')


	field = 'Temperature'
	pp = yt.ProjectionPlot(ds, 'z', field, weight_field ='number_density', 
		               data_source = reg_plane, width = width)
	pp.set_cmap(field,'RdYlBu_r')
	pp.set_zlim(field, 10.0,1.0E7)
	#pp.show()

	xpos = (reg_plane['particle_position_x'] - ds.domain_center[0]).to('pc').value
	ypos = (reg_plane['particle_position_y'] - ds.domain_center[1]).to('pc').value
	m    = (reg_plane['birth_mass'].value)

	pp.set_buff_size(800)
	pp.annotate_text(xyannotate, 'Temperature', 
		         text_args=text_args,
		         coord_system='plot',inset_box_args=inset_box_args)
	pp[field].image.axes.scatter( xpos[m>8.0], ypos[m>8.0], s = 80, 
		                      color = 'black', marker = '*')
	pp.hide_axes()
	pp.hide_colorbar()
	pp.save("multipanel/" + dsname + '_' + field + '_particles_face_on_square.png',mpl_kwargs={'dpi':100})
	#pp[field].image.write_png(field + '_particles_face_on_square.png')


	#
	#
	#
	#
	#
	#

	field = 'O_over_H'
	pp = yt.ProjectionPlot(ds, 'z', field, weight_field ='number_density', 
		               data_source = reg_plane, width = width)
	pp.set_cmap(field,'plasma')
	pp.set_log(field, False)
	pp.set_zlim(field, -6, 1)
	#pp.annotate_particles(0.9)
	pp.annotate_text(xyannotate, '[O/H]', 
		         text_args=text_args,
		         coord_system='plot',inset_box_args=inset_box_args)
	#pp.annotate_particles(0.9)
	#pp.show()
	pp.hide_axes()
	pp.hide_colorbar()
	pp.save("multipanel/" + dsname + '_' + field + '_face_on_square.png',mpl_kwargs={'dpi':100})
	#
	#
	#
	#
	field = 'N_over_O'
	pp = yt.ProjectionPlot(ds, 'z', field, weight_field ='number_density', 
		               data_source = reg_plane, width = width)
	pp.set_cmap(field,'PRGn')
	pp.set_log(field, False)
	pp.set_zlim(field, -2., 2.)
	#pp.annotate_particles(0.9)
	pp.annotate_text(xyannotate, '[N/O]',
		         text_args=text_args,
		         coord_system='plot',inset_box_args=inset_box_args)
	#pp.annotate_particles(0.9)
	#pp.show()
	pp.hide_axes()
	pp.hide_colorbar()
	pp.save("multipanel/" + dsname + '_' + field + '_face_on_square.png',mpl_kwargs={'dpi':100})
	#
	#
	#
	#
	field = 'Mg_over_O'
	pp = yt.ProjectionPlot(ds, 'z', field, weight_field ='number_density', 
		               data_source = reg_plane, width = width)
	pp.set_cmap(field,'PuOr')
	pp.set_log(field, False)
	pp.set_zlim(field, -2,2)
	#pp.annotate_particles(0.9)
	pp.annotate_text(xyannotate, '[Mg/O]',
		         text_args=text_args,
		         coord_system='plot',inset_box_args=inset_box_args)
	#pp.annotate_particles(0.9)
	#pp.show()
	pp.hide_axes()
	pp.hide_colorbar()
	pp.save("multipanel/" + dsname + '_' + field + '_face_on_square.png',mpl_kwargs={'dpi':100})
	#
	#
	#
	#
	field = 'O_over_Fe'
	pp = yt.ProjectionPlot(ds, 'z', field, weight_field ='number_density', 
		               data_source = reg_plane, width = width)
	pp.set_cmap(field,'BrBG')
	pp.set_log(field, False)
	pp.set_zlim(field, -2,2)
	#pp.annotate_particles(0.9)
	pp.annotate_text(xyannotate, '[O/Fe]',
		         text_args=text_args,
		         coord_system='plot',inset_box_args=inset_box_args)
	#pp.annotate_particles(0.9)
	#pp.show()
	pp.hide_axes()
	pp.hide_colorbar()
	pp.save("multipanel/" + dsname + '_' + field + '_face_on_square.png',mpl_kwargs={'dpi':100})
	#
	#
	#
	#
	field = 'a_rad_over_a_grav'
	pp = yt.ProjectionPlot(ds, 'z', field, weight_field ='number_density', 
		               data_source = reg_plane, width = width)
	pp.set_cmap(field,'RdBu_r')
	#pp.set_log(field, False)
	pp.set_zlim(field, 1.0E-3, 1.0E3)
	#pp.annotate_particles(0.9)
	pp.annotate_text(xyannotate, 'Rad. Accel / Self-gravity',
		         text_args=text_args,
		         coord_system='plot',inset_box_args=inset_box_args)
	#pp.annotate_particles(0.9)
	#pp.show()
	pp.hide_axes()
	pp.hide_colorbar()
	pp.save("multipanel/" + dsname + '_' + field + '_face_on_square.png',mpl_kwargs={'dpi':100})
	#
	#
	#
	field = 'H2_p0_number_density'
	pp = yt.ProjectionPlot(ds, 'z', field, weight_field = None, 
		               data_source = disk, width = width)
	pp.set_cmap(field,'plasma')
	pp.set_unit(field, 'cm**(-2)')
	pp.set_colorbar_label(field, r'H$_2$ Column Density')
	pp.set_buff_size(800)
	pp.set_zlim(field, 1.0E14, 1.0E19)
	#pp.set_zlim(field, 1.0, 3.0E7)
	pp.annotate_text(xyannotate, 'H$_2$ Column Density', 
		         text_args=text_args,
		         coord_system='plot',inset_box_args=inset_box_args)
	#pp.annotate_particles(0.9)
	#pp.show()
	pp.hide_axes()
	pp.hide_colorbar()
	pp.save("multipanel/" + dsname + '_' + field + '_face_on_square.png',mpl_kwargs={'dpi':100})
	#
	#
	#
	field = 'H_p0_number_density'
	pp = yt.ProjectionPlot(ds, 'z', field, weight_field = None, 
		               data_source = disk, width = width)
	pp.set_cmap(field,'Greys')
	pp.set_unit(field, 'cm**(-2)')
	pp.set_colorbar_label(field, 'HI Column Density')
	pp.set_buff_size(800)
	pp.set_zlim(field, 1.0E16, 1.0E22)
	pp.annotate_text(xyannotate, 'HI Column Density', 
		         text_args=text_args,
		         coord_system='plot',inset_box_args=inset_box_args)
	#pp.annotate_particles(0.9)
	#pp.show()
	pp.hide_axes()
	pp.hide_colorbar()
	pp.save("multipanel/" + dsname + '_' + field + '_face_on_square.png',mpl_kwargs={'dpi':100})

        del(pp)
        del(gal)

