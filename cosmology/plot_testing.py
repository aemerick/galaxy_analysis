import matplotlib
matplotlib.use('agg')

import yt
import numpy as np
import glob
import sys
from galaxy_analysis.yt_fields import field_generators as fg



yt.enable_parallelism()

#if len(sys.argv) > 1:
#    imin = int(sys.argv[1])
#    imax = int(sys.argv[2])
#else:
imin = 0
imax = -1


dsfiles = np.sort(glob.glob("DD????/DD????"))


fields = [('gas','number_density'),'Temperature',
          ('gas','photo_gamma'),
          ('gas','metallicity'),
          ('gas','C_over_N'),
          ('gas','C_over_H'),
          ('gas','O_over_H'),
          ('gas','Fe_over_H')]

#
# fields = [('gas','number_density')]
#

def format_plot(ds, sp, fields):

    if ('gas', 'C_over_N') in fields:
        sp.set_log(('gas','C_over_N'),False)
        sp.set_zlim(('gas','C_over_N'),-3,3)
        sp.set_cmap(('gas','C_over_N'),"PRGn")

    if ('gas','O_over_H') in fields:
        sp.set_log(('gas','O_over_H'),False)
        sp.set_zlim(('gas','O_over_H'), -8,0)

    if ('gas','C_over_H') in fields:
        sp.set_log(('gas','C_over_H'),False)
        sp.set_zlim(('gas','C_over_H'), -8,0)

    if ('gas','Fe_over_H') in fields:
        sp.set_log(('gas','Fe_over_H'),False)
        sp.set_zlim(('gas','Fe_over_H'),-8,0)

    if 'HM_kph' in fields:
        sp.set_unit("HM_kph","1/s")
        sp.set_zlim("HM_kph",1.0E-17,1.0E-7)

    if ('enzo','HI_kph') in fields:
        sp.set_unit("HI_kph","1/s")
        sp.set_zlim("HI_kph",1.0E-16,1.0E-6)

    if ('gas','photo_gamma') in fields:
        sp.set_unit(('gas','photo_gamma'),'erg/s')
        sp.set_zlim(('gas','photo_gamma'),1.0E-26,1.0E-20)
        sp.set_cmap(('gas','photo_gamma'),'cubehelix')

    if 'FUV_FluxDensity' in fields:
        sp.set_zlim("FUV_FluxDensity",1.0E-6,0.06)

    if 'Temperature' in fields:
        sp.set_zlim('Temperature',10.0,1.0E7)
        sp.set_cmap('Temperature', 'RdYlBu_r')

    if ('gas','number_density') in fields:
        sp.set_zlim(('gas','number_density'),1.0E-5, 1.0E3)
        sp.set_unit(('gas','number_density'),'(cm)**(-3)')
        sp.set_cmap(('gas','number_density'),'viridis')

    if ('gas','metallicity') in fields:
        sp.set_zlim(('gas','metallicity'), 1.0E-7,1.0)

#    sp.set_buff_size(512)

    sp.annotate_particles(0.9,ptype='all_stars')


    sp.annotate_title("t = %.1f Myr - z = %.2f"%( ds.current_time.to('Myr'), ds.current_redshift))


    return


weight_field = 'metallicity'

for dsname in dsfiles[imin:imax]:
    ds = yt.load(dsname)
    fg.generate_derived_fields(ds)
    ds = yt.load(dsname)
    fg.generate_particle_filters(ds)

    data = ds.all_data()

    px = data['particle_position_x'].to('kpc').value
    py = data['particle_position_y'].to('kpc').value
    pz = data['particle_position_z'].to('kpc').value
    pid = data['particle_index']

    select = (pid == 100)

    pcenter = np.array([px[select][0], py[select][0], pz[select][0]]) * yt.units.kpc

    width= (data.right_edge[0].to('kpc') - data.left_edge[0].to('kpc'))
    sp = yt.ProjectionPlot(ds, 'z', fields, center=data.center, width=width,
                           weight_field=weight_field, data_source=data)
    format_plot(ds,sp,fields)
    if yt.is_root():
        sp.save("./proj/")
    del(sp)

    #
    # now make new regions
    #
    for factor in [4,8,16,32]:
        region = ds.region( pcenter, pcenter - width/(factor*2.0), pcenter + width / (factor*2.0))
        sp = yt.ProjectionPlot(ds, 'z', fields, center = region.center, width=width/(1.0*factor),
                               weight_field=weight_field, data_source=region)

        format_plot(ds,sp,fields)

        if yt.is_root():
            sp.save("./proj/zoomx%i/"%(factor))
        del(sp)

    del(ds)
