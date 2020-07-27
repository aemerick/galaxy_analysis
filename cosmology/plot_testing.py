import matplotlib
matplotlib.use('agg')

import yt
import numpy as np
import glob
import sys, os
from galaxy_analysis.yt_fields import field_generators as fg

import h5py
import deepdish as dd


def get_output_H5(outdata_path):
    if os.path.isfile(outdata_path):
        outH5 = dd.io.load(outdata_path)
    else:
        outH5 = {}

    return outH5

def save_to_file(outdata_path, sp, data, center):

    outH5 = get_output_H5(outdata_path)

    for f in sp.frb.data.keys():
        if not f[1] in outH5.keys():
            outH5[f[1]] = {}
        outH5[f[1]][weight_field] = sp.frb.data[f]

    if not 'particle_positions' in outH5.keys():
        outH5['particle_positions'] = {}

    for i,ax in enumerate(['x','y','z']):
        outH5['particle_positions'][ax] = (data[('all_stars','particle_position_'+ax)].to('kpc') - center.to('kpc')[i])/sp.width[0].to('kpc')

    dd.io.save(outdata_path, outH5)

    return outH5



save_frb = True

yt.enable_parallelism()

weight_field = "number_density"

imin,imax = None, None
if len(sys.argv) > 1:
    imin = int(sys.argv[1])
    imax = int(sys.argv[2])
if len(sys.argv) > 3:
    if str(sys.argv[3]) != "&":
        weight_field = str(sys.argv[3])

if imin is None:
    dsfiles = [np.sort(glob.glob("DD????/DD????"))[-1]]
else:
    dsfiles = ["DD%0004i/DD%0004i"%(i,i) for i in np.arange(imin,imax+1,1)]

if yt.is_root():
    print(dsfiles)


fields = [('gas','number_density'),'Temperature',
          ('gas','photo_gamma'),
          ('gas','metallicity'),
          ('gas','C_over_N'),
          ('gas','C_over_H'),
          ('gas','O_over_H'),
          ('gas','Fe_over_H')]

fields = [('gas','number_density'),'Temperature',('gas','O_over_H'),('gas','photo_gamma')]

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

for dsname in dsfiles:

    ds = yt.load(dsname)
    fg.generate_derived_fields(ds)
    ds = yt.load(dsname)
    fg.generate_particle_filters(ds)

    data = ds.all_data()

#    px = data['particle_position_x'].to('kpc').value
#    py = data['particle_position_y'].to('kpc').value
#    pz = data['particle_position_z'].to('kpc').value
#    pid = data['particle_index']

#    select = (pid == 16777372)
#
#    if np.size(pid[select]) == 0:
#        select = (pid == 16777218)

#    pcenter = np.array([px[select][0], py[select][0], pz[select][0]]) * yt.units.kpc

    pcenter = np.array([44.04561761 ,13.03435496, 10.58975735]) * yt.units.kpc


    width= (data.right_edge[0].to('kpc') - data.left_edge[0].to('kpc'))
    sp = yt.ProjectionPlot(ds, 'z', fields, center=data.center, width=width,
                           weight_field=weight_field, data_source=data)
    format_plot(ds,sp,fields)
    if yt.is_root():
        sp.save("./proj/")

    if yt.is_root() and save_frb:
        outdata_path = './proj/' + str(ds) + '_Proj_data.h5'
        save_to_file(outdata_path, sp, data, data.center)

    del(sp)

    #
    # now make new regions
    #
    for factor in [2,4]:
        region = ds.region( pcenter, pcenter - width/(factor*2.0), pcenter + width / (factor*2.0))
        sp = yt.ProjectionPlot(ds, 'z', fields, center = region.center, width=width/(1.0*factor),
                               weight_field=weight_field, data_source=region)

        format_plot(ds,sp,fields)

        if yt.is_root():
            sp.save("./proj/zoomx%i/"%(factor))

        if yt.is_root() and save_frb:
            outdata_path = "./proj/zoomx%i/"%(factor) + str(ds) + '_Proj_data_zoomx%i.h5'%(factor)
            save_to_file(outdata_path, sp, region, region.center)

        del(sp)

    del(ds)
