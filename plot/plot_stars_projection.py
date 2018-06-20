
import yt
from yt import derived_field

import numpy as np
import matplotlib.pyplot as plt
from galaxy_analysis import Galaxy
import glob
import sys, os

#
# As Copied from stack overflow with some minor renaming edits:
#
import matplotlib.colors as mpl_colors
class MidpointNormalize(mpl_colors.Normalize):
    """
    Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

    e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        mpl_colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))
# ------------------
# end copy
# ------------------

#
# Some globals for convenience
#
colors     = {'massive_star_winds' : 'white', 'AGB_winds' : 'C1', 'SN' : 'black', 'other_stars' : 'black'}
markers    = {'massive_star_winds' : '*',     'AGB_winds' : 'D', 'SN' : '*', 'other_stars' : '.'}
ps         = {'massive_star_winds' :  75, 'AGB_winds' : 300, 'SN' : 300, 'other_stars' : 15}
all_fields = ['number_density','Temperature','N_over_O_filtered','logNO_filtered','O_over_H_filtered']

tol      = 1.0E-25

def make_filtered_field(ds, fieldname, filter_fields = [], tolerance = tol):
    """
    Creates a derived field that filters 'fieldname' if any field in
    'filter_fields' has a value below 'tolerance'. Filtered by setting
    cells where this is true to np.nan (i.e. they will show up as
    white / blank in projection / slice plots).

    Derived field will have the name  : ('gas',fieldname + "_filtered")
    Assumes that the original field is: ('gas',fieldname)
    """
    def _filtered_field(field, data):
        x = data[('gas',fieldname)]

        select = data[filter_fields[0]] < 0
        for f in filter_fields:
            select = select + (data[f] < tolerance)
        x[select] = np.nan

        return x

    ds.add_field(('gas',fieldname + '_filtered'), function = _filtered_field, units = "")
    return

def plot(dsname, wdir = './', width = 500.0, dt = 5.0*yt.units.Myr, fields = all_fields,
         thickness = 20.0, outdir = './enrichment_plots'):
    """
    Given a datasetname, compute thin projections with different types of stars
    using different points. Meant to illustrate evolution of metal enrichment from various
    stars.
    """

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    gal = Galaxy(dsname, wdir = wdir)
    data = gal.df

    @derived_field(name="logNO", units="")
    def _logNO(field, data):
        return  np.log10(data['N_Abundance'] / data['O_Abundance'])
    gal.ds.add_field(("gas", "logNO"), function=_logNO, units="")

    make_filtered_field(gal.ds, 'logNO', ['O_Fraction','N_Fraction'])
    make_filtered_field(gal.ds, 'O_over_H', ['O_Fraction'])
    make_filtered_field(gal.ds, 'N_over_O', ['O_Fraction','N_Fraction'])
#    def _logNO_filtered(field,data):
#        x = data[('gas','logNO')]
#
#        f1 = data[('gas','O_Fraction')]
#        f2 = data[('gas','N_Fraction')]
#
#        x[ (f1 < tol) + (f2 < tol)] = np.nan
#
#        return x
#    gal.ds.add_field(('gas','logNO_filtered'), function = _logNO_filtered, units = "")

    M            = data['birth_mass']
    t_o          = data['creation_time'].convert_to_units('Myr')
    MS_lifetime  = data[('io','particle_model_lifetime')].to('Myr')
    MS_death     = t_o + MS_lifetime
    px           = (data['particle_position_x'] - gal.ds.domain_center[0]).to('pc')
    py           = (data['particle_position_y'] - gal.ds.domain_center[1]).to('pc')
    pz           = (data['particle_position_z'] - gal.ds.domain_center[2]).to('pc')

    recent_death = (MS_death > gal.ds.current_time - dt) * (MS_death <= gal.ds.current_time + 0.001*yt.units.Myr)
    alive        = MS_death > gal.ds.current_time + 0.001*yt.units.Myr

    AGB           = M < 8.0
    massive_star  = (M > 8.0) * (M < 25.0)

    boxdim = np.array([width*1.25,width*1.25,thickness])*yt.units.pc
    region = gal.ds.box(gal.ds.domain_center - boxdim*0.5, gal.ds.domain_center + boxdim*0.5)

    proj = yt.ProjectionPlot(gal.ds, 'z', fields,
                        weight_field = 'number_density', data_source = region, width = (width,'pc'))

    if 'number_density' in fields:
        proj.set_unit('number_density','cm**(-3)')
        proj.set_cmap('number_density','viridis')
        proj.set_zlim('number_density',1.0E-4,200.0)

    if 'O_over_H_filtered' in fields:
        proj.set_cmap('O_over_H_filtered','cubehelix')
        proj.set_log('O_over_H_filtered', False)
        proj.set_zlim('O_over_H_filtered', -5, 1)
        proj.set_colorbar_label('O_over_H_filtered', r'[O/H]')

    if 'N_over_O_filtered' in fields:
        proj.set_cmap('N_over_O_filtered','PRGn')
        proj.set_log('N_over_O_filtered',False)
        proj.set_zlim('N_over_O_filtered',-2,2)
        proj.set_colorbar_label('N_over_O_filtered', r'[N/O]')

    if 'logNO' in fields:
        proj.set_cmap('logNO','PRGn')
        proj.set_log('logNO',False)
        proj.set_zlim('logNO',-2,0.5)
        proj.set_colorbar_label('logNO', r'log( N / O )')

    if 'logNO_filtered' in fields:
        proj.set_cmap('logNO_filtered','PRGn')
        proj.set_log('logNO_filtered',False)
        proj.set_zlim('logNO_filtered',-2,0.5)
        proj.set_colorbar_label('logNO_filtered', r'log( N / O )')

    if 'Temperature' in fields:
        proj.set_cmap('Temperature', 'RdYlBu_r')
        proj.set_log('Temperature',True)
        proj.set_zlim('Temperature',10.0, 1.0E7)
        proj.set_colorbar_label('Temperature', r'Temperature (K)')

    proj.save(outdir + '/') # necessary


    dt = 5.0 * yt.units.Myr
    in_image     = (np.abs(pz) <= boxdim[2]*0.5) * (np.abs(px) <= width*0.5) * (np.abs(py) <= width*0.5)

    pp = {}
    pp['massive_star_winds'] = in_image * alive * massive_star
    pp['AGB_winds']          = in_image * recent_death * AGB
    pp['SN']                 = in_image * recent_death * massive_star
    #pp['other_stars']        = in_image * alive * (np.logical_not(pp['massive_star_winds']))

    for k in proj.plots.keys():
        image = proj.plots[k]

    #
    # Now select and annotate the points we want
    #
        for s in pp.keys():
            if np.size(px[pp[s]].value) > 0:
                print np.size(px[pp[s]]), 'Particles in ', s, px[pp[s]], py[pp[s]]
                image.axes.scatter(px[pp[s]].value,py[pp[s]].value, s = ps[s], marker = markers[s], color = colors[s])
            else:
                print 'No particles in ', s

#    proj.refresh()
    proj.save(outdir + '/') # necessary

    if 'N_over_O' in fields:
        vmin,vmax = -2,2
        x = proj.plots['N_over_O']
        x.image.set_norm( MidpointNormalize(midpoint= 0.5*(vmin+vmax), vmin=vmin,vmax=vmax))
        x.cb.set_norm(MidpointNormalize(midpoint=0.5*(vmin+vmax),vmin=vmin,vmax=vmax))
        x.cb.update_normal(x.image)
        x.save(outdir + '/' + str(gal.ds) + '_Projection_z_N_over_O_number_density.png')

    if 'logNO' in fields:
        vmin, vmax = -2, 0.25
        x = proj.plots['logNO']
        x.image.set_norm( MidpointNormalize(midpoint= 0.0, vmin=vmin,vmax=vmax))
        x.cb.set_norm(MidpointNormalize(midpoint=0.0, vmin=vmin,vmax=vmax))
        x.cb.update_normal(x.image)
        x.save(outdir + '/' + str(gal.ds) + '_Projection_z_logNO_number_density.png')

    del(proj)
    del(gal)

    return


if __name__ == "__main__":

    all_ds = np.sort(glob.glob('DD????/DD????'))

    if len(sys.argv) > 1:
        start = int(sys.argv[1])
        fin   = int(sys.argv[2])

        if len(sys.argv) == 4:
            di = int(sys.argv[3])
        else:
            di = 1

        for dsi in np.arange(start, fin, di):
            dsname = "DD%0004i"%(dsi)
            print "Plotting for ", dsname
            plot(dsname)
    else:
        plot(all_ds[0].split('/')[0])

