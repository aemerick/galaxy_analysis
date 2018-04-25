
import yt
from yt import derived_field

import numpy as np
import matplotlib.pyplot as plt
from galaxy_analysis import Galaxy
import glob
import sys

######################
# set the colormap and centre the colorbar
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





#######################


# In[79]:


colors   = {'massive_star_winds' : 'black', 'AGB_winds' : 'C1', 'SN' : 'C4', 'other_stars' : 'white'}
markers  = {'massive_star_winds' : '*',     'AGB_winds' : 'D', 'SN' : '*', 'other_stars' : '.'}
ps       = {'massive_star_winds' : 50, 'AGB_winds' : 200, 'SN' : 200, 'other_stars' : 15}


# In[26]:


#ds = yt.load('./../example_data/DD0401/DD0401')
#data = ds.all_data()

def plot(dsname, wdir = './', width = 500.0, dt = 5.0*yt.units.Myr):

    gal = Galaxy(dsname, wdir = wdir)
    data = gal.df

    @derived_field(name="logNO", units="")
    def _logNO(field, data):
        return  np.log10(data['N_Abundance'] / data['O_Abundance'])
    gal.ds.add_field(("gas", "logNO"), function=_logNO, units="")


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



    boxdim = np.array([width*1.25,width*1.25,30.0])*yt.units.pc
    region = gal.ds.box(gal.ds.domain_center - boxdim*0.5, gal.ds.domain_center + boxdim*0.5)

    proj = yt.ProjectionPlot(gal.ds, 'z', ['number_density','N_over_O','O_over_H','logNO'],
                        weight_field = 'number_density', data_source = region, width = (width,'pc'))
    proj.set_unit('number_density','cm**(-3)')
    proj.set_cmap('number_density','viridis')
    proj.set_zlim('number_density',1.0E-4,200.0)

    proj.set_cmap('O_over_H','cubehelix')
    proj.set_log('O_over_H', False)
    proj.set_zlim('O_over_H',-5,0)
    proj.set_colorbar_label('O_over_H', r'[O/H]')

    proj.set_cmap('N_over_O','PRGn')
    proj.set_log('N_over_O',False)
    proj.set_zlim('N_over_O',-2,2)
    proj.set_colorbar_label('N_over_O', r'[N/O]')


    proj.set_cmap('logNO','PRGn')
    proj.set_log('logNO',False)
    proj.set_zlim('logNO',-2,0.5)
    proj.set_colorbar_label('logNO', r'log( N / O )')
    proj.save() # necessary


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
    proj.save() # necessary

    vmin,vmax = -2,2
    x = proj.plots['N_over_O']
    x.image.set_norm( MidpointNormalize(midpoint= 0.5*(vmin+vmax), vmin=vmin,vmax=vmax))
    x.cb.set_norm(MidpointNormalize(midpoint=0.5*(vmin+vmax),vmin=vmin,vmax=vmax))
    x.cb.update_normal(x.image)
    x.save(str(gal.ds) + '_Projection_z_N_over_O_number_density.png')

    vmin, vmax = -2, 0.25
    x = proj.plots['logNO']
    x.image.set_norm( MidpointNormalize(midpoint= 0.0, vmin=vmin,vmax=vmax))
    x.cb.set_norm(MidpointNormalize(midpoint=0.0, vmin=vmin,vmax=vmax))
    x.cb.update_normal(x.image)
    x.save(str(gal.ds) + '_Projection_z_logNO_number_density.png')

    del(proj)
    del(gal)

    return


if __name__ == "__main__":

    all_ds = np.sort(glob.glob('DD????/DD????'))

    if len(all_ds) > 1:
        for dsname in all_ds[int(sys.argv[1]):int(sys.argv[2])]:
            print dsname
            plot(dsname.split('/')[0])
    else:
        plot(all_ds[0].split('/')[0])

