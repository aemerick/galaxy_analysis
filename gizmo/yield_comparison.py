import numpy as np
import matplotlib.pyplot as plt
import yt

from galaxy_analysis.gizmo import yield_model
from galaxy_analysis.utilities import cy_convert_abundances as ca
from galaxy_analysis.plot.plot_styles import *
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


# need to loop through files and compute:
#    1) Total mass in each element + error
#    2) X/H for all elements in stars. Save
#       distribution for eac
#
class GizmoData:

    def __init__(self, snapshot_path, yield_table = None,
                       yield_model_Z = 1.0E-3, yield_model_FIRE_Z_scaling = True):

        self.work_dir = snapshot_path.strip("output/snapshot_")[0]
        self.snapshot_name = snapshot_path.strip("output/")[-1]

        if not (os.isfile(self.work_dir 'output/snapshot_000.hdf5')):
            print("Cannot find initial snapshot file")
            raise RuntimeError

        self.age_bins = yield_model.get_bins(config_file = self.work_dir + "/gizmo.out",
                                             binfile     = self.work_dir + "/age_bins.txt")

        self.yield_model_Z = None
        self.yield_model_FIRE_Z_scaling = None
        self.yield_model_age_fraction = None
        if yield_table is None:
            self.generate_yield_table(yield_model_Z, yield_model_FIRE_Z_scaling)
        else:
            self.yield_table = yield_table

        # load datasets and define fields
        self.ds0      = yt.load(self.work_dir + 'output/snapshot_000.hdf5')
        self.data0    = self.ds0.all_data()

        yield_model.generate_metal_fields(self.ds0, _agebins=self.age_bins,
                                                    _yields=self.yield_table,
                                                    _age_is_fraction=True)

        self.metals   = np.unique([x[1] for x in self.ds0.field_list if ((x[0]=='PartType0')and('Metal' in x[1]))])

        self.initial_abundance = np.zeros(15)
        for i in np.arange(np.size(self.initial_abundance)):
            z = self.data0[('PartType0','Metallicity_%02i'%(i))]
            initial_abundance[i] = np.average(z).value # should all be identical
        self._logH = np.log10(self.ds0.hubble_constant)
        self._littleh = 1.0

        # load actual dataset
        self.ds       = yt.load(snapshot_path)
        self.data     = self.ds.all_data()

        yield_model.generate_metal_fields(self.ds,_agebins=self.age_bins,
                                         _yields=self.yield_table,age_is_fraction=True)
        yield_model._generate_star_metal_fields(self.ds, _agebins = self.age_bins,
                                         _yields = self.yield_table,age_is_fraction=True)


        return

    def generate_yield_table(self, yield_model_Z = 1.0E-3, yield_model_FIRE_Z_scaling = True):
        self.yield_model_Z              = yield_model_Z
        self.yield_model_FIRE_Z_scaling = yield_model_FIRE_scaling
        self.yield_table = yield_model.construct_yields(self.age_bins/1000.0, # pass bins as Gyr, Z = Z,
                                    Z = self.yield_model_Z, yieldtype = 'total'
                                    FIRE_Z_scaling= self.yield_model_FIRE_Z_scaling)
        return

    def get_abund(self,e,ptype='star'):
        if ptype == 'star':
            ptype = "PartType4"
        elif ptype == 'gas':
            ptype = "PartType0"
        if e == "H":
            return 1.0 - self.data[(ptype,"Metallicity_00")] - self.data[(ptype,"Metallicity_01")]
        else:
            ei = yield_model.elements.index(e)

        return self.data[(ptype,"Metallicity_%02i"%(ei))]

    def get_age_abund(self,e,ptype='star'):
        if ptype == 'star':
            ptype = "PartType4"
        elif ptype == 'gas':
            ptype = "PartType0"
        if e == "H":
            # H_frac = 1.0 - self.data[(ptype,"Metallicity_00")] - self.data[(ptype,"Metallicity_01")]
            M = self.data[(ptype,'particle_mass')].to('Msun')
            H_frac = M*(1.0-initial_abundance[0]-initial_abundance[1])-\
                       self.data[('all',ptype+'_Total_mass')] / self._littleh - self.data[('all',ptype+'_He_mass')] / self._littleh
            H_frac = H_frac / self.data[(ptype,'particle_mass')].to('Msun')

            return H_frac
        else:
            ei = yield_model.elements.index(e)
            # very bad!!!
            mass = self.data[('all', ptype + '_' + e + '_mass')].to('Msun') / self._littleh
            norm = self.data[(ptype,'particle_mass')].to('Msun')
           # M_norm # (16752.063237698454*yt.units.Msun)
            Z = mass / norm
            Z = Z + initial_abundance[yield_model.elements.index(e)]
            return Z # mass/norm #+ initial_abundance[yield_model.elements.index(e)]
        
    def get_ratio(self,e1,e2,age=True):
        if age:
            vals1 = get_age_abund(e1,dat)
            vals2 = get_age_abund(e2,dat)
        else:
            vals1 = get_abund(e1,dat)
            vals2 = get_abund(e2,dat)
        return ca.abundance_ratio_array(e1,vals1,e2,vals2,
                                      input_type="mass")

    def MDF(x1e,x2e,rmin=None,rmax=None,
            dbin=0.25, age=True, ptype='star', diff = False, absval=False):
        """
        Return MDF
        """

        if (absval) and (not diff):
            print("Are you sure you want to take the absolute value of hte abundance if it is not a diff?")
            raise ValueError

        if diff:
            x1_vals_age = get_age_abund(x1e,self.data,ptype=ptype)
            x2_vals_age = get_age_abund(x2e,self.data,ptype=ptype)        

            x1_vals = get_abund(x1e,self.data,ptype=ptype)
            x2_vals = get_abund(x2e,self.data,ptype=ptype)

            abund_age  = ca.abundance_ratio_array(x1e, x1_vals_age, x2e, x2_vals_age,
                                          input_type="mass")
            abund  = ca.abundance_ratio_array(x1e, x1_vals, x2e, x2_vals,
                                          input_type="mass")
            cutvals1 = get_abund('O',self.data,ptype=ptype)
            cutvals2 = get_abund('H',self.data,ptype=ptype)       
            H_cut = ca.abundance_ratio_array('O',cutvals1,'H',cutvals2,input_type="mass")
            #abund     = abund[ (H_cut > -2.6)]        
            #abund_age = abund_age[ (H_cut > -2.6)]

            if absval:
                abund = np.abs(abund - abund_age) # diff
            else:
                abund = abund-abund_age

        else:
            if age:
                x1_vals = get_age_abund(x1e,self.data,ptype=ptype)
                x2_vals = get_age_abund(x2e,self.data,ptype=ptype)
            else:
                x1_vals = get_abund(x1e,self.data,ptype=ptype)
                x2_vals = get_abund(x2e,self.data,ptype=ptype)

            abund  = ca.abundance_ratio_array(x1e, x1_vals, x2e, x2_vals,
                                      input_type="mass")        
            if ptype == 'gas':
                cutvals1 = get_abund('O',self.data,ptype=ptype)
                cutvals2 = get_abund('H',self.data,ptype=ptype)       

                H_cut = ca.abundance_ratio_array('O',cutvals1,'H',cutvals2,input_type="mass")

                abund     = abund[ (H_cut > -2.6)]        

        if rmin is None:
            rmin = np.min(abund)
        if rmax is None:
            rmax = np.max(abund)

        nbins = int((rmax - rmin)/dbin)
        hist, bins = np.histogram(abund, bins = nbins, range = (rmin,rmax))
        hist2 = np.ones(np.size(hist)+1)
        hist2[:-1] = hist
        hist2[-1] = hist2[-2]

        stats = {'median' : np.median(abund), 'mean' : np.average(abund),
                 'Q1'     : np.quantile(abund,0.25), 'Q3' : np.quantile(abund,0.75),
                 'IQR'    : np.quantile(abund,0.75) - np.quantile(abund,0.25),
                 'std'    : np.std(abund)}

        # compute fraction < a given offset
        if diff:
            stats['0.2dex']   = np.size( abund[ np.abs(abund) < 0.2  ]) / (1.0*np.size(abund))
            stats['0.1dex']   = np.size( abund[ np.abs(abund) < 0.1  ]) / (1.0*np.size(abund))
            stats['0.05dex']  = np.size( abund[ np.abs(abund) < 0.05 ]) / (1.0*np.size(abund))
            stats['0.02dex']  = np.size( abund[ np.abs(abund) < 0.02 ]) / (1.0*np.size(abund))
            stats['0.01dex']  = np.size( abund[ np.abs(abund) < 0.01 ]) / (1.0*np.size(abund))
            stats['0.005dex'] = np.size( abund[ np.abs(abund) < 0.005]) / (1.0*np.size(abund))
        if diff:
            return bins,hist2,stats
        else:
            return bins, hist2

#
def logbin_comparison(snapshot):

    all_data = {'log8'  : GizmoData("./m12q_res5700_logbin/output/" + snapshot)
                'log16' : GizmoData("./m12q_res5700_logbin/output/" + snapshot)
                'log32' : GizmoData("./m12q_res5700_logbin/output/" + snapshot) }

    fig,axes = plt.subplots(1,2,sharey=True,sharex=True)
    fig.set_size_inches(12,6)
    fig.sublots_adjust(wspace=0,hspace=0)

    xy=(0.05,0.20)

    def plot_ax(ax, e1, e2, gizmo_data, db = 0.1,
                            amin=-3, amax=3, age=True,
                            color = None, label = None):
        if age:
            if color is None:
                color = 'C0'
            if label is None:
                label = "Post-process"

            bins, hist1 = gizmo_data.MDF(e1,e2,amin,amax,age=True,dbin=db)

        else: # plot simulation abundances
            if color is None:
                color = 'black'
            if label is None:
                label = 'Simulation'

            bins, hist1 = gizmo_data.MDF(e1,e2,amin,amax,age=False,dbin=db)

        ax.step(bins,hist1/(1.0*np.sum(hist1)),where='post',
                         color=color, label = label,lw=3)

        return


    for i,k in enumerate(all_data.keys()):
        plot_ax(axes[0], 'O','H', all_data[k], age = True,
                color = 'C%i'%(i), label = k)
        plot_ax(axes[1],'Fe','H',all_data[k], age = True,
                color = 'C%i'%(i), label = k)
    plot_ax(axes[0],'O','H',all_data['log16'],age=False,color='black',label='FIRE')
    plot_ax(axes[1],'Fe','H',all_data['log16'],age=False,color='black',label='FIRE')

    for a in axes:
        a.set_ylim(0.0,0.2)
        a.set_xlim(-4,2)
        a.set_xlabel('[X/H] [dex]')
    axes[0].set_ylabel('Fraction of Stars')

    axes[0].annotate('[O/H] ', xy=xy,xycoords='axes fraction', size = 30)
    axes[1].annotate('[Fe/H]', xy=xy,xycoords='axes fraction', size = 30)

    plt.minorticks_on()

    ax[0].legend(loc='upper left')
    fig.savefig("logbin_comparison.png")

    return


if __name__ == "__main__":

    # test this out:
    logbin_comparison(str(sys.argv[1]))
