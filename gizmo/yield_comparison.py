import matplotlib
matplotlib.use('agg')

import numpy as np
import matplotlib.pyplot as plt
import yt

from galaxy_analysis.gizmo import yield_model
from galaxy_analysis.utilities import cy_convert_abundances as ca
from galaxy_analysis.plot.plot_styles import *
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import os, sys

# need to loop through files and compute:
#    1) Total mass in each element + error
#    2) X/H for all elements in stars. Save
#       distribution for eac
#
class GizmoData:

    def __init__(self, snapshot_path, yield_table = None,
                       yield_model_Z = 0.002, yield_model_FIRE_Z_scaling = True):

        self.work_dir = snapshot_path.split("output/snapshot_")[0]
        self.snapshot_name = snapshot_path.split("output/")[-1]

        if not (os.path.isfile(self.work_dir + 'output/snapshot_000.hdf5')):
            print("Cannot find initial snapshot file at " + self.work_dir + 'output/snapshot_000.hdf5')
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
                                                    age_is_fraction=True)

        self.metals   = np.unique([x[1] for x in self.ds0.field_list if ((x[0]=='PartType0')and('Metal' in x[1]))])

        self.initial_abundance = np.zeros(15)
        for i in np.arange(np.size(self.initial_abundance)):
            z = self.data0[('PartType0','Metallicity_%02i'%(i))]
            self.initial_abundance[i] = np.average(z).value # should all be identical
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
        self.yield_model_FIRE_Z_scaling = yield_model_FIRE_Z_scaling
        self.yield_table = yield_model.construct_yields(self.age_bins/1000.0, # pass bins as Gyr, Z = Z,
                                    Z = self.yield_model_Z, yieldtype = 'total',
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
            H_frac = M*(1.0-self.initial_abundance[0]-self.initial_abundance[1])-\
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
            Z = Z + self.initial_abundance[yield_model.elements.index(e)]
            return Z # mass/norm #+ self.initial_abundance[yield_model.elements.index(e)]
        
    def get_ratio(self,e1,e2,age=True):
        if age:
            vals1 = self.get_age_abund(e1)
            vals2 = self.get_age_abund(e2)
        else:
            vals1 = self.get_abund(e1)
            vals2 = self.get_abund(e2)
        return ca.abundance_ratio_array(e1,vals1,e2,vals2,
                                      input_type="mass")

    def get_difference(self,e1,e2, absval=True, ptype = 'star'):
        x1_vals_age = self.get_age_abund(e1,ptype=ptype)
        x2_vals_age = self.get_age_abund(e2,ptype=ptype)

        x1_vals = self.get_abund(e1,ptype=ptype)
        x2_vals = self.get_abund(e2,ptype=ptype)

        abund_age  = ca.abundance_ratio_array(e1, x1_vals_age, e2, x2_vals_age,
                                      input_type="mass")
        abund  = ca.abundance_ratio_array(e1, x1_vals, e2, x2_vals,
                                      input_type="mass")
        cutvals1 = self.get_abund('O',ptype=ptype)
        cutvals2 = self.get_abund('H',ptype=ptype)
        H_cut = ca.abundance_ratio_array('O',cutvals1,'H',cutvals2,input_type="mass")
        #abund     = abund[ (H_cut > -2.6)]
        #abund_age = abund_age[ (H_cut > -2.6)]

        if absval:
            abund = np.abs(abund - abund_age) # diff
        else:
            abund = abund-abund_age

        return abund

    def MDF(self, x1e,x2e,rmin=None,rmax=None,
            dbin=0.25, age=True, ptype='star', diff = False, absval=False):
        """
        Return MDF
        """

        if (absval) and (not diff):
            print("Are you sure you want to take the absolute value of hte abundance if it is not a diff?")
            raise ValueError

        if diff:
            x1_vals_age = self.get_age_abund(x1e,ptype=ptype)
            x2_vals_age = self.get_age_abund(x2e,ptype=ptype)        

            x1_vals = self.get_abund(x1e,ptype=ptype)
            x2_vals = self.get_abund(x2e,ptype=ptype)

            abund_age  = ca.abundance_ratio_array(x1e, x1_vals_age, x2e, x2_vals_age,
                                          input_type="mass")
            abund  = ca.abundance_ratio_array(x1e, x1_vals, x2e, x2_vals,
                                          input_type="mass")
            cutvals1 = self.get_abund('O',ptype=ptype)
            cutvals2 = self.get_abund('H',ptype=ptype)       
            H_cut = ca.abundance_ratio_array('O',cutvals1,'H',cutvals2,input_type="mass")
            #abund     = abund[ (H_cut > -2.6)]        
            #abund_age = abund_age[ (H_cut > -2.6)]

            if absval:
                abund = np.abs(abund - abund_age) # diff
            else:
                abund = abund-abund_age

        else:
            if age:
                x1_vals = self.get_age_abund(x1e,ptype=ptype)
                x2_vals = self.get_age_abund(x2e,ptype=ptype)
            else:
                x1_vals = self.get_abund(x1e,ptype=ptype)
                x2_vals = self.get_abund(x2e,ptype=ptype)

            abund  = ca.abundance_ratio_array(x1e, x1_vals, x2e, x2_vals,
                                      input_type="mass")        
            if ptype == 'gas':
                cutvals1 = self.get_abund('O',ptype=ptype)
                cutvals2 = self.get_abund('H',ptype=ptype)       

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
                 'std'    : np.std(abund), 'D9' : np.quantile(abund,0.9), 'D1' : np.quantile(abund,0.1) }

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


def identify_correlations(snapshot):

    amin = 0
    amax = 3
    db   = 0.001
    basename = snapshot.strip('.hdf5')

    gizmo_data = GizmoData("./m12q_res5700_logage/output/" + snapshot, yield_model_Z = 0.002)

    # plot correlations between different things

    fields = {'metallicity' : ('PartType4','Metallicity_00'),
              'age' : ('PartType4', 'StellarFormationTime') }

    field_units = {'metallicity' : None, 'age' : None }

    log_field = {'metallicity' : True, 'age' : False}

    elements = ['O','Fe','N']

    nrow = len(elements)
    ncol = len(list(fields.keys()))

    fig, ax = plt.subplots( nrow, ncol )
    fig.set_size_inches( ncol * 6, nrow * 6)


    for i,e in enumerate(elements):
        diff = gizmo_data.get_difference(e,'H',absval=True)

        for j,f in enumerate(fields.keys()):
            axindex = (i,j)
            field_val = gizmo_data.data[ fields[f] ]

            if not (field_units[f] is None):
                field_val = field_val.to( field_units[f] ).value

            if log_field[f]:
                field_val = np.log10( field_val )
            
            
            ax[axindex].scatter( diff, field_val)

            ax[axindex].set_xlabel("["+e+"/H] difference [dex]")
            ax[axindex].set_ylabel(f)
            ax[axindex].set_xlim(0.0001, 1.0)
            ax[axindex].semilogx()

            print(f, "correlation", np.corrcoef(diff,field_val)[1,0])


    plt.tight_layout()    
    fig.savefig(basename + "_correlations.png")

    return


def test_metallicity(snapshot):

    amin = 0
    amax = 3
    db = 0.001

    basename = snapshot.strip(".hdf5")

    all_stats = {}
    all_stats['O'] = {}
    all_stats['N'] = {}
    for Z in [0.0002, 0.001, 0.002, 0.01, 0.02, 0.04]:
        zstr = "%6.4f"%(Z)


        gizmo_data = GizmoData("./m12q_res5700/output/" + snapshot, yield_model_Z = Z)

        e1 = "O"
        e2 = "H"
        bins, hist1, all_stats[e1][zstr] = gizmo_data.MDF(e1,e2,amin,amax, diff=True,
                                             absval=True, dbin=db)

        e1 = "N"
        e2 = "H"
        bins, hist1, all_stats[e1][zstr] = gizmo_data.MDF(e1,e2,amin,amax, diff=True,
                                            absval=True, dbin=db)

    for e in all_stats:
        outname = basename + "_" + e + "_Z_comparison_stats.dat"
        outfile = open(outname, 'w')

        outfile.write("# Name ")
        for s in all_stats[e][ list(all_stats[e].keys())[0] ].keys():
            outfile.write(s + " ")
        outfile.write("\n")

        for k in all_stats[e].keys():
            outfile.write("%12s"%(k) + " ")

            for s in all_stats[e][k].keys():
                outfile.write("%8.2E "%(all_stats[e][k][s]))
            outfile.write("\n")

        outfile.close()



    return 




#
def logbin_comparison(snapshot, ptype = 'star'):

#    all_data = {'match16' : GizmoData("./m12q_res5700/output/" + snapshot),
#                'log4'  : GizmoData("./m12q_res5700_logage4/output/" + snapshot),
#                'log8'  : GizmoData("./m12q_res5700_logage8/output/" + snapshot),
#                'log16' : GizmoData("./m12q_res5700_logage/output/" + snapshot),
#                'log32' : GizmoData("./m12q_res5700_logage32/output/" + snapshot) }

    filepaths = {'match16': './m12q_res5700/', 'log4' : './m12q_res5700_logage4/',
                 'log8' : './m12q_res5700_logage8/', 'log16' : './m12q_res5700_logage/',
                 'log32' : './m12q_res5700_logage32/'}

    all_data = {}

    for k in filepaths:
        if os.path.isfile( filepaths[k] + 'output/' + snapshot):

            all_data[k] = GizmoData( filepaths[k] + 'output/' + snapshot)

    if len( list(all_data.keys())) == 0:
        print("No data files found in any of the filepaths for snapshot " + snapshot)
        print(filepaths.values())
        raise RuntimeError


    basename = snapshot.strip(".hdf5")

    #
    # Plot the MDFs
    #

    fig,axes = plt.subplots(1,3,sharey=True,sharex=True)
    fig.set_size_inches(18,6)
    fig.subplots_adjust(wspace=0,hspace=0)

    xy=(0.05,0.20)

    def plot_ax(ax, e1, e2, gizmo_data, db = 0.1,
                            amin=-3, amax=3, age=True,
                            color = None, label = None):
        if age:
            if color is None:
                color = 'C0'
            if label is None:
                label = "Post-process"

            bins, hist1 = gizmo_data.MDF(e1,e2,amin,amax, dbin=db, age=True, ptype = ptype)

        else: # plot simulation abundances
            if color is None:
                color = 'black'
            if label is None:
                label = 'Simulation'

            bins, hist1 = gizmo_data.MDF(e1,e2,amin,amax,dbin=db, age=False, ptype = ptype)

        ax.step(bins,hist1/(1.0*np.sum(hist1)),where='post',
                         color=color, label = label,lw=3)

        return


    for i,k in enumerate(all_data.keys()):
        plot_ax(axes[0], 'O','H', all_data[k], age = True,
                color = 'C%i'%(i), label = k)
        plot_ax(axes[1],'Fe','H',all_data[k], age = True,
                color = 'C%i'%(i), label = k)
        plot_ax(axes[2], 'N', 'H',all_data[k],age=True,color='C%i'%(i),label=k)

    plot_ax(axes[0],'O','H',all_data['log16'],age=False,color='black',label='FIRE')
    plot_ax(axes[1],'Fe','H',all_data['log16'],age=False,color='black',label='FIRE')

    for a in axes:
        a.set_ylim(0.0,0.2)
        a.set_xlim(-4,2)
        a.set_xlabel('[X/H] [dex]')
    axes[0].set_ylabel('Fraction of Stars')

    axes[0].annotate('[O/H] ', xy=xy,xycoords='axes fraction', size = 30)
    axes[1].annotate('[Fe/H]', xy=xy,xycoords='axes fraction', size = 30)
    axes[2].annotate('[N/H]',xy=xy,xycoords='axes fraction',size=30)

    plt.minorticks_on()

    axes[0].legend(loc='upper left')
    fig.savefig(basename + "_logbin_comparison_" + ptype + ".png")

    plt.close()

    #
    #
    # Plot the offsets
    #
    #

    fig,axes = plt.subplots(1,3,sharey=True,sharex=True)
    fig.set_size_inches(18,6)
    fig.subplots_adjust(wspace=0,hspace=0)

    xy=(0.7, 0.90)

    def plot_ax_2(ax, e1, e2, gizmo_data, db = 0.001,
                            amin=0, amax=3,
                            color = None, label = None, yshift=0):
        if color is None:
            color = 'C0'
        if label is None:
            label = "Post-process"

        bins, hist1, stats = gizmo_data.MDF(e1,e2,amin,amax, diff=True, absval=True, dbin=db, ptype=ptype)

        ax.step(bins,np.cumsum(hist1/(1.0*np.sum(hist1))),where='post',
                         color=color, label = label,lw=3)

        yoff = 0.6 + 0.05*yshift
        ax.annotate(label + ' < 0.05 = %0.3f dex'%stats['0.05dex'],
                    xy=(xy[0]-0.32,xy[1]-0.2-yoff),xycoords='axes fraction',
                    size = 15)

        return stats


    all_stats = {'O' : {}, 'Fe' : {}, 'N' : {}}
    for i,k in enumerate(all_data.keys()):
        all_stats['O'][k] = plot_ax_2(axes[0], 'O','H', all_data[k],
                color = 'C%i'%(i), label = k, yshift=i)
        all_stats['Fe'][k] = plot_ax_2(axes[1],'Fe','H',all_data[k],
                color = 'C%i'%(i), label = k, yshift=i)
        all_stats['N'][k] = plot_ax_2(axes[2], 'N','H',all_data[k],
                color = 'C%i'%(i), label = k, yshift=i)

    axes[0].set_ylabel("Fraction of Stars")
    for i in np.arange(3):
        axes[i].set_ylim(0.01,1.0)
        axes[i].semilogx()
        axes[i].set_xlabel("Abundance Difference [dex]")
        axes[i].set_xticks([0.001,0.01,0.1,1.0])
        axes[i].set_xticklabels(["0.001","0.01","0.1","1.0"])
        axes[i].yaxis.set_minor_locator(AutoMinorLocator())

    fig.savefig(basename + "_logbin_comparison_diff_" + ptype + ".png")

    for e in all_stats:
        outname = basename + "_" + e + "_comparison_stats_" + ptype + ".dat"
        outfile = open(outname, 'w')

        outfile.write("# Name ")
        for s in all_stats[e][ list(all_stats[e].keys())[0] ].keys():
            outfile.write(s + " ")
        outfile.write("\n")

        for k in all_stats[e].keys():
            outfile.write("%8s"%(k) + " ")

            for s in all_stats[e][k].keys():
                outfile.write("%8.2E "%(all_stats[e][k][s]))
            outfile.write("\n")

        outfile.close()

    return


if __name__ == "__main__":


    # 
#    identify_correlations(   str(sys.argv[1]))

    # test this out:
    logbin_comparison(str(sys.argv[1]))
#    logbin_comparison(str(sys.argv[1]), ptype = 'gas')

    # test for best metallicity
    
#    test_metallicity( str(sys.argv[1]) )

