import numpy as np
from scipy import integrate
import yt
import os
import matplotlib.pyplot as plt
import glob

# AE: Comment out below import unless you feel like
#     installing a bunch of stuff:
# from galaxy_analysis.plot.plot_styles import *



# globals since they are ifdefs in code and I'm lazy:
AGE_BIN_START      = 1.0     # Myr
AGE_BIN_END        = 14000.0 # Myr

elements = ['Total','He','C','N','O','Ne','Mg','Si','S','Ca','Fe']
element_num = {}
i = 0
for e in elements:
    element_num[e] = i
    i = i + 1

def generate_metal_fields(ds, _agebins=None,
                              _elements=elements,
                              _yields=None,
                              ptype='PartType0'):
    """
    Generate derived fields mapping the age tracers to 
    actual elemental abundances using the given set of
    yields. yields must be a NxM array with N = the
    number of age bins and M = the number of elements
    (M = 11 for default FIRE). Each value here should be the
    yield (in solar masses) per solar mass of star formation
    in each age bin for each element.

    The derived fields computing elemental mass loop through
    all the age bin tracer fields.

    Derived fields will be of form:
        (ptype,"ELEMENTNAME_mass")
        (ptype,"ELEMENTNAME_fraction")
        (ptype,"ELEMENTNAME_actual_mass")

    where 'ptype' is the passed particle type ("PartType0" or 
    "PartType4" probably) and the "ELEMENTNAME_actual_mass" field
    is the actual mass of that element in the yields followed in
    the simulation (e.g. something like:
       (ptype,"Metallicity_X")*(ptype,"particle_mass"), where X
       is the metallicity number corresponding to that given element).
    """
    
    offset = 15 # Metallicity_15 is the first age tracer field
                #   -- 11 elements + 4 r-process... change if this changes

    def _metal_mass_test(_ptype, _ei):
        # test metals in bin zero 
    
        def temp(field,data):
            mass_p = np.zeros(np.shape( data[(_ptype,'particle_mass')]))
            for i in np.arange(np.size(_agebins)-1):
                fname    = 'Metallicity_%02i'%(offset + i)
                mass_p += data[(ptype,fname)].value * _yields[i,_ei] #(agebinnum,  elementnum)
            
            mass_p = mass_p * yt.units.Msun
          
            return mass_p
        
        return temp
    
    def _metal_fraction_test(_ptype, _e):
        def temp(field,data):
            Mp    = data[(_ptype,'particle_mass')].to('Msun')
            abund = data[('all',_ptype + '_' + _e + '_mass')].to('Msun') / Mp
            abund[Mp==0.0] = 0.0
            return abund
        return temp
    
    def _metal_mass_actual_test(_ptype,_ei):
        def temp(field,data):
            Mp    = data[(_ptype,'particle_mass')].to('Msun')
            abund = data[(_ptype,"Metallicity_%02i"%(_ei))]
        
            
            return abund*Mp
        return temp
            
    
    
    for ei,e in enumerate(_elements):

        ds.add_field( ('all', ptype + '_' + e + '_mass'), sampling_type='particle', 
                      function=_metal_mass_test(ptype, ei),
                      units = 'Msun', force_override=True)
        ds.add_field(  ('all',ptype + '_' + e + '_fraction'), sampling_type='particle',
                      function=_metal_fraction_test(ptype,e),
                      units='', force_override=True)
        
        ds.add_field( ('all',ptype + '_' + e + '_actual_mass'), sampling_type='particle',
                       function=_metal_mass_actual_test(ptype,ei),
                       units='Msun', force_override=True)
        
    return    

def _generate_star_metal_fields(ds,
                                _agebins=None,
                                _elements=elements,
                                _yields=None,
                                ptype='PartType4'):

    """
    See above function. Computes surface abundances
    of given star as derived from the age tracers
    """
    
    offset = 15
    
    def _star_metal_mass_test(_ptype, _ei):
        # test metals in bin zero 
    
        def temp(field,data):
            mass_p = np.zeros(np.shape( data[(_ptype,'particle_mass')]))
            for i in np.arange(np.size(_agebins)-1):
                fname    = 'Metallicity_%02i'%(offset + i)
                mass_p += data[(ptype,fname)].value * _yields[i,_ei] #(agebinnum,  elementnum)
            
            mass_p = mass_p * yt.units.Msun
          
            return mass_p
        
        return temp
    
        
    for ei,e in enumerate(_elements):

        ds.add_field( ('all', ptype + '_' + e + '_mass'), sampling_type='particle', 
                      function=_star_metal_mass_test(ptype, ei),
                      units = 'Msun', force_override=True)
        
    return


#
# Extracted FIRE yield model:
#   - Model is as-is from the code, but does not include any
#     metallicity dependence (wind rates are fixed to a chosen
#     metallicity, default is solar)
#
#

def sn_rate(t):
    """
    CCSNE rate

    SN / Gyr  per solar msas of star formation
        Changed output to /Gyr to keep same units as input t
    """
    agemin = 0.003401   # Gyr
    agebrk = 0.010370   # Gyr
    agemax = 0.03753    # Gyr
    
    RSNE = 0.0
    if (t>agemin):
        if (t <= agebrk):
            RSNE = 5.408E-4
        elif (t<=agemax):
            RSNE=2.516E-4
        
        if (t > agemax):
            #RSNE=5.3E-8+1.6*np.exp(-0.5*((t-0.05)/0.01)*((t-0.05)/0.01)) # This is JUST SNIa
            RSNE=0.0 # set to zero for CCSNE
            
    return RSNE * 1000.0

def snIa_rate(t):
    """
    SNIa rate (SN/Gyr)  - t in Gyr
    """

    agemin = 0.003401   # Gyr
    agebrk = 0.010370   # Gyr
    agemax = 0.03753    # Gyr
    
    RSNE = 0.0
    if (t > agemax):
        RSNE=5.3E-8+1.6E-5*np.exp(-0.5*((t-0.05)/0.01)*((t-0.05)/0.01)) # This is JUST SNIa
            
    return RSNE * 1000.0

def wind_yields(i,element=None):
    """
    Yields (in fraction) per element with winds
    """
    yields = [0.0, 0.36,0.016,0.0041,0.0118] + [0.0]*6
    yields[0] = np.sum(yields) # total yield

    # if element passed, use that - otherwise use yield indeces
    if not (element is None):
        if element == 'all':
            return yields

    return yields[i]

def wind_rate(t, Z = 1.0, GasReturnFraction = 1.0):
    """
    Mass loss rate from stellar winds. Z is in solar.
    """    
    p = 0.0
    if (t <= 0.001):
        p = 11.6846
    else:
        if (t <=0.0035):
            logZ=np.log10(Z)
            p=11.6846*Z*10.0**(1.838*(0.79+logZ)*(np.log10(t)-(-3.00)))
        else:
            if (t<=0.1):
                p=72.1215*(t/0.0035)**(-3.25)+0.0103
            else:
                p=1.03*t**(-1.1) / (12.9-np.log10(t))
    
    #if (t < 0.1):
    #    p = p * 1.0
    
    # assuming wind_rate is in Msun / Myr per solar mass of SF
    
    rate = p * GasReturnFraction * 1.4 * 0.291175
    
    return rate # might already be / Gyr

def snIa_yields(i, element = None):
    yields = [1.4,0.0,0.049,1.2E-6,0.143,0.0045,0.0086,0.156,0.087,0.012,0.743]
    
    # if element passed, use that - otherwise use yield indeces
    if not (element is None):
        if element == 'all':
            return yields

    return yields[i]
    
def snII_yields(i, element = None):
    # if element passed, use that - otherwise use yield indeces
    
    yields = [2.0,3.87,0.133,0.0479,1.17,0.30,0.0987,0.0933,0.0397,0.00458,0.0741]

    if not (element is None):
        if element == 'all':
            return yields
        
        i = element_num[element]
    

    return yields[i]    

def construct_yields(agebins, yieldtype = 'total'):

    points = [0.003401, 0.010370, 0.03753, 0.001, 0.05, 0.10]
    
    yields = np.zeros( (np.size(agebins)-1 , np.size(elements))) 
    
    if yieldtype == 'snII' or yieldtype == 'total' or yieldtype == 'sn_only':
        numsnII = np.zeros(np.size(agebins)-1)
        
        for i in np.arange(np.size(agebins) - 1):
            
            if i == 0:
                mint = 0.0
            else:
                mint = agebins[i]
            maxt = agebins[i+1]

            numsnII[i] = integrate.quad( sn_rate, mint, maxt, points = points)[0]
            
        yields += np.outer(numsnII, snII_yields(-1, element = 'all'))
    
    if yieldtype == 'snIa' or yieldtype == 'total' or yieldtype == 'sn_only':    
        numsnIa = np.zeros(np.size(agebins)-1)
        
        for i in np.arange(np.size(agebins) - 1 ):
            
            if i == 0:
                mint = 0.0
            else:
                mint = agebins[i]
            maxt = agebins[i+1]
            numsnIa[i] = integrate.quad( snIa_rate, mint, maxt, points = points)[0]

        yields += np.outer(numsnIa, snIa_yields(-1, element = 'all'))        
 
    if yieldtype == 'winds' or yieldtype == 'total':    
        wind_mass  = np.zeros(np.size(agebins)-1)
        windpoints = [0.001, 0.0035, 0.1]
        
        for i in np.arange(np.size(agebins) - 1 ):
            
            if i == 0:
                mint = 0.0
            else:
                mint = agebins[i]
            maxt         = agebins[i+1]
            wind_mass[i] = integrate.quad(wind_rate, mint, maxt, points = windpoints)[0]

        yields += np.outer(wind_mass, wind_yields(-1, element = 'all'))  
    return yields


def get_bins(infile = "./gizmo.out", binfile = "age_bins.txt"):
    """
    Assuming gizmo.out is included in directory, generate
    age bins. Requires age bin file if custom bins used
    """
    count   = 0
    logbins = True
    for line in open(infile,'r'):
        if "GALSF_FB_FIRE_AGE_TRACERS=" in line:
            num_tracers = int(line.split("=")[1])

        if "GALSF_FB_FIRE_AGE_TRACERS_CUSTOM" in line:
            logbins = False

        if count > 100:
            break
        count = count + 1


    if logbins:
        NBINS    = num_tracers + 1
        binstart = np.log10(AGE_BIN_START)
        binend   = np.log10(AGE_BIN_END)
        bins     = np.logspace(binstart, binend, NBINS)

    else:
        # read bins from file
        try:
            bins     = np.genfromtxt(binfile)
        except:
            print("Custom bins. Problem loading binfile " + binfile)
            raise ValueError


    return bins

def compute_error(outfile = 'error.dat', overwrite=False, limit_input=False):
    """
    Compute the error in the given age tracer model, defined
    for each element (i) as:

          (M_i_age - M_i_sim) / M_i_sim

    where M_i_age is the mass of that element as derived from
    convolving age tracers with extracted yield model from FIRE
    and M_i_sim is the actual mass of that element in the simulation

    """
    
    bins = get_bins()


    total_yields = construct_yields(bins/1000.0, # pass bins as Gyr
                                    yieldtype = 'total')


    if limit_input:
        ds_list = np.sort(glob.glob('./output/snapshot_*0.hdf5'))
    else:
        ds_list = np.sort(glob.glob('./output/snapshot_*.hdf5'))

    #
    nlines = 0
    open_mode = 'w'

    if os.path.exists(outfile) and (not overwrite):
        data = np.genfromtxt(outfile)
        nlines = np.size(data[:,0])

        if nlines == (np.size(ds_list)):
            print("Not overwriting output")

            if nlines > np.size(ds_list):
                print("This file seems to exist for more data than available. Double check that you want to compute")

            return
        else:
            print("Only computing for files that don't already exist in error output - assuming contiguous")
            open_mode = 'a'

    f = open(outfile, open_mode)

    ds0    = yt.load(ds_list[0])
    data0   = ds0.all_data()
    generate_metal_fields(ds0, _agebins = bins, _yields = total_yields)

    mtrue_initial = {}
    for e in elements:
        m = data0[('all','PartType0_'+e+'_actual_mass')]
        mtrue_initial[e] = np.sum( m ).to('Msun')

    dsstart = 0
    if nlines > 0:
        ds_list = ds_list[nlines:]
        dsstart = nlines

    for dsi, dsname in enumerate(ds_list):
        dsi = dsi + dsstart

        ds = yt.load(dsname)
        data = ds.all_data()
        fields = ds.field_list
        generate_metal_fields(ds, _agebins=bins, _yields=total_yields)

        ptypes = np.unique([x[0] for x in ds.field_list])
        metals = np.unique([x[1] for x in ds.field_list if ((x[0] == 'PartType0') and ('Metal' in x[1]))])

        M_new_stars = 0.0 * yt.units.Msun
        if 'PartType4' in ptypes:
            print("File %003i:  Number of new stars    %00005i"%(dsi, np.size(data[('PartType4','Metallicity_00')])))
            M_new_stars = data[('PartType4','particle_mass')].to('Msun')
            _generate_star_metal_fields(ds, _agebins = bins, _yields = total_yields)
        else:
            print("File %003i:  No star particles    "%(dsi))

        if ds.cosmological_simulation:
            f.write("%03i %3.3f "%(dsi, ds.current_redshift))
        else:
            f.write("%03i %3.3f "%(dsi, 0.0000))

        tolerance = 1.0E-6 * yt.units.Msun

        for ei,e in enumerate(elements):
#            if e=='Total':
#                continue

            m = data[('all','PartType0_' + e + '_mass')]
            mstar = 0.0 * yt.units.Msun
            mtrue_mstar = 0.0 * yt.units.Msun
            if 'PartType4' in ptypes:
                mstar = data[('all','PartType4_' + e + '_mass')]
                mtrue_mstar = data[('PartType4','Metallicity_%02i'%(ei))].value * data[('PartType4','particle_mass')].to('Msun')

            mtrue = data[('all','PartType0_'+e+'_actual_mass')]
            mtrue_total = np.max( [(np.sum(mtrue) + np.sum(mtrue_mstar) - mtrue_initial[e]), 0.0])
            m_total     = np.sum(m) + np.sum(mstar)

#            if mtrue_total == 0.0:
#                error = 0.00
#            else:
#                error       = (m_total.value - mtrue_total) / mtrue_total

            f.write("%5.5E %5.5E %5.5E %5.5E "%(m_total, mtrue_total, np.sum(mstar), np.sum(mtrue_mstar) ))

        f.write("\n")
        f.flush()
        os.fsync(f.fileno())
        del(ds)   # just in case
        del(data) # just in case


    f.close()


    return


def plot_error(infile = 'error.dat'):

    data = np.genfromtxt(infile)

    fig, ax = plt.subplots()
    fig.set_size_inches(6,6)

    time     = data[:,0] # Myr
    redshift = data[:,1] # redshift

    cosmological = False
    if np.size(redshift[redshift>0]) > 0:
        plot_x = redshift + 1
        cosmological = True
        HubbleParam  = 0.702
    else:
        plot_x = time
        HubbleParam  = 1.0


    i = 0
    for ei, e in enumerate(elements):
        if e == 'Total':
            continue

        if (not(e == 'C')) and (not(e == 'O')) and (not(e == 'N')) and (not(e == 'Fe')) and (not(e=='Mg')):
            continue

        mass   = data[:,2 + 4*ei]
        mtrue  = data[:,3 + 4*ei] * HubbleParam
        error  = np.zeros(np.size(mtrue))

        error[mtrue>0]  = (mass[mtrue>0]-mtrue[mtrue>0])/mtrue[mtrue>0]

        ax.plot(plot_x, np.abs(error), lw = 3, ls = '-', label = e, color = 'C%i'%(i))
#        ax.plot(plot_x[error<0], np.abs(error[error<0]), lw = 3, ls = ':', color = 'C%i'%(i))
        i = i + 1

    if cosmological:
        ax.set_xlim(15.0+1, 0.0 + 1)
        ax.set_xlabel("1 + z")
        ax.semilogx()
    else:
        ax.set_xlim(0.0, np.max([np.max(time),50.0]))
        ax.set_xlabel("Time (Myr)")

    ax.set_ylim(1.0E-4,1.0)
    ax.semilogy()

    ax.set_ylabel("Fractional Error")
    ax.legend(loc="best",ncol=3)

    plt.minorticks_on()

    plt.tight_layout()
    fig.savefig("error.png")

    return


if __name__ == "__main__":

    compute_error()
    plot_error()
