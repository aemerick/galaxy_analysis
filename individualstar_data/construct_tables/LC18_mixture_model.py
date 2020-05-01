"""
    Generate a mixture model of the Limongi+2018 (LC18) yields using the
    v_rot population fractions of Prantzos+2018. Output yields in same format
    as the input.

    Also, allows generation of re-sampled data in metallicity using linear
    interpolation in log space.

"""


import numpy as np
import h5py
import copy

available_yields = {'FRUITY' : 'fruity_elem_stable.npy',
                    'LC18_R_0' : 'LC18_R_0_elem_stable.npy',
                    'LC18_R_150': 'LC18_R_150_elem_stable.npy',
                    'LC18_R_300': 'LC18_R_300_elem_stable.npy',
                    'LC18_winds_R_0' : 'LC18_winds_R_0_elem_stable.npy',
                    'LC18_winds_R_150': 'LC18_winds_R_150_elem_stable.npy',
                    'LC18_winds_R_300': 'LC18_winds_R_300_elem_stable.npy',
                    'Monash' : 'monash_elem_stable.npy',
                    'Nomoto' : 'nomoto_elem_stable.npy',
                    'Nomoto-HNe' : 'nomoto_HNe_elem_stable.npy',
                    'PopIII' : 'popIII_yields.in'}



# not best practice, but keep as global to load only once
prantzos_data = np.genfromtxt("./../../physics_data/Prantzos2018/prantzos_fig4_fractions.dat",names=True)

def prantzos_mixture(Z, vtype='v0'):
    """
    Given a metallicity, returns fration of stars at a given velocity
    (vtype = 'v0','v150', OR 'v300').
    """

    # zsun = 0.01345
    logZ = np.log10(Z)

    m = (-2.0 - (-1.0))/(np.log10(3.236E-4) - np.log10(3.236E-3))
    FeH  = m * logZ + (-2 - m*np.log10(3.236E-4))

    return np.interp(FeH, prantzos_data['FeH'], prantzos_data[vtype])


def make_mixture(wdir = "./", resample=True, nsamples=3):
    """
    Generate the mixture model with optional re-sampling and interpolation
    in metallicity.

    If resampling is True, `nsamples` provides the number of new metallicity
    points in between each current Z value (default = 2 gives 10 total
    metallicity points: 4 original + 2 in between each).
    """


    if resample:
        # need to generate a resampling of the yields

        tempyields = np.load( wdir + available_yields['LC18_R_0'], allow_pickle=True).item()
        oldZ   = tempyields['Z_list']
        newZ   = np.append(
                  np.concatenate(
                     [np.logspace(np.log10(oldZ[i]),np.log10(oldZ[i+1]),2+nsamples)[:-1] for i in np.arange(np.size(oldZ)-1)]),
                  oldZ[-1]
                )

        oldZkeys = list(tempyields['yields'].keys())
        del tempyields
        newformat="Z=%6.4e"
        reform_oldZkeys = [newformat%(z) for z in oldZ] # reformat for consistency

        for typestr in ['_','_winds_']:
            v0yields = np.load( wdir + available_yields['LC18' + typestr + 'R_0'], allow_pickle=True).item()
            v150yields = np.load( wdir + available_yields['LC18' + typestr + 'R_150'], allow_pickle=True).item()
            v300yields = np.load( wdir + available_yields['LC18' + typestr + 'R_300'], allow_pickle=True).item()

            newset = copy.deepcopy(v0yields)

            newset['Z_list'] = newZ


            ziold = -1 # Only works if znew includes same grid points as zold
            # otherwise this iterator does not work
            for zinew,znew in enumerate(newset['Z_list']):

                Zkeynew = newformat%(znew)

                #print(zinew,znew)

                f0   = prantzos_mixture( znew, vtype='v0')
                f150 = prantzos_mixture( znew, vtype='v150')
                f300 = prantzos_mixture( znew, vtype='v300')

                if Zkeynew in reform_oldZkeys:
                    #print("new key same as old ", Zkeynew, ziold, ziold+1)
                    ziold = ziold + 1 # see note on this iterator above

                    Zkeyold = oldZkeys[ziold]

                    del newset['yields'][Zkeyold] # to make consistent with new key
                    newset['yields'][Zkeynew] = {} # new dict

                    #newset[Zkeynew] = {}
                    # just average here , this Z already exists

                    for mi,Mkey in enumerate(list(newset['yields'][ list(newset['yields'].keys())[0] ].keys())):
                        newset['yields'][Zkeynew][Mkey] = {}

                        for ei,e in enumerate(newset['elem_list']):
                            newset['yields'][Zkeynew][Mkey][e] = v0yields['yields'][Zkeyold][Mkey][e] * f0 +\
                                                                 v150yields['yields'][Zkeyold][Mkey][e] * f150 +\
                                                                 v300yields['yields'][Zkeyold][Mkey][e] * f300


                else: # interpolate
                #    print("new key ", Zkeynew, ziold)
                    Zkeyold      = oldZkeys[ziold]
                    Zkeyold_next = oldZkeys[ziold+1]

                    newset['yields'][Zkeynew] = {}

                    # Z's are in descending order
                    if (znew > oldZ[ziold]) or (znew < oldZ[ziold+1]):
                        print(znew, oldZ[ziold], oldZ[ziold+1], ziold)
                        print("New Z is out of bounds")
                        raise RuntimeError

                    for mi,Mkey in enumerate(list(newset['yields'][ list(newset['yields'].keys())[0] ].keys())):
                        newset['yields'][Zkeynew][Mkey] = {}

                        for ei,e in enumerate(newset['elem_list']):

                            # see note on ziold iterator above
                            # fractional distance between points i and i + 1
                            #    interpval = t * array[i+1] + (1.0-t) * (array[i])
                            #    Z's are stored in descending order... so make negative
                            # t = (np.log10(znew) - np.log10(oldZ[ziold])) / (np.log10(oldZ[ziold+1]) - np.log10(oldZ[ziold]))
                            t = (znew - oldZ[ziold]) / (oldZ[ziold+1] - oldZ[ziold])


                            newset['yields'][Zkeynew][Mkey][e] =\
                              (t*v0yields['yields'][Zkeyold_next][Mkey][e] + (1.0-t)*v0yields['yields'][Zkeyold][Mkey][e]) * f0 +\
                              (t*v150yields['yields'][Zkeyold_next][Mkey][e] + (1.0-t)*v150yields['yields'][Zkeyold][Mkey][e]) * f150 +\
                              (t*v300yields['yields'][Zkeyold_next][Mkey][e] + (1.0-t)*v300yields['yields'][Zkeyold][Mkey][e]) * f300


            np.save("LC18" + typestr + "resampled_mixture_elem_stable.npy", newset, allow_pickle=True)

    else: # we do not need to resample. This is now easy

        for typestr in ['_','_winds_']: # do it for total and wind yields
            v0yields = np.load( wdir + available_yields['LC18' + typestr + 'R_0'], allow_pickle=True).item()
            v150yields = np.load( wdir + available_yields['LC18' + typestr + 'R_150'], allow_pickle=True).item()
            v300yields = np.load( wdir + available_yields['LC18' + typestr + 'R_300'], allow_pickle=True).item()

            newset = copy.deepcopy(v0yields)

            for zi,Zkey in enumerate(list(newset['yields'].keys())):
                f0   = prantzos_mixture( v0yields['Z_list'][zi], vtype='v0')
                f150 = prantzos_mixture( v0yields['Z_list'][zi], vtype='v150')
                f300 = prantzos_mixture( v0yields['Z_list'][zi], vtype='v300')

                #print(f0,f150,f300,f0+f150+f300)

                for mi,Mkey in enumerate(list(newset['yields'][ list(newset['yields'].keys())[0] ].keys())):
                    for ei,e in enumerate(newset['elem_list']):
                        newset['yields'][Zkey][Mkey][e] = v0yields['yields'][Zkey][Mkey][e] * f0 +\
                                                          v150yields['yields'][Zkey][Mkey][e] * f150 +\
                                                          v300yields['yields'][Zkey][Mkey][e] * f300

            np.save("LC18" + typestr + "mixture_elem_stable.npy", newset, allow_pickle=True)
    return

if __name__ == "__main__":

    make_mixture()
    make_mixture(resample=False)

    #test = np.load( 'test.npy', allow_pickle=True).item()

    #print(test.keys())
