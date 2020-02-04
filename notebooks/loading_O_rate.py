import numpy as np
import glob
import deepdish as dd

from galaxy_analysis.utilities.utilities import simple_rebin

def compute_orate(directory):
    
    files     = np.sort(glob.glob(directory + '/*galaxy_data*.h5'))
    all_times  = np.zeros(np.size(files))

    for i,f in enumerate(files):
        meta_data = dd.io.load(f,'/meta_data')
        all_times[i] = meta_data['Time']

    Omass = np.zeros(np.size(all_times))
    for i,t in enumerate(all_times):
        Omass[i] = dd.io.load(files[i],
                             '/gas_meta_data/masses/FullBox')['O'] +\
                   dd.io.load(files[i],
                              '/gas_meta_data/masses/OutsideBox')['O']

    Orate = np.zeros(np.size(all_times))
    Orate[:-1] = (Omass[1:] - Omass[:-1]) / (all_times[1:] - all_times[:-1]) / 1.0E6 # want Msun / yr
    Orate[-1]  = Orate[-2]

    # now rebin
#    temp_all_times, temp_Orate = simple_rebin(np.round(all_times,0), Orate[:-1], 10.0, 'average', force = True)

#    if len(temp_Orate) <= 2:
#        temp_all_times, temp_Orate = simple_rebin(np.round(all_times,0), Orate[:-1], 25.0, 'average', force = True)
#    all_times = 1.0 * temp_all_times
#    Orate = 1.0 * temp_Orate

#    print(all_times, Orate)

    f = open("orate.dat","w")
    f.write("#time O_rate\n")

    Orate_out = np.zeros(np.size(Orate))
    for i in np.arange(np.size(Orate)):

        Orate_out[i] = np.average( Orate[((all_times > all_times[i] - 50.0)) * (all_times <all_times[i])])

#        if Orate[i] < 0:
#            Orate[i] = np.average([ Orate[i-1], np.max([Orate[i+1],0.0])])
    Orate_out[0] = 1.0 * Orate_out[1]

    for i in np.arange(np.size(Orate)):
        f.write("%5.5E %5.5E\n"%(all_times[i], Orate_out[i]))
    f.close()

    return

if __name__ == "__main__":

    compute_orate("./")
