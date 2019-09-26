import yt
import glob.glob


from galaxy_analysis.analysis import Galaxy


ds_list = np.sort(glob.glob('DD????_galaxy_data.h5'))

for ds in ds_list:
