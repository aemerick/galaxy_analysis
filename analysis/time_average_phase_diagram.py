from galaxy_analysis import Galaxy
import numpy as np
import yt
import glob
import os

parallel_on = yt.enable_parallelism()


def _create_region(ds, region_type, prop):

    if (region_type is None) or (region_type == 'FullBox') or (region_type == 'Fullbox'):
        if len(prop.keys()) > 0:
            print "Following keys do not do anything as full box was selected"
            print prop.keys()
            print "Proceeding anyway"

        region = ds.all_data()

    elif region_type == 'Disk' or region_type == 'disk':
        if not 'center' in prop.keys():
            prop['center'] = ds.domain_center
        if not 'normal' in prop.keys():
            prop['normal'] = [0,0,1]
        if (not 'radius' in prop.keys()) or (not 'height' in prop.keys()):
            print "If using disk region, must include radius and height"
            raise ValueError

        region = ds.disk(prop['center'], prop['normal'],
                         prop['radius'], prop['height'])

    elif (region_type == 'sphere' or region_type == 'Sphere'):
        if not 'center' in prop.keys():
            prop['center'] = ds.domain_center
        if not 'radius' in prop.keys():
            print "If using sphere, must include radius"
            raise ValueError

        region = ds.sphere(prop['center'], prop['radius'])

    elif ((region_type == 'Box') or (region_type == 'Region') or\
          (region_type == 'box') or (region_type == 'region')):

        region = ds.region(prop['center'], prop['left_edge'], prop['right_edge'])


    return region

def time_average_phase_diagram(tmin, tmax, wdir = './',
                               ds_list = None, xfield = 'number_density',yfield = 'temperature',
                               zfield = 'cell_mass', x_unit = None, y_unit = None,
                               weight_field = None,
                               x_bins=None, y_bins=None, outname = None,
                               zunit = None, zlim = None,
                               region_type = 'FullBox', region_kwargs = {},
                               xlog = True, ylog = True, zlog=False, cmap = 'cube_helix'):

    if (zfield == 'cell_mass') and (zunit is None):
        zunit = 'Msun'
    # filter based on times, unless ds list is already provided
    # this is really not fun when parallel computation is on, as it does
    # this for every processor. Preference is to use already filtered dataset
    # list. (could just do this on root and communicate, but I'm a lazy coder)
    if outname is None:
        outname = wdir + xfield + '_' + yfield + '_time_avg_PD.png'

    if ds_list is None:
        ds_list = np.sort(glob.glob(wdir + '/DD????/DD????'))
        times = np.zeros(len(ds_list))
        i = 0
        dataset_series = [None]*len(ds_list)
        for name in ds_list:
            dataset_series[i] = yt.load(name)
            times[i] = dataset_series[i].current_time.convert_to_units('Myr')
            i = i + 1
        dataset_series = np.array(dataset_series)
        dataset_series = dataset_series[ 1 * ((times >= tmin) * (times < tmax))]
    else:
        dataset_series = [None] * len(ds_list)
        for i,name in enumerate(ds_list):
            dataset_series[i] = yt.load(name)

    # if no xbins, make profile plot with most recent data set and use defaults
    if x_bins is None or y_bins is None:
        pd = yt.PhasePlot(dataset_series[-1].all_data(), xfield, yfield, zfield, weight_field=weight_field)

        if x_bins is None:
            if not x_unit is None:
                pd.profile.set_x_unit(x_unit)
            x_bins = pd.profile.x_bins

        if y_bins is None:
            if not y_unit is None:
                pd.profile.set_y_unit(y_unit)
            y_bins = pd.profile.y_bins

    # set units if not specified
    if x_unit is None:
        x_unit = x_bins.units
    if y_unit is None:
        y_unit = y_bins.units

    # set ranges
    ymin,ymax = np.min(y_bins), np.max(y_bins)
    xmin,xmax = np.min(x_bins), np.max(x_bins)
    x_bins = np.size(x_bins)
    y_bins = np.size(y_bins)

    def _set_axis_lim(plot):
        plot.set_log(xfield, xlog)
        plot.set_log(yfield, ylog)
        plot.set_unit(xfield, x_unit)
        plot.set_unit(yfield, y_unit)
        plot.set_xlim(xmin,xmax)
        plot.set_ylim(ymin,ymax)
        return


    i = 0
    phase_data = {}

    if parallel_on:
        num_dat = len(dataset_series)
        # iterate through all datasets except the first one, which will be used
        # later as the base plot. For whatever reason including the first dataset here
        # and then re-computing it later on (see below) causes a hang on the root processor
        # this is a bit hacky, but whatever... it works...
        DS = yt.DatasetSeries(dataset_series[1:], parallel = True)

        for sto, dataset in DS.piter(storage=phase_data):
            region = _create_region(dataset, region_type, region_kwargs)

            pd = yt.PhasePlot(region, xfield, yfield, zfield,
                                 weight_field=weight_field, x_bins = x_bins, y_bins = y_bins)
            _set_axis_lim(pd)

            sto.result = pd.profile.field_data
            sto.result_id = str(dataset)


    else:
        # loop through all data sets in a consistent way as the above
        for dataset in dataset_series[1:]:
            pd = yt.PhasePlot(dataset.all_data(), xfield, yfield, zfield,
                                 weight_field=weight_field, x_bins = x_bins, y_bins = y_bins)
            _set_axis_lim(pd)
            phase_data[str(dataset)] = pd.profile.field_data
            del(dataset)
        num_dat = len(dataset_series)

    # now save the data and plot
    # for whatever reason, this fails if re-making a PD from a dataset that
    # was already used when parallelism is ON (not sure why)
    region = _create_region(dataset_series[0], region_type, region_kwargs)

    main_pd = yt.PhasePlot(region, xfield, yfield, zfield,
                             weight_field=weight_field, x_bins = x_bins, y_bins = y_bins)
    _set_axis_lim(main_pd)

    # set the color field properties and limits:
    main_pd.set_cmap(zfield, cmap)
    if ('gas',zfield) in dataset_series[0].derived_field_list:
        zf = ('gas',zfield)
    else:
        zf = zfield

    for index in np.arange(1, num_dat):
        main_pd.profile.field_data[zf] +=\
            phase_data[str(dataset_series[index])][zf]
    main_pd.profile.field_data[zf] /= (1.0 * num_dat)
    if not (zunit is None):
        main_pd.set_unit(zfield, zunit)
    main_pd.set_log(zfield, zlog)
    if not (zlim is None):
        main_pd.set_zlim(zfield, zlim[0], zlim[1])

    # only write one image
    if yt.is_root():
        main_pd.save(outname)

    return main_pd


def spatial_example():

    imin, imax = 124, 126
    ds_list = [None]*len(np.arange(imin,imax,1))
    j = 0
    for x in np.arange(imin,imax,1):
        name = 'DD%0004i/DD%0004i'%(x,x)
        if os.path.isfile(name):
            ds_list[j] = name
            j = j + 1

    ds_list = ds_list[:j]


    # use galaxy class to get consistent region parameters
    region_type = 'disk'
    gal = Galaxy( (ds_list[0]).split('/')[0] )
    region_kwargs = {}
    for k in ['center','normal','radius','height']:
        region_kwargs[k] = gal.disk_region[k]
    nbin = 100
    r = gal.disk_region['radius'].convert_to_units('pc').value
    z = gal.disk_region['height'].convert_to_units('pc').value * 0.5
    xbins = np.linspace(0.0, r, nbin)*yt.units.pc
    ybins = np.linspace(-z, z, nbin*(2.0*z)/r)*yt.units.pc

    del(gal) # don't need anymore

    pd = time_average_phase_diagram(120.0, 130.0, ds_list = ds_list,
                                    xfield = 'radius', yfield = 'cylindrical_z',
                                    x_bins = xbins, y_bins = ybins,
                                    zfield = 'Fe_over_H', zunit = None, zlim = [-7,1.5], cmap = 'algae',
                                    region_type   = region_type, weight_field = 'cell_mass',
                                    region_kwargs = region_kwargs, xlog=False, ylog=False, zlog=False)


    return

def n_T_example():

    imin, imax = 101, 282
    ds_list = [None]*len(np.arange(imin,imax,1))
    j = 0
    for x in np.arange(imin,imax,1):
        name = 'DD%0004i/DD%0004i'%(x,x)

        if os.path.isfile(name):
            ds_list[j] = name
            j = j + 1

    ds_list = ds_list[:j]

    pd = time_average_phase_diagram(290.0, 300.0, ds_list = ds_list,
                                    x_bins = np.logspace(-4,3,128)*yt.units.cm**(-3),
                                    y_bins = np.logspace(0,7,128)*yt.units.K,
                                    region_type = 'disk',
                                    region_kwargs = {'center' : [0.5,0.5,0.5],
                                                     'normal' : [0,0,1],
                                                     'radius' : (4,'kpc'),
                                                     'height' : (1,'kpc')  })
    return

if __name__ == "__main__":


#    n_T_example()
    spatial_example()
