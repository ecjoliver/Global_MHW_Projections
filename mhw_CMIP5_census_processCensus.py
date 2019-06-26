'''

  MHW time series from CMIP5 GCMs

'''

import numpy as np
import scipy.signal as sig
from scipy import linalg
from scipy import stats
from scipy import io
from scipy import interpolate as interp

import ecoliver as ecj

from matplotlib import pyplot as plt
import mpl_toolkits.basemap as bm


models = ['ACCESS1-3', 'CSIRO-Mk3-6-0', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'CanESM2', 'CNRM-CM5']
#models = ['IPSL-CM5A-MR']
Ly = {'ACCESS1-3': 365.25, 'CSIRO-Mk3-6-0':365, 'HadGEM2-ES': 360, 'IPSL-CM5A-LR': 365, 'IPSL-CM5A-MR': 365, 'CanESM2': 365, 'CNRM-CM5': 365.25}
header = '/home/oliver/data/MHWs/Projections/'
dl = 2 # 0.5 # Horizontal resolution of regridded data

for model in models:
    # Load in data
    #infile = 'mhw_census_' + model + '.npz'
    #infile = 'mhw_census_' + model + '_climPeriod_1982_2005.npz'
    infile = 'mhw_census_' + model + '_climPeriod_1861_1890.npz'
    infile = 'mhw_census_' + model + '_climPeriod_1861_1890_natClim.npz'
    infile = 'mhw_census_' + model + '_climPeriod_1982_2005_histClim.npz'
    data = np.load(header + infile)
    lon = data['lon']
    lat = data['lat']
    years = data['years'].item()
    X = data['X']
    Y = data['Y']
    Ty = data['Ty'].item()
    NENS = data['NENS'].item()
    NEXPS = data['NEXPS']
    models = data['models']
    experiments = data['experiments']
    MHW_ts_orig = data['MHW_ts'].item()
    if (model == 'HadGEM2-ES') + (model == 'CanESM2') + (model == 'CSIRO-Mk3-6-0'):
        lon, lat = np.meshgrid(lon, lat)
    # Variables needed for the calculates
    #mhw_keys = ['count', 'duration', 'intensity_max_max', 'intensity_mean', 'intensity_cumulative', 'total_days', 'temp_mean']
    mhw_keys = ['count', 'count_new', 'duration', 'duration_new', 'intensity_max_max', 'intensity_max_max_new', 'total_days', 'moderate_days', 'strong_days', 'severe_days', 'extreme_days', 'temp_mean', 'temp_max']
    meanYears = {}
    meanYears['hist'] = [1982, 2005] #[1961, 1990]
    meanYears['histNat'] = [1982, 2005] #[1961, 1990]
    meanYears['rcp45'] = [2031, 2060]
    meanYears['rcp85'] = [2031, 2060]
    # New regular grid, shared across all processed models
    lon_regrid = np.arange(0, 360, dl)
    lat_regrid = np.arange(-90, 90+0.5, dl)
    lon_regrid, lat_regrid = np.meshgrid(lon_regrid, lat_regrid)
    Y_regrid, X_regrid = lon_regrid.shape
    # Load iceMean map for making mask
    matobj = io.loadmat(header + '../Trends/NOAAOISST_iceMean.mat')
    ice_longestRun = matobj['ice_longestRun']
    lon_ice, lat_ice = np.meshgrid(matobj['lon'][:,0], matobj['lat'][:,0])
    # Make data mask based on land and ice
    ice_longestRun = interp.griddata((lon_ice.flatten(), lat_ice.flatten()), ice_longestRun.flatten() , (lon_regrid, lat_regrid), method='nearest')
    icemask = np.ones(ice_longestRun.shape)
    icemask[ice_longestRun>=6.] = np.nan # mask where ice had runs of 6 days or longer
    icemask[lat_regrid<=-65.] = np.nan
    # 
    # Processed blocked census data, produce regridded output
    #
    MHW_ts = {}
    MHW_mean = {}
    for metric in mhw_keys:
        MHW_ts[metric] = {}
        MHW_mean[metric] = {}
        if (metric == 'total_days') + (metric == 'moderate_days') + (metric == 'strong_days') + (metric == 'severe_days') + (metric == 'extreme_days'): # Permanent MHW state
            MHW_ts[metric + '_permMHW'] = {}
            MHW_mean[metric + '_permMHW'] = {}
        # Generate [lat,lon,time] array for mhw metric of interest
        MHW_ts_regrid = {}
        for exp in experiments:
            #MHW_ts_full = np.zeros((Y, X, Ty[exp], NENS[exp]))
            MHW_ts_regrid[exp] = np.zeros((Y_regrid, X_regrid, Ty[exp], NENS[exp]))
            icemask_ens = np.swapaxes(np.swapaxes(np.tile(icemask, (NENS[exp],1,1)), 0, 1), 1, 2)
            icemask_full = np.swapaxes(np.swapaxes(np.tile(icemask_ens, (Ty[exp],1,1,1)), 0, 1), 1, 2)
            print model, metric, exp
            for ens in range(NENS[exp]):
        # Re-grid
                for tt in range(MHW_ts_orig[metric][exp].shape[2]):
                    MHW_ts_regrid[exp][:,:,tt,ens] = interp.griddata((lon.flatten(), lat.flatten()), MHW_ts_orig[metric][exp][:,:,tt,ens].flatten() , (lon_regrid, lat_regrid), method='nearest')
        # Time-mean maps
            tt = (years[exp] >= meanYears[exp][0]) * (years[exp] <= meanYears[exp][1])
            MHW_mean[metric][exp] = {}
            MHW_mean[metric][exp]['full'] = np.nanmean(MHW_ts_regrid[exp][:,:,tt,:], axis=2)
            MHW_mean[metric][exp]['icemask'] = np.nanmean(MHW_ts_regrid[exp][:,:,tt,:], axis=2)*icemask_ens
        # Calculate global averaged time series
            MHW_ts[metric][exp] = {}
            MHW_ts[metric][exp]['full'] = np.nanmean(np.nanmean(MHW_ts_regrid[exp], axis=0), axis=0)
            MHW_ts[metric][exp]['icemask'] = np.nanmean(np.nanmean(MHW_ts_regrid[exp]*icemask_full, axis=0), axis=0)
        # Permanent MHW state
            if (metric == 'total_days') + (metric == 'moderate_days') + (metric == 'strong_days') + (metric == 'severe_days') + (metric == 'extreme_days'):
                mask = np.ones(MHW_ts_regrid[exp].shape)
                mask[np.isnan(MHW_ts_regrid[exp])] = np.nan
                permMHW = (MHW_ts_regrid[exp] >= np.floor(Ly[model]-1))*mask
                MHW_ts[metric + '_permMHW'][exp] = {} # Proportion of globe in permanent MHW state
                MHW_ts[metric + '_permMHW'][exp]['full'] = np.nansum(np.nansum(permMHW, axis=0), axis=0) / np.nansum(np.nansum(mask, axis=0), axis=0)
                MHW_ts[metric + '_permMHW'][exp]['icemask'] = np.nansum(np.nansum(permMHW*icemask_full, axis=0), axis=0) / np.nansum(np.nansum(mask*icemask_full, axis=0), axis=0)
                MHW_mean[metric + '_permMHW'][exp] = {} # First year of permanent MHW state
                MHW_mean[metric + '_permMHW'][exp]['full'] = np.nan*np.zeros(permMHW[:,:,0,:].shape)
                for i in range(permMHW.shape[1]):
                    for j in range(permMHW.shape[0]):
                        for ens in range(NENS[exp]):
                            tt = np.where(permMHW[j,i,:,ens]==1)[0]
                            if len(tt) >= 1:
                                MHW_mean[metric + '_permMHW'][exp]['full'][j,i,ens] = years[exp][tt[0]]
                MHW_mean[metric + '_permMHW'][exp]['icemask'] = MHW_mean[metric + '_permMHW'][exp]['full']*icemask_ens

    #
    # Save
    #
    #outfile = 'mhw_census_' + model + '.regrid'
    #outfile = 'mhw_census_' + model + '_climPeriod_1982_2005.regrid'
    outfile = 'mhw_census_' + model + '_climPeriod_1861_1890.regrid'
    outfile = 'mhw_census_' + model + '_climPeriod_1861_1890_natClim.regrid'
    outfile = 'mhw_census_' + model + '_climPeriod_1982_2005_histClim.regrid'
    outfile = 'mhw_census_' + model + '_climPeriod_1982_2005_histClim.regrid.2deg'
    np.savez(header + outfile, lon_regrid=lon_regrid, lat_regrid=lat_regrid, X_regrid=X_regrid, Y_regrid=Y_regrid, years=years, MHW_mean=MHW_mean, MHW_ts=MHW_ts, NENS=NENS, models=models, experiments=experiments, icemask=icemask, meanYears=meanYears)


