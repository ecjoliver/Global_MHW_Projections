'''

  Generate MHW census from CMIP5 daily 'tos' data

'''

import numpy as np
import scipy as sp
from scipy import stats
from scipy import signal
from scipy import io
import scipy.optimize as opt
from datetime import date
import os
import sys
import ecoliver as ecj
import deseason as ds

from matplotlib import pyplot as plt
import mpl_toolkits.basemap as bm

import marineHeatWaves as mhw

#
# Load data
#

# Which climate models
models = ['ACCESS1-3', 'CSIRO-Mk3-6-0', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'CanESM2', 'CNRM-CM5']
models = ['IPSL-CM5A-MR', 'IPSL-CM5A-LR', 'CanESM2', 'CSIRO-Mk3-6-0', 'ACCESS1-3', 'HadGEM2-ES', 'CNRM-CM5'] # Re-order with fastest first
experiments = ['hist', 'rcp45', 'rcp85', 'histNat']
#experiments = ['hist', 'rcp85']

# Load in meta-data about CMIP5 time series
header_CMIP5ts = '/media/ecoliver/Pentagram/data/CMIP5/'
header_CMIP5ts = '/data/home/oliver/data/CMIP5/'
outheader = '/data/home/oliver/data/MHW/Projections/'
metafile = header_CMIP5ts + 'meta/timeSeries_metaData.npz'
meta = np.load(metafile)
Ly = meta['Ly'].item()
t = meta['t'].item()
years = meta['years'].item()
Ty = meta['Ty'].item()
NEXPS = meta['NEXPS'].item()
NENS = meta['NENS'].item()
modelPaths = meta['modelPaths'].item()
meta.close()

T = {}
year = {}
for model in models:
    T[model] = {}
    year[model] = {}
    for exp in experiments:
        T[model][exp] = len(t[model][exp])
        year[model][exp] = [date.fromordinal(t[model][exp][tt]).year for tt in range(T[model][exp])]

#
# Loop over models, lons, and calculate MHW stats
#

mhw_keys = ['count', 'count_new', 'duration', 'duration_new', 'intensity_max_max', 'intensity_max_max_new', 'total_days', 'moderate_days', 'strong_days', 'severe_days', 'extreme_days', 'temp_mean', 'temp_max']
#climPeriod = [1861,1890]
#climPeriod = [1961,1990]
climPeriod = [1982,2005]

sst = {}
# Loop over models
for model in models:
    filename = 'sst_hist_ens_00_lon_0001.mat'
    matobj = io.loadmat(header_CMIP5ts + modelPaths[model] + 'timeseries/' + filename)
    X = matobj['X'][0][0]
    Y = matobj['Y'][0][0]
    lon = matobj['lon']
    lat = matobj['lat']
    MHW_ts = {}
    # Initialise MHW time series variables
    for key in mhw_keys:
        MHW_ts[key] = {}
        for exp in experiments:
            MHW_ts[key][exp] = np.NaN*np.zeros((Y, X, Ty[model][exp], NENS[model][exp]))
    # Loop over lon
    for i in range(X):
        print model, ',', i, '/', X
        # Loop over experiments, ensemble members
        for exp in experiments:
            sst[exp] = np.zeros((T[model][exp], Y, NENS[model][exp]))
            for ens in range(NENS[model][exp]):
                # SSTs
                filename = 'sst_' + exp + '_ens_' + str(ens).zfill(2) + '_lon_' + str(i+1).zfill(4) + '.mat'
                matobj = io.loadmat(header_CMIP5ts + modelPaths[model] + 'timeseries/' + filename)
                sst[exp][:,:,ens] = matobj['sst']
                # Loop over lats and do the calculations
                for j in range(Y):
                    which_pre2100 = np.array(year[model][exp]) < 2100
                    if np.logical_not(np.isfinite(sst[exp][which_pre2100,j,ens].sum())) + ((sst[exp][which_pre2100,j,ens]<-0).sum()>0): # check for land, ice
                        continue
                    # MHW detection
                    if (exp == 'histNat') + (exp == 'rcp45') + (exp == 'rcp85'):
                        if ens > NENS[model]['hist']-1:
                            ens_clim = NENS[model]['hist']-1
                        else:
                            ens_clim = ens
                    else:
                        ens_clim = ens
                    mhws, clim = mhw.detect(t[model][exp][which_pre2100], sst[exp][which_pre2100,j,ens], climatologyPeriod=climPeriod, alternateClimatology=[t[model]['hist'], sst['hist'][:,j,ens_clim]])
                    mhwBlock = mhw.blockAverage(t[model][exp][which_pre2100], mhws, temp=sst[exp][which_pre2100,j,ens], clim=clim)
                    #else:
                    #    mhws, clim = mhw.detect(t[model][exp], sst[exp][:,j,ens], climatologyPeriod=climPeriod)
                    #    mhwBlock = mhw.blockAverage(t[model][exp], mhws, temp=sst[exp][:,j,ens], clim=clim)
                    # Calculate new measure of annual duration, intensity
                    mhwBlock['count_new'] = np.zeros(mhwBlock['count'].shape)
                    mhwBlock['duration_new'] = np.zeros(mhwBlock['duration'].shape)
                    mhwBlock['intensity_max_max_new'] = np.zeros(mhwBlock['intensity_max_max'].shape)
                    for ev in range(mhws['n_events']):
                        # Block index for year of each MHW (MHW year defined by start year)
                        iBlock = np.where((mhwBlock['years_start'] >= mhws['date_start'][ev].year) * (mhwBlock['years_end'] <= mhws['date_end'][ev].year))[0]
                        # Add MHW properties to block count
                        mhwBlock['count_new'][iBlock] += 1
                        mhwBlock['duration_new'][iBlock] += mhws['duration'][ev]
                        for iB in iBlock:
                            mhwBlock['intensity_max_max_new'][iB] = np.max([mhwBlock['intensity_max_max_new'][iB], mhws['intensity_max'][ev]])
                    count = 1.*mhwBlock['count_new']
                    count[count==0] = np.nan
                    mhwBlock['duration_new'] = mhwBlock['duration_new'] / count
                    mhwBlock['intensity_max_max_new'][np.isnan(mhwBlock['intensity_max'])] = np.nan
                    # Save some statistics
                    TBlock = len(mhwBlock['years_centre'])
                    for key in mhw_keys:
                        MHW_ts[key][exp][j,i,:TBlock,ens] = mhwBlock[key]
    # Save data
    outfile = outheader + '/mhw_census_' + model + '_climPeriod_' + str(climPeriod[0]) + '_' + str(climPeriod[1]) + '_histClim'
    np.savez(outfile, lon=lon, lat=lat, MHW_ts=MHW_ts, years=years[model], X=X, Y=Y, Ty=Ty[model], NEXPS=NEXPS, NENS=NENS[model], models=models, experiments=experiments)



