'''

  Software which uses the MHW definition
  of Hobday et al. (2015) applied to 
  select SST time series around the globe

'''

# Load required modules

import numpy as np
from scipy import io
from scipy import linalg
from scipy import stats
from scipy import interpolate
from scipy import signal
from datetime import date
from netCDF4 import Dataset

from matplotlib import pyplot as plt
import mpl_toolkits.basemap as bm

import marineHeatWaves as mhw

import ecoliver as ecj

#
# observations
#

pathroot = '/data/home/oliver/'
header = pathroot+'data/sst/noaa_oi_v2/avhrr/'
file0 = header + '1982/avhrr-only-v2.19820101.nc'
t, dates, T, year, month, day, doy = ecj.timevector([1982, 1, 1], [2017, 12, 31])
tt_1982_2005 = year <= 2005
t, dates, T, year, month, day, doy = ecj.timevector([1982, 1, 1], [2005, 12, 31])

#
# lat and lons of obs
#

fileobj = Dataset(file0, mode='r')
lon = fileobj.variables['lon'][:].astype(float)
lat = fileobj.variables['lat'][:].astype(float)
fill_value = fileobj.variables['sst']._FillValue.astype(float)
scale = fileobj.variables['sst'].scale_factor.astype(float)
offset = fileobj.variables['sst'].add_offset.astype(float)
fileobj.close()

#
# Size of mhwBlock variable
#

matobj = io.loadmat(header + 'timeseries/avhrr-only-v2.ts.' + str(300).zfill(4) + '.mat')
mhws, clim = mhw.detect(t, matobj['sst_ts'][300,tt_1982_2005])
mhwBlock = mhw.blockAverage(t, mhws)
years = mhwBlock['years_centre']
NB = len(years)

#
# initialize some variables
#

pctile = 90 # Percentile for calculation of MHWs
X = len(lon)
Y = len(lat)
i_which = range(0,X)
j_which = range(0,Y)
DIM = (len(j_which), len(i_which))
SST_mean = np.NaN*np.zeros(DIM)
MHW_total = np.NaN*np.zeros(DIM)
MHW_cnt = np.NaN*np.zeros(DIM)
MHW_dur = np.NaN*np.zeros(DIM)
MHW_max = np.NaN*np.zeros(DIM)
MHW_mean = np.NaN*np.zeros(DIM)
MHW_cum = np.NaN*np.zeros(DIM)
MHW_var = np.NaN*np.zeros(DIM)
MHW_td = np.NaN*np.zeros(DIM)
MHW_tc = np.NaN*np.zeros(DIM)
N_ts = np.zeros((len(j_which), len(i_which), NB))
SST_ts = np.zeros((len(j_which), len(i_which), NB))
MHW_cnt_ts = np.zeros((len(j_which), len(i_which), NB))
MHW_dur_ts = np.zeros((len(j_which), len(i_which), NB))
MHW_max_ts = np.zeros((len(j_which), len(i_which), NB))
MHW_mean_ts = np.zeros((len(j_which), len(i_which), NB))
MHW_cum_ts = np.zeros((len(j_which), len(i_which), NB))
MHW_var_ts = np.zeros((len(j_which), len(i_which), NB))
MHW_td_ts = np.zeros((len(j_which), len(i_which), NB))
MHW_tc_ts = np.zeros((len(j_which), len(i_which), NB))
lon_map =  np.NaN*np.zeros(len(i_which))
lat_map =  np.NaN*np.zeros(len(j_which))

#
# loop through locations
#

# Tropical Pacific: i = 800, j = 360
# WA: i = 450, j = 260
icnt = 0
for i in i_which:
    print i, 'of', len(lon)-1
#   load SST
    matobj = io.loadmat(header + 'timeseries/avhrr-only-v2.ts.' + str(i+1).zfill(4) + '.mat')
    sst_ts = matobj['sst_ts']
    lon_map[icnt] = lon[i]
#   loop over j
    jcnt = 0
    for j in j_which:
        lat_map[jcnt] = lat[j]
        sst = sst_ts[j,tt_1982_2005].copy()
        #if np.logical_not(np.isfinite(sst.sum())): # check for land
        if np.logical_not(np.isfinite(sst.sum())) + ((sst<-1).sum()>0): # check for land, ice
            jcnt += 1
            continue
        # MHW detection
        mhws, clim = mhw.detect(t, sst, pctile=pctile)
        mhwBlock = mhw.blockAverage(t, mhws, temp=sst)
        # Total count
        MHW_total[jcnt,icnt] = mhwBlock['count'].sum()
        # Mean and trend
        mean, trend, dtrend = mhw.meanTrend(mhwBlock)
        # Mean and trend
        MHW_cnt[jcnt,icnt], tmp, tmp = mean['count'], trend['count'], dtrend['count']
        MHW_dur[jcnt,icnt], tmp, tmp = mean['duration'], trend['duration'], dtrend['duration']
        MHW_max[jcnt,icnt], tmp, tmp = mean['intensity_max_max'], trend['intensity_max_max'], dtrend['intensity_max_max']
        MHW_mean[jcnt,icnt], tmp, tmp = mean['intensity_mean'], trend['intensity_mean'], dtrend['intensity_mean']
        MHW_cum[jcnt,icnt], tmp, tmp = mean['intensity_cumulative'], trend['intensity_cumulative'], dtrend['intensity_cumulative']
        MHW_var[jcnt,icnt], tmp, tmp = mean['intensity_var'], trend['intensity_var'], dtrend['intensity_var']
        MHW_td[jcnt,icnt], tmp, tmp = mean['total_days'], trend['total_days'], dtrend['total_days']
        MHW_tc[jcnt,icnt], tmp, tmp = mean['total_icum'], trend['total_icum'], dtrend['total_icum']
        SST_mean[jcnt,icnt], tmp, tmp = mean['temp_mean'], trend['temp_mean'], dtrend['temp_mean']
        # Time series
        MHW_cnt_ts[jcnt,icnt,:] += mhwBlock['count']
        MHW_dur_ts[jcnt,icnt,np.where(~np.isnan(mhwBlock['duration']))[0]] = mhwBlock['duration'][np.where(~np.isnan(mhwBlock['duration']))[0]]
        MHW_max_ts[jcnt,icnt,np.where(~np.isnan(mhwBlock['intensity_max_max']))[0]] = mhwBlock['intensity_max_max'][np.where(~np.isnan(mhwBlock['intensity_max_max']))[0]]
        MHW_mean_ts[jcnt,icnt,np.where(~np.isnan(mhwBlock['intensity_mean']))[0]] = mhwBlock['intensity_mean'][np.where(~np.isnan(mhwBlock['intensity_mean']))[0]]
        MHW_cum_ts[jcnt,icnt,np.where(~np.isnan(mhwBlock['intensity_cumulative']))[0]] = mhwBlock['intensity_cumulative'][np.where(~np.isnan(mhwBlock['intensity_cumulative']))[0]]
        MHW_var_ts[jcnt,icnt,np.where(~np.isnan(mhwBlock['intensity_var']))[0]] = mhwBlock['intensity_var'][np.where(~np.isnan(mhwBlock['intensity_var']))[0]]
        MHW_td_ts[jcnt,icnt,np.where(~np.isnan(mhwBlock['total_days']))[0]] = mhwBlock['total_days'][np.where(~np.isnan(mhwBlock['total_days']))[0]]
        MHW_tc_ts[jcnt,icnt,np.where(~np.isnan(mhwBlock['total_icum']))[0]] = mhwBlock['total_icum'][np.where(~np.isnan(mhwBlock['total_icum']))[0]]
        N_ts[jcnt,icnt,:] += (~np.isnan(mhwBlock['duration'])).astype(int)
        SST_ts[jcnt,icnt,:] = mhwBlock['temp_mean']
        # Up counts
        jcnt += 1
    icnt += 1
    # Save data so far
    outfile = '/data/home/oliver/data/MHWs/Projections/mhw_census.NOAAOISST.1982_2005'
    if (i % 100) + (i == i_which[-1]):
        np.savez(outfile, lon_map=lon_map, lat_map=lat_map, SST_mean=SST_mean, MHW_total=MHW_total, MHW_cnt=MHW_cnt, MHW_dur=MHW_dur, MHW_max=MHW_max, MHW_mean=MHW_mean, MHW_cum=MHW_cum, MHW_var=MHW_var, MHW_td=MHW_td, MHW_tc=MHW_tc, SST_ts=SST_ts, MHW_cnt_ts=MHW_cnt_ts, MHW_dur_ts=MHW_dur_ts, MHW_max_ts=MHW_max_ts, MHW_mean_ts=MHW_mean_ts, MHW_cum_ts=MHW_cum_ts, MHW_var_ts=MHW_var_ts, MHW_td_ts=MHW_td_ts, MHW_tc_ts=MHW_tc_ts, N_ts=N_ts, years=years)


