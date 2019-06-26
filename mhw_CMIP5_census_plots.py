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

import matplotlib as mpl
from matplotlib import pyplot as plt
import mpl_toolkits.basemap as bm


#
# Load data and make plots
#

header = '/home/oliver/data/MHWs/Projections/'
figfolder = 'figures/'

#model = 'ACCESS1-3'
#model = 'IPSL-CM5A-LR'
#model = 'IPSL-CM5A-MR'
#model = 'HadGEM2-ES'
#model = 'CanESM2'
#model = 'CNRM-CM5'
#model = 'CSIRO-Mk3-6-0'
#model = 'ALL'
models = ['ACCESS1-3', 'CSIRO-Mk3-6-0', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'CanESM2'] #, 'CNRM-CM5']
modelDir = 'ALL'

#basePeriod = [1961, 1990]
#basePeriod = [1861, 1890]
basePeriod = [1982, 2005]
rcp = 'rcp85'
natClim = False
histClim = True

years = {}
NENS = {}
MHW_ts = {}
MHW_mean = {}
for model in models:
    print model
    outfile = 'mhw_census_' + model
    #data = np.load(header + outfile + '.regrid.npz')
    if natClim:
        data = np.load(header + outfile + '_climPeriod_' + str(basePeriod[0]) + '_' + str(basePeriod[1]) + '_natClim.regrid.npz')
    elif histClim:
        #data = np.load(header + outfile + '_climPeriod_' + str(basePeriod[0]) + '_' + str(basePeriod[1]) + '_histClim.regrid.npz')
        data = np.load(header + outfile + '_climPeriod_' + str(basePeriod[0]) + '_' + str(basePeriod[1]) + '_histClim.regrid.2deg.npz')
    else:
        data = np.load(header + outfile + '_climPeriod_' + str(basePeriod[0]) + '_' + str(basePeriod[1]) + '.regrid.npz')
    lon = data['lon_regrid']
    lat = data['lat_regrid']
    years[model] = data['years'].item()
    meanYears = data['meanYears'].item()
    NENS[model] = data['NENS'].item()
    experiments = data['experiments']
    MHW_ts[model] = data['MHW_ts'].item()
    MHW_mean[model] = data['MHW_mean'].item()
    mhw_metrics = MHW_mean[model].keys()
    # Model-specific fixes
    if model == 'HadGEM2-ES':
        for metric in mhw_metrics:
            exp = 'hist'
            for icemask in MHW_ts[model][metric][exp].keys():
                MHW_ts[model][metric][exp][icemask] = np.delete(MHW_ts[model][metric][exp][icemask], [0, 3], axis=1)
                MHW_mean[model][metric][exp][icemask] = np.delete(MHW_mean[model][metric][exp][icemask], [0, 3], axis=2)
                NENS[model][exp] = MHW_ts[model][metric][exp][icemask].shape[1]
            exp = 'histNat'
            for icemask in MHW_ts[model][metric][exp].keys():
                MHW_ts[model][metric][exp][icemask] = np.delete(MHW_ts[model][metric][exp][icemask], [0, 2], axis=1)
                MHW_mean[model][metric][exp][icemask] = np.delete(MHW_mean[model][metric][exp][icemask], [0, 2], axis=2)
                NENS[model][exp] = MHW_ts[model][metric][exp][icemask].shape[1]
    if model == 'CNRM-CM5':
        for metric in mhw_metrics:
            exp = 'hist'
            for icemask in MHW_ts[model][metric][exp].keys():
                MHW_ts[model][metric][exp][icemask] = np.delete(MHW_ts[model][metric][exp][icemask], [0,1,2,3,4,5,6,8,9], axis=1)
                MHW_mean[model][metric][exp][icemask] = np.delete(MHW_mean[model][metric][exp][icemask], [1,2,3,4,5,6,8,9], axis=2)
                NENS[model][exp] = MHW_ts[model][metric][exp][icemask].shape[1]
            exp = 'histNat'
            for icemask in MHW_ts[model][metric][exp].keys():
                MHW_ts[model][metric][exp][icemask] = np.delete(MHW_ts[model][metric][exp][icemask], [0], axis=1)
                MHW_mean[model][metric][exp][icemask] = np.delete(MHW_mean[model][metric][exp][icemask], [0], axis=2)
                NENS[model][exp] = MHW_ts[model][metric][exp][icemask].shape[1]
            exp = 'rcp85'
            for icemask in MHW_ts[model][metric][exp].keys():
                MHW_ts[model][metric][exp][icemask] = np.delete(MHW_ts[model][metric][exp][icemask], [1,2,3,4], axis=1)
                MHW_mean[model][metric][exp][icemask] = np.delete(MHW_mean[model][metric][exp][icemask], [1,2,3,4], axis=2)
                NENS[model][exp] = MHW_ts[model][metric][exp][icemask].shape[1]
    if model == 'CSIRO-Mk3-6-0':
        for metric in mhw_metrics:
            exp = 'histNat'
            for icemask in MHW_ts[model][metric][exp].keys():
                MHW_ts[model][metric][exp][icemask] = np.delete(MHW_ts[model][metric][exp][icemask], [8], axis=1)
                MHW_mean[model][metric][exp][icemask] = np.delete(MHW_mean[model][metric][exp][icemask], [8], axis=2)
                NENS[model][exp] = MHW_ts[model][metric][exp][icemask].shape[1]
    ## Centre temp_mean about basePeriod mean
    exp = 'hist'
    tt = (years[model][exp] >= basePeriod[0]) * (years[model][exp] <= basePeriod[1])
    climMean = {}
    climMean['icemask'] = np.zeros(NENS[model][exp])
    climMean['full'] = np.zeros(NENS[model][exp])
    for ens in range(NENS[model][exp]):
        # Calculate mean over base period
        climMean['icemask'][ens] = np.mean(MHW_ts[model]['temp_mean'][exp]['icemask'][tt,ens])
        climMean['full'][ens] = np.mean(MHW_ts[model]['temp_mean'][exp]['full'][tt,ens])
        # Correct time series
        #MHW_ts[model]['temp_mean'][exp]['icemask'][:,ens] = MHW_ts[model]['temp_mean'][exp]['icemask'][:,ens] - climMean['icemask'][ens]
        #MHW_ts[model]['temp_mean'][exp]['full'][:,ens] = MHW_ts[model]['temp_mean'][exp]['full'][:,ens] - climMean['full'][ens]
    # Now correct all runs with same offset
    #exp = rcp
    for exp in experiments:
        for ens in range(NENS[model][exp]):
            # Sometimes there are more rcp85 ens members than hist
            if ens > NENS[model]['hist']-1:
                ens_clim = NENS[model]['hist']-1
            else:
                ens_clim = ens
            # Correct time series
            MHW_ts[model]['temp_mean'][exp]['icemask'][:,ens] = MHW_ts[model]['temp_mean'][exp]['icemask'][:,ens] - climMean['icemask'][ens_clim]
            MHW_ts[model]['temp_mean'][exp]['full'][:,ens] = MHW_ts[model]['temp_mean'][exp]['full'][:,ens] - climMean['full'][ens_clim]

# Generate ensemble-mean results by model
NENSNat = 0
for model in models:
    NENSNat += MHW_ts[model]['count']['histNat']['full'].shape[1]
ices = MHW_mean[model]['count'][exp].keys()
MHW_meanEM = {}
MHW_meanMM = {}
MHW_tsEM = {}
MHW_tsMM = {}
MHW_tsNat = {}
MHW_tsNatWeights = {}
MHW_diffEM = {}
MHW_diffMM = {}
MHW_agrMM = {}
years_ts = {}
years_ts['hist'] = np.arange(1850, 2005+1)
years_ts['histNat'] = np.arange(1850, 2005+1)
years_ts['rcp45'] = np.arange(2006, 2100+1)
years_ts['rcp85'] = np.arange(2006, 2100+1)
for metric in mhw_metrics:
    # Mean
    MHW_meanEM[metric] = {}
    MHW_meanMM[metric] = {}
    MHW_tsEM[metric] = {}
    MHW_tsMM[metric] = {}
    MHW_tsNat[metric] = {}
    MHW_tsNatWeights[metric] = {}
    for exp in experiments:
        MHW_meanEM[metric][exp] = {}
        MHW_meanMM[metric][exp] = {}
        MHW_tsEM[metric][exp] = {}
        MHW_tsMM[metric][exp] = {}
        MHW_tsNat[metric][exp] = {}
        MHW_tsNatWeights[metric][exp] = {}
        for ice in ices:
            MHW_meanEM[metric][exp][ice] = np.nan*np.zeros((lat.shape[0], lat.shape[1], len(models)))
            MHW_tsEM[metric][exp][ice] = np.nan*np.zeros((len(years_ts[exp]), len(models)))
            if exp == 'histNat':
                MHW_tsNat[metric][exp][ice] = np.nan*np.zeros((len(years_ts[exp]), NENSNat))
                MHW_tsNatWeights[metric][exp][ice] = np.nan*np.zeros((len(years_ts[exp]), NENSNat))
            cnt = 0
            ensCnt = 0
            for model in models:
                MHW_meanEM[metric][exp][ice][:,:,cnt] = np.mean(MHW_mean[model][metric][exp][ice], axis=2)
                tt1 = np.in1d(years_ts[exp], years[model][exp])
                tt2 = np.in1d(years[model][exp], years_ts[exp])
                MHW_tsEM[metric][exp][ice][tt1,cnt] = np.mean(MHW_ts[model][metric][exp][ice][tt2,:], axis=1)
                cnt += 1
                if exp == 'histNat':
                    for ens in range(MHW_ts[model][metric][exp][ice].shape[1]):
                        MHW_tsNat[metric][exp][ice][tt1,ensCnt] = MHW_ts[model][metric][exp][ice][tt2,ens]
                        MHW_tsNatWeights[metric][exp][ice][tt1,ensCnt] = 1./NENS[model]['histNat']
                        ensCnt += 1
            # Multi-model mean
            MHW_meanMM[metric][exp][ice] = np.nanmean(MHW_meanEM[metric][exp][ice], axis=2)
            MHW_tsMM[metric][exp][ice] = np.nanmean(MHW_tsEM[metric][exp][ice], axis=1)
    # Difference of RCP w/ hist
    MHW_diffEM[metric] = {}
    MHW_diffMM[metric] = {}
    MHW_agrMM[metric] = {}
    for exp in ['rcp45', 'rcp85']:
        MHW_diffEM[metric][exp] = {}
        MHW_diffMM[metric][exp] = {}
        MHW_agrMM[metric][rcp] = {}
        for ice in ices:
            MHW_diffEM[metric][exp][ice] = np.nan*np.zeros((lat.shape[0], lat.shape[1], len(models)))
            cnt = 0
            for model in models:
                MHW_diffEM[metric][exp][ice][:,:,cnt] = MHW_meanEM[metric][exp][ice][:,:,cnt] - MHW_meanEM[metric]['hist'][ice][:,:,cnt]
                cnt += 1
            # Multi-model mean and level of agreement
            MHW_diffMM[metric][exp][ice] = MHW_meanMM[metric][exp][ice] - MHW_meanMM[metric]['hist'][ice]
            MHW_agrMM[metric][rcp][ice] = np.zeros((lat.shape[0], lat.shape[1]))
            cnt = 0
            for model in models:
                MHW_agrMM[metric][rcp][ice] += (np.sign(MHW_diffEM[metric][exp][ice][:,:,cnt]) == np.sign(MHW_diffMM[metric][exp][ice])).astype(int)*1.
                cnt += 1

# Precalculate some histNat percentiles
exp = 'histNat'
MHW_pNat = {}
for metric in mhw_metrics:
    MHW_pNat[metric] = {}
    for ice in ices:
        MHW_pNat[metric][ice] = {}
        MHW_pNat[metric][ice]['min'] = np.min(ecj.nonans(MHW_tsNat[metric][exp][ice].flatten()))
        MHW_pNat[metric][ice]['max'] = np.max(ecj.nonans(MHW_tsNat[metric][exp][ice].flatten()))
        MHW_pNat[metric][ice]['2.5'] = ecj.weighted_quantile(ecj.nonans(MHW_tsNat[metric][exp][ice].flatten()), 0.025, sample_weight=ecj.nonans(MHW_tsNatWeights[metric][exp][ice].flatten()))
        MHW_pNat[metric][ice]['16.5'] = ecj.weighted_quantile(ecj.nonans(MHW_tsNat[metric][exp][ice].flatten()), 0.165, sample_weight=ecj.nonans(MHW_tsNatWeights[metric][exp][ice].flatten()))
        MHW_pNat[metric][ice]['83.5'] = ecj.weighted_quantile(ecj.nonans(MHW_tsNat[metric][exp][ice].flatten()), 0.835, sample_weight=ecj.nonans(MHW_tsNatWeights[metric][exp][ice].flatten()))
        MHW_pNat[metric][ice]['97.5'] = ecj.weighted_quantile(ecj.nonans(MHW_tsNat[metric][exp][ice].flatten()), 0.975, sample_weight=ecj.nonans(MHW_tsNatWeights[metric][exp][ice].flatten()))

# FAR / Probability ratios
#metric = 'total_days'
#metric = 'extreme_days'
daysMetrics = ['moderate_days', 'strong_days', 'severe_days', 'extreme_days']
cnt = 0
FAR = {}
for metric in daysMetrics:
    p_MHW_Nat_num = 0.
    p_MHW_Nat_den = 0.
    for k in range(cnt, 4):
        p_MHW_Nat_num += ecj.nonans(MHW_tsEM[daysMetrics[k]]['histNat'][ice].flatten()).sum() # Number of MHW / Category days
        p_MHW_Nat_den = len(ecj.nonans(MHW_tsEM[daysMetrics[k]]['histNat'][ice].flatten()))*365. # Total number of days
    #p_MHW_Nat = ecj.nonans(MHW_tsEM[metric]['histNat'][ice].flatten()).sum() / (len(ecj.nonans(MHW_tsEM[metric]['histNat'][ice].flatten()))*365.)
    p_MHW_Nat = p_MHW_Nat_num / p_MHW_Nat_den
    FAR[metric] = {}
    for exp in ['hist', 'rcp45', 'rcp85']:
        FAR[metric][exp] = {}
        for ice in ices:
            FAR[metric][exp][ice] = {}
            p_MHW = 0
            for k in range(cnt, 4):
                p_MHW += MHW_tsMM[daysMetrics[k]][exp][ice] / 365.
            FAR[metric][exp][ice]['mean'] = 1. - p_MHW_Nat / p_MHW
            p_MHW = np.zeros(MHW_tsEM[metric][exp][ice].shape)
            for k in range(cnt, 4):
                p_MHW += MHW_tsEM[daysMetrics[k]][exp][ice] / 365.
            p_MHW = np.min(p_MHW, axis=1)
            #p_MHW = np.min(MHW_tsEM[metric][exp][ice], axis=1) / 365.
            FAR[metric][exp][ice]['min'] = 1. - p_MHW_Nat / p_MHW
            p_MHW = np.zeros(MHW_tsEM[metric][exp][ice].shape)
            for k in range(cnt, 4):
                p_MHW += MHW_tsEM[daysMetrics[k]][exp][ice] / 365.
            p_MHW = np.max(p_MHW, axis=1)
            #p_MHW = np.max(MHW_tsEM[metric][exp][ice], axis=1) / 365.
            FAR[metric][exp][ice]['max'] = 1. - p_MHW_Nat / p_MHW
    cnt += 1

#
# Fix up some ranges
#

# Remove set all durations >durationMax days equal to durationMax days
durationMax = 1000
for metric in ['duration', 'duration_new']:
    for exp in experiments:
        for ice in ices:
            for model in models:
                MHW_mean[model][metric][exp][ice][MHW_mean[model][metric][exp][ice]>durationMax] = durationMax
            MHW_meanEM[metric][exp][ice][MHW_meanEM[metric][exp][ice]>durationMax] = durationMax
            MHW_meanMM[metric][exp][ice][MHW_meanMM[metric][exp][ice]>durationMax] = durationMax

#metric = 'intensity_max_max'
#iMax = 2.5
#for exp in experiments:
#    for ice in ices:
#        for model in models:
#            MHW_mean[model][metric][exp][ice][MHW_mean[model][metric][exp][ice]>iMax] = iMax
#        MHW_meanEM[metric][exp][ice][MHW_meanEM[metric][exp][ice]>iMax] = iMax
#        MHW_meanMM[metric][exp][ice][MHW_meanMM[metric][exp][ice]>iMax] = iMax

#
# Make some plots
#

# Plotting information
domain = [-65, 0, 70, 360]
domain_draw = [-60, 40, 60, 360]
dlat = 30
dlon = 60
bg_col = '0.6'
cont_col = '1.0'
mod_col = ['#ffb3ba', '#ffdfba', '#c3f6ce', '#bae1ff', '#c9c9ff', '#f1cbff', '#dcafa2']
#nat_col = {'max': '0.9', '95': '0.8', '67': '0.7'} # Greys
nat_col = {'max': (0.9,0.9,1.0), '95': (0.8,0.8,1.0), '67': (0.7,0.7,1.0)} # Blue-ish
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
lonproj, latproj = proj(lon, lat)
ice = 'full' #'icemask' # 'full'
hatch = '////'
Nagr = len(models)

plt.figure(figsize=(22,9))

# Code below will map ens-mean for each model, and multi-model mean
# For maps of individual models including ens. members, use the script: mhw_CMIP5_census_plots.indiv.py

for metric in ['count', 'count_new']:
    plt.clf()
    plt.subplot(2,2,1, axisbg=bg_col)
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True])
    plt.contourf(lonproj, latproj, MHW_meanMM[metric]['hist'][ice], levels=[0,1,1.5,2,2.5,3,3.5,5], cmap=plt.cm.YlOrRd)
    H = plt.colorbar()
    H.set_label('[count]')
    plt.clim(0.75,3.25)
    plt.title('Count (historical ' + str(meanYears['hist'][0]) + '-' + str(meanYears['hist'][1]) + ')')
    plt.subplot(2,2,2, axisbg=bg_col)
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True])
    plt.contourf(lonproj, latproj, MHW_meanMM[metric][rcp][ice], levels=[0,1,1.5,2,2.5,3,3.5,5], cmap=plt.cm.YlOrRd)
    H = plt.colorbar()
    H.set_label('[count]')
    plt.clim(0.75,3.25)
    plt.title('Count (RCP8.5 ' + str(meanYears[rcp][0]) + '-' + str(meanYears[rcp][1]) + ')')
    plt.subplot(2,2,3, axisbg=bg_col)
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True])
    plt.contourf(lonproj, latproj, MHW_diffMM[metric][rcp][ice], levels=np.arange(-3,3+0.5,0.5), cmap=plt.cm.RdBu_r)
    H = plt.colorbar()
    H.set_label('[count]')
    plt.clim(-2, 2)
    plt.contourf(lonproj, latproj, (MHW_agrMM[metric][rcp][ice] >= Nagr).astype(float), hatches=['', hatch], levels=[0., 0.5, 1.0], colors='none')
    plt.title('Count (Difference)')
    plt.subplot(2,2,4)
    plt.title('Global mean')
    plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['min'], MHW_pNat[metric][ice]['max'], color=nat_col['max'])
    plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['2.5'], MHW_pNat[metric][ice]['97.5'], color=nat_col['95'])
    plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['16.5'], MHW_pNat[metric][ice]['83.5'], color=nat_col['67'])
    for i in range(len(models)):
        plt.plot(years_ts['hist'], MHW_tsEM[metric]['hist'][ice][:,i], '-', color=mod_col[i])
    plt.plot(years_ts['hist'], MHW_tsMM[metric]['hist'][ice], 'k-', linewidth=2)
    for i in range(len(models)):
        plt.plot(years_ts[rcp], MHW_tsEM[metric][rcp][ice][:,i], '-', color=mod_col[i])
    plt.plot(years_ts[rcp], MHW_tsMM[metric][rcp][ice], 'r-', linewidth=2)
    plt.xlim(1900, 2060)
    #plt.ylim(0, 3.5)
    plt.ylabel('Count')
    #
    plt.savefig(figfolder + modelDir + '/MHW_' + metric + '.png', bbox_inches='tight', dpi=150)

for metric in ['duration', 'duration_new']:
    plt.clf()
    plt.subplot(2,2,1, axisbg=bg_col)
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True])
    plt.contourf(lonproj, latproj, MHW_meanMM[metric]['hist'][ice], levels=[5,15,30,60,125,250,1000], cmap=plt.cm.gist_heat_r)
    H = plt.colorbar()
    H.set_label('[days]')
    plt.clim(10, 300)
    plt.title('Duration (historical ' + str(meanYears['hist'][0]) + '-' + str(meanYears['hist'][1]) + ')')
    plt.subplot(2,2,2, axisbg=bg_col)
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True])
    plt.contourf(lonproj, latproj, MHW_meanMM[metric][rcp][ice], levels=[5,15,30,60,125,250,1000], cmap=plt.cm.gist_heat_r)
    H = plt.colorbar()
    H.set_label('[days]')
    plt.clim(10, 300)
    plt.title('Duration (RCP8.5 ' + str(meanYears[rcp][0]) + '-' + str(meanYears[rcp][1]) + ')')
    plt.subplot(2,2,3, axisbg=bg_col)
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True])
    plt.contourf(lonproj, latproj, MHW_diffMM[metric][rcp][ice], levels=[-10000,-1000,-750,-500,-250,-100,100,250,500,750,1000,10000], cmap=plt.cm.RdBu_r)
    H = plt.colorbar()
    H.set_label('[days]')
    plt.clim(-750, 750)
    plt.contourf(lonproj, latproj, (MHW_agrMM[metric][rcp][ice] >= Nagr).astype(float), hatches=['', hatch], levels=[0., 0.5, 1.0], colors='none')
    plt.title('Duration (Difference)')
    plt.subplot(2,2,4)
    plt.title('Global mean')
    plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['min'], MHW_pNat[metric][ice]['max'], color=nat_col['max'])
    plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['2.5'], MHW_pNat[metric][ice]['97.5'], color=nat_col['95'])
    plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['16.5'], MHW_pNat[metric][ice]['83.5'], color=nat_col['67'])
    for i in range(len(models)):
        plt.semilogy(years_ts['hist'], MHW_tsEM[metric]['hist'][ice][:,i], '-', color=mod_col[i])
    plt.semilogy(years_ts['hist'], MHW_tsMM[metric]['hist'][ice], 'k-', linewidth=2)
    for i in range(len(models)):
        plt.semilogy(years_ts[rcp], MHW_tsEM[metric][rcp][ice][:,i], '-', color=mod_col[i])
    plt.semilogy(years_ts[rcp], MHW_tsMM[metric][rcp][ice], 'r-', linewidth=2)
    plt.xlim(1900, 2060)
    #plt.ylim(15, 180)
    plt.ylabel('Duration [days]')
    #
    plt.savefig(figfolder + modelDir + '/MHW_' + metric + '.png', bbox_inches='tight', dpi=150)

for metric in ['intensity_max_max', 'intensity_max_max_new']:
    plt.clf()
    plt.subplot(2,2,1, axisbg=bg_col)
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True])
    plt.contourf(lonproj, latproj, MHW_meanMM[metric]['hist'][ice], levels=np.arange(0,6+0.5,0.5), cmap=plt.cm.gist_heat_r)
    H = plt.colorbar()
    H.set_label(r'[$^\circ$C]')
    plt.clim(0.5,5.75)
    plt.title('Max intensity (historical ' + str(meanYears['hist'][0]) + '-' + str(meanYears['hist'][1]) + ')')
    plt.subplot(2,2,2, axisbg=bg_col)
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True])
    plt.contourf(lonproj, latproj, MHW_meanMM[metric][rcp][ice], levels=np.arange(0,6+0.5,0.5), cmap=plt.cm.gist_heat_r)
    H = plt.colorbar()
    H.set_label(r'[$^\circ$C]')
    plt.clim(0.5,5.75)
    plt.title('Max intensity (RCP8.5 ' + str(meanYears[rcp][0]) + '-' + str(meanYears[rcp][1]) + ')')
    plt.subplot(2,2,3, axisbg=bg_col)
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True])
    plt.contourf(lonproj, latproj, MHW_diffMM[metric][rcp][ice], levels=[-10,-5,-4,-3,-2,-1,1,2,3,4,5,10], cmap=plt.cm.RdBu_r)
    H = plt.colorbar()
    H.set_label(r'[$^\circ$C]')
    plt.clim(-5, 5)
    plt.contourf(lonproj, latproj, (MHW_agrMM[metric][rcp][ice] >= Nagr).astype(float), hatches=['', hatch], levels=[0., 0.5, 1.0], colors='none')
    plt.title('Max intensity (Difference)')
    plt.subplot(2,2,4)
    plt.title('Global mean')
    plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['min'], MHW_pNat[metric][ice]['max'], color=nat_col['max'])
    plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['2.5'], MHW_pNat[metric][ice]['97.5'], color=nat_col['95'])
    plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['16.5'], MHW_pNat[metric][ice]['83.5'], color=nat_col['67'])
    for i in range(len(models)):
        plt.plot(years_ts['hist'], MHW_tsEM[metric]['hist'][ice][:,i], '-', color=mod_col[i])
    plt.plot(years_ts['hist'], MHW_tsMM[metric]['hist'][ice], 'k-', linewidth=2)
    for i in range(len(models)):
        plt.plot(years_ts[rcp], MHW_tsEM[metric][rcp][ice][:,i], '-', color=mod_col[i])
    plt.plot(years_ts[rcp], MHW_tsMM[metric][rcp][ice], 'r-', linewidth=2)
    plt.xlim(1900, 2060)
    #plt.ylim(0.5, 3.5)
    plt.ylabel(r'Intensity [$^\circ$C]')
    #
    plt.savefig(figfolder + modelDir + '/MHW_' + metric + '.png', bbox_inches='tight', dpi=150)

metric = 'total_days'
plt.clf()
plt.subplot(2,2,1, axisbg=bg_col)
proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False])
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True])
plt.contourf(lonproj, latproj, MHW_meanMM[metric]['hist'][ice], levels=[20,25,30,35,40,45,50,365], cmap=plt.cm.gist_heat_r)
H = plt.colorbar()
H.set_label('[days]')
plt.clim(25, 50)
plt.title('Total days (historical ' + str(meanYears['hist'][0]) + '-' + str(meanYears['hist'][1]) + ')')
plt.subplot(2,2,2, axisbg=bg_col)
proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False])
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True])
plt.contourf(lonproj, latproj, MHW_meanMM[metric][rcp][ice], levels=np.arange(0,400+1,40), cmap=plt.cm.gist_heat_r)
H = plt.colorbar()
H.set_label('[days]')
plt.clim(0, 450)
plt.title('Total days (RCP8.5 ' + str(meanYears[rcp][0]) + '-' + str(meanYears[rcp][1]) + ')')
plt.subplot(2,2,3, axisbg=bg_col)
proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False])
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True])
plt.contourf(lonproj, latproj, MHW_diffMM[metric][rcp][ice], levels=np.arange(-380,380+1,40), cmap=plt.cm.RdBu_r)
H = plt.colorbar()
H.set_label('[days]')
plt.clim(-365, 365)
plt.contourf(lonproj, latproj, (MHW_agrMM[metric][rcp][ice] >= Nagr).astype(float), hatches=['', hatch], levels=[0., 0.5, 1.0], colors='none')
plt.title('Total days (Difference)')
plt.subplot(2,2,4)
plt.title('Global mean')
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['min'], MHW_pNat[metric][ice]['max'], color=nat_col['max'])
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['2.5'], MHW_pNat[metric][ice]['97.5'], color=nat_col['95'])
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['16.5'], MHW_pNat[metric][ice]['83.5'], color=nat_col['67'])
for i in range(len(models)):
    plt.plot(years_ts['hist'], MHW_tsEM[metric]['hist'][ice][:,i], '-', color=mod_col[i])
plt.plot(years_ts['hist'], MHW_tsMM[metric]['hist'][ice], 'k-', linewidth=2)
for i in range(len(models)):
    plt.plot(years_ts[rcp], MHW_tsEM[metric][rcp][ice][:,i], '-', color=mod_col[i])
plt.plot(years_ts[rcp], MHW_tsMM[metric][rcp][ice], 'r-', linewidth=2)
plt.xlim(1900, 2060)
plt.ylim(5, 370)
plt.ylabel(r'Total Days [days]')
#
plt.savefig(figfolder + modelDir + '/MHW_' + metric + '.png', bbox_inches='tight', dpi=150)

for metric in ['temp_mean', 'temp_max']:
    plt.clf()
    plt.subplot(2,2,1, axisbg=bg_col)
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True])
    plt.contourf(lonproj, latproj, MHW_meanMM[metric]['hist'][ice], levels=np.arange(0,36+1,2))
    H = plt.colorbar()
    H.set_label(r'[$^\circ$C]')
    plt.clim(2, 34)
    plt.title('SST (historical ' + str(meanYears['hist'][0]) + '-' + str(meanYears['hist'][1]) + ')')
    plt.subplot(2,2,2, axisbg=bg_col)
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True])
    plt.contourf(lonproj, latproj, MHW_meanMM[metric][rcp][ice], levels=np.arange(0,36+1,2))
    H = plt.colorbar()
    H.set_label(r'[$^\circ$C]')
    plt.clim(2, 34)
    plt.title('SST (RCP8.5 ' + str(meanYears[rcp][0]) + '-' + str(meanYears[rcp][1]) + ')')
    plt.subplot(2,2,3, axisbg=bg_col)
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True])
    plt.contourf(lonproj, latproj, MHW_diffMM[metric][rcp][ice], levels=[-10,-5,-4,-3,-2,-1,1,2,3,4,5,10], cmap=plt.cm.RdBu_r)
    H = plt.colorbar()
    H.set_label(r'[$^\circ$C]')
    plt.clim(-5, 5)
    plt.contourf(lonproj, latproj, (MHW_agrMM[metric][rcp][ice] >= Nagr).astype(float), hatches=['', hatch], levels=[0., 0.5, 1.0], colors='none')
    plt.title('SST (Difference)')
    plt.subplot(2,2,4)
    plt.title('Global mean')
    plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['min'], MHW_pNat[metric][ice]['max'], color=nat_col['max'])
    plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['2.5'], MHW_pNat[metric][ice]['97.5'], color=nat_col['95'])
    plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['16.5'], MHW_pNat[metric][ice]['83.5'], color=nat_col['67'])
    for i in range(len(models)):
        plt.plot(years_ts['hist'], MHW_tsEM[metric]['hist'][ice][:,i], '-', color=mod_col[i])
    plt.plot(years_ts['hist'], MHW_tsMM[metric]['hist'][ice], 'k-', linewidth=2)
    for i in range(len(models)):
        plt.plot(years_ts[rcp], MHW_tsEM[metric][rcp][ice][:,i], '-', color=mod_col[i])
    plt.plot(years_ts[rcp], MHW_tsMM[metric][rcp][ice], 'r-', linewidth=2)
    plt.xlim(1900, 2060)
    #plt.ylim(np.nanmin(MHW_ts[metric]['hist'][ice])-0.25, np.nanmax(MHW_ts[metric][rcp][ice])+0.25)
    #plt.ylim(-1, 3.5)
    plt.ylabel(r'SST [$^\circ$C]')
    #
    plt.savefig(figfolder + modelDir + '/MHW_' + metric + '.png', bbox_inches='tight', dpi=150)

# Main big figure for paper

#year_15 = {}
#year_20 = {}
#for exp in ['rcp45', 'rcp85']:
#    tt = np.where(np.mean(MHW_ts['temp_mean'][exp][ice], axis=1) >= 1.5)[0]
#    if len(tt) > 0:
#        year_15[exp] = years[exp][tt[0]]
#    else:
#        year_15[exp] = 9999
#    tt = np.where(np.mean(MHW_ts['temp_mean'][exp][ice], axis=1) >= 2.0)[0]
#    if len(tt) > 0:
#        year_20[exp] = years[exp][tt[0]]
#    else:
#        year_20[exp] = 9999

plt.figure(figsize=(19,11))
plt.clf()

metric = 'count_new'
plt.subplot(4,3,1, axisbg=bg_col)
proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], dashes=[3,900])
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True], dashes=[3,900])
plt.contourf(lonproj, latproj, MHW_meanMM[metric]['hist'][ice], levels=[0,0.5,1,1.5,2,2.5,3], cmap=plt.cm.gist_heat_r)
H = plt.colorbar()
H.set_label('[count]')
plt.clim(0.75,5)
plt.title('Count (historical ' + str(meanYears['hist'][0]) + '-' + str(meanYears['hist'][1]) + ')')
plt.subplot(4,3,2, axisbg=bg_col)
proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], dashes=[3,900])
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True], dashes=[3,900])
plt.contourf(lonproj, latproj, MHW_diffMM[metric][rcp][ice], levels=np.arange(-4,4+0.5,1), cmap=plt.cm.RdBu_r)
H = plt.colorbar()
H.set_label('[count]')
plt.clim(-3.5, 3.5)
plt.contourf(lonproj, latproj, (MHW_agrMM[metric][rcp][ice] >= Nagr).astype(float), hatches=['', hatch], levels=[0., 0.5, 1.0], colors='none')
plt.title('(historical ' + str(meanYears['hist'][0]) + '-' + str(meanYears['hist'][1]) + ') - (RCP8.5 ' + str(meanYears['rcp85'][0]) + '-' + str(meanYears['rcp85'][1]) + ')')
plt.subplot(4,3,3)
plt.title('Global mean')
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['min'], MHW_pNat[metric][ice]['max'], color=nat_col['max'])
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['2.5'], MHW_pNat[metric][ice]['97.5'], color=nat_col['95'])
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['16.5'], MHW_pNat[metric][ice]['83.5'], color=nat_col['67'])
plt.fill_between(years_ts['hist'], np.nanmin(MHW_tsEM[metric]['hist'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['hist'][ice], axis=1), color='0.5', alpha=0.8)
plt.plot(years_ts['hist'], MHW_tsMM[metric]['hist'][ice], 'k-', linewidth=2)
plt.fill_between(years_ts['rcp45'], np.nanmin(MHW_tsEM[metric]['rcp45'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['rcp45'][ice], axis=1), color=(0.6,0.75,0.6), alpha=0.8)
plt.plot(years_ts[rcp], MHW_tsMM[metric]['rcp45'][ice], '-', linewidth=2, color=(0,0.6,0))
plt.fill_between(years_ts['rcp85'], np.nanmin(MHW_tsEM[metric]['rcp85'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['rcp85'][ice], axis=1), color=(0.9,0.5,0.5), alpha=0.8)
plt.plot(years_ts[rcp], MHW_tsMM[metric]['rcp85'][ice], '-', linewidth=2, color=(0.9,0,0))
plt.xlim(1850, 2100)
#plt.ylim(0, 3.5)
plt.ylabel('Count')

metric = 'intensity_max_max_new'
plt.subplot(4,3,4, axisbg=bg_col)
proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], dashes=[3,900])
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True], dashes=[3,900])
plt.contourf(lonproj, latproj, MHW_meanMM[metric]['hist'][ice], levels=np.append(np.arange(0,2.5+0.25,0.25),5), cmap=plt.cm.gist_heat_r)
H = plt.colorbar()
H.set_label(r'[$^\circ$C]')
plt.clim(0.5,2.75)
plt.title('Max intensity')
plt.subplot(4,3,5, axisbg=bg_col)
proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], dashes=[3,900])
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True], dashes=[3,900])
plt.contourf(lonproj, latproj, MHW_diffMM[metric][rcp][ice], levels=[-10,-5,-4,-3,-2,-1,1,2,3,4,5,10], cmap=plt.cm.RdBu_r)
H = plt.colorbar()
H.set_label(r'[$^\circ$C]')
plt.clim(-5, 5)
plt.contourf(lonproj, latproj, (MHW_agrMM[metric][rcp][ice] >= Nagr).astype(float), hatches=['', hatch], levels=[0., 0.5, 1.0], colors='none')
plt.subplot(4,3,6)
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['min'], MHW_pNat[metric][ice]['max'], color=nat_col['max'])
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['2.5'], MHW_pNat[metric][ice]['97.5'], color=nat_col['95'])
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['16.5'], MHW_pNat[metric][ice]['83.5'], color=nat_col['67'])
plt.fill_between(years_ts['hist'], np.nanmin(MHW_tsEM[metric]['hist'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['hist'][ice], axis=1), color='0.5', alpha=0.8)
plt.plot(years_ts['hist'], MHW_tsMM[metric]['hist'][ice], 'k-', linewidth=2)
plt.fill_between(years_ts['rcp45'], np.nanmin(MHW_tsEM[metric]['rcp45'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['rcp45'][ice], axis=1), color=(0.6,0.75,0.6), alpha=0.8)
plt.plot(years_ts[rcp], MHW_tsMM[metric]['rcp45'][ice], '-', linewidth=2, color=(0,0.6,0))
plt.fill_between(years_ts['rcp85'], np.nanmin(MHW_tsEM[metric]['rcp85'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['rcp85'][ice], axis=1), color=(0.9,0.5,0.5), alpha=0.8)
plt.plot(years_ts[rcp], MHW_tsMM[metric]['rcp85'][ice], '-', linewidth=2, color=(0.9,0,0))
plt.xlim(1850, 2100)
#plt.ylim(0.25, 5.25)
plt.ylabel(r'Intensity [$^\circ$C]')

metric = 'duration_new'
plt.subplot(4,3,7, axisbg=bg_col)
proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], dashes=[3,900])
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True], dashes=[3,900])
plt.contourf(lonproj, latproj, MHW_meanMM[metric]['hist'][ice], levels=[5,10,15,30,45,60,120,240], cmap=plt.cm.gist_heat_r)
H = plt.colorbar()
H.set_label('[days]')
#plt.clim(10, 300)
plt.title('Duration')
plt.subplot(4,3,8, axisbg=bg_col)
proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], dashes=[3,900])
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True], dashes=[3,900])
plt.contourf(lonproj, latproj, MHW_diffMM[metric][rcp][ice], levels=[-35000,-1000,-750,-500,-250,-100,100,250,500,750,1000,35000], cmap=plt.cm.RdBu_r)
H = plt.colorbar()
H.set_label('[days]')
plt.clim(-1000, 1000)
plt.contourf(lonproj, latproj, (MHW_agrMM[metric][rcp][ice] >= Nagr).astype(float), hatches=['', hatch], levels=[0., 0.5, 1.0], colors='none')
plt.subplot(4,3,9)
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['min'], MHW_pNat[metric][ice]['max'], color=nat_col['max'])
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['2.5'], MHW_pNat[metric][ice]['97.5'], color=nat_col['95'])
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['16.5'], MHW_pNat[metric][ice]['83.5'], color=nat_col['67'])
plt.fill_between(years_ts['hist'], np.nanmin(MHW_tsEM[metric]['hist'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['hist'][ice], axis=1), color='0.5', alpha=0.8)
plt.semilogy(years_ts['hist'], MHW_tsMM[metric]['hist'][ice], 'k-', linewidth=2)
plt.fill_between(years_ts['rcp45'], np.nanmin(MHW_tsEM[metric]['rcp45'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['rcp45'][ice], axis=1), color=(0.6,0.75,0.6), alpha=0.8)
plt.semilogy(years_ts[rcp], MHW_tsMM[metric]['rcp45'][ice], '-', linewidth=2, color=(0,0.6,0))
plt.fill_between(years_ts['rcp85'], np.nanmin(MHW_tsEM[metric]['rcp85'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['rcp85'][ice], axis=1), color=(0.9,0.5,0.5), alpha=0.8)
plt.semilogy(years_ts[rcp], MHW_tsMM[metric]['rcp85'][ice], '-', linewidth=2, color=(0.9,0,0))
plt.xlim(1850, 2100)
#plt.ylim(15, 180)
plt.ylabel('Duration [days]')

metric = 'total_days'
plt.subplot(4,3,10, axisbg=bg_col)
proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], dashes=[3,900])
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True], dashes=[3,900])
plt.contourf(lonproj, latproj, MHW_meanMM[metric]['hist'][ice], levels=np.arange(25,40+2.5,2.5), cmap=plt.cm.gist_heat_r)
H = plt.colorbar()
H.set_label('[days]')
plt.clim(27, 49)
plt.title('Total days')
plt.subplot(4,3,11, axisbg=bg_col)
proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], dashes=[3,900])
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True], dashes=[3,900])
plt.contourf(lonproj, latproj, MHW_diffMM[metric][rcp][ice], levels=np.arange(-400,400+1,80), cmap=plt.cm.RdBu_r)
H = plt.colorbar()
H.set_label('[days]')
#plt.clim(-365, 365)
plt.contourf(lonproj, latproj, (MHW_agrMM[metric][rcp][ice] >= Nagr).astype(float), hatches=['', hatch], levels=[0., 0.5, 1.0], colors='none')
plt.subplot(4,3,12)
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['min'], MHW_pNat[metric][ice]['max'], color=nat_col['max'])
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['2.5'], MHW_pNat[metric][ice]['97.5'], color=nat_col['95'])
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['16.5'], MHW_pNat[metric][ice]['83.5'], color=nat_col['67'])
plt.plot([1850, 2100], [365, 365], '--', color='0.4')
plt.fill_between(years_ts['hist'], np.nanmin(MHW_tsEM[metric]['hist'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['hist'][ice], axis=1), color='0.5', alpha=0.8)
plt.plot(years_ts['hist'], MHW_tsMM[metric]['hist'][ice], 'k-', linewidth=2)
plt.fill_between(years_ts['rcp45'], np.nanmin(MHW_tsEM[metric]['rcp45'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['rcp45'][ice], axis=1), color=(0.6,0.75,0.6), alpha=0.8)
plt.plot(years_ts[rcp], MHW_tsMM[metric]['rcp45'][ice], '-', linewidth=2, color=(0,0.6,0))
plt.fill_between(years_ts['rcp85'], np.nanmin(MHW_tsEM[metric]['rcp85'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['rcp85'][ice], axis=1), color=(0.9,0.5,0.5), alpha=0.8)
plt.plot(years_ts[rcp], MHW_tsMM[metric]['rcp85'][ice], '-', linewidth=2, color=(0.9,0,0))
plt.xlim(1850, 2100)
plt.ylim(0, 400)
plt.ylabel(r'Total Days [days]')

# plt.savefig(figfolder + modelDir + '/MHW_meanDiffTimeseries.png', bbox_inches='tight', dpi=150)

plt.figure(figsize=(12,8))
plt.clf()

metric = 'intensity_max_max_new'
plt.subplot(3,2,1, axisbg=bg_col)
proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], dashes=[3,900])
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,False], dashes=[3,900])
plt.contourf(lonproj, latproj, MHW_meanMM[metric]['hist'][ice], levels=np.arange(0, 5+0.5, 0.5), cmap=plt.cm.gist_heat_r)
H = plt.colorbar()
H.set_label(r'[$^\circ$C]')
plt.clim(0.5, 5)
plt.title('(A) Maximum intensity (historical ' + str(meanYears['hist'][0]) + '-' + str(meanYears['hist'][1]) + ')')
plt.subplot(3,2,3, axisbg=bg_col)
proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], dashes=[3,900])
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True], dashes=[3,900])
plt.contourf(lonproj, latproj, MHW_diffMM[metric][rcp][ice], levels=[-4,-2,-1.5,-1,-0.5,0.5,1,1.5,2,4], cmap=plt.cm.RdBu_r)
H = plt.colorbar()
H.set_label(r'[$^\circ$C]')
plt.clim(-3.5, 3.5)
plt.contourf(lonproj, latproj, (MHW_agrMM[metric][rcp][ice] >= Nagr).astype(float), hatches=['', hatch], levels=[0., 0.5, 1.0], colors='none')
plt.title('(B) (RCP8.5 ' + str(meanYears[rcp][0]) + '-' + str(meanYears[rcp][1]) + ') - (historical ' + str(meanYears['hist'][0]) + '-' + str(meanYears['hist'][1]) + ')')
#plt.subplot(3,2,5)
plt.subplot2grid((3,10), (2,0), rowspan=1, colspan=4)
plt.plot(years_ts['hist'], MHW_tsMM[metric]['hist'][ice], 'k-', linewidth=2, zorder=30)
plt.plot(years_ts[rcp], MHW_tsMM[metric]['rcp45'][ice], '-', linewidth=2, color=(0.6,0.3,0), zorder=30)
plt.plot(years_ts[rcp], MHW_tsMM[metric]['rcp85'][ice], '-', linewidth=2, color=(0.9,0,0), zorder=30)
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['min'], MHW_pNat[metric][ice]['max'], color=nat_col['max'], zorder=10)
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['2.5'], MHW_pNat[metric][ice]['97.5'], color=nat_col['95'], zorder=10)
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['16.5'], MHW_pNat[metric][ice]['83.5'], color=nat_col['67'], zorder=10)
plt.fill_between(years_ts['hist'], np.nanmin(MHW_tsEM[metric]['hist'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['hist'][ice], axis=1), color='0.5', alpha=0.8, zorder=20)
plt.fill_between(years_ts['rcp45'], np.nanmin(MHW_tsEM[metric]['rcp45'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['rcp45'][ice], axis=1), color=(1.0,0.7,0.3), alpha=0.8, zorder=20)
plt.fill_between(years_ts['rcp85'], np.nanmin(MHW_tsEM[metric]['rcp85'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['rcp85'][ice], axis=1), color=(0.9,0.5,0.5), alpha=0.8, zorder=20)
plt.xlim(1850, 2100)
#plt.ylim(0.25, 5.25)
plt.ylabel(r'Intensity [$^\circ$C]')
plt.title('(C) Global mean time series')
plt.legend(['Historical', 'RCP4.5', 'RCP8.5', 'Nat. (min/max)', 'Nat. (95% CI)', 'Nat. (66% CI)', 'Hist. (min/max)', 'RCP4.5 (min/max)', 'RCP8.5 (min/max)'], loc='upper left', ncol=2, fontsize=8)

metric = 'total_days'
plt.subplot(3,2,2, axisbg=bg_col)
proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], dashes=[3,900])
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,False], dashes=[3,900])
plt.contourf(lonproj, latproj, MHW_meanMM[metric]['hist'][ice], levels=np.arange(23, 37+2, 2), cmap=plt.cm.gist_heat_r)
H = plt.colorbar()
H.set_label('[days]')
plt.clim(25, 45)
plt.title('(D) Total Annual MHW Days (historical ' + str(meanYears['hist'][0]) + '-' + str(meanYears['hist'][1]) + ')')
plt.subplot(3,2,4, axisbg=bg_col)
proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], dashes=[3,900])
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True], dashes=[3,900])
plt.contourf(lonproj, latproj, MHW_diffMM[metric][rcp][ice], levels=np.arange(-400,400+1,80), cmap=plt.cm.RdBu_r)
H = plt.colorbar()
H.set_label('[days]')
plt.clim(-365, 365)
plt.contourf(lonproj, latproj, (MHW_agrMM[metric][rcp][ice] >= Nagr).astype(float), hatches=['', hatch], levels=[0., 0.5, 1.0], colors='none')
plt.title('(E) (RCP8.5 ' + str(meanYears[rcp][0]) + '-' + str(meanYears[rcp][1]) + ') - (historical ' + str(meanYears['hist'][0]) + '-' + str(meanYears['hist'][1]) + ')')
#plt.subplot(3,2,6)
plt.subplot2grid((3,20), (2,11), rowspan=1, colspan=8)
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['min'], MHW_pNat[metric][ice]['max'], color=nat_col['max'])
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['2.5'], MHW_pNat[metric][ice]['97.5'], color=nat_col['95'])
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['16.5'], MHW_pNat[metric][ice]['83.5'], color=nat_col['67'])
plt.plot([1850, 2100], [365, 365], '--', color='0.4')
plt.fill_between(years_ts['hist'], np.nanmin(MHW_tsEM[metric]['hist'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['hist'][ice], axis=1), color='0.5', alpha=0.8)
plt.plot(years_ts['hist'], MHW_tsMM[metric]['hist'][ice], 'k-', linewidth=2)
plt.fill_between(years_ts['rcp45'], np.nanmin(MHW_tsEM[metric]['rcp45'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['rcp45'][ice], axis=1), color=(1,0.7,0.3), alpha=0.8)
plt.plot(years_ts[rcp], MHW_tsMM[metric]['rcp45'][ice], '-', linewidth=2, color=(0.6,0.3,0))
plt.fill_between(years_ts['rcp85'], np.nanmin(MHW_tsEM[metric]['rcp85'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['rcp85'][ice], axis=1), color=(0.9,0.5,0.5), alpha=0.8)
plt.plot(years_ts[rcp], MHW_tsMM[metric]['rcp85'][ice], '-', linewidth=2, color=(0.9,0,0))
plt.xlim(1850, 2100)
plt.ylim(0, 400)
plt.ylabel(r'Total Annual MHW Days [days]')
plt.title('(F) Global mean time series')

# plt.savefig(figfolder + modelDir + '/MHW_meanDiffTimeseries_Vert.png', bbox_inches='tight', dpi=300)

# Annual Max temp
plt.figure(figsize=(12,8))
plt.clf()

metric = 'temp_max'
plt.subplot(3,2,1, axisbg=bg_col)
proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], dashes=[3,900])
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,False], dashes=[3,900])
plt.contourf(lonproj, latproj, MHW_meanMM[metric]['hist'][ice], levels=np.arange(0, 5+0.5, 0.5), cmap=plt.cm.gist_heat_r)
H = plt.colorbar()
H.set_label(r'[$^\circ$C]')
plt.clim(0.5, 5)
plt.title('Maximum intensity (historical ' + str(meanYears['hist'][0]) + '-' + str(meanYears['hist'][1]) + ')')
plt.subplot(3,2,3, axisbg=bg_col)
proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], dashes=[3,900])
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True], dashes=[3,900])
plt.contourf(lonproj, latproj, MHW_diffMM[metric][rcp][ice], levels=[-4,-2,-1.5,-1,-0.5,0.5,1,1.5,2,4], cmap=plt.cm.RdBu_r)
H = plt.colorbar()
H.set_label(r'[$^\circ$C]')
plt.clim(-3.5, 3.5)
plt.contourf(lonproj, latproj, (MHW_agrMM[metric][rcp][ice] >= Nagr).astype(float), hatches=['', hatch], levels=[0., 0.5, 1.0], colors='none')
plt.title('(RCP8.5 ' + str(meanYears[rcp][0]) + '-' + str(meanYears[rcp][1]) + ') - (historical ' + str(meanYears['hist'][0]) + '-' + str(meanYears['hist'][1]) + ')')
#plt.subplot(3,2,5)
plt.subplot2grid((3,10), (2,0), rowspan=1, colspan=4)
plt.plot(years_ts['hist'], MHW_tsMM[metric]['hist'][ice], 'k-', linewidth=2, zorder=30)
plt.plot(years_ts[rcp], MHW_tsMM[metric]['rcp45'][ice], '-', linewidth=2, color=(0,0.6,0), zorder=30)
plt.plot(years_ts[rcp], MHW_tsMM[metric]['rcp85'][ice], '-', linewidth=2, color=(0.9,0,0), zorder=30)
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['min'], MHW_pNat[metric][ice]['max'], color=nat_col['max'], zorder=10)
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['2.5'], MHW_pNat[metric][ice]['97.5'], color=nat_col['95'], zorder=10)
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['16.5'], MHW_pNat[metric][ice]['83.5'], color=nat_col['67'], zorder=10)
plt.fill_between(years_ts['hist'], np.nanmin(MHW_tsEM[metric]['hist'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['hist'][ice], axis=1), color='0.5', alpha=0.8, zorder=20)
plt.fill_between(years_ts['rcp45'], np.nanmin(MHW_tsEM[metric]['rcp45'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['rcp45'][ice], axis=1), color=(0.6,0.75,0.6), alpha=0.8, zorder=20)
plt.fill_between(years_ts['rcp85'], np.nanmin(MHW_tsEM[metric]['rcp85'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['rcp85'][ice], axis=1), color=(0.9,0.5,0.5), alpha=0.8, zorder=20)
plt.xlim(1850, 2100)
#plt.ylim(0.25, 5.25)
plt.ylabel(r'Intensity [$^\circ$C]')
plt.title('Global mean time series')
plt.legend(['Historical', 'RCP4.5', 'RCP8.5', 'Nat. (min/max)', 'Nat. (95% CI)', 'Nat. (66% CI)', 'Hist. (min/max)', 'RCP4.5 (min/max)', 'RCP8.5 (min/max)'], loc='upper left', ncol=2, fontsize=8)

# Emergence years
for metric in ['intensity_max_max_new', 'total_days']:
    for exp in ['rcp45', 'rcp85']:
        tmp = np.where(np.nanmax(MHW_tsEM[metric][exp][ice], axis=1) > MHW_pNat[metric][ice]['max'])[0]
        tt1 = tmp[0]
        tmp = np.where(np.nanmin(MHW_tsEM[metric][exp][ice], axis=1) > MHW_pNat[metric][ice]['max'])[0]
        tt2 = tmp[0]
        print metric, exp, years_ts[exp][tt1], years_ts[exp][tt2]

metric = 'total_days'
exp = 'hist'
tmp = np.where(np.nanmax(MHW_tsEM[metric][exp][ice], axis=1) > MHW_pNat[metric][ice]['max'])[0]
tt1 = tmp[0]
#tmp = np.where(np.nanmin(MHW_tsEM[metric][exp][ice], axis=1) > MHW_pNat[metric][ice]['max'])[0]
#tt2 = tmp[0]
print metric, exp, years_ts[exp][tt1], years_ts[exp][tt2]


# Latitudinal structure

plt.clf()
plt.subplot(1,3,1)
metric = 'intensity_max_max_new'
#plt.fill_between([MHW_ts[metric]['histNat'][ice].min(), MHW_ts[metric]['histNat'][ice].max()], -60, 60, color=nat_col['max'])
#plt.fill_between(np.percentile(MHW_ts[metric]['histNat'][ice], [2.5, 97.5]), -60, 60, color=nat_col['95'])
#plt.fill_between(np.percentile(MHW_ts[metric]['histNat'][ice], [16.5, 83.5]), -60, 60, color=nat_col['67'])
plt.plot(np.nanmean(np.mean(MHW_mean[metric]['hist'][ice], axis=2), axis=1), lat[:,0], 'k-', linewidth=2)
plt.plot(np.nanmean(np.mean(MHW_mean[metric]['rcp45'][ice], axis=2), axis=1), lat[:,0], '-', linewidth=2, color=(0,0.6,0))
plt.plot(np.nanmean(np.mean(MHW_mean[metric]['rcp85'][ice], axis=2), axis=1), lat[:,0], 'r-', linewidth=2)
plt.plot(np.nanmean(np.mean(MHW_mean[metric]['rcp45'][ice], axis=2), axis=1) - np.nanmean(np.mean(MHW_mean[metric]['hist'][ice], axis=2), axis=1), lat[:,0], '--', linewidth=2, color=(0,0.6,0))
plt.plot(np.nanmean(np.mean(MHW_mean[metric]['rcp85'][ice], axis=2), axis=1) - np.nanmean(np.mean(MHW_mean[metric]['hist'][ice], axis=2), axis=1), lat[:,0], 'r--', linewidth=2)
plt.ylim(-60, 60)
plt.ylabel('Latitude')
plt.xlabel('Maximum intensity [$^\circ$C]')
plt.legend(['historical (' + str(meanYears['hist'][0]) + '-' + str(meanYears['hist'][1]) + ')', 'RCP45 (' + str(meanYears[rcp][0]) + '-' + str(meanYears[rcp][1]) + ')', 'RCP85 (' + str(meanYears[rcp][0]) + '-' + str(meanYears[rcp][1]) + ')', 'Difference (RCP45)', 'Difference (RCP85)'], loc='lower right')
plt.subplot(1,3,2)
metric = 'duration'
plt.semilogx(np.nanmean(np.mean(MHW_mean[metric]['hist'][ice], axis=2), axis=1), lat[:,0], 'k-', linewidth=2)
plt.semilogx(np.nanmean(np.mean(MHW_mean[metric]['rcp45'][ice], axis=2), axis=1), lat[:,0], '-', linewidth=2, color=(0,0.6,0))
plt.semilogx(np.nanmean(np.mean(MHW_mean[metric]['rcp85'][ice], axis=2), axis=1), lat[:,0], 'r-', linewidth=2)
plt.semilogx(np.nanmean(np.mean(MHW_mean[metric]['rcp45'][ice], axis=2), axis=1) - np.nanmean(np.mean(MHW_mean[metric]['hist'][ice], axis=2), axis=1), lat[:,0], '--', linewidth=2, color=(0,0.6,0))
plt.semilogx(np.nanmean(np.mean(MHW_mean[metric]['rcp85'][ice], axis=2), axis=1) - np.nanmean(np.mean(MHW_mean[metric]['hist'][ice], axis=2), axis=1), lat[:,0], 'r--', linewidth=2)
plt.ylim(-60, 60)
plt.xlabel('Duration [days]')
plt.subplot(1,3,3)
metric = 'total_days'
plt.plot(np.nanmean(np.mean(MHW_mean[metric]['hist'][ice], axis=2), axis=1), lat[:,0], 'k-', linewidth=2)
plt.plot(np.nanmean(np.mean(MHW_mean[metric]['rcp45'][ice], axis=2), axis=1), lat[:,0], '-', linewidth=2, color=(0,0.6,0))
plt.plot(np.nanmean(np.mean(MHW_mean[metric]['rcp85'][ice], axis=2), axis=1), lat[:,0], 'r-', linewidth=2)
plt.plot(np.nanmean(np.mean(MHW_mean[metric]['rcp45'][ice], axis=2), axis=1) - np.nanmean(np.mean(MHW_mean[metric]['hist'][ice], axis=2), axis=1), lat[:,0], '--', linewidth=2, color=(0,0.6,0))
plt.plot(np.nanmean(np.mean(MHW_mean[metric]['rcp85'][ice], axis=2), axis=1) - np.nanmean(np.mean(MHW_mean[metric]['hist'][ice], axis=2), axis=1), lat[:,0], 'r--', linewidth=2)
plt.ylim(-60, 60)
plt.xlabel('Annual MHW days [days]')
#
# plt.savefig(figfolder + modelDir + '/MHW_zonalMean.png', bbox_inches='tight', dpi=150)

# Year of Permanent MHW state
metric = 'total_days_permMHW'
#cmap1 = plt.get_cmap('viridis_r')
#cmap1 = mpl.colors.ListedColormap(np.append(cmap1(np.floor(np.linspace(10, len(cmap1.colors)-1, 5)).astype(int)), np.array([[1,1,1,1]]), axis=0))
plt.figure()
plt.clf()
#plt.subplot(1,3,1, axisbg=bg_col)
PMHW45 = MHW_meanMM[metric]['rcp45'][ice].copy()
never = np.isnan(PMHW45) * ~np.isnan(MHW_meanMM['count']['rcp45'][ice])
PMHW45[never] = 2105
#proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
#proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], dashes=[3,900])
#proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True], dashes=[3,900])
##plt.contourf(lonproj, latproj, MHW_meanMM[metric]['rcp45'][ice], levels=np.arange(2000, 2120+1, 20), cmap=cmap1) #plt.cm.viridis_r)
#plt.contourf(lonproj, latproj, PMHW45, levels=np.arange(2000, 2120+1, 20), cmap=plt.cm.viridis_r)
#H = plt.colorbar()
#plt.clim(2000, 2100)
#H.set_label(r'Year')
#plt.title('(A) Date of first permanent MHW (RCP4.5)')
plt.subplot(1,2,1, axisbg=bg_col)
PMHW85 = MHW_meanMM[metric]['rcp85'][ice].copy()
never = np.isnan(PMHW85) * ~np.isnan(MHW_meanMM['count']['rcp85'][ice])
PMHW85[never] = 2105
proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], dashes=[3,900])
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True], dashes=[3,900])
#plt.contourf(lonproj, latproj, MHW_meanMM[metric]['rcp85'][ice], levels=np.arange(2000, 2120+1, 20), cmap=plt.cm.viridis_r)
plt.contourf(lonproj, latproj, PMHW85, levels=np.arange(2000, 2120+1, 20), cmap=plt.cm.viridis_r)
H = plt.colorbar()
plt.clim(2000, 2100)
#H.set_label(r'Year')
H.set_ticks([2000, 2020, 2040, 2060, 2080, 2100, 2120])
H.set_ticklabels([2000, 2020, 2040, 2060, 2080, 2100, '>2100'])
plt.title('(A) Date of first permanent MHW (RCP8.5)')

plt.subplot2grid((1,20), (0,10), colspan=7, rowspan=1)
#plt.plot(years_ts['hist'], MHW_tsMM[metric]['hist'][ice], 'k-', linewidth=2)
plt.plot(years_ts[rcp], MHW_tsMM[metric]['rcp45'][ice], '-', linewidth=2, color=(0.6,0.3,0), zorder=30)
plt.plot(years_ts[rcp], MHW_tsMM[metric]['rcp85'][ice], '-', linewidth=2, color=(0.9,0,0), zorder=30)
plt.fill_between(years_ts['rcp45'], np.nanmin(MHW_tsEM[metric]['rcp45'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['rcp45'][ice], axis=1), color=(1,0.7,0.3), alpha=0.8, zorder=20)
plt.fill_between(years_ts['rcp85'], np.nanmin(MHW_tsEM[metric]['rcp85'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['rcp85'][ice], axis=1), color=(0.9,0.5,0.5), alpha=0.8, zorder=20)
plt.legend(['RCP4.5', 'RCP8.5', 'RCP4.5 (min/max)', 'RCP8.5 (min/max)'])
plt.xlim(2005, 2100)
plt.ylim(0, 1)
plt.ylabel(r'Proportion of global ocean')
plt.title('(B) Proportion of global ocean with permanent MHW')

# plt.savefig(figfolder + modelDir + '/MHW_yearPermMHW.png', bbox_inches='tight', dpi=300)

# Partition by MHW category
mhw_colours = [[1., 200./255, 102./255], [1., 105./255, 0.], [158./255, 0., 0.], [45./255, 0., 0.]]
plt.figure()
plt.clf()
plt.subplot(1,2,1)
for exp in ['hist', 'rcp45']:
    plt.fill_between(years_ts[exp], 0,MHW_tsMM['moderate_days'][exp][ice] + MHW_tsMM['strong_days'][exp][ice] + MHW_tsMM['severe_days'][exp][ice] + MHW_tsMM['extreme_days'][exp][ice], color=mhw_colours[3])
    plt.fill_between(years_ts[exp], 0, MHW_tsMM['moderate_days'][exp][ice] + MHW_tsMM['strong_days'][exp][ice] + MHW_tsMM['severe_days'][exp][ice], color=mhw_colours[2])
    plt.fill_between(years_ts[exp], 0, MHW_tsMM['moderate_days'][exp][ice] + MHW_tsMM['strong_days'][exp][ice], color=mhw_colours[1])
    plt.fill_between(years_ts[exp], 0, MHW_tsMM['moderate_days'][exp][ice], color=mhw_colours[0])
plt.legend(['Extreme', 'Severe', 'Strong', 'Moderate'], loc='center left')
plt.plot([1850, 2100], [365, 365], '--', color='0.4')
plt.xlim(1850, 2100)
plt.ylim(0, 400)
plt.ylabel('Global mean days per year')
plt.title('(A) RCP4.5')
plt.subplot(1,2,2)
for exp in ['hist', 'rcp85']:
    plt.fill_between(years_ts[exp], 0, MHW_tsMM['moderate_days'][exp][ice] + MHW_tsMM['strong_days'][exp][ice] + MHW_tsMM['severe_days'][exp][ice] + MHW_tsMM['extreme_days'][exp][ice], color=mhw_colours[3])
    plt.fill_between(years_ts[exp], 0, MHW_tsMM['moderate_days'][exp][ice] + MHW_tsMM['strong_days'][exp][ice] + MHW_tsMM['severe_days'][exp][ice], color=mhw_colours[2])
    plt.fill_between(years_ts[exp], 0, MHW_tsMM['moderate_days'][exp][ice] + MHW_tsMM['strong_days'][exp][ice], color=mhw_colours[1])
    plt.fill_between(years_ts[exp], 0, MHW_tsMM['moderate_days'][exp][ice], color=mhw_colours[0])
plt.plot([1850, 2100], [365, 365], '--', color='0.4')
plt.xlim(1850, 2100)
plt.ylim(0, 400)
plt.ylabel('Global mean days per year')
plt.title('(B) RCP8.5')

print MHW_tsMM['moderate_days']['rcp45'][ice][-2]/MHW_tsMM['total_days']['rcp45'][ice][-2]
print MHW_tsMM['strong_days']['rcp45'][ice][-2]/MHW_tsMM['total_days']['rcp45'][ice][-2]
print MHW_tsMM['severe_days']['rcp45'][ice][-2]/MHW_tsMM['total_days']['rcp45'][ice][-2]
print MHW_tsMM['extreme_days']['rcp45'][ice][-2]/MHW_tsMM['total_days']['rcp45'][ice][-2]

print MHW_tsMM['moderate_days']['rcp85'][ice][-2]/MHW_tsMM['total_days']['rcp85'][ice][-2]
print MHW_tsMM['extreme_days']['rcp85'][ice][-2]/MHW_tsMM['total_days']['rcp85'][ice][-2]

# plt.savefig(figfolder + modelDir + '/MHW_CategoryDays.png', bbox_inches='tight', dpi=300)

plt.figure()
plt.clf()
cnt = 0
titles = ['(A) $\geq$Moderate', '(B) $\geq$Strong', '(C) $\geq$Severe', '(D) Extreme']
for metric in ['moderate_days', 'strong_days', 'severe_days', 'extreme_days']:
    cnt += 1
    AX = plt.subplot(2,2,cnt)
    plt.plot(years_ts['hist'], 1./(1 - FAR[metric]['hist'][ice]['mean']), 'k-', linewidth=2, zorder=30)
    plt.plot(years_ts['rcp45'], 1./(1 - FAR[metric]['rcp45'][ice]['mean']), '-', linewidth=2, color=(0.6,0.3,0), zorder=30)
    plt.plot(years_ts['rcp85'], 1./(1 - FAR[metric]['rcp85'][ice]['mean']), '-', linewidth=2, color=(0.9,0,0), zorder=30)
    plt.fill_between(years_ts['hist'], 1./(1 - FAR[metric]['hist'][ice]['min']), 1./(1 - FAR[metric]['hist'][ice]['max']), color='0.5', alpha=0.8, zorder=20)
    plt.fill_between(years_ts['rcp45'], 1./(1 - FAR[metric]['rcp45'][ice]['min']), 1./(1 - FAR[metric]['rcp45'][ice]['max']), color=(1,0.7,0.3), alpha=0.8, zorder=20)
    plt.fill_between(years_ts['rcp85'], 1./(1 - FAR[metric]['rcp85'][ice]['min']), 1./(1 - FAR[metric]['rcp85'][ice]['max']), color=(0.9,0.5,0.5), alpha=0.8, zorder=20)
    print metric, 'hist', years_ts['hist'][-1], 1./(1 - FAR[metric]['hist'][ice]['mean'][-1]), 1./(1 - FAR[metric]['hist'][ice]['min'][-1]), 1./(1 - FAR[metric]['hist'][ice]['max'][-1])
    print metric, 'rcp45', years_ts['rcp45'][-2], 1./(1 - FAR[metric]['rcp45'][ice]['mean'][-2]), 1./(1 - FAR[metric]['rcp45'][ice]['min'][-2]), 1./(1 - FAR[metric]['rcp45'][ice]['max'][-2])
    print metric, 'rcp85', years_ts['rcp85'][-2], 1./(1 - FAR[metric]['rcp85'][ice]['mean'][-2]), 1./(1 - FAR[metric]['rcp85'][ice]['min'][-2]), 1./(1 - FAR[metric]['rcp85'][ice]['max'][-2])
    if cnt == 4:
        plt.legend(['Historical', 'RCP4.5', 'RCP8.5', 'Hist. (min/max)', 'RCP4.5 (min/max)', 'RCP8.5 (min/max)'], loc='upper left', ncol=1, fontsize=10)
    plt.plot([1850, 2100], [1, 1], 'k--', color='0.4', zorder=10)
    plt.xlim(1850, 2100)
    #plt.ylim(0, 400)
    plt.title(titles[cnt-1] + ' MHW days')
    if np.mod(cnt,2) == 1:
        plt.ylabel(r'Probability Ratio')
    if cnt <= 2:
        AX.set_xticklabels([])

# plt.savefig(figfolder + modelDir + '/MHW_ProbRatio.png', bbox_inches='tight', dpi=300)



#
# I-D Phase space
#

# PDF of histNat data
X, Y = np.mgrid[0:2000:1, 0:4:0.05]
positions = np.vstack([X.ravel(), Y.ravel()])
kernel = stats.gaussian_kde(np.array([ecj.nonans(MHW_tsNat['duration_new']['histNat'][ice].flatten()), ecj.nonans(MHW_tsNat['intensity_max_max_new']['histNat'][ice].flatten())]))
pdf = np.reshape(kernel(positions), X.shape)

# Theoretical biological impact
impact = np.zeros(X.shape)
a = 0.005
pmort = [0, 0.6, 0.8, 0.9, 0.95, 0.99, 1]
for y in pmort[1:-1]:
    I = np.log(y/(1-y))/(X*a)
    I = np.sqrt(I)
    impact[Y>=I] = y

plt.figure()
plt.clf()

plt.contourf(X, Y, impact, levels=pmort, cmap=plt.cm.afmhot_r) #, zorder =5)
plt.clim(0.6, 1.3)
plt.xlim(10, 1000)
plt.ylim(1, 2.7)
H = plt.colorbar()
H.set_label('Probability of mortality')

plt.semilogx(MHW_tsNat['duration_new']['histNat'][ice].flatten(), MHW_tsNat['intensity_max_max_new']['histNat'][ice].flatten(), 'o', linewidth=2, markerfacecolor=(0.5,0.5,1), markeredgecolor=(0.5,0.5,1), zorder=10)

plt.contour(X, Y, pdf, levels=[0.05], color='b', linewidth=2, zorder=20)
plt.contour(X, Y, pdf, levels=[0.01], color='b', linewidth=2, zorder=20)
plt.contour(X, Y, pdf, levels=[0.001], color='b', linewidth=2, zorder=20)
tt = years_ts['hist'] >= 1850
plt.semilogx(MHW_tsMM['duration_new']['hist'][ice][tt], MHW_tsMM['intensity_max_max_new']['hist'][ice][tt], '-', color=(0.3,0.3,0.3), linewidth=2, zorder=30)
plt.semilogx(MHW_tsMM['duration_new'][rcp][ice], MHW_tsMM['intensity_max_max_new'][rcp][ice], '-', color=(1,0.,0.), linewidth=2, zorder=30)
for year in [1980, 1990, 2000]:
    tt = years_ts['hist'] == year
    plt.semilogx(MHW_tsMM['duration_new']['hist'][ice][tt], MHW_tsMM['intensity_max_max_new']['hist'][ice][tt], 'ko', markersize=10, zorder=40)
    plt.text(MHW_tsMM['duration_new']['hist'][ice][tt], MHW_tsMM['intensity_max_max_new']['hist'][ice][tt] + 0.08, str(year), fontsize=12, zorder=50)
for year in [2010, 2020, 2030]:
    tt = years_ts[rcp] == year
    plt.semilogx(MHW_tsMM['duration_new'][rcp][ice][tt], MHW_tsMM['intensity_max_max_new'][rcp][ice][tt], 'ko', markersize=10, zorder=40)
    plt.text(MHW_tsMM['duration_new'][rcp][ice][tt], MHW_tsMM['intensity_max_max_new'][rcp][ice][tt] + 0.08, str(year), fontsize=12, zorder=50)
plt.xlim(9, 2000)
plt.ylim(0.1, 2.8)
plt.legend(['Natural world', 'Historical world', 'Future world'], loc='lower right')

plt.xlabel('Duration [days]')
plt.ylabel('Intensity [$^\circ$C]')

# plt.savefig(figfolder + modelDir + '/phaseSpace_ID_emergence.png', bbox_inches='tight', dpi=300)
# plt.savefig(figfolder + modelDir + '/phaseSpace_ID_emergence.pdf', bbox_inches='tight', dpi=150)


# Moving envelope version

def gaussianKernelSum(X, Y, Z, std):
    '''
    Estimates the 2D probability distribution of scattered data Z on a
    regular grid (X, Y). Does so by summing Gaussian pdfs centred on
    each location Z with variance std[0] in the X dimension and STD[1]
    in the Y dimention.

     X, Y    Matrices of coordinates, in np.meshgrid format
     Z       (2 x N) array of data points, where N is the number of
             points. The first index of dimension zero corresponds
             to X, while the second index corresponds to Y.
     std     Standard deviation of Gaussians summed to generate pdf.
             List of length 2.
    '''
    # Initialize pdf to zero
    pdf = np.zeros(X.shape)
    # Add to it a gaussian centred on each of the data points, with the appropriate variance
    for i in range(Z.shape[1]):
        pdf += (1./(2*np.pi*std[0]*std[1]))*np.exp(-0.5*(X - Z[0,i])**2/std[0]**2 - 0.5*(Y - Z[1,i])**2/std[1]**2)
    #
    return pdf

def findConfdidenceContour(pdf, alpha):
    '''
    Find contour level that encloses (1-alpha)*100% of a 2D PDF.
    Assume underlying dimensions are regularly spaced.
    Does not assume integral of PDF is one.
    '''
    # Normalise pdf by integral
    pdf /= pdf.sum()
    # 
    dp = 0.01
    pmax = pdf.max()
    for lev in np.arange(dp*pmax, pmax, dp*pmax):
        pdfSum = pdf[pdf >= lev].sum()
        if pdfSum >= 1-alpha:
            break
    #
    return lev
    

dX = 1
dY = 0.05
X, Y = np.mgrid[0:3500:dX, 0:10:dY]
positions = np.vstack([X.ravel(), Y.ravel()])

# Theoretical biological impact
impact = np.zeros(X.shape)
a = 0.005
pmort = [0, 0.6, 0.8, 0.9, 0.95, 0.99, 1]
pmort_all = np.arange(0, 1, 0.001)
for y in pmort_all[1:-1]:
    I = np.log(y/(1-y))/(X*a)
    I = np.sqrt(I)
    impact[Y>=I] = y

# PDF of histNat data
pdf = {}
# histNat
kernel = stats.gaussian_kde(np.array([ecj.nonans(MHW_tsNat['duration_new']['histNat'][ice].flatten()), ecj.nonans(MHW_tsNat['intensity_max_max_new']['histNat'][ice].flatten())]))
pdf['histNat'] = np.reshape(kernel(positions), X.shape)

# hist, 1900-1950
tmp = np.zeros((2,0))
for model in models:
    tt = (years[model]['hist'] >= 1900) * (years[model]['hist'] <= 1950)
    tmp = np.append(tmp, np.array([ecj.nonans(MHW_ts[model]['duration']['hist'][ice][tt,:].flatten()), ecj.nonans(MHW_ts[model]['intensity_max_max']['hist'][ice][tt,:].flatten())]), axis=1)
#kernel = stats.gaussian_kde(tmp)
pdf['1900_1950'] = gaussianKernelSum(X, Y, tmp, [5, 0.1]) #np.reshape(kernel(positions), X.shape)

# hist, 1970-2000
tmp = np.zeros((2,0))
for model in models:
    tt = (years[model]['hist'] >= 1970) * (years[model]['hist'] <= 2000)
    tmp = np.append(tmp, np.array([ecj.nonans(MHW_ts[model]['duration']['hist'][ice][tt,:].flatten()), ecj.nonans(MHW_ts[model]['intensity_max_max']['hist'][ice][tt,:].flatten())]), axis=1)
#kernel = stats.gaussian_kde(tmp)
pdf['1970_2000'] = gaussianKernelSum(X, Y, tmp, [5, 0.1]) #np.reshape(kernel(positions), X.shape)

# rcp85, 2010-2040
tmp = np.zeros((2,0))
for model in models:
    tt = (years[model]['rcp85'] >= 2005) * (years[model]['rcp85'] <= 2035)
    tmp = np.append(tmp, np.array([ecj.nonans(MHW_ts[model]['duration']['rcp85'][ice][tt,:].flatten()), ecj.nonans(MHW_ts[model]['intensity_max_max']['rcp85'][ice][tt,:].flatten())]), axis=1)
#kernel = stats.gaussian_kde(tmp)
pdf['2005_2035'] = gaussianKernelSum(X, Y, tmp, [10, 0.1]) #np.reshape(kernel(positions), X.shape)

# rcp85, 2050-2080
tmp = np.zeros((2,0))
for model in models:
    tt = (years[model]['rcp85'] >= 2050) * (years[model]['rcp85'] <= 2080)
    tmp = np.append(tmp, np.array([ecj.nonans(MHW_ts[model]['duration']['rcp85'][ice][tt,:].flatten()), ecj.nonans(MHW_ts[model]['intensity_max_max']['rcp85'][ice][tt,:].flatten())]), axis=1)
#kernel = stats.gaussian_kde(tmp)
pdf['2050_2080'] = gaussianKernelSum(X, Y, tmp, [100, 0.15]) #np.reshape(kernel(positions), X.shape)




plt.figure()
plt.clf()

plt.subplot(1,2,1)

plt.contourf(X, Y, impact, levels=pmort, cmap=plt.cm.afmhot_r) #, zorder =5)
plt.clim(0.6, 1.3)
plt.xlim(1, 3.5e3) #10, 2000)
plt.ylim(0, 7) #10)
H = plt.colorbar()
H.set_label('Probability of ecological impact')

plt.contour(X, Y, pdf['histNat'], levels=[findConfdidenceContour(pdf['histNat'], 0.05)], colors='b', linewidth=2, zorder=20)
plt.contour(X, Y, pdf['1900_1950'], levels=[findConfdidenceContour(pdf['1900_1950'], 0.05)], colors='0.5', linewidth=2, zorder=30)
plt.contour(X, Y, pdf['1970_2000'], levels=[findConfdidenceContour(pdf['1970_2000'], 0.05)], colors='k', linewidth=2, zorder=30)
plt.contour(X, Y, pdf['2005_2035'], levels=[findConfdidenceContour(pdf['2005_2035'], 0.05)], colors='r', linewidth=2, zorder=30)
plt.contour(X, Y, pdf['2050_2080'], levels=[findConfdidenceContour(pdf['2050_2080'], 0.05)], colors='m', linewidth=2, zorder=30)

plt.plot([-100, -10], [1, 2], 'b', linewidth=2)
plt.plot([-100, -10], [1, 2], color='0.5', linewidth=2)
plt.plot([-100, -10], [1, 2], 'k', linewidth=2)
plt.plot([-100, -10], [1, 2], 'r', linewidth=2)
plt.plot([-100, -10], [1, 2], 'm', linewidth=2)

plt.legend(['Natural World', 'Historical: 1900-1950', 'Historical: 1970-2000', 'Present (RCP8.5): 2005-2035', 'Future (RCP8.5): 2050-2080'], loc='upper left')

plt.semilogx(MHW_tsNat['duration_new']['histNat'][ice].flatten(), MHW_tsNat['intensity_max_max_new']['histNat'][ice].flatten(), 'o', linewidth=2, markerfacecolor=(0.5,0.5,1), markeredgecolor=(0.5,0.5,1), zorder=10, markersize=0)

plt.semilogx(MHW_tsNat['duration_new']['histNat'][ice].flatten(), MHW_tsNat['intensity_max_max_new']['histNat'][ice].flatten(), 'o', linewidth=2, markerfacecolor=(0.5,0.5,1), markeredgecolor=(0.5,0.5,1))

plt.xlabel('MHW Duration [days]')
plt.ylabel('MHW Intensity [$^\circ$C]')

plt.subplot(1,2,2)

for p in pmort[1:-1]:
    plt.semilogx(np.arange(1,3500,1), np.sqrt(np.log(p/(1-p))/(np.arange(1,3500,1)*a)), '-', color='0.25', linewidth=0.75)
plt.xlim(1, 3.5e3) #10, 2000)
plt.ylim(0, 7) #10)

plt.xlabel('MHW Duration [days]')
plt.ylabel('MHW Intensity [$^\circ$C]')

plt.colorbar() # To squish the right panel to match the left panel, will result in an error but that's OK

# plt.savefig(figfolder + modelDir + '/phaseSpace_ID_emergence_PDFs.png', bbox_inches='tight', dpi=300)

# Alternate, based on barnacles and seagrass shit
# Theoretical lines
impact_DHWs = np.zeros(X.shape)
for DHWs in np.arange(0.05, 20, 0.05):
    impact_DHWs[X*Y >= DHWs*7] = DHWs

plt.clf()

plt.subplot(1,2,1)

#plt.contourf(X, Y, impact_DHWs, levels=[0, 4, 8, 12], cmap=plt.cm.afmhot_r) #, zorder =5)
plt.contourf(X, Y, impact_DHWs, levels=[0, 4, 6, 8, 10, 12, 14, 16, 18, 20], cmap=plt.cm.afmhot_r) #, zorder =5)
plt.clim(2, 33)
plt.xlim(1, 3.5e3) #10, 2000)
plt.ylim(0, 7) #10)
H = plt.colorbar()
H.set_label('Increasing probability of ecological impact')
H.set_ticklabels([])
H.set_ticks([])

plt.contour(X, Y, pdf['histNat'], levels=[findConfdidenceContour(pdf['histNat'], 0.05)], colors='b', linewidth=2, zorder=20)
plt.contour(X, Y, pdf['1900_1950'], levels=[findConfdidenceContour(pdf['1900_1950'], 0.05)], colors='0.5', linewidth=2, zorder=30)
plt.contour(X, Y, pdf['1970_2000'], levels=[findConfdidenceContour(pdf['1970_2000'], 0.05)], colors='k', linewidth=2, zorder=30)
plt.contour(X, Y, pdf['2005_2035'], levels=[findConfdidenceContour(pdf['2005_2035'], 0.05)], colors='r', linewidth=2, zorder=30)
plt.contour(X, Y, pdf['2050_2080'], levels=[findConfdidenceContour(pdf['2050_2080'], 0.05)], colors='m', linewidth=2, zorder=30)

plt.plot([-100, -10], [1, 2], 'b', linewidth=2)
plt.plot([-100, -10], [1, 2], color='0.5', linewidth=2)
plt.plot([-100, -10], [1, 2], 'k', linewidth=2)
plt.plot([-100, -10], [1, 2], 'r', linewidth=2)
plt.plot([-100, -10], [1, 2], 'm', linewidth=2)

plt.legend(['Natural World', 'Historical: 1900-1950', 'Historical: 1970-2000', 'Present (RCP8.5): 2005-2035', 'Future (RCP8.5): 2050-2080'], loc='upper left')

plt.semilogx(MHW_tsNat['duration_new']['histNat'][ice].flatten(), MHW_tsNat['intensity_max_max_new']['histNat'][ice].flatten(), 'o', linewidth=2, markerfacecolor=(0.5,0.5,1), markeredgecolor=(0.5,0.5,1), zorder=10, markersize=0)

plt.semilogx(MHW_tsNat['duration_new']['histNat'][ice].flatten(), MHW_tsNat['intensity_max_max_new']['histNat'][ice].flatten(), 'o', linewidth=2, markerfacecolor=(0.5,0.5,1), markeredgecolor=(0.5,0.5,1))

plt.xlabel('MHW Duration [days]')
plt.ylabel('MHW Intensity [$^\circ$C]')

# plt.savefig(figfolder + modelDir + '/phaseSpace_ID_emergence_PDFs_DHWs.png', bbox_inches='tight', dpi=300)


#for model in models:
#    for ens in range(MHW_ts[model]['intensity_max_max']['hist'][ice].shape[1]):
#        plt.semilogx(MHW_ts[model]['duration']['hist'][ice][:,ens], MHW_ts[model]['intensity_max_max']['hist'][ice][:,ens], 'ko')

#for model in models:
#    for ens in range(MHW_ts[model]['intensity_max_max']['rcp85'][ice].shape[1]):
#        plt.semilogx(MHW_ts[model]['duration']['rcp85'][ice][:,ens], MHW_ts[model]['intensity_max_max']['rcp85'][ice][:,ens], 'ro')

#
# Load in NOAA OI SST for validation/comparison
#

outfile = 'mhw_census.NOAAOISST.1982_2005.npz'
data = np.load(header + outfile)
years_obs = data['years']
lon_obs = data['lon_map']
lat_obs = data['lat_map']
llon_obs, llat_obs = np.meshgrid(lon_obs, lat_obs)
SST_mean = data['SST_mean']
MHW_cnt = data['MHW_cnt']
MHW_dur = data['MHW_dur']
MHW_max = data['MHW_max']
MHW_td = data['MHW_td']
SST_ts = data['SST_ts']
MHW_cnt_ts = data['MHW_cnt_ts']
MHW_dur_ts = data['MHW_dur_ts']
MHW_max_ts = data['MHW_max_ts']
MHW_td_ts = data['MHW_td_ts']
# Fix some stuff...
SST_ts[SST_ts==0] = np.nan
MHW_dur_ts[MHW_cnt_ts==0] = np.nan
MHW_max_ts[MHW_cnt_ts==0] = np.nan
MHW_td_ts[MHW_cnt_ts==0] = np.nan
# Global mean time series
SST_ts_glob = np.zeros(SST_ts.shape[2])
MHW_cnt_ts_glob = np.zeros(MHW_cnt_ts.shape[2])
MHW_dur_ts_glob = np.zeros(MHW_dur_ts.shape[2])
MHW_max_ts_glob = np.zeros(MHW_max_ts.shape[2])
MHW_td_ts_glob = np.zeros(MHW_td_ts.shape[2])
scaling = np.cos(llat_obs*np.pi/180)
for tt in range(MHW_cnt_ts.shape[2]):
    # Create mask
    mask = np.ones(llat_obs.shape)
    mask[np.isnan(SST_ts[:,:,tt])] = np.nan
    # SST
    SST_ts_glob[tt] = np.average(SST_ts[:,:,tt][~np.isnan(mask)], weights=scaling[~np.isnan(mask)])
    # Create mask
    mask = np.ones(llat_obs.shape)
    mask[np.isnan(MHW_dur_ts[:,:,tt])] = np.nan
    # Count, Duration, Maximum intensity, Total MHW days
    MHW_cnt_ts_glob[tt] = np.average(MHW_cnt_ts[:,:,tt][~np.isnan(mask)], weights=scaling[~np.isnan(mask)])
    MHW_dur_ts_glob[tt] = np.average(MHW_dur_ts[:,:,tt][~np.isnan(mask)], weights=scaling[~np.isnan(mask)])
    MHW_max_ts_glob[tt] = np.average(MHW_max_ts[:,:,tt][~np.isnan(mask)], weights=scaling[~np.isnan(mask)])
    MHW_td_ts_glob[tt] = np.average(MHW_td_ts[:,:,tt][~np.isnan(mask)], weights=scaling[~np.isnan(mask)])

# Regrid to same grid as models
MHW_cnt_regrid = interp.griddata((llon_obs.flatten(), llat_obs.flatten()), MHW_cnt.flatten() , (lon, lat), method='nearest')
MHW_dur_regrid = interp.griddata((llon_obs.flatten(), llat_obs.flatten()), MHW_dur.flatten() , (lon, lat), method='nearest')
MHW_max_regrid = interp.griddata((llon_obs.flatten(), llat_obs.flatten()), MHW_max.flatten() , (lon, lat), method='nearest')
MHW_td_regrid = interp.griddata((llon_obs.flatten(), llat_obs.flatten()), MHW_td.flatten() , (lon, lat), method='nearest')
# Spatial average to same grid as models
MHW_cnt_regrid = np.nan*np.zeros(lon.shape)
MHW_dur_regrid = np.nan*np.zeros(lon.shape)
MHW_max_regrid = np.nan*np.zeros(lon.shape)
MHW_td_regrid = np.nan*np.zeros(lon.shape)
lon_vec = lon[0,:]
lat_vec = lat[:,0]
dl = 2
for i in range(len(lon_vec)):
    for j in range(len(lat_vec)):
        ii = np.where((lon_obs >= lon_vec[i]-dl/2.) * (lon_obs <= lon_vec[i] + dl/2.))[0]
        jj = np.where((lat_obs >= lat_vec[j]-dl/2.) * (lat_obs <= lat_vec[j] + dl/2.))[0]
        MHW_cnt_regrid[j,i] = np.nanmean(MHW_cnt[jj,:][:,ii].flatten())
        MHW_dur_regrid[j,i] = np.nanmean(MHW_dur[jj,:][:,ii].flatten())
        MHW_max_regrid[j,i] = np.nanmean(MHW_max[jj,:][:,ii].flatten())
        MHW_td_regrid[j,i] = np.nanmean(MHW_td[jj,:][:,ii].flatten())



plt.clf()
exp = 'hist'
plt.subplot(4,3,1)
plt.contourf(lon, lat, MHW_cnt_regrid)
plt.colorbar()
plt.subplot(4,3,2)
plt.contourf(lon, lat, MHW_meanMM['count_new']['hist'][ice] - MHW_cnt_regrid, cmap=plt.cm.RdBu_r)
plt.colorbar()
plt.subplot(4,3,3)
plt.plot(years_obs, MHW_cnt_ts_glob, 'k-', linewidth=2)
for model in models:
    for ens in range(NENS[model][exp]):
        plt.plot(years[model][exp], MHW_ts[model]['count_new'][exp][ice][:,ens], 'r-')
plt.xlim(1982, 2005)

plt.subplot(4,3,4)
plt.contourf(lon, lat, MHW_max_regrid)
plt.colorbar()
plt.subplot(4,3,5)
plt.contourf(lon, lat, MHW_meanMM['intensity_max_max_new']['hist'][ice] - MHW_max_regrid, cmap=plt.cm.RdBu_r)
plt.colorbar()
plt.subplot(4,3,6)
plt.plot(years_obs, MHW_max_ts_glob, 'k-', linewidth=2)
for model in models:
    for ens in range(NENS[model][exp]):
        plt.plot(years[model][exp], MHW_ts[model]['intensity_max_max_new'][exp][ice][:,ens], 'r-')
plt.xlim(1982, 2005)

plt.subplot(4,3,7)
plt.contourf(lon, lat, MHW_dur_regrid)
plt.colorbar()
plt.subplot(4,3,8)
plt.contourf(lon, lat, MHW_meanMM['duration_new']['hist'][ice] - MHW_dur_regrid, cmap=plt.cm.RdBu_r)
plt.colorbar()
plt.subplot(4,3,9)
plt.plot(years_obs, MHW_dur_ts_glob, 'k-', linewidth=2)
for model in models:
    for ens in range(NENS[model][exp]):
        plt.plot(years[model][exp], MHW_ts[model]['duration_new'][exp][ice][:,ens], 'r-')
plt.xlim(1982, 2005)

plt.subplot(4,3,10)
plt.contourf(lon, lat, MHW_td_regrid)
plt.colorbar()
plt.subplot(4,3,11)
plt.contourf(lon, lat, MHW_meanMM['total_days']['hist'][ice] - MHW_td_regrid, cmap=plt.cm.RdBu_r)
plt.colorbar()
plt.subplot(4,3,12)
plt.plot(years_obs, MHW_td_ts_glob, 'k-', linewidth=2)
for model in models:
    for ens in range(NENS[model][exp]):
        plt.plot(years[model][exp], MHW_ts[model]['total_days'][exp][ice][:,ens], 'r-')
plt.xlim(1982, 2005)

# Just max intensity and total days
plt.figure(figsize=(13,8))
plt.clf()

exp = 'hist'
plt.subplot(3,2,1, axisbg=bg_col)
proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], dashes=[3,900])
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True], dashes=[3,900])
plt.contourf(lonproj, latproj, MHW_max_regrid, levels=np.arange(0, 5+0.5, 0.5), cmap=plt.cm.gist_heat_r)
H = plt.colorbar()
H.set_label(r'[$^\circ$C]')
plt.clim(0.5, 5)
plt.title('(A) Maximum intensity (Obs)')

plt.subplot(3,2,3, axisbg=bg_col)
proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], dashes=[3,900])
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True], dashes=[3,900])
plt.contourf(lonproj, latproj, MHW_meanMM['intensity_max_max_new']['hist'][ice] - MHW_max_regrid, levels=[-3,-2,-1.5,-1,-0.5,0.5,1,1.5,2,3], cmap=plt.cm.RdBu_r)
H = plt.colorbar()
H.set_label(r'[$^\circ$C]')
#plt.clim(0.5, 5)
plt.title('(B) Model - Obs')

plt.subplot(3,2,5)
for model in models:
    for ens in range(NENS[model][exp]):
        plt.plot(years[model][exp], MHW_ts[model]['intensity_max_max_new'][exp][ice][:,ens], '-', color=(1,0.5,0.5))
plt.plot(years_obs, MHW_max_ts_glob, 'k-', linewidth=2)
plt.xlim(1982, 2005)
plt.ylabel('Intensity [$^\circ$C]')
plt.title('(C) Global mean time series')

plt.subplot(3,2,2, axisbg=bg_col)
proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], dashes=[3,900])
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True], dashes=[3,900])
plt.contourf(lonproj, latproj, MHW_td_regrid, levels=np.arange(10, 38+4, 4), cmap=plt.cm.gist_heat_r)
H = plt.colorbar()
H.set_label('[days]')
plt.clim(14, 45)
plt.title('(D) Total annual MHW days (Obs)')

plt.subplot(3,2,4, axisbg=bg_col)
proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], dashes=[3,900])
proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True], dashes=[3,900])
plt.contourf(lonproj, latproj, MHW_meanMM['total_days']['hist'][ice] - MHW_td_regrid, levels=[-25,-20,-15,-10,-5,5,10,15,20,25], cmap=plt.cm.RdBu_r)
H = plt.colorbar()
H.set_label('[days]')
#plt.clim(25, 45)
plt.title('(E) Model - Obs')

plt.subplot(3,2,6)
for model in models:
    for ens in range(NENS[model][exp]):
        lineMod, = plt.plot(years[model][exp], MHW_ts[model]['total_days'][exp][ice][:,ens], '-', color=(1,0.5,0.5))
lineObs, = plt.plot(years_obs, MHW_td_ts_glob, 'k-', linewidth=2)
plt.xlim(1982, 2005)
plt.ylabel('Total annual MHW days [days]')
plt.legend([lineObs, lineMod], ['NOAA OI SST (Obs.)', 'CMIP5 (indiv. ens. members)'])
plt.title('(F) Global mean time series')

# plt.savefig(figfolder + modelDir + '/validation.png', bbox_inches='tight', dpi=300)

#
# Generate plots for individual models
#

plt.figure(figsize=(12,8))
cnt = 0
for model in models:
    print model
    plt.clf()
    #
    metric = 'intensity_max_max_new'
    plt.subplot(3,2,1, axisbg=bg_col)
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], dashes=[3,900])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,False], dashes=[3,900])
    plt.contourf(lonproj, latproj, MHW_meanEM[metric]['hist'][ice][:,:,cnt], levels=np.arange(0, 5+0.5, 0.5), cmap=plt.cm.gist_heat_r)
    H = plt.colorbar()
    H.set_label(r'[$^\circ$C]')
    plt.clim(0.5, 5)
    plt.title('(A) Maximum intensity (historical ' + str(meanYears['hist'][0]) + '-' + str(meanYears['hist'][1]) + ')')
    plt.subplot(3,2,3, axisbg=bg_col)
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], dashes=[3,900])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True], dashes=[3,900])
    plt.contourf(lonproj, latproj, MHW_diffEM[metric][rcp][ice][:,:,cnt], levels=[-4,-2,-1.5,-1,-0.5,0.5,1,1.5,2,4], cmap=plt.cm.RdBu_r)
    H = plt.colorbar()
    H.set_label(r'[$^\circ$C]')
    plt.clim(-3.5, 3.5)
    #plt.contourf(lonproj, latproj, (MHW_agrMM[metric][rcp][ice] >= Nagr).astype(float), hatches=['', hatch], levels=[0., 0.5, 1.0], colors='none')
    plt.title('(B) (RCP8.5 ' + str(meanYears[rcp][0]) + '-' + str(meanYears[rcp][1]) + ') - (historical ' + str(meanYears['hist'][0]) + '-' + str(meanYears['hist'][1]) + ')')
    #plt.subplot(3,2,5)
    plt.subplot2grid((3,10), (2,0), rowspan=1, colspan=4)
    plt.plot(years_ts['hist'], MHW_tsEM[metric]['hist'][ice][:,cnt], 'k-', linewidth=2, zorder=30)
    plt.plot(years_ts[rcp], MHW_tsEM[metric]['rcp45'][ice][:,cnt], '-', linewidth=2, color=(0.6,0.3,0), zorder=30)
    plt.plot(years_ts[rcp], MHW_tsEM[metric]['rcp85'][ice][:,cnt], '-', linewidth=2, color=(0.9,0,0), zorder=30)
    #
    tsNat = MHW_ts[model][metric]['histNat'][ice][np.in1d(years[model]['histNat'], years_ts['histNat']),:].flatten()
    plt.fill_between([1850, 2100], np.min(tsNat), np.max(tsNat), color=nat_col['max'], zorder=10)
    plt.fill_between([1850, 2100], np.percentile(tsNat, 2.5), np.percentile(tsNat, 97.5), color=nat_col['95'], zorder=10)
    plt.fill_between([1850, 2100], np.percentile(tsNat, 16.5), np.percentile(tsNat, 83.5), color=nat_col['67'], zorder=10)
    plt.fill_between(years[model]['hist'], np.nanmin(MHW_ts[model][metric]['hist'][ice], axis=1), np.nanmax(MHW_ts[model][metric]['hist'][ice], axis=1), color='0.5', alpha=0.8, zorder=20)
    plt.fill_between(years[model]['rcp45'], np.nanmin(MHW_ts[model][metric]['rcp45'][ice], axis=1), np.nanmax(MHW_ts[model][metric]['rcp45'][ice], axis=1), color=(1,0.7,0.3), alpha=0.8, zorder=20)
    plt.fill_between(years[model]['rcp85'], np.nanmin(MHW_ts[model][metric]['rcp85'][ice], axis=1), np.nanmax(MHW_ts[model][metric]['rcp85'][ice], axis=1), color=(0.9,0.5,0.5), alpha=0.8, zorder=20)
    plt.xlim(1850, 2100)
    plt.ylabel(r'Intensity [$^\circ$C]')
    plt.title('(C) Global mean time series')
    plt.legend(['Historical', 'RCP4.5', 'RCP8.5', 'Nat. (min/max)', 'Nat. (95% CI)', 'Nat. (66% CI)', 'Hist. (min/max)', 'RCP4.5 (min/max)', 'RCP8.5 (min/max)'], loc='upper left', ncol=2, fontsize=8)
    #
    metric = 'total_days'
    plt.subplot(3,2,2, axisbg=bg_col)
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], dashes=[3,900])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,False], dashes=[3,900])
    plt.contourf(lonproj, latproj, MHW_meanEM[metric]['hist'][ice][:,:,cnt], levels=np.arange(23, 37+2, 2), cmap=plt.cm.gist_heat_r)
    H = plt.colorbar()
    H.set_label('[days]')
    plt.clim(25, 45)
    plt.title('(D) Total Annual MHW Days (historical ' + str(meanYears['hist'][0]) + '-' + str(meanYears['hist'][1]) + ')')
    plt.subplot(3,2,4, axisbg=bg_col)
    proj.fillcontinents(color=(0,0,0), lake_color=None, ax=None, zorder=None, alpha=None)
    proj.drawparallels(range(domain_draw[0],domain_draw[2]+1,dlat), labels=[True,False,False,False], dashes=[3,900])
    proj.drawmeridians(range(domain_draw[1],domain_draw[3]+1,dlon), labels=[False,False,False,True], dashes=[3,900])
    plt.contourf(lonproj, latproj, MHW_diffEM[metric][rcp][ice][:,:,cnt], levels=np.arange(-400,400+1,80), cmap=plt.cm.RdBu_r)
    H = plt.colorbar()
    H.set_label('[days]')
    plt.clim(-365, 365)
    #plt.contourf(lonproj, latproj, (MHW_agrMM[metric][rcp][ice] >= Nagr).astype(float), hatches=['', hatch], levels=[0., 0.5, 1.0], colors='none')
    plt.title('(E) (RCP8.5 ' + str(meanYears[rcp][0]) + '-' + str(meanYears[rcp][1]) + ') - (historical ' + str(meanYears['hist'][0]) + '-' + str(meanYears['hist'][1]) + ')')
    #plt.subplot(3,2,6)
    plt.subplot2grid((3,20), (2,11), rowspan=1, colspan=8)
    #
    tsNat = MHW_ts[model][metric]['histNat'][ice][np.in1d(years[model]['histNat'], years_ts['histNat']),:].flatten()
    plt.fill_between([1850, 2100], np.min(tsNat), np.max(tsNat), color=nat_col['max'])
    plt.fill_between([1850, 2100], np.percentile(tsNat, 2.5), np.percentile(tsNat, 97.5), color=nat_col['95'])
    plt.fill_between([1850, 2100], np.percentile(tsNat, 16.5), np.percentile(tsNat, 83.5), color=nat_col['67'])
    plt.plot([1850, 2100], [365, 365], '--', color='0.4')
    plt.fill_between(years[model]['hist'], np.nanmin(MHW_ts[model][metric]['hist'][ice], axis=1), np.nanmax(MHW_ts[model][metric]['hist'][ice], axis=1), color='0.5', alpha=0.8)
    plt.plot(years_ts['hist'], MHW_tsEM[metric]['hist'][ice][:,cnt], 'k-', linewidth=2)
    plt.fill_between(years[model]['rcp45'], np.nanmin(MHW_ts[model][metric]['rcp45'][ice], axis=1), np.nanmax(MHW_ts[model][metric]['rcp45'][ice], axis=1), color=(1,0.7,0.3), alpha=0.8)
    plt.plot(years_ts[rcp], MHW_tsEM[metric]['rcp45'][ice][:,cnt], '-', linewidth=2, color=(0.6,0.3,0))
    plt.fill_between(years[model]['rcp85'], np.nanmin(MHW_ts[model][metric]['rcp85'][ice], axis=1), np.nanmax(MHW_ts[model][metric]['rcp85'][ice], axis=1), color=(0.9,0.5,0.5), alpha=0.8)
    plt.plot(years_ts[rcp], MHW_tsEM[metric]['rcp85'][ice][:,cnt], '-', linewidth=2, color=(0.9,0,0))
    plt.xlim(1850, 2100)
    plt.ylim(0, 400)
    plt.ylabel(r'Total Annual MHW Days [days]')
    plt.title('(F) Global mean time series')
    #
    cnt += 1
    #
    plt.savefig(figfolder + modelDir + '/MHW_meanDiffTimeseries_Vert_' + model + '.png', bbox_inches='tight', dpi=300)

#
# Include observed time series
#

# Load obs
data = np.load('/home/oliver/data/MHWs/Trends/mhw_proxies_GLM.1900.2016.ALL_ts.npz')
years_obs = data['years_data']
MHW_f_ts_glob_aggAll = data['MHW_f_ts_glob_aggAll']
MHW_d_ts_glob_aggAll = data['MHW_d_ts_glob_aggAll']
MHW_td_ts_glob_aggAll = data['MHW_td_ts_glob_aggAll']
MHW_f_sigMod_agg = data['MHW_f_sigMod_agg']
MHW_d_sigMod_agg = data['MHW_d_sigMod_agg']
MHW_td_sigMod_agg = data['MHW_td_sigMod_agg']


metric = 'total_days'
#plt.subplot2grid((3,20), (2,11), rowspan=1, colspan=8)
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['min'], MHW_pNat[metric][ice]['max'], color=nat_col['max'])
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['2.5'], MHW_pNat[metric][ice]['97.5'], color=nat_col['95'])
plt.fill_between([1850, 2100], MHW_pNat[metric][ice]['16.5'], MHW_pNat[metric][ice]['83.5'], color=nat_col['67'])
plt.plot([1850, 2100], [365, 365], '--', color='0.4')
plt.fill_between(years_ts['hist'], np.nanmin(MHW_tsEM[metric]['hist'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['hist'][ice], axis=1), color='0.5', alpha=0.8)
plt.plot(years_ts['hist'], MHW_tsMM[metric]['hist'][ice], 'k-', linewidth=2)
plt.fill_between(years_ts['rcp45'], np.nanmin(MHW_tsEM[metric]['rcp45'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['rcp45'][ice], axis=1), color=(0.6,0.75,0.6), alpha=0.8)
plt.plot(years_ts[rcp], MHW_tsMM[metric]['rcp45'][ice], '-', linewidth=2, color=(0,0.6,0))
plt.fill_between(years_ts['rcp85'], np.nanmin(MHW_tsEM[metric]['rcp85'][ice], axis=1), np.nanmax(MHW_tsEM[metric]['rcp85'][ice], axis=1), color=(0.9,0.5,0.5), alpha=0.8)
plt.plot(years_ts[rcp], MHW_tsMM[metric]['rcp85'][ice], '-', linewidth=2, color=(0.9,0,0))
plt.xlim(1850, 2100)
plt.ylim(0, 400)
plt.ylabel(r'Total Annual MHW Days [days]')
plt.title('Global mean time series')

plt.plot(years_obs, MHW_td_ts_glob_aggAll)
