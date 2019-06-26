'''

  Generate time series data from CMIP5 daily 'tos' data

'''

import numpy as np
import scipy as sp
from scipy import stats
from scipy import signal
from scipy import io
import scipy.optimize as opt
from datetime import date
from netCDF4 import Dataset
import os
import ecoliver as ecj

from matplotlib import pyplot as plt
import mpl_toolkits.basemap as bm

#
# GLobal info on runs and ensemble members
#

# mokami
#pathroot = '/mnt/insect/'
# sverdrup
pathroot = '/home/ecoliver/'
# pakasai
pathroot = '/media/ecoliver/Pentagram/'
# tikoraluk
pathroot = '/data/home/oliver/'

models = ['ACCESS1-3', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'HadGEM2-ES', 'CanESM2', 'CSIRO-Mk3-6-0', 'CNRM-CM5']
#models = ['IPSL-CM5A-LR']
experiments = ['histNat', 'hist', 'rcp45', 'rcp85']
header = {}
Ly = {}
t = {}
T = {}
year = {}
month = {}
day = {}
doy = {}
lon = {}
lat = {}
X = {}
Y = {}
years = {}
Ty = {}
KelvinFix = {}
for model in models:
    header[model] = {}
    t[model] = {}
    T[model] = {}
    year[model] = {}
    month[model] = {}
    day[model] = {}
    doy[model] = {}
    years[model] = {}
    Ty[model] = {}
    KelvinFix[model] = {}

#
# Setup variables for each model
#

# ACCESS1-3
model = 'ACCESS1-3'
header[model]['hist'] = [pathroot+'data/CMIP5/CSIRO-BOM/ACCESS1-3/historical/day/ocean/day/r1i1p1/files/tos_20121112/',
                         pathroot+'data/CMIP5/CSIRO-BOM/ACCESS1-3/historical/day/ocean/day/r2i1p1/files/tos_20121112/',
                         pathroot+'data/CMIP5/CSIRO-BOM/ACCESS1-3/historical/day/ocean/day/r3i1p1/files/tos_20121112/']
header[model]['histNat'] = [pathroot+'data/CMIP5/CSIRO-BOM/ACCESS1-3/historicalNat/day/ocean/day/r1i1p1/files/tos_20130906/',
                            pathroot+'data/CMIP5/CSIRO-BOM/ACCESS1-3/historicalNat/day/ocean/day/r2i1p1/files/tos_20130912/',
                            pathroot+'data/CMIP5/CSIRO-BOM/ACCESS1-3/historicalNat/day/ocean/day/r3i1p1/files/tos_20140501/']
header[model]['rcp85'] = [pathroot+'data/CMIP5/CSIRO-BOM/ACCESS1-3/rcp85/day/ocean/day/r1i1p1/files/tos_20120413/']
header[model]['rcp45'] = [pathroot+'data/CMIP5/CSIRO-BOM/ACCESS1-3/rcp45/day/ocean/day/r1i1p1/v20120413/tos/']
# Time and date vectors
Ly[model] = 365.25
t[model]['hist'], tmp, T[model]['hist'], year[model]['hist'], month[model]['hist'], day[model]['hist'], doy[model]['hist'] = ecj.timevector([1850,1,1], [2005,12,31])
t[model]['histNat'], tmp, T[model]['histNat'], year[model]['histNat'], month[model]['histNat'], day[model]['histNat'], doy[model]['histNat'] = ecj.timevector([1850,1,1], [2005,12,31])
t[model]['rcp85'], tmp, T[model]['rcp85'], year[model]['rcp85'], month[model]['rcp85'], day[model]['rcp85'], doy[model]['rcp85'] = ecj.timevector([2006,1,1], [2100,12,31])
t[model]['rcp45'], tmp, T[model]['rcp45'], year[model]['rcp45'], month[model]['rcp45'], day[model]['rcp45'], doy[model]['rcp45'] = ecj.timevector([2006,1,1], [2100,12,31])
for exp in experiments:
    years[model][exp] = np.unique(year[model][exp])
    Ty[model][exp] = len(years[model][exp])
# Kelvin -> Celcius (only not required for ACCESS hist ens > 0 and non RCP)
# ALL OTHER MODELS WILL NEED A -273.15 OFFSET FOR K->C
KelvinFix[model]['hist'] = [-273.15, 0., 0.]
KelvinFix[model]['histNat'] = [0., 0., 0.]
KelvinFix[model]['rcp85'] = [-273.15]
KelvinFix[model]['rcp45'] = [-273.15]

# IPSL-CM5A-LR
model = 'IPSL-CM5A-LR'
header[model]['hist'] = [pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/historical/day/ocean/day/r1i1p1/v20111010/tos/',
                         pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/historical/day/ocean/day/r2i1p1/v20111010/tos/',
                         pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/historical/day/ocean/day/r3i1p1/v20111010/tos/',
                         pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/historical/day/ocean/day/r4i1p1/v20111010/tos/',
                         pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/historical/day/ocean/day/r5i1p1/v20111119/tos/',
                         pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/historical/day/ocean/day/r6i1p1/v20120526/tos/']
header[model]['histNat'] = [pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/historicalNat/day/ocean/day/r1i1p1/v20130506/tos/',
                            pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/historicalNat/day/ocean/day/r2i1p1/v20130506/tos/',
                            pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/historicalNat/day/ocean/day/r3i1p1/v20130506/tos/']
header[model]['rcp85'] = [pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/rcp85/day/ocean/day/r1i1p1/v20111103/tos/',
                          pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/rcp85/day/ocean/day/r2i1p1/v20110901/tos/',
                          pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/rcp85/day/ocean/day/r3i1p1/v20110901/tos/',
                          pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/rcp85/day/ocean/day/r4i1p1/v20110901/tos/']
header[model]['rcp45'] = [pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/rcp45/day/ocean/day/r1i1p1/',
                          pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/rcp45/day/ocean/day/r2i1p1/',
                          pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/rcp45/day/ocean/day/r3i1p1/',
                          pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/rcp45/day/ocean/day/r4i1p1/']
# Time and date vectors
Ly[model] = 365
t[model]['hist'], tmp, T[model]['hist'], year[model]['hist'], month[model]['hist'], day[model]['hist'], doy[model]['hist'] = ecj.timevector([1850,1,1], [2005,12,31])
t[model]['histNat'], tmp, T[model]['histNat'], year[model]['histNat'], month[model]['histNat'], day[model]['histNat'], doy[model]['histNat'] = ecj.timevector([1850,1,1], [2012,12,31])
t[model]['rcp85'], tmp, T[model]['rcp85'], year[model]['rcp85'], month[model]['rcp85'], day[model]['rcp85'], doy[model]['rcp85'] = ecj.timevector([2006,1,1], [2300,12,31])
t[model]['rcp45'], tmp, T[model]['rcp45'], year[model]['rcp45'], month[model]['rcp45'], day[model]['rcp45'], doy[model]['rcp45'] = ecj.timevector([2006,1,1], [2300,12,31])
for exp in experiments:
    feb29s = ~((month[model][exp]==2) * (day[model][exp]==29))
    t[model][exp] = t[model][exp][feb29s]
    year[model][exp] = year[model][exp][feb29s]
    month[model][exp] = month[model][exp][feb29s]
    day[model][exp] = day[model][exp][feb29s]
    doy[model][exp] = doy[model][exp][feb29s]
    T[model][exp] = len(t[model][exp])
    years[model][exp] = np.unique(year[model][exp])
    Ty[model][exp] = len(years[model][exp])
# Kelvin -> Celcius
KelvinFix[model]['hist'] = [-273.15, -273.15, -273.15, -273.15, -273.15, -273.15]
KelvinFix[model]['histNat'] = [-273.15, -273.15, -273.15]
KelvinFix[model]['rcp85'] = [-273.15, -273.15, -273.15, -273.15]
KelvinFix[model]['rcp45'] = [-273.15, -273.15, -273.15, -273.15]

# IPSL-CM5A-MR
model = 'IPSL-CM5A-MR'
header[model]['hist'] = [pathroot+'data/CMIP5/IPSL/IPSL-CM5A-MR/historical/day/ocean/day/r1i1p1/v20111119/tos/',
                         pathroot+'data/CMIP5/IPSL/IPSL-CM5A-MR/historical/day/ocean/day/r2i1p1/v20120430/tos/',
                         pathroot+'data/CMIP5/IPSL/IPSL-CM5A-MR/historical/day/ocean/day/r3i1p1/v20120804/tos/']

header[model]['histNat'] = [pathroot+'data/CMIP5/IPSL/IPSL-CM5A-MR/historicalNat/day/ocean/day/r1i1p1/v20120804/tos/',
                            pathroot+'data/CMIP5/IPSL/IPSL-CM5A-MR/historicalNat/day/ocean/day/r2i1p1/v20120804/tos/',
                            pathroot+'data/CMIP5/IPSL/IPSL-CM5A-MR/historicalNat/day/ocean/day/r3i1p1/v20120804/tos/']
header[model]['rcp85'] = [pathroot+'data/CMIP5/IPSL/IPSL-CM5A-MR/rcp85/day/ocean/day/r1i1p1/v20111119/tos/']
header[model]['rcp45'] = [pathroot+'data/CMIP5/IPSL/IPSL-CM5A-MR/rcp45/day/ocean/day/r1i1p1/']
# Time and date vectors
Ly[model] = 365
t[model]['hist'], tmp, T[model]['hist'], year[model]['hist'], month[model]['hist'], day[model]['hist'], doy[model]['hist'] = ecj.timevector([1850,1,1], [2005,12,31])
t[model]['histNat'], tmp, T[model]['histNat'], year[model]['histNat'], month[model]['histNat'], day[model]['histNat'], doy[model]['histNat'] = ecj.timevector([1850,1,1], [2012,12,31])
t[model]['rcp85'], tmp, T[model]['rcp85'], year[model]['rcp85'], month[model]['rcp85'], day[model]['rcp85'], doy[model]['rcp85'] = ecj.timevector([2006,1,1], [2300,12,31])
t[model]['rcp45'], tmp, T[model]['rcp45'], year[model]['rcp45'], month[model]['rcp45'], day[model]['rcp45'], doy[model]['rcp45'] = ecj.timevector([2006,1,1], [2300,12,31])
for exp in experiments:
    feb29s = ~((month[model][exp]==2) * (day[model][exp]==29))
    t[model][exp] = t[model][exp][feb29s]
    year[model][exp] = year[model][exp][feb29s]
    month[model][exp] = month[model][exp][feb29s]
    day[model][exp] = day[model][exp][feb29s]
    doy[model][exp] = doy[model][exp][feb29s]
    T[model][exp] = len(t[model][exp])
    years[model][exp] = np.unique(year[model][exp])
    Ty[model][exp] = len(years[model][exp])
# Kelvin -> Celcius
KelvinFix[model]['hist'] = [-273.15, -273.15, -273.15]
KelvinFix[model]['histNat'] = [-273.15, -273.15, -273.15]
KelvinFix[model]['rcp85'] = [-273.15]
KelvinFix[model]['rcp45'] = [-273.15]

# HadGEM2-ES
model = 'HadGEM2-ES'
header[model]['hist'] = [pathroot+'data/CMIP5/MOHC/HadGEM2-ES/historical/day/ocean/day/r1i1p1/v20110131/tos/',
                         pathroot+'data/CMIP5/MOHC/HadGEM2-ES/historical/day/ocean/day/r2i1p1/v20110418/tos/',
                         pathroot+'data/CMIP5/MOHC/HadGEM2-ES/historical/day/ocean/day/r3i1p1/v20110418/tos/',
                         pathroot+'data/CMIP5/MOHC/HadGEM2-ES/historical/day/ocean/day/r4i1p1/v20110418/tos/']
header[model]['histNat'] = [pathroot+'data/CMIP5/MOHC/HadGEM2-ES/historicalNat/day/ocean/day/r1i1p1/v20110728/tos/',
                            pathroot+'data/CMIP5/MOHC/HadGEM2-ES/historicalNat/day/ocean/day/r2i1p1/v20110609/tos/',
                            pathroot+'data/CMIP5/MOHC/HadGEM2-ES/historicalNat/day/ocean/day/r3i1p1/v20110609/tos/',
                            pathroot+'data/CMIP5/MOHC/HadGEM2-ES/historicalNat/day/ocean/day/r4i1p1/v20110609/tos/']
header[model]['rcp85'] = [pathroot+'data/CMIP5/MOHC/HadGEM2-ES/rcp85/day/ocean/day/r1i1p1/v20111128/tos/',
                          pathroot+'data/CMIP5/MOHC/HadGEM2-ES/rcp85/day/ocean/day/r2i1p1/v20120114/tos/',
                          pathroot+'data/CMIP5/MOHC/HadGEM2-ES/rcp85/day/ocean/day/r3i1p1/v20120114/tos/',
                          pathroot+'data/CMIP5/MOHC/HadGEM2-ES/rcp85/day/ocean/day/r4i1p1/v20120114/tos/']
header[model]['rcp45'] = [pathroot+'data/CMIP5/MOHC/HadGEM2-ES/rcp45/day/ocean/day/r1i1p1/',
                          pathroot+'data/CMIP5/MOHC/HadGEM2-ES/rcp45/day/ocean/day/r2i1p1/',
                          pathroot+'data/CMIP5/MOHC/HadGEM2-ES/rcp45/day/ocean/day/r3i1p1/',
                          pathroot+'data/CMIP5/MOHC/HadGEM2-ES/rcp45/day/ocean/day/r4i1p1/']
# Time and date vectors
Ly[model] = 360.
t[model]['hist'], tmp, T[model]['hist'], year[model]['hist'], month[model]['hist'], day[model]['hist'], doy[model]['hist'] = ecj.timevector([1859,12,1], [2005,12,30])
t[model]['histNat'], tmp, T[model]['histNat'], year[model]['histNat'], month[model]['histNat'], day[model]['histNat'], doy[model]['histNat'] = ecj.timevector([1859,12,1], [2019,12,30])
t[model]['rcp85'], tmp, T[model]['rcp85'], year[model]['rcp85'], month[model]['rcp85'], day[model]['rcp85'], doy[model]['rcp85'] = ecj.timevector([2005,12,1], [2300,1,30])
t[model]['rcp45'], tmp, T[model]['rcp45'], year[model]['rcp45'], month[model]['rcp45'], day[model]['rcp45'], doy[model]['rcp45'] = ecj.timevector([2005,12,1], [2300,1,30])
for exp in experiments:
    feb29s = ~((month[model][exp]==2) * (day[model][exp]==29))
    t[model][exp] = t[model][exp][feb29s]
    year[model][exp] = year[model][exp][feb29s]
    month[model][exp] = month[model][exp][feb29s]
    day[model][exp] = day[model][exp][feb29s]
    doy[model][exp] = doy[model][exp][feb29s]
    y360 = ~((day[model][exp] > 30) * (month[model][exp] > 3))
    t[model][exp] = t[model][exp][y360]
    year[model][exp] = year[model][exp][y360]
    month[model][exp] = month[model][exp][y360]
    day[model][exp] = day[model][exp][y360]
    doy[model][exp] = doy[model][exp][y360]
    T[model][exp] = len(t[model][exp])
    years[model][exp] = np.unique(year[model][exp])
    Ty[model][exp] = len(years[model][exp])
# Kelvin -> Celcius
KelvinFix[model]['hist'] = [-273.15, -273.15, -273.15, -273.15]
KelvinFix[model]['histNat'] = [-273.15, -273.15, -273.15, -273.15]
KelvinFix[model]['rcp85'] = [-273.15, -273.15, -273.15, -273.15]
KelvinFix[model]['rcp45'] = [-273.15, -273.15, -273.15, -273.15]

# CanESM2
model = 'CanESM2'
header[model]['hist'] = [pathroot+'data/CMIP5/CCCma/CanESM2/historical/day/ocean/day/r1i1p1/tos/1/']
header[model]['histNat'] = [pathroot+'data/CMIP5/CCCma/CanESM2/historicalNat/day/ocean/day/r1i1p1/tos/1/',
                            pathroot+'data/CMIP5/CCCma/CanESM2/historicalNat/day/ocean/day/r3i1p1/tos/1/',
                            pathroot+'data/CMIP5/CCCma/CanESM2/historicalNat/day/ocean/day/r5i1p1/tos/1/']
header[model]['rcp85'] = [pathroot+'data/CMIP5/CCCma/CanESM2/rcp85/day/ocean/day/r1i1p1/tos/1/',
                          pathroot+'data/CMIP5/CCCma/CanESM2/rcp85/day/ocean/day/r2i1p1/tos/1/',
                          pathroot+'data/CMIP5/CCCma/CanESM2/rcp85/day/ocean/day/r3i1p1/tos/1/',
                          pathroot+'data/CMIP5/CCCma/CanESM2/rcp85/day/ocean/day/r4i1p1/tos/1/',
                          pathroot+'data/CMIP5/CCCma/CanESM2/rcp85/day/ocean/day/r5i1p1/tos/1/']
header[model]['rcp45'] = [pathroot+'data/CMIP5/CCCma/CanESM2/rcp45/day/ocean/day/r1i1p1/',
                          pathroot+'data/CMIP5/CCCma/CanESM2/rcp45/day/ocean/day/r2i1p1/',
                          pathroot+'data/CMIP5/CCCma/CanESM2/rcp45/day/ocean/day/r3i1p1/',
                          pathroot+'data/CMIP5/CCCma/CanESM2/rcp45/day/ocean/day/r4i1p1/',
                          pathroot+'data/CMIP5/CCCma/CanESM2/rcp45/day/ocean/day/r5i1p1/']
# Time and date vectors
Ly[model] = 365
t[model]['hist'], tmp, T[model]['hist'], year[model]['hist'], month[model]['hist'], day[model]['hist'], doy[model]['hist'] = ecj.timevector([1850,1,1], [2005,12,31])
t[model]['histNat'], tmp, T[model]['histNat'], year[model]['histNat'], month[model]['histNat'], day[model]['histNat'], doy[model]['histNat'] = ecj.timevector([1850,1,1], [2012,12,31])
t[model]['rcp45'], tmp, T[model]['rcp45'], year[model]['rcp45'], month[model]['rcp45'], day[model]['rcp45'], doy[model]['rcp45'] = ecj.timevector([2006,1,1], [2300,12,31])
t[model]['rcp85'], tmp, T[model]['rcp85'], year[model]['rcp85'], month[model]['rcp85'], day[model]['rcp85'], doy[model]['rcp85'] = ecj.timevector([2006,1,1], [2100,12,31])
for exp in experiments:
    feb29s = ~((month[model][exp]==2) * (day[model][exp]==29))
    t[model][exp] = t[model][exp][feb29s]
    year[model][exp] = year[model][exp][feb29s]
    month[model][exp] = month[model][exp][feb29s]
    day[model][exp] = day[model][exp][feb29s]
    doy[model][exp] = doy[model][exp][feb29s]
    T[model][exp] = len(t[model][exp])
    years[model][exp] = np.unique(year[model][exp])
    Ty[model][exp] = len(years[model][exp])
# Kelvin -> Celcius
KelvinFix[model]['hist'] = [-273.15]
KelvinFix[model]['histNat'] = [-273.15, -273.15, -273.15]
KelvinFix[model]['rcp85'] = [-273.15, -273.15, -273.15, -273.15, -273.15]

# CSIRO-Mk3-6-0
model = 'CSIRO-Mk3-6-0'
header[model]['hist'] = [pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historical/day/ocean/day/r1i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historical/day/ocean/day/r2i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historical/day/ocean/day/r3i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historical/day/ocean/day/r4i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historical/day/ocean/day/r5i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historical/day/ocean/day/r6i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historical/day/ocean/day/r7i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historical/day/ocean/day/r8i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historical/day/ocean/day/r9i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historical/day/ocean/day/r10i1p1/v20111222/tos/']
header[model]['histNat'] = [pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historicalNat/day/ocean/day/r1i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historicalNat/day/ocean/day/r2i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historicalNat/day/ocean/day/r3i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historicalNat/day/ocean/day/r4i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historicalNat/day/ocean/day/r5i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historicalNat/day/ocean/day/r6i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historicalNat/day/ocean/day/r7i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historicalNat/day/ocean/day/r8i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historicalNat/day/ocean/day/r9i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historicalNat/day/ocean/day/r10i1p1/v20111222/tos/']
header[model]['rcp85'] = [pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/rcp85/day/ocean/day/r1i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/rcp85/day/ocean/day/r2i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/rcp85/day/ocean/day/r3i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/rcp85/day/ocean/day/r4i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/rcp85/day/ocean/day/r5i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/rcp85/day/ocean/day/r6i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/rcp85/day/ocean/day/r7i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/rcp85/day/ocean/day/r8i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/rcp85/day/ocean/day/r9i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/rcp85/day/ocean/day/r10i1p1/v20111222/tos/']
header[model]['rcp45'] = [pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/rcp45/day/ocean/day/r1i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/rcp45/day/ocean/day/r2i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/rcp45/day/ocean/day/r3i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/rcp45/day/ocean/day/r4i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/rcp45/day/ocean/day/r5i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/rcp45/day/ocean/day/r6i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/rcp45/day/ocean/day/r7i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/rcp45/day/ocean/day/r8i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/rcp45/day/ocean/day/r9i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/rcp45/day/ocean/day/r10i1p1/v20111222/tos/']
# Time and date vectors
Ly[model] = 365
t[model]['hist'], tmp, T[model]['hist'], year[model]['hist'], month[model]['hist'], day[model]['hist'], doy[model]['hist'] = ecj.timevector([1850,1,1], [2005,12,31])
t[model]['histNat'], tmp, T[model]['histNat'], year[model]['histNat'], month[model]['histNat'], day[model]['histNat'], doy[model]['histNat'] = ecj.timevector([1850,1,1], [2012,12,31])
t[model]['rcp85'], tmp, T[model]['rcp85'], year[model]['rcp85'], month[model]['rcp85'], day[model]['rcp85'], doy[model]['rcp85'] = ecj.timevector([2006,1,1], [2300,12,31])
t[model]['rcp45'], tmp, T[model]['rcp45'], year[model]['rcp45'], month[model]['rcp45'], day[model]['rcp45'], doy[model]['rcp45'] = ecj.timevector([2006,1,1], [2300,12,31])
for exp in experiments:
    feb29s = ~((month[model][exp]==2) * (day[model][exp]==29))
    t[model][exp] = t[model][exp][feb29s]
    year[model][exp] = year[model][exp][feb29s]
    month[model][exp] = month[model][exp][feb29s]
    day[model][exp] = day[model][exp][feb29s]
    doy[model][exp] = doy[model][exp][feb29s]
    T[model][exp] = len(t[model][exp])
    years[model][exp] = np.unique(year[model][exp])
    Ty[model][exp] = len(years[model][exp])
# Kelvin -> Celcius
KelvinFix[model]['hist'] = [-273.15, -273.15, -273.15, -273.15, -273.15, -273.15, -273.15, -273.15, -273.15, -273.15]
KelvinFix[model]['histNat'] = [-273.15, -273.15, -273.15, -273.15, -273.15, -273.15, -273.15, -273.15, -273.15, -273.15]
KelvinFix[model]['rcp85'] = [-273.15, -273.15, -273.15, -273.15, -273.15, -273.15, -273.15, -273.15, -273.15, -273.15]
KelvinFix[model]['rcp45'] = [-273.15, -273.15, -273.15, -273.15, -273.15, -273.15, -273.15, -273.15, -273.15, -273.15]

# CNRM-CM5
model = 'CNRM-CM5'
header[model]['hist'] = [pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historical/day/ocean/day/r1i1p1/v20111213/tos/',
                         pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historical/day/ocean/day/r2i1p1/v20111114/tos/',
                         pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historical/day/ocean/day/r3i1p1/v20111114/tos/',
                         pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historical/day/ocean/day/r4i1p1/v20111115/tos/',
                         pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historical/day/ocean/day/r5i1p1/v20111115/tos/',
                         pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historical/day/ocean/day/r6i1p1/v20111115/tos/',
                         pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historical/day/ocean/day/r7i1p1/v20111125/tos/',
                         pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historical/day/ocean/day/r8i1p1/v20111117/tos/',
                         pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historical/day/ocean/day/r9i1p1/v20111115/tos/',
                         pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historical/day/ocean/day/r10i1p1/v20111122/tos/']
header[model]['histNat'] = [pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historicalNat/day/ocean/day/r1i1p1/v20120702/tos/',
                            pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historicalNat/day/ocean/day/r2i1p1/v20111213/tos/',
                            pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historicalNat/day/ocean/day/r3i1p1/v20111213/tos/',
                            pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historicalNat/day/ocean/day/r4i1p1/v20111213/tos/',
                            pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historicalNat/day/ocean/day/r5i1p1/v20111213/tos/',
                            pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historicalNat/day/ocean/day/r8i1p1/v20111213/tos/']
header[model]['rcp85'] = [pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/rcp85/day/ocean/day/r1i1p1/v20111025/tos/',
                          pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/rcp85/day/ocean/day/r2i1p1/v20120614/tos/',
                          pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/rcp85/day/ocean/day/r4i1p1/v20111025/tos/',
                          pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/rcp85/day/ocean/day/r6i1p1/v20111025/tos/',
                          pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/rcp85/day/ocean/day/r10i1p1/v20120120/tos/']
header[model]['rcp45'] = [pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/rcp45/day/ocean/day/r1i1p1/']
# Time and date vectors
Ly[model] = 365.25
t[model]['hist'], tmp, T[model]['hist'], year[model]['hist'], month[model]['hist'], day[model]['hist'], doy[model]['hist'] = ecj.timevector([1850,1,1], [2005,12,31])
t[model]['histNat'], tmp, T[model]['histNat'], year[model]['histNat'], month[model]['histNat'], day[model]['histNat'], doy[model]['histNat'] = ecj.timevector([1850,1,1], [2012,12,31])
t[model]['rcp45'], tmp, T[model]['rcp45'], year[model]['rcp45'], month[model]['rcp45'], day[model]['rcp45'], doy[model]['rcp45'] = ecj.timevector([2006,1,1], [2300,12,31])
t[model]['rcp85'], tmp, T[model]['rcp85'], year[model]['rcp85'], month[model]['rcp85'], day[model]['rcp85'], doy[model]['rcp85'] = ecj.timevector([2006,1,1], [2300,12,31])
for exp in experiments:
    years[model][exp] = np.unique(year[model][exp])
    Ty[model][exp] = len(years[model][exp])
# Kelvin -> Celcius (only not required for ACCESS hist ens > 0 and non RCP)
# ALL OTHER MODELS WILL NEED A -273.15 OFFSET FOR K->C
KelvinFix[model]['hist'] = [-273.15, -273.15, -273.15, -273.15, -273.15, -273.15, -273.15, -273.15, -273.15, -273.15]
KelvinFix[model]['histNat'] = [-273.15, -273.15, -273.15, -273.15, -273.15, -273.15]
KelvinFix[model]['rcp85'] = [-273.15, -273.15, -273.15, -273.15, -273.15]
KelvinFix[model]['rcp45'] = [-273.15]

#
# Calculate number of runs, time experiments, ensemble members, etc...
#

NEXPS = len(experiments)

NENS = {}
for model in models:
    NENS[model] = {}
    for exp in experiments:
        NENS[model][exp] = len(header[model][exp])


modelPaths = {}
modelPaths['ACCESS1-3'] = 'CSIRO-BOM/ACCESS1-3/'
modelPaths['IPSL-CM5A-LR'] = 'IPSL/IPSL-CM5A-LR/'
modelPaths['IPSL-CM5A-MR'] = 'IPSL/IPSL-CM5A-MR/'
modelPaths['HadGEM2-ES'] = 'MOHC/HadGEM2-ES/'
modelPaths['CanESM2'] = 'CCCma/CanESM2/'
modelPaths['CSIRO-Mk3-6-0'] = 'CSIRO-QCCCE/CSIRO-Mk3-6-0/'
modelPaths['CNRM-CM5'] = 'CNRM-CERFACS/CNRM-CM5/'

# Save meta-data
np.savez(pathroot + 'data/CMIP5/meta/timeSeries_metaData.npz', experiments=experiments, Ly=Ly, t=t, years=years, Ty=Ty, KelvinFix=KelvinFix, NEXPS=NEXPS, NENS=NENS, modelPaths=modelPaths)

#
# Load data from each run/ensemble member
#

# LOOP OVER MODELS
X0 = 18 # Size of lon-blocks to load separately

# LOOP OVER MODELS
for model in ['CNRM-CM5']:
    # Setup lats and lons
    files = []
    os.chdir(header[model]['hist'][0]) # Assume all files have same lat/lon
    file0 = []
    for file in os.listdir('.'):
        if file.endswith('.nc'):
            file0.append(header[model]['hist'][0] + file)
            break
    fileobj = Dataset(file0[0], mode='r')
    lon[model] = fileobj.variables['lon'][:].astype(float)
    lat[model] = fileobj.variables['lat'][:].astype(float)
    fill_value = fileobj.variables['tos']._FillValue.astype(float)
    landmask = fileobj.variables['tos'][0,:,:].data
    fileobj.close()
    if (model == 'HadGEM2-ES') + (model == 'CanESM2') + (model == 'CSIRO-Mk3-6-0'):
        Y[model] = len(lat[model])
        X[model] = len(lon[model])
    else:
        Y[model], X[model] = lon[model].shape
    NBLOCKS = np.ceil(1.*X[model]/X0).astype(int) # Number of lon blocks to loop over
    # Generate landmask
    landmask[landmask==fill_value] = np.nan
    landmask[~np.isnan(landmask)] = 1.
    # Load data in lon-blocks
    for n in range(NBLOCKS):
        # Go through each experiment
        for exp in ['histNat', 'hist', 'rcp45', 'rcp85']: #experiments:
            for ens in range(NENS[model][exp]):
                print 'Model:', model, ', Experiment:', exp, ', Ensemble member:', ens+1, '/', NENS[model][exp], ', Block:', n+1, '/', NBLOCKS
                # Files
                files = []
                os.chdir(header[model][exp][ens])
                for file in os.listdir('.'):
                    if file.endswith('.nc'):
                        files.append(header[model][exp][ens] + file)
                N_files = len(files)
                file0 = files[0]
                fileobj = Dataset(file0, mode='r')
                fill_value = fileobj.variables['tos']._FillValue.astype(float)
                fileobj.close()
                #   load SST
                sst_ts = np.nan*np.ones((T[model][exp], Y[model], X0))
                t_ts = np.nan*np.ones((T[model][exp],))
                t0 = 0
                for file in files:
                    fileobj = Dataset(file, mode='r')
                    t_file = fileobj.variables['time'][:]
                    T_file = len(t_file)
                    if (n+1)*X0 > X[model]: # Last block is a partial-block
                        sst_file = fileobj.variables['tos'][:,:,n*X0:]
                    else: # Blocks fit integrally into X
                        sst_file = fileobj.variables['tos'][:,:,n*X0:(n+1)*X0]
                    sst_file = sst_file.data
                    sst_file[sst_file==fill_value] = np.nan
                    sst_ts[t0:(t0+T_file),:,:(sst_file.shape[2])] = sst_file
                    t_ts[range(t0, t0 + T_file)] = t_file
                    t0 = t0 + T_file
                    fileobj.close()
                #   Sort based on time-vector and add offset (if required)
                tt = np.argsort(t_ts)
                sst_ts = sst_ts[tt,:,:] + KelvinFix[model][exp][ens]
                # loop over i, j
                if (n+1)*X0 > X[model]: # Last block is a partial-block
                    X00 = X[model] - n*X0
                else:
                    X00 = X0 + 0
                ii = range(n*X0, n*X0 + X00) # Actual i-index of data
                for i in range(X00):
                    outfile = pathroot + 'data/CMIP5/' + modelPaths[model] + 'timeseries/sst_' + exp + '_ens_' + str(ens).zfill(2) + '_lon_' + str(ii[i]+1).zfill(4) + '.mat'
                    io.savemat(outfile, {'t': t[model][exp], 'sst': sst_ts[:,:,i], 'lon': lon[model], 'lat': lat[model], 'X': X[model], 'Y': Y[model], 'landmask': landmask})


