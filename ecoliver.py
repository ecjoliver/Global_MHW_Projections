import numpy as np
import scipy as sp
from scipy import stats
from scipy import linalg
import scipy.ndimage as ndimage
from datetime import date

def find_nearest(array, value):
    idx=(np.abs(array-value)).argmin()
    return array[idx], idx

def findxy(x, y, loc):
    '''
    Find (i,j) coordinates of loc = (x0,y0) in 2D irregularly spaced
    coordinate matrices (2D numpy arrays) x and y.
    '''
    # Dimensions
    Y, X = x.shape
    # Make matrix of indices
    iijj = np.meshgrid(range(X), range(Y))
    ii = iijj[0].flatten()
    jj = iijj[1].flatten()
    # Calculate distance-squared
    dist2 = (x-loc[0])**2 + (y-loc[1])**2
    k = np.argmin(dist2)
    return ii[k], jj[k]

def nanmean(array, axis=None):
    return np.mean(np.ma.masked_array(array, np.isnan(array)), axis)

def nanvar(array, axis=None):
    return np.var(np.ma.masked_array(array, np.isnan(array)), axis)

def nanskew(array, axis=None):
    # only woks for 1D data
    return stats.skew(array[np.logical_not(np.isnan(array))])

def nanmax(array, axis=None):
    maxes = np.max(np.ma.masked_array(array, np.isnan(array)), axis)
    data = maxes.data
    mask = maxes.mask
    return data[np.logical_not(mask)]

def nonans(array):
    '''
    Return input array [1D numpy array] with
    all nan values removed
    '''
    return array[~np.isnan(array)]

def nozeros(array):
    '''
    Return input array [1D numpy array] with
    all zeros removed
    '''
    return array[~(array==0)]


def latlon2km(lon1, lat1, lon2, lat2):
    EARTH_RADIUS = 6378.1
    c = np.sin(np.radians(lat1)) * np.sin(np.radians(lat2)) + np.cos(np.radians(lon1-lon2)) * np.cos(np.radians(lat1)) * np.cos(np.radians(lat2))
    d = EARTH_RADIUS * np.arccos(c)
    return d

def latlonArea(lon1, lat1, lon2, lat2):
    '''
    Surface area (in km2^) of a lat/lon "rectangle" included between specified longitudes and latitudes
    '''
    EARTH_RADIUS = 6378.1
    return (EARTH_RADIUS**2) * np.abs(lon1*np.pi/180. - lon2*np.pi/180.) * np.abs(np.sin(lat1*np.pi/180.) - np.sin(lat2*np.pi/180.))

def dxdy(lon, lat):
    '''
    Takes M+1 length lat and N+1 length lon vectors
    and returns MxN 2D arrays of distances across cells
    in x and y directions
    '''
    X = len(lon)-1
    Y = len(lat)-1
    dx = np.zeros((Y,X))
    dy = np.zeros((Y,X))
    for j in range(dx.shape[0]):
        for i in range(dx.shape[1]):
            dx[j,i] = 1e3 * latlon2km(lon[i+1], lat[j], lon[i], lat[j])
            dy[j,i] = 1e3 * latlon2km(lon[i], lat[j+1], lon[i], lat[j])
    return dx, dy

def gradient(field, dx, dy):
    '''
    Performs the gradient of input field given
    dx and dy fields (in metres)
    '''
    field_y, field_x = np.gradient(field)
    field_x = field_x / dx
    field_y = field_y / dy
    return field_x, field_y

def acf(x):
    result = np.correlate(x, x, mode = 'full')
    maxcorr = np.argmax(result)
    result = result / result[maxcorr]     # <=== normalization
    return result[result.size/2:]

def ccf(x, y):
    '''
    Cross-Correlation Function
    +ve lags mean x leads y
    -ve lags mean x lags  y
    '''

    x = (x-np.mean(x))/(np.std(x)*len(x))
    y = (y-np.mean(y))/np.std(y)
    ccf = np.correlate(x, y, mode='full')
    lags = np.arange(len(ccf)) - (len(x)-1)
    return lags, ccf

def ttest_serialcorr(x, y):
    '''
    Calculates the t-test for the means of two samples under an assumption of serial
    correlation, following the technique of Zwiers and von Storch (Journal of Climate, 1995)
    '''
    # Valid (non-Nan) data, and return NaN if insufficient valid data
    validx = ~np.isnan(x)
    validy = ~np.isnan(y)
    if (validx.sum() <= 1) + (validy.sum() <= 1):
        return np.nan, np.nan
    else:
        # Sample lengths
        nx = len(x[validx])
        ny = len(y[validy])
        # Autocorrelation Function (pad NaN values for an approximation)
        rhox = acf(pad(x - np.nanmean(x)))
        rhoy = acf(pad(y - np.nanmean(y)))
        # Equivalent sample lengths
        nx = nx / (1 + ((1-np.arange(1, int(nx))/nx)*rhox[validx][:-1]).sum())
        ny = ny / (1 + ((1-np.arange(1, int(ny))/ny)*rhoy[validy][:-1]).sum())
        #if (nx < 30) or (ny < 30):
        #    print 'Effective sample size(s) are less than 30: distribution of t statistics will deviate significantly from the t-distribution'
        # Sample standard deviations
        sx = np.sqrt(x[validx].var())
        sy = np.sqrt(y[validy].var())
        s = np.sqrt(sx**2/nx + sy**2/ny)
        # t-statistic
        t = (np.nanmean(x) - np.nanmean(y))/s
        # Degrees of freedom
        df = (sx**2/nx + sy**2/ny)**2 / ((sx**2/nx)**2/(nx-1) + (sy**2/ny)**2/(ny-1))
        # p-value
        p = 1 - stats.t.cdf(t, df)
        return t, p

def pad(data, maxPadLength=False):
    '''

    Linearly interpolate over missing data (NaNs) in a time series.

    Inputs:

      data           Time series [1D numpy array]
      maxPadLength   Specifies the maximum length over which to interpolate,
                     i.e., any consecutive blocks of NaNs with length greater
                     than maxPadLength will be left as NaN. Set as an integer.
                     maxPadLength=False (default) interpolates over all NaNs.

    Written by Eric Oliver, Institue for Marine and Antarctic Studies, University of Tasmania, Jun 2015

    '''
    data_padded = data.copy()
    if len(data) == np.isnan(data).sum():
        return np.nan*data_padded
    else:
        bad_indexes = np.isnan(data)
        good_indexes = np.logical_not(bad_indexes)
        good_data = data[good_indexes]
        interpolated = np.interp(bad_indexes.nonzero()[0], good_indexes.nonzero()[0], good_data)
        data_padded[bad_indexes] = interpolated
        if maxPadLength:
            blocks, n_blocks = ndimage.label(np.isnan(data))
            for bl in range(1, n_blocks+1):
                if (blocks==bl).sum() > maxPadLength:
                    data_padded[blocks==bl] = np.nan
        return data_padded

def runavg_periodic(ts, w):
    '''
    Perform running average of ts (1D numpy array) using uniform window
    of width w (w must be odd). Assumes periodicity of ts.
    '''
    N = len(ts)
    ts = np.append(ts, np.append(ts, ts))
    ts_smooth = np.convolve(ts, np.ones(w)/w, mode='same')
    ts = ts_smooth[N:2*N]
    return ts

def runavg(ts, w, mode='same'):
    '''
    Perform running average of ts (1D numpy array) using uniform window
    of width w (w must be odd). Pads with NaNs outside of valid range.
    Option 'mode' specifies if output should be defined over
    '''
    if mode == 'same':
        ts_smooth = np.convolve(ts, np.ones(w)/w, mode=mode)
    elif mode == 'valid':
        ts_smooth = np.append(np.append(np.nan*np.ones((w-1)/2), np.convolve(ts, np.ones(w)/w, mode=mode)), np.nan*np.ones((w-1)/2))
    return ts_smooth

def timevector(date_start, date_end):
    '''
    Generated daily time vector, along with year, month, day, day-of-year,
    and full date information, given start and and date. Format is a 3-element
    list so that a start date of 3 May 2005 is specified date_start = [2005,5,3]
    Note that day-of year (doy) is [0 to 59, 61 to 366] for non-leap years and [0 to 366]
    for leap years.
    returns: t, dates, T, year, month, day, doy
    '''
    # Time vector
    t = np.arange(date(date_start[0],date_start[1],date_start[2]).toordinal(),date(date_end[0],date_end[1],date_end[2]).toordinal()+1)
    T = len(t)
    # Date list
    dates = [date.fromordinal(tt.astype(int)) for tt in t]
    # Vectors for year, month, day-of-month
    year = np.zeros((T))
    month = np.zeros((T))
    day = np.zeros((T))
    for tt in range(T):
        year[tt] = date.fromordinal(t[tt]).year
        month[tt] = date.fromordinal(t[tt]).month
        day[tt] = date.fromordinal(t[tt]).day
    year = year.astype(int)
    month = month.astype(int)
    day = day.astype(int)
    # Leap-year baseline for defining day-of-year values
    year_leapYear = 2012 # This year was a leap-year and therefore doy in range of 1 to 366
    t_leapYear = np.arange(date(year_leapYear, 1, 1).toordinal(),date(year_leapYear, 12, 31).toordinal()+1)
    dates_leapYear = [date.fromordinal(tt.astype(int)) for tt in t_leapYear]
    month_leapYear = np.zeros((len(t_leapYear)))
    day_leapYear = np.zeros((len(t_leapYear)))
    doy_leapYear = np.zeros((len(t_leapYear)))
    for tt in range(len(t_leapYear)):
        month_leapYear[tt] = date.fromordinal(t_leapYear[tt]).month
        day_leapYear[tt] = date.fromordinal(t_leapYear[tt]).day
        doy_leapYear[tt] = t_leapYear[tt] - date(date.fromordinal(t_leapYear[tt]).year,1,1).toordinal() + 1
    # Calculate day-of-year values
    doy = np.zeros((T))
    for tt in range(T):
        doy[tt] = doy_leapYear[(month_leapYear == month[tt]) * (day_leapYear == day[tt])]
    doy = doy.astype(int)

    return t, dates, T, year, month, day, doy

def spatial_filter(field, res, cut_lon, cut_lat):
    '''
    Performs a spatial filter, removing all features with
    wavelenth scales larger than cut_lon in longitude and
    cut_lat in latitude from field. Field has spatial
    resolution of res and land identified by np.nan's
    '''

    field_filt = np.zeros(field.shape)

    # see Chelton et al, Prog. Ocean., 2011 for explanation of factor of 1/5
    sig_lon = (cut_lon/5.) / res
    sig_lat = (cut_lat/5.) / res

    land = np.isnan(field)
    field[land] = nanmean(field)

    field_filt = field - ndimage.gaussian_filter(field, [sig_lat, sig_lon])

    field_filt[land] = np.nan

    return field_filt

def trend(x, y, alpha=0.05):
    '''
    Calculates the trend of y given the linear
    independent variable x. Outputs the mean,
    trend, and alpha-level (e.g., 0.05 for 95%)
    confidence limit on the trend.
    returns mean, trend, dtrend_95
    '''
    valid = ~np.isnan(y)
    if valid.sum() <= 1:
        return np.nan, np.nan, np.nan
    else:
        X = np.array([np.ones(len(x)), x-x.mean()])
        beta = linalg.lstsq(X[:,valid].T, y[valid])[0]
        yhat = np.sum(beta*X.T, axis=1)
        t_stat = stats.t.isf(alpha/2, len(x[valid])-2)
        s = np.sqrt(np.sum((y[valid] - yhat[valid])**2) / (len(x[valid])-2))
        Sxx = np.sum(X[1,valid]**2) - (np.sum(X[1,valid])**2)/len(x[valid]) # np.var(X, axis=1)[1]
        return beta[0], beta[1], t_stat * s / np.sqrt(Sxx)

def trend_TheilSen(x, y, alpha=0.05):
    '''
    Calculates the trend of y given the linear
    independent variable x. Outputs the mean,
    trend, and alpha-level (e.g., 0.05 for 95%)
    confidence limit on the trend. Estimate of
    the trend uses a Theil-Sen estimator.
    returns mean, trend, dtrend_95
    '''
    # Construct matrix of predictors, first column is all ones to estimate the mean,
    # second column is the time vector, equal to zero at mid-point.
    X = x-x.mean()
#
    # Predictand (MHW property of interest)
    valid = ~np.isnan(y) # non-NaN indices
#
    # Perform linear regression over valid indices
    if np.sum(~np.isnan(y)) > 0: # If at least one non-NaN value
        slope, y0, beta_lr, beta_up = stats.mstats.theilslopes(y[valid], X[valid], alpha=1-alpha)
        beta = np.array([y0, slope])
    else:
        beta_lr, beta_up = [np.nan, np.nan]
        beta = [np.nan, np.nan]
#
    return beta[0], beta[1], [beta_lr, beta_up]


def acf(x):
    result = np.correlate(x, x, mode = 'full')
    maxcorr = np.argmax(result)
    result = result / result[maxcorr]     # <=== normalization
    return result[result.size/2:]

def ttest_unequalvar(x, y):
    '''
    Calculates the t-test for the means of two samples under an assumption of no serial
    correlation (but different variances), following the technique of Zwiers and von Storch
    (Journal of Climate, 1995)
    '''
    # Sample lengths
    nx = len(x)
    ny = len(y)
    # Sample standard deviations
    sx = np.sqrt(x.var())
    sy = np.sqrt(y.var())
    s = np.sqrt(sx**2/nx + sy**2/ny)
    # t-statistic
    t = np.abs(x.mean() - y.mean())/s
    # Degrees of freedom
    df = (sx**2/nx + sy**2/ny)**2 / ((sx**2/nx)**2/(nx-1) + (sy**2/ny)**2/(ny-1))
    # p-value
    p = 1 - stats.t.cdf(t, df)
    return t, p

def pattern_correlation(x1, x2, centred=True):
    '''
    Calculates the pattern correlation of x1 and x2.
    Assumes and x1 and x2 are 2D numpy arrays. Can
    handle missing values, even if missing values are
    distributed differently in x1 and x2. By default
    calculated the centred pattern correlation (centred
    =True) in which the spatial means of x1 and x2 are
    removed prior to calculation. Can calculated uncentred
    pattern correlation (centred=False) in which these
    means are not removed.
    .
    Written by Eric Oliver, IMAS/UTAS, Nov 2015
    '''
    # Flatten 2D arrays and find shared valid (non-nan) indices
    X1 = x1.flatten()
    X2 = x2.flatten()
    valid = ~(np.isnan(X1) + np.isnan(X2))
    # Create Nx2 array of valid data
    X = np.zeros((valid.sum(), 2))
    X[:,0] = X1[valid]
    X[:,1] = X2[valid]
    # Centre data if desired
    if centred:
        X[:,0] = X[:,0] - np.mean(X[:,0])
        X[:,1] = X[:,1] - np.mean(X[:,1])
    #
    # Calculate pattern correlation
    pcorr = np.corrcoef(X.T)[0,1]
    return pcorr

def polyArea(x, y):
    '''
    Area of a simple polygon, defined by vertices (x,y)
    in the plane. Assumes x and y and numpy arrays of the
    same length. Algorithm is based on the Shoelace Formula
    https://en.wikipedia.org/wiki/Shoelace_formula
    http://stackoverflow.com/questions/24467972/calculate-area-of-polygon-given-x-y-coordinates
    '''
    return 0.5*np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))

def point_inside_polygon(x, y, poly):
    '''
    Determine if a point is inside a given polygon or not
    Polygon is a list of (x,y) pairs.
    http://www.ariel.com.au/a/python-point-int-poly.html
    '''
    n = len(poly)
    inside =False
    #
    p1x, p1y = poly[0]
    for i in range(n+1):
        p2x, p2y = poly[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xinters = (y - p1y)*(p2x - p1x)/(p2y - p1y)+p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x, p1y = p2x, p2y
    #
    return inside

def pAgree(k, n, p=0.5):
    '''
    Returns the significance level of agreement of k results across n datasets, based on a binomial distribution.
    Assumes a "fair coin toss" (p=0.5) but this can be optionally changed.
    For example, if k=9 out of n=10 datasets show a trend of the same sign, this returns the significance level for
     such a result, assuming independence of the datasets (ha!) and a probability of a trend with that sign of p=0.5.
    
    '''
    return stats.binom.pmf(np.arange(k,n+1), n, p).sum()

def weighted_quantile(values, quantiles, sample_weight=None, values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of initial array
    :param old_style: if True, will correct output to be consistent with numpy.percentile.
    :return: numpy.array with computed quantiles.

    Authored by Alleo, Apr 16 '15 at 14:22
    Source: https://stackoverflow.com/questions/21844024/weighted-percentile-using-numpy
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), 'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)

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

