#!/usr/bin/python


import sys
import numpy as np
import scipy.optimize, scipy.ndimage, scipy.stats, scipy.odr
import sklearn.preprocessing
import statsmodels.api as sm

import grid_and_interpolation as gi

######################################################################
######################################################################
######################################################################
#####                          #######################################
##### STUFF FOR STATISTICS     #######################################
#####                          #######################################
######################################################################
######################################################################
######################################################################


######################################################################
######################################################################

# 
# ERROR HERE: OLS regression does not allow for variations in 
#             independent variable!
#
# def get_fit_range(func, x, y,
#                   perc = [16, 50, 84], 
#                   xerr = 0, 
#                   yerr = 0, 
#                   Nsample = 100):

    
#     ''' 
#     Calculates percentile values of possible fitting curve using a bootstrap
#     approach.

#     (x,y) are interpreted as mean values and error bars as standard deviations of
#     2d Gaussian distribution per point and the fit is repeated for random draws
#     out of these. After Nsampe fits, the range of possible fits is analyzed for their
#     percentiles.


#     USAGE:
#     =====
#     pfit = get_fit_range(func, x, y)


#     INPUT:
#     =====
#     func: fitting function func(x, *args), args: parameters of the function
#     x: independent values
#     y: dependent values


#     OUPUT:
#     =====
#     pfit: percentiles of possible fit curves obtained from bootstrapping

#     '''
    

#     # get dimension ..................................................
#     Nsize = x.shape[0]



#     # init output ....................................................
#     xout = np.linspace(x.min(), x.max(), Nsize)
#     yfit = []


#     # loop over random realizations ----------------------------------
#     for n in range(Nsample):
        
#         # randomize values ...........................................
#         xnew = x + np.abs(xerr) * np.random.randn(Nsize)
#         ynew = y + np.abs(yerr) * np.random.randn(Nsize)


#         # do the fit with scipy ......................................
#         popt, pcov = scipy.optimize.curve_fit(func, xnew, ynew)

        
#         # collect fitted values ......................................
#         yfit.append( func(xout, *popt) ) 


#     # make numpy array ...............................................
#     yfit = np.row_stack(yfit)
#     # ================================================================


#     # return percentiles of fitted functions
#     return xout, np.percentile(yfit, perc, axis = 0)


######################################################################
######################################################################

def symmetric_linear_fit(x, y, 
                         reg_type = 'robust',
                         xerr = 1.,
                         yerr = 1.,
                         add_constant = True,
                         Nsample = 1000,
                         output_all = False):

    # select either y = ax OR y = ax + b
    if add_constant:
        X = sm.add_constant(x)
        Y = sm.add_constant(y)
    else:
        X = x
        Y = y


    if reg_type == 'robust':
        # do fit minimizing y-residuals
        yfit = sm.RLM(y, X).fit()

        # do fit minimizing y-residuals
        xfit = sm.RLM(x, Y).fit()

        # get parameter covariance
        pcovy =  yfit.bcov_scaled
        pcovx =  xfit.bcov_scaled

    elif reg_type == 'weighted':

        # assume that error is standard deviation
        xvar = xerr**2
        yvar = yerr**2

        # do fit minimizing y-residuals
        yfit = sm.WLS(y, X, weights = 1./ yvar).fit()

        # do fit minimizing y-residuals
        xfit = sm.WLS(x, Y, weights = 1./ xvar).fit()

        # get parameter covariance
        pcovy =  yfit.cov_params()
        pcovx =  xfit.cov_params()
        
        
    pmy = yfit.params
    pmx = xfit.params

    
    # generate parameter set
    ypars = np.random.multivariate_normal(pmy, pcovy, Nsample)
    xpars = np.random.multivariate_normal(pmx, pcovx, Nsample)
 

    # calculate arithmetic average for linear fit curve
    a1 = ypars[:,-1]
    alpha2 = xpars[:,-1]
    a2 = np.ma.divide(1, alpha2)
    a = 0.5 * (a1 + a2)

    
    if add_constant:
        b1 = ypars[:, 0]
        beta2 = xpars[:,0]


        b2 = -beta2 * a2
        b = 0.5 * (b1 + b2)


        pm = np.array([b.mean(), a.mean()] ) 
        pcov = np.ma.cov(b, a)

    else:
        pm = np.array([a.mean()])
        pcov = np.array([a.var()]).reshape(1,1)

    if output_all:
        return pm, pcov, yfit.params, yfit.bcov_scaled, xfit.params, xfit.bcov_scaled
    else:
        return pm, pcov

    
######################################################################
######################################################################



def para_bootstrapp_func_range(func, x, pm, pcov, 
                               perc = [5, 50, 95],
                               Nsample = 1000 ):


    '''
    Calculates percentile values of a function func(params, x) where the parameter set p
    is varied as multivariate normal distribution.
    

    USAGE:
    =====
    frange =  para_bootstrapp_func_range(func, x, pm, pcov)


    INPUT:
    =====
    func: fitting function func(params, x), params: parameters of the function ! ARGS HAVE BEEN TURNED !
    x: function variables
    pm: mean value of function parameters
    pcov: covaraince matrix of function parameters


    OUPUT:
    =====
    frange: percentiles of possible function curves obtained from bootstrapping
    
    '''


    # generate random sample the parameter values ....................
    p_sample = np.random.multivariate_normal(pm, pcov, Nsample)


    yset = []

    # loop over random realizations ----------------------------------
    for n in range(Nsample):
        

        # select parameters from sample set ..........................
        p = p_sample[n]

        # do fitting .................................................
        y = func(p, x)

        # collect fitted values ......................................
        yset.append( y )


    # make numpy array ...............................................
    yset = np.row_stack(yset)
    # ================================================================


    # return percentiles for function
    return  np.percentile(yset, perc, axis = 0)


######################################################################
######################################################################



def get_fit_range(func, x, y,
                  perc = [5, 50, 95], 
                  xerr = 1, 
                  yerr = 1, 
                  p0 = None,
                  reg_type = 'ortho_dist',
                  outlier_thresh = 3.,
                  Nsample = 1000):

    
    ''' 
    Calculates percentile values of possible fitting curve using a orthogonal
    distance regression.

    (x,y) are interpreted as mean values and error bars as standard deviations of
    2d Gaussian distribution per point. The range of possible fits is analyzed for their
    percentiles.


    USAGE:
    =====
    pfit = get_fit_range(func, x, y)


    INPUT:
    =====
    func: fitting function func(args, x), args: parameters of the function ! ARGS HAVE BEEN TURNED !
    x: independent values
    y: dependent values


    OUPUT:
    =====
    pfit: percentiles of possible fit curves obtained from bootstrapping

    '''
    


    # make the (linear) fit ------------------------------------------
    if reg_type == 'ortho_dist':
        myres = odrfit_with_outlier_removal(func, x, y,
                                            xerr = xerr, 
                                            yerr = yerr, 
                                            outlier_thresh = outlier_thresh,
                                            p0 = p0)


        # average parameters
        pm = myres.beta 
        pcov = myres.cov_beta


    print
    print 'fit parameters'
    print '=============='
    for i, bm in enumerate(pm):
        print 'para %d = %f +- %f' % (i, bm, myres.sd_beta[i])
    # ================================================================



    # get fit range via parameter bootstrapping ---------------------- 

    # get dimension ..................................................
    Nsize = x.shape[0]

    # init output ....................................................
    xout = np.linspace(x.min(), x.max(), Nsize)


    # 
    yfit = para_bootstrapp_func_range(func, xout, pm, pcov, 
                                   perc = perc,
                                   Nsample = Nsample)
    # ================================================================



    return  xout, yfit


######################################################################
######################################################################

def linear_func(p, x):

    n, m = p

    return m*x + n


######################################################################
######################################################################


def linear_func_no_const(p, x):

    return p*x


######################################################################
######################################################################

def get_linear_fit_range(x, y,
                         perc = [5, 50, 95], 
                         xerr = 1, 
                         yerr = 1, 
                         reg_type = 'robust',
                         add_constant = True,
                         Nsample = 1000):

    
    ''' 
    Calculates percentile values of possible linear fitting curve using either robust
    or weight least-squares regression.

    (x,y) are interpreted as mean values and error bars as standard deviations of
    2d Gaussian distribution per point. The range of possible fits is analyzed for their
    percentiles.


    USAGE:
    =====
    yfit = get_linear_fit_range(x, y)


    INPUT:
    =====
    x: independent values
    y: dependent values


    OUPUT:
    =====
    yfit: percentiles of possible fit curves obtained from bootstrapping

    '''
    


    # make the (linear) fit ------------------------------------------
    pm, pcov =  symmetric_linear_fit(x, y, 
                                     reg_type = reg_type, 
                                     add_constant = add_constant,
                                     xerr = xerr,
                                     yerr = yerr,
                                     Nsample = Nsample)
                                     

    print
    print 'fit parameters'
    print '=============='
    for i, bm in enumerate(pm):
        print 'para %d = %f +- %f' % (i, bm, np.sqrt(np.diag(pcov))[i])
    # ================================================================



    # get fit range via parameter bootstrapping ---------------------- 

    # get dimension ..................................................
    Nsize = x.shape[0]

    # init output ....................................................
    xout = np.linspace(x.min(), x.max(), Nsize)


    # 
    if add_constant:
        f = linear_func
    else:
        f = linear_func_no_const

    yfit = para_bootstrapp_func_range(f, xout, pm, pcov, 
                                   perc = perc,
                                   Nsample = Nsample)
    # ================================================================



    return  xout, yfit


######################################################################
######################################################################


def get_outlier_from_residuals(e, thresh = 3.):

    '''
    Calculate outliers from residuals using percentile-based
    standardization.
    
    INPUT
    =====
    e: input residuals

    OUTPUT
    ======
    m: outlier mask
    '''

    # get percentiles 
    ep = np.percentile(e, [16, 50, 84])


    # standardize
    es = 2. * ( e  - ep[1] ) / ( ep[2] - ep[0] )


    return np.abs(es) > thresh
    

######################################################################
######################################################################


def odrfit_with_outlier_removal(func, x, y,
                                xerr = 0, 
                                yerr = 0, 
                                n_iter = 10,
                                outlier_thresh = 3.,
                                p0 = None):
    
    ''' 
    Calculates fitting parameter object using a orthogonal distance regression.

    (x,y) are interpreted as mean values and error bars as standard deviations of
    2d Gaussian distribution per point. The range of possible fits is analyzed for their
    percentiles.


    USAGE:
    =====
    fit_result = odrfit_with_outlier_removal(func, x, y)

    INPUT:
    =====
    func: fitting function func(args, x), args: parameters of the function ! ARGS HAVE BEEN TURNED !
    x: independent values
    y: dependent values


    OUPUT:
    =====
    fit_result: result of the ODR fit

    '''
    

    # extend xerr and yerr to array if needed
    xerr = xerr * np.ones_like(x)
    yerr = yerr * np.ones_like(x)


    # set fit model for orthogonal distance regression
    mod = scipy.odr.Model(func)

    # starting with trivial mask
    m = np.ones_like(x).astype(np.bool)
    n = 0

    while n < n_iter: 

        print 'number of iteration', n
        mydata = scipy.odr.RealData(x[m], y[m], sx = xerr[m], sy = yerr[m])
        myodr = scipy.odr.ODR(mydata, mod, beta0 = p0)
        myres = myodr.run()

        # sum of squared errors
        sqerror = myres.delta**2 + myres.eps**2

        # get outliers
        m_out =  get_outlier_from_residuals(sqerror, thresh = outlier_thresh)
        m_in = np.logical_not(m_out)


        if (m_in != True).sum() == 0:
            break
        else:
            m[m] = m_in[:]

        n += 1
        
    return myres


######################################################################
######################################################################

def correlation_bootstrap(x, y,
                            perc = [2.5, 50, 97.5], 
                            xerr = 0, 
                            yerr = 0, 
                            Nsample = 1000,
                            corr = scipy.stats.pearsonr):

    
    ''' 
    Calculates percentile values of possible correlation using a bootstrap
    approach.

    (x,y) are interpreted as mean values and error bars as standard deviations of
    2d Gaussian distribution per point and the fit is repeated for random draws
    out of these. After Nsampe fits, the range of possible fits is analyzed for their
    percentiles.


    USAGE:
    =====
    corr = correlation_bootstrap(x, y)


    INPUT:
    =====
    x: independent values
    y: dependent values


    OUPUT:
    =====
    corr: percentiles of possible correlation obtained from bootstrapping

    '''
    

    # get dimension ..................................................
    Nsize = x.shape[0]


    cvec = []
    # loop over random realizations ----------------------------------
    for n in range(Nsample):
        
        # randomize values ...........................................
        xnew = x + np.abs(xerr) * np.random.randn(Nsize)
        ynew = y + np.abs(yerr) * np.random.randn(Nsize)


        # do the fit with scipy ......................................
        c, pval = corr(xnew, ynew)

        
        # collect fitted values ......................................
        cvec.append( c )


    # make numpy array ...............................................
    cvec = np.row_stack(cvec)
    # ================================================================


    # return percentiles of fitted functions
    return np.percentile(cvec, perc, axis = 0)

######################################################################
######################################################################

def cond_perc_simple(v1, v2, x1, p = [10, 25, 50, 75, 90], sig = 0, medsize = 0):

    '''
    This is a simple routine to emulate the Rosenfeld T-Re plots.

    Input are v1 (T) and v2 (Re) and percentiles of v2 conditioned 
    on intervals of v1 given by x1 are calculated.
    '''

    # number of intervals 
    Ninterval = len(x1) - 1

    v2perc = [] 

    # collect the percentiles by selective masking
    for i in range(Ninterval):

        m = np.logical_and(v1 > x1[i], v1 < x1[i + 1])

        try:
            vp = np.percentile(v2[m], p)
        except:
            vp = np.nan * np.ones(len(p))
#        print vp

        v2perc.append( vp )

    v2perc = np.column_stack(v2perc)


    # do smoothing if wanted
    if sig != 0:
        v2p = scipy.ndimage.gaussian_filter(v2perc, (0, sig))
    elif medsize != 0:
        v2p = scipy.ndimage.median_filter(v2perc, size = (1, medsize))
    else:
        v2p = v2perc


    return v2p 


######################################################################
######################################################################

def autocorr(x, Nmax = None):

    '''
    Calculates autocorrelation function of masked array.

    INPUT
    =====
    x: masked input array
    Nmax: optional, maximum lag 

    OUTPUT
    ======
    c: autocorrelation function for lag 0 ... Nmax
    '''
    
    N = len(x)
    c = []

    if Nmax == None:
        Nmax = N - 1
        

    for n in range(Nmax):

        x1 = x[n:]
        x2 = x[:N-n]
        c.append( np.ma.corrcoef(x1, x2)[0,1] )

    c = np.ma.array(c)

    return c

######################################################################
######################################################################


def KStest(x, y, Nx_eff = None, Ny_eff = None):

    '''
    Calculates Kolmogorov Smirnov test for the distributions of two samples
    in a standard way, but allows for input of effective degrees of freedom.

    Nullhypothesis that the two samples share the same distribution is rejected
    if p-value is smaller than a threshold (typically 0.05 or 0.01).


    INPUT
    =====
    x: 1st sample array
    y: 2nd sample array
    Nmax: optional, maximum lag 

    OUTPUT
    ======
    D: KS test statistic
    pval: p-value
    '''

    # remove invalid values first ....................................
    xma = np.ma.masked_invalid(x).compressed()
    yma = np.ma.masked_invalid(y).compressed()



    # get (effective) dimensions of arrays ...........................
    Nx = xma.size
    Ny = yma.size

    if Nx_eff == None:
        Nx_eff = Nx

    if Ny_eff == None:
        Ny_eff = Ny

    # calculate weighted dimension number ............................
    N = Nx * Ny / np.float( Nx + Ny )
    Neff = Nx_eff * Ny_eff / np.float(Nx_eff + Ny_eff)


 
    # perform standard KS test .......................................
    D, pval = scipy.stats.ks_2samp(xma, yma)


    # perform weighting ..............................................    
    pval_eff = 2 *  (pval / 2)**(Neff / N) 


    return D, pval_eff

######################################################################
######################################################################


def crosscorr(x, y, Nmax = None):

    '''
    Calculates time-lagged crosscorrelation function of masked array.

    INPUT
    =====
    x: 1st masked input array, evaluated at t + tau
    y: 2nd masked input array
    Nmax: optional, maximum lag 

    OUTPUT
    ======
    c: autocorrelation function for lag 0 ... Nmax, 

    '''
    
    N = len(x)
    c = []

    if Nmax == None:
        Nmax = N - 1
        

    for n in range(Nmax):

        x1 = x[n:]
        x2 = y[:N-n]
        c.append( np.ma.corrcoef(x1, x2)[0,1] )

    c = np.ma.array(c)

    return c




######################################################################
######################################################################


def rank_transformation(f, normalize = True, gamma = 1.):

    '''
    The routine performs rank transformation of a field, meaning that 
    the ranks of the individual field values are returned.


    INPUT
    =====
    f: field
    normalize: optional, option if the rank field should be normalized to one
    gamma: optional, gamma correction factor


    OUTPUT
    ======
    f: field ranks
    '''


    # (i) make vector out of field matrix............................
    fvec = f.flatten()

    # (ii) sort the field ............................................
    fsort = np.sort( fvec )

    # (iii) calculate the ranks of the field .........................
    frank = np.searchsorted( fsort, fvec ).reshape( f.shape )


    # optional things ................................................
    if normalize:
        frank = (1. * frank ) / frank.max()

    if gamma != 1:
        frank = frank**gamma

    return frank

######################################################################
######################################################################


def fdistrib_mapping(fin, fmap, keep_range = False):

    '''
    Performs a transformation of field fin to have the same distribution
    as fmap.

    INPUT
    =====
    fin: input field
    fmap: field from which distribution is taken


    OUTPUT
    ======
    fout: transformed field fin
    '''


    # (i) make vector out of field matrix............................
    fvec = fin.flatten()
    gvec = fmap.flatten()

    # (ii) sort the field ............................................
    fsort = np.sort( fvec )
    gsort = np.sort( gvec )

    # (iii) calculate the ranks of the field .........................
    frank = np.searchsorted( fsort, fvec )

    fout = gsort[ frank ].reshape( fin.shape ) 

    if keep_range:

        # transform output field
        fout_trans = ( fout - fout.min() ) / (fout.max() - fout.min())

        # use range of input field
        fout = (fin.max() - fin.min()) * fout_trans  +  fin.min()

    return fout


######################################################################
######################################################################


def cumsum_data_fraction(h, divide_by_total_sum = True):

    '''
    The function uses iso-lines of equal density and maps the fraction
    of data enclosed by these lines onto it.

    INPUT
    =====
    h : histogram / density


    OUTPUT
    ======
    f : data fraction field
    '''


    # (i) make vector out of field matrix............................
    hvec = h.flatten()

    # (ii) sort the field ............................................
    hsort = np.sort( hvec )

    # (iii) calculate the ranks of the field .........................
    hrank = np.searchsorted( hsort, hvec )

    # (iv) cummulative sum of hsort ..................................
    hs = hsort[::-1]  # turn order
    
    if divide_by_total_sum:
        hc = 100. * hs.cumsum() / hs.sum() # get cummulative sum / in percent
    else:
        hc = hs.cumsum()
        
    hc = hc[::-1] # turn order back again

    # (v) map it onto the field structure ............................
    f = hc[ hrank ].reshape( h.shape )

    return f 


######################################################################
######################################################################


def draw_from_empirical_dist(bins, hist, Nsamp = 100, discrete_version = False):
    
    '''
    Draw random number from an empirical distribution. Implemented for 2d histograms.
    
    INPUT
    ======
    bins: list of bins edges (2dim)
    hist: histogram values (either absolute frequencies or relative)
    
    OUTPUT
    ======
    rvalues: random values (2d) that are distributes like hist
    '''
    
    
    #LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    # (i) draw discrete bin values from given distribution
    #TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

    x_bins, y_bins = bins
    
    # get bin mid points
    x_bin_midpoints = gi.lmean(x_bins)
    y_bin_midpoints = gi.lmean(y_bins)
    
    # resort the empirical distribution and calculate cdf
    cdf = np.cumsum(hist.ravel())
    cdf = cdf / np.float( cdf[-1] )
    
    # draw the a random realization of cdf values
    cvalues = np.random.rand( Nsamp )
    
    # get bin index by inversion of cdf values
    value_bins = np.searchsorted(cdf, cvalues)
    
    # and find the correspoding tuble
    ix, iy = np.unravel_index(value_bins, hist.shape)
    rvalues = [x_bin_midpoints[ix], y_bin_midpoints[iy]]
            
    
    #LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    # (ii) add random noise in the size of the grid box
    #TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    
    if not discrete_version:
        # grid box sizes
        dx = np.diff(x_bins)
        dy = np.diff(y_bins)
    
    
        # uniformly-distributed random numbers (centered around zero !)
        rx_uniform = np.random.rand( Nsamp ) - 0.5
        ry_uniform = np.random.rand( Nsamp ) - 0.5
    
        # update the random draws
        rvalues[0] += rx_uniform * dx[ix]
        rvalues[1] += ry_uniform * dy[iy]
    

    return rvalues

######################################################################
######################################################################


def draw_from_empirical_1ddist(bins, hist, Nsamp = 100):

    '''
    Draw random number from an empirical distribution. Implemented for 1d histograms.
    
    
    INPUT
    ======
    bins: list of bins edges (1dim)
    hist: histogram values (either absolute frequencies or relative)
    
    OUTPUT
    ======
    rvalues: random values (2d) that are distributes like hist
    '''


    # get bin mid points
    bin_midpoints = gi.lmean(bins)

    # resort the empirical distribution and calculate cdf
    cdf = np.cumsum(hist.ravel())
    cdf = cdf / np.float(cdf[-1])
    
    # draw the a random realization of cdf values
    cvalues = np.random.rand( Nsamp )

    # get bin index by inversion of cdf values
    value_bins = np.searchsorted(cdf, cvalues)

    # and find the correspoding values
    rvalues = bin_midpoints[value_bins]


    return rvalues
    

######################################################################
######################################################################

def normalize_field( f,
                     vmin = None, 
                     vmax = None ):

    '''
    Normalization of a field is applied.


    Parameters
    ----------
    f : numpy array, 2dim
         2dim field 

    vmin : float, optional, default = None
        lower bound of fields for rescaling

    vmax : float, optional, default = None
        upper bound of fields for rescaling


    Returns
    --------
    f_norm  : numpy array, 2dim
         2dim field 
    '''


    # get min/max if not specified
    if vmin is None:
        vmin = np.ma.masked_invalid( f ).min()

    if vmax is None:
        vmax = np.ma.masked_invalid( f ).max()


    # field normalization
    f_norm = (1. * f - vmin) / (vmax - vmin)
    f_norm = np.clip( f_norm, 0., 1.)

    return f_norm

######################################################################
######################################################################



if __name__ == '__main__':


    T = np.linspace(210, 300, 500)
    Re1 = 20* (T - 300) / (210 - 300) + 10 + 0.5*np.random.randn(500)
    Re2 = 40* (T - 300) / (210 - 300) + 10 + 0.5*np.random.randn(500)

    Re = Re1
    Re[T < 250] = Re2[T < 250]


    Tout = np.arange(230, 290)
    Tm = 0.5 * ( Tout[1:] + Tout[:-1])

    Rperc = cond_perc_simple(T, Re, Tout, medsize = 5)

    import pylab as pl
    cs = ['b', 'g', 'r', 'g', 'b']

    pl.plot(Re, T, 'o', alpha = 0.6)
    for i in range(5):
        pl.plot( Rperc[i] , Tm, c = cs[i], lw = 3)
    pl.ylim(300, 230)
