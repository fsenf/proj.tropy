#!/bin/env python

import numpy as np
import pylab as pl
import scipy.ndimage

import analysis_tools.grid_and_interpolation as gi

######################################################################
######################################################################

def hist2d_scatter( x, y, bins=200, bubbles=True, **kwargs ):

    '''
    A 2d histogram is constructed and displayed as
    scatter plot.

    USAGE
    =====
    shist2d_scatter( x, y, bins=200, bubbles=True )

    
    INPUT
    =====
    x, y : numpy arrays which define two input quantities for the histogram
    
    OUTPUT
    ======
    hxy, xs, ys: frequence of occurence and x-y grid

    '''


    if kwargs.has_key('cmap'):
        cmap = kwargs['cmap']
        del kwargs['cmap']
    else:
        cmap = pl.cm.jet

    hxy, xs, ys = np.histogram2d( x, y, bins=bins, **kwargs )

    yy,xx = np.meshgrid(ys,xs)


    hxy = hxy / hxy.max()

    ran = np.linspace(0,1,10)

    for i, r in enumerate(ran[1:]):
        m = np.ma.masked_inside(hxy, ran[i-1], r).mask
        
        pl.scatter(xx[m], yy[m], edgecolors='none',
                   vmin = 0, vmax = 1.,
                   marker='o',c= hxy[m], s=40 * hxy[m], cmap = cmap)

    return hxy, xs, ys

######################################################################
######################################################################



def ave_sig_from_hist(xe, ye, h):

    '''
    Use output from the numpy 2d histogram routine 
    to calculate y-mean and standrad deviations
    along y direction.


    INPUT
    =====
    xe: edge values of x variable
    ye: edge values of y variable
    h: occurence rates (relative or absolute)

    OUTPUT
    ======
    yave, ysig: mean and std of y weigthed by h along y-direction
    '''


#    y, x = np.meshgrid(gi.lmean(ye), gi.lmean(xe))
    y, x = np.meshgrid(ye[:-1], xe[:-1])

    # calculate weighted mean value
    yave = np.where(h.sum(axis = 1) == 0, 0, (h*y).sum(axis=1) / h.sum(axis=1))

    # calculate weighted standard deviation
    o, ym = np.meshgrid(np.ones(ye.shape[0]-1), yave)
    
    yvar = (h * (y - ym)**2).sum(axis=1) / h.sum(axis=1)
    ysig = np.sqrt(yvar)

    return yave, ysig
    

######################################################################
######################################################################

def max_from_hist(xe, ye, h):

    '''
    Use output from the numpy 2d histogram routine 
    to calculate y-positions where maximum occurences are located.


    INPUT
    =====
    xe: edge values of x variable
    ye: edge values of y variable
    h: occurence rates (relative or absolute)

    OUTPUT
    ======
    ymax: y-position where h is maximal along y-direction
    '''

    norm1, dum  = np.meshgrid(h.sum(axis=1), np.ones(h.shape[1]))
    norm1 = norm1.transpose()

    h1 = np.where(norm1 == 0, 0, h / norm1)

    
    imax = h1.argmax(axis=1)

    return ye[imax+1]
    

######################################################################
######################################################################


def percentiles_from_hist(xe, ye, h, p = [25, 50, 75], axis = 0, sig  = 0):

    '''
    Use output from the numpy 2d histogram routine 
    to calculate percentiles of the y variable
    based on relative occurence rates.


    INPUT
    =====
    xe: edge values of x variable
    ye: edge values of y variable
    h: occurence rates (relative or absolute)
    p: list of percentiles values (from 0 to 100)

    OUTPUT
    ======
    yperc: list of arrays containing the percentiles
    '''
    
    # normalization
    hn = conditioned_hist(1. * h, axis = axis)


    # get the right coordinate
    if axis == 0:
        z = xe[1:]
        hn = hn.transpose()

    elif axis == 1:
        z = gi.lmean(ye)
        z = ye[1:]

        


    # what follows assumes axis = 1 LLLLLLLLLLLLLLLLLLLLLLL


    # add edge 
    z0 = z[0] - ( z[1] - z[0] ) / 2.
    zN = z[-1] + ( z[-1] - z[-2] ) / 2.

    z = np.hstack([z0, z, zN])

    h0 = np.zeros_like(hn[:,0])
    hn = np.column_stack([h0, hn,h0])

    

    sigvec = (0, sig)
    hn = scipy.ndimage.gaussian_filter(hn, sigvec)
    
    hcum = hn.cumsum(axis = 1) * 100.


    zperc = []

    for pval in p:
        

        # percentile crossing
        hr = hcum - pval

        # this is the point closest to crossing
        imax = np.abs(hr).argmin(axis = 1)


        # neighboring indices
        iminus = np.where(imax - 1 < 0, 0, imax - 1)
        iplus = np.where(imax + 1 < hn.shape[1], imax + 1, imax)


        # indexing stuff!!!!
        
        # for axis = 1
        i0 = range(imax.shape[0])

        hmask = hr[i0, imax] < 0
        
        # lower bound
        dh0 = np.where(hmask,  hr[i0, imax],  hr[i0, iminus])
        z0 = np.where(hmask,  z[imax], z[iminus])

        # upper bound
        dh1 = np.where(hmask,  hr[i0, iplus],  hr[i0, imax])
        z1 = np.where(hmask,  z[iplus], z[imax])
        
        zp = (dh1 * z0 - dh0 * z1) / (dh1 - dh0)
        

        zperc.append(zp)

        

    return zperc
    

######################################################################
######################################################################

def plot_cond_hist(xe, ye, h, ax, 
                   axis = 0, logscale = True, **kwargs):

    '''
    Use output from the numpy 2d histogram routine 
    to calculate percentiles of the y variable
    based on relative occurence rates.


    INPUT
    =====
    xe: edge values of x variable
    ye: edge values of y variable
    h: occurence rates (relative or absolute)
    ax: axis object where to plot the hist

    optional:
    ========
    axis: the axis used for normalizeation (default: 0)
    logscale: if h is plotted in logarithmic scale (default: True)
    **kwargs: pcolormesh keywords



    OUTPUT
    ======
    pcm pcolormesh object
    '''

    # calculate norm
    norm0, dum  = np.meshgrid(h.sum(axis=0), np.ones(h.shape[0]))
    norm1, dum  = np.meshgrid(h.sum(axis=1), np.ones(h.shape[1]))

    norm1 = norm1.transpose()


    # do normalization
    if axis == 0:
        hn = np.where(norm0 == 0, 0, h / norm0)
    elif axis == 1:
        hn = np.where(norm1 == 0, 0, h / norm1)


    # transform to logscale if wanted
    if logscale:
        v = np.ma.log10(hn)
    else:
        v = hn

    # set the grid
    y, x = np.meshgrid(ye[:-1], xe[:-1])
    

    # do plotting
    pcm = ax.pcolormesh(x, y, v, **kwargs)

    return pcm


######################################################################
######################################################################


def plot_hist_median_and_intervals(xe, ye, h, ax):

    '''
    Use output from the numpy 2d histogram routine 
    to calculate percentiles and maximum, etc.
    and plot it onto on axes object.


    INPUT
    =====
    xe: edge values of x variable
    ye: edge values of y variable
    h: occurence rates (relative or absolute)
    ax: axis object where to plot the hist

    optional:
    ========
    **kwargs: plot keywords



    OUTPUT
    ======
    ax: axes object
    '''

 
    yave, ysig = ave_sig_from_hist(xe, ye, h)
    ymax =  max_from_hist(xe, ye, h)

    y25, y50, y75 = percentiles_from_hist(xe, ye, h, [25, 50, 75], axis = 1)

    x = 0.5*(xe[1:] + xe[:-1])
    #x = xe[:-1]

    sig = 2
   # ymax = scipy.ndimage.gaussian_filter1d(ymax, sig)
    y25 = scipy.ndimage.gaussian_filter1d(y25, sig)
    y50 = scipy.ndimage.gaussian_filter1d(y50, sig)
    y75 = scipy.ndimage.gaussian_filter1d(y75, sig)

#    ax.plot(x, ymax, 'w-', lw = 3)

    ax.plot(x, y25, 'k--', lw = 3)
    ax.plot(x, y50, 'k-', lw = 3)
    ax.plot(x, y75, 'k--', lw = 3)

    return ax

######################################################################
######################################################################


def hist3d(v1, v2, v3, bins):

    '''
    A wrapper for 3d histogram binning. 

    It additionally generates the average bin values with the same 3d shape 
    as histogram itself.

    '''

    h, bs = np.histogramdd([v1, v2, v3], bins)

    # generate average bin value
    bave = []
    for b in bs:
        bave.append(gi.lmean(b))


    bins3d = np.meshgrid(*bave, indexing = 'ij')

    return bins3d, h


######################################################################
######################################################################


def axis_average_from_hist3d(bins3d, h, axis = 0):


    y = bins3d[axis]

    ave =  np.where(h.sum(axis = axis) == 0, 
                    0, (h * y).sum(axis = axis) / h.sum(axis = axis))

    return ave

######################################################################
######################################################################


def conditioned_hist(h, axis = 0):

    '''
    Make conditioned histogram. 


    Prob per bin (not PDF!!!).

    '''


    
    hsum0 = h.sum(axis = 0)        
    hsum1 = h.sum(axis = 1)        

    hs0, hs1 = np.meshgrid(hsum0, hsum1)

    if axis == 0:
        hsum = hs0
    elif axis == 1:
        hsum = hs1

    return np.where( hsum == 0, 0, h / hsum) 

######################################################################
######################################################################



if __name__ == '__main__':

    print 'test'

    rstate = np.random.get_state()
    np.random.set_state(rstate)
    
    x = 2* np.random.randn(1000)
    y = x**2 + np.random.randn(1000)
#    y = x + np.random.randn(1000)


    pl.figure()
    pl.plot(y,x, 'o')

    h, xe, ye = np.histogram2d(y, x, (18,20), normed = True)

    hn = conditioned_hist(h, axis = 0)

    pl.pcolormesh(xe, ye, hn.transpose(), alpha = 0.5)
    
    x25, x50, x75 = percentiles_from_hist(xe, ye, h) 
    pl.plot(x50, gi.lmean(ye), 'r-', lw = 3)
    pl.plot(x25, gi.lmean(ye), 'g-', lw = 3)
    pl.plot(x75, gi.lmean(ye), 'g-', lw = 3)

    pl.figure()
    pl.plot(gi.lmean(ye),x50)
    pl.plot(x,x**2,'o')

# 
#    pl.figure()
#    pl.plot(x,y, 'o')
# 
#    h, xe, ye = np.histogram2d(x, y, (20,8), normed = True)
# 
#    hn = conditioned_hist(h, axis = 1)
# 
#    pl.pcolormesh(xe, ye, hn.transpose(), alpha = 0.5)
# 
#    y25, y50, y75 = percentiles_from_hist(xe, ye, h, axis  = 1) 
#    pl.plot(gi.lmean(xe), y50, 'r-', lw = 3)
#    pl.plot(gi.lmean(xe), y25, 'g-', lw = 3)
#    pl.plot(gi.lmean(xe), y75, 'g-', lw = 3)
# 
#    pl.figure()
#    pl.plot(gi.lmean(xe), y50)
#    pl.plot(x,x,'o')
