#!/usr/bin/python

"""
A module to make shaded plots.

It contains a function for non-linear colormaps.
"""

import pylab as py
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.mlab import stineman_interp

######################################################################
######################################################################

class nlcmap(LinearSegmentedColormap):
    
    """
    A nonlinear colormap class.
    
    Parameters
    ----------
    cmap : matplotlib.cmap
        a matplotlib colormap, e.g. matplotlib.cm.jet
    
    levels : list or np.array
        list of values where different colors should appear


    Notes
    -----
    Derived from the matplotlib.colors.LinearSegmentedColormap 
    
    Copyright (c) 2006-2007, Robert Hetland <hetland@tamu.edu>
    Release under MIT license.
    """
    
    name = 'nlcmap'
    
    def __init__(self, cmap, levels):

        '''
        Setting up a non-linear colormap.
        '''


        self.cmap = cmap
        self.monochrome = self.cmap.monochrome
        self.levels = np.asarray(levels, dtype='float64')
 
        lmin = self.levels.min()
        lmax = self.levels.max()

        self._x = (self.levels - lmin) / (lmax - lmin)
        self._y = np.linspace(0.0, 1.0, len(self.levels))

    # ==================================================

    def __call__(self, xi, **kwargs): 

        '''
        Applying the color mapping.

        
        Parameters
        ----------
        xi : np.array
            values mapped into color space
        **kwargs : dict
            optional arguments passed to the cmap instance
        
        '''

        yi = stineman_interp(xi, self._x, self._y)
        return self.cmap(yi, **kwargs)

######################################################################
######################################################################


def set_levs(nmax, depth, largest = 5, sym = True, sign = True, zmax = 'None'):

    '''
    Makes a non-linear level set based on pre-defined numbers.

    
    Parameters
    ----------
    nmax : int
        exponent of maximum number

    depth : int
        number of iterations used down to scales of 10**(nmax - depth)

    sym : {True, False}, optional
        switch if levels are symmetric around origin

    sign : {True, False}, optional
        switches sign if negative

    zmax : float, optional, default = 'none'
        limiter for levels, |levs| > zmax are not allowed

    
    Returns
    --------
    levs : np.array
        set of non-linear levels
    
    '''

    # nmax: 
    # depth: 

    # base
    if largest == 5:
        b = np.array([1,2,3,5], dtype = 'float')
    elif largest == 8:
        b = np.array([1,2,3,5,8], dtype = 'float')
    
    # exponents
    nmin = nmax - depth + 1
    ex = np.arange(nmin, nmax + 1, dtype = 'float')

    levs = np.array([])

    # calculation of levels ------------------------------------------
    # range below one
    for n in ex:
        levs = np.append(levs, b * 10**n)

    # ================================================================

    if not zmax == 'None':
       levs = levs[np.logical_and(levs>=-zmax,levs<=zmax)] 

    if sym:
        return np.append(-levs[::-1],np.append(0,levs))
    elif not sign:
        return np.append(-levs[::-1],0)
    else:
        return np.append(0,levs)
 
######################################################################
######################################################################


def shaded(x, y, z, *args, **kwargs):
    '''
    The function shaded is a wrapper for the pylab function contourf.


    In addition to contourf, shaded can plot filled contours for which
    a non-linear colormap is used.


    Parameters
    ----------
    x : np.array 
        x-values passed to `plt.contourf`

    y : np.array 
        y-values passed to `plt.contourf`

    z : np.array 
        z-values passed to `plt.contourf` (color value)

    *args : list
        other positional arguments passed to `plt.contourf`

    **kwargs : dict
        other optional arguments passed to `plt.contourf`

        special keywords:
 
        * 'levels' : numpy array of color / contour levels

        * 'lev_depth' : gives the depth of an automatically generated level set
          i.e. an initial base (e.g. [1,2,3,5] ) that contains
          the maximum value of the field (e.g. 4.3) is downscaled
          by 10**-1, ..., 10**-lev_depth.

          Example:
          Given lev_depth = 2, for a positive field with maximum 4.3
          a level set [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5]
          is generated.
      
    Returns
    --------
    pl.contourf instance
    '''

    # get level set from arguments if needed -------------------------
    if args:
#        print 'args not empty'

        # assume that the first argument are levels !!!
        kwargs['levels'] = args[0]
        DO_NONLIN_CMAP = True
    # ================================================================



    # check distribution of z ----------------------------------------
    if z.min() > 0:
        IS_SYMMETRIC = False
        IS_POSITIVE = True

    elif z.max() < 0:
        IS_SYMMETRIC = False
        IS_POSITIVE = False

    else:
        IS_SYMMETRIC = True
        IS_POSITIVE = False
    # ================================================================



    # get levels for the non-linear colormap -------------------------
    # from keyword ...................................................
    if 'levels' in kwargs and 'lev_depth' not in kwargs:
        clev = kwargs['levels']


    # check if lev to be setted ......................................

    # NOTE: 
    # =====
    # if the lev_depth keyword is set it will overwrite any
    # level information given additionally

    else:


        if 'lev_depth' in kwargs:
            depth = kwargs['lev_depth']
        else:
            # default value 
            depth = 3

        zmax = np.abs(z).max()
        largest_in_base = 5   # only 5 or 8 can be chosen !

        imax = int(0.9999999999 + (np.log(zmax) - \
                                       np.log(largest_in_base)) / np.log(10))


        kwargs['levels'] = set_levs(imax, depth, \
                                        sign = IS_POSITIVE, \
                                        largest = largest_in_base, \
                                        sym = IS_SYMMETRIC,\
                                        zmax = zmax)
    # ================================================================
            
        
    # get the current colormap ---------------------------------------
    if 'cmap' in kwargs:
        cmap =  kwargs['cmap']
    else:
        # default is set to jet !
        cmap = py.cm.jet
    # ================================================================


    # make a non-linear colormap -------------------------------------
    clev = kwargs['levels']
    kwargs['cmap'] = nlcmap(cmap, clev)
    # ================================================================

    # range of colors ------------------------------------------------
    if 'extend' not in kwargs:
        kwargs['extend'] = 'both'
    # ================================================================
    
        
    return py.contourf(x, y, z, *args, **kwargs)
    
######################################################################
######################################################################

if __name__ == '__main__':
    
    y, x = np.mgrid[0.0:3.0:100j, 0.0:5.0:100j]
    H = 50.0 * np.exp( -(x**2 + y**2) / 4.0 )

    clev = set_levs(0,1,sym=True)
    
    py.subplot(2,2,1)
    shaded(x, y, H)
    py.colorbar()

    py.subplot(2,2,2)
    shaded(x, y, H, clev, cmap=py.cm.jet_r)
    py.colorbar()

    py.subplot(2,2,3)
    shaded(x, y, H, lev_depth = 5)
    py.colorbar()

    py.subplot(2,2,4)
    shaded(x, y, -H)
    py.colorbar()


    py.show()
