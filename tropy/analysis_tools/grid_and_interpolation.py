#!/usr/bin/env python


# Copyright (c) 2012- 2016

# TROPOS,
# Permoserstr. 15
# 04318 Leipzig, Germany. 

# Author:
# ====== 
# Fabian Senf <senf@tropos.de>


# This program is free software; you can redistribute it and/or 
# modify it under the terms of the GNU General Public License as 
# published by the Free Software Foundation; either version 3 of 
# the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU General Public License for more details.


import sys, copy

import numpy as np
from scipy.interpolate import griddata
import scipy.interpolate

try:
    import SimpleITK
except:
    print('SimpleITK not available.')


######################################################################
######################################################################
######################################################################
#####                          #######################################
##### GRID MANIPULATIONS, ETC. #######################################
#####                          #######################################
######################################################################
######################################################################
######################################################################



def COSMO_mean_vertical_levels(Nlev):

    ''' 
    Outputs the mean vertical levels of COSMO-DE taken from the 
    model description of COSMO (see reference below).


    Parameters
    ----------
    Nlev : int
        number of vertical levels, i.e.
        
        Nlev = 51 for intermediate half-levels (e.g. for w)
        Nlev = 50 for main levels (e.g. for u)


    Returns
    --------
    lev : numpy array (1-dim)
        mean vertical COSMO-DE levels 

           
    References
    ----------
    M. Baldauf, J. Foerstner, S. Klink, T. Reinhardt, C. Schraff, A. Seifert und K. Stephan
    "Kurze Beschreibung des Lokal-Modells Kuerzestfrist COSMO-DE (LMK) und seiner Datenbanken
    auf dem Datenserver des DWD", Version 1.6, Stand: 31. 03. 2011
    
    Tabelle 1, Seite 29

    '''



    zm = np.array([22000.0000,21000.0000,20028.5703,19085.3594,18170.0000,17282.1406,16421.4297,\
                   15587.5000,14780.0000,13998.5703,13242.8594,12512.5000,11807.1367,11126.4297,\
                   10470.0000,9837.5000,9228.5703,8642.8594,8080.0000,7539.6367,7021.4297,\
                   6525.0000,6050.0000,5596.0664,5162.8594,4750.0000,4357.1367,3983.9299,\
                   3630.0000,3295.0000,2978.5701,2680.3601,2400.0000,2137.1399,1891.4299,\
                   1662.5000,1450.0000,1253.5698,1072.8599,907.5000,757.1399,621.4299,\
                   500.0000,392.5000,298.5698,217.8600,150.0000,94.6400,51.4300,20.0000,0.0000])

    if Nlev == 51:
        return zm / 1e3
    elif Nlev == 50:
        return 0.5e-3 * (zm[1:] + zm[:-1])


######################################################################
######################################################################

def ll2xy(lon, lat, lon0 = 10., lat0 = 50., R = 6380.):

    ''' 
    Transformation between longitude and latitude and local Cartesian coordinates
    in west-east direction (x) and south-north direction (y). The local projection 
    is called "sinusoidal".


    Parameters
    ----------
    lon : numpy array
        longitude

    lat : numpy array
        latitude

    lon0 : float, optional, default = 10.
        arbitrary initial longitude for x = 0

    lat0 : float, optional, default = 50.
        arbitrary initial latitude for y = 0

    R : float, optional, default = 6380.
        approximation for the radius of Earth
    

    Returns
    --------
    x : numpy array
        Cartesian coordinate in west-east direction (increasing towards the East)

    y : numpy array
        Cartesian coordinate in south-north direction (increasing towards the North)

    
    Notes 
    -----
    The unit of Earth's radius determines the units of x and y. 
    Default (km).


    References
    ----------
    See https://en.wikipedia.org/wiki/Sinusoidal_projection for further info.

    '''


    # arbritary initial value for longitude and latitude
    lam0 = np.deg2rad(lon0)
    phi0 = np.deg2rad(lat0)


    # conversion from degree to radiant
    lam = np.deg2rad(lon)
    phi = np.deg2rad(lat)
    
    # get x and y from lam and phi
    x = R  * np.cos(phi - phi0) * (lam - lam0)
    y = R * (phi - phi0)

    return x,y 

######################################################################
######################################################################

def ll2xyc(lon, lat, mlon = None, mlat = None, lon0 = 0, lat0 = 0.):

    '''
    Applies a centered version of ll2xy. The centering is around the
    mean lon/lat value and lon0, lat0 only define the x/y offset.


    Parameters
    ----------
    lon : numpy array
        longitude
    
    lat : numpy array
        latitude
    
    mlon : float, optional, default = None
        longitude center point of projection

        if ``None`` then average of ``lon`` is used

    mlat : float, optional, defaulft = None
        latitutde center point of projection

        if ``None`` then average of ``lat`` is used

    lon0 : float, optional, default = 0.
        zero longitude (for x-offset)

    lat0 : float, optional, default = 0.
        zero latitude (for y-offset)
        

    Returns
    -------
    x : numpy array
        x-coordinate

    y : numpy array
        y-coordinate
    '''

    # define the center point
    if mlon is None:
        mlon = lon.mean()

    if mlat is None:
        mlat = lat.mean()

    # apply sinusoidal transformation
    x, y = ll2xy( lon, lat, lon0 = mlon, lat0 = mlat )

    # and finally shift x/y with offset
    x0, y0 = ll2xy( lon0, lat0, lon0 = mlon, lat0 = mlat)
    x -= x0
    y -= y0

    return x, y


######################################################################
######################################################################




def xy2ll(x, y, lon0=10., lat0=50., R = 6380):

    ''' 
    Transformation between local Cartesian coordinates in west-east direction (x) 
    and south-north direction (y) and longitude and latitude. This assumes a 
    sinusoidal projection. Inversion of the function ``ll2xy``.



    Parameters
    ----------
    x : numpy array
        Cartesian coordinate in west-east direction (increasing towards the East)

    y : numpy array
        Cartesian coordinate in south-north direction (increasing towards the North)

    lon0 : float, optional, default = 10.
        arbitrary initial longitude for x = 0

    lat0 : float, optional, default = 50.
        arbitrary initial latitude for y = 0

    R : float, optional, default = 6380.
        approximation for the radius of Earth
    

    Returns
    --------
    lon : numpy array
        longitude

    lat : numpy array
        latitude

    
    Notes 
    -----
    The unit of Earth's radius determines the units of x and y. 
    Default (km).


    References
    ----------
    See https://en.wikipedia.org/wiki/Sinusoidal_projection for further info.

    '''


    # arbritary initial value for longitude and latitude
    lam0 = np.deg2rad(lon0)
    phi0 = np.deg2rad(lat0)

    # get phi and lam from x and y 
    phi = y / R + phi0
    lam = x / ( R * np.cos(phi - phi0) ) + lam0
    
    # conversion from radiant to degree
    lon = np.rad2deg(lam)
    lat = np.rad2deg(phi)
    

    return lon, lat


######################################################################
######################################################################


def mid2edge_gridvector(x_mid, x_edge0):
    
    '''
    Converts an mid-point based grid vector into an edge-based grid vector.
    

    Parameters
    ----------
    x_mid : numpy array (1-dim)
        mid-point-based grid vector (length N)

    x_edge0 : float
        left edge
    

    Returns
    --------
    x_edge : numpy array (1-dim)
        edge-based grid vector (length N + 1)
    '''
    
    # get length of grid
    N = len(x_mid) + 1
    
    # init edge-based vector
    x_edge = np.zeros(N)
    x_edge[0] = x_edge0
    
    # loop over all positions
    for i in range(1, N):
        x_edge[i] = 2 * x_mid[i - 1] - x_edge[i - 1]
        
    return x_edge


######################################################################
######################################################################
######################################################################
#####                          #######################################
##### HRV Upscaling            #######################################
#####                          #######################################
######################################################################
######################################################################
######################################################################


def make_hrv_upscaling(v):
    
    '''
    A 2d field is upscaled to 3-fold higher resolution using linear 
    interpolation. Edge values are determined by linear extrapolation.
   

    Parameters
    ----------
    v : numpy array (2-dim)
        2d variable field (Nrow, Ncols)


    Returns
    --------
    vhigh : numpy array (2-dim)
        variable field with 3-times higher resolution (3*Nrows, 3*Ncols)
    
    '''

    var = make_add_edge(v)

    # x shift first ..................................................
    dvar_x = var[1:,:] - var[:-1,:]
    var0 = var[:-1,:]
    var1 = var0 + 1./3. * dvar_x
    var2 = var0 + 2./3. * dvar_x

    Nx,Ny = var0.shape

    varx = np.zeros((3*Nx,Ny))
    
    varx[0::3,:] = var0
    varx[1::3,:] = var1
    varx[2::3,:] = var2
    
    dvar_y = varx[:,1:] - varx[:,:-1]
    var0 = varx[:,:-1]
    var1 = var0 + 1./3. * dvar_y
    var2 = var0 + 2./3. * dvar_y

    Nx,Ny = var0.shape

    vary = np.zeros((Nx,3*Ny))
    
    vary[:,0::3] = var0
    vary[:,1::3] = var1
    vary[:,2::3] = var2

    vhigh =  vary[2:-1,2:-1]
    

    return vhigh



######################################################################
######################################################################


def make_add_edge(var):

    '''
    Adds egde region to 2d data field. Values in added left and right column, 
    and lower and upper row are linearily extrapolated from their neighbors.
   

    Parameters
    ----------
    var : numpy array (2-dim)
        2d variable field (Nrow, Ncols)

    
    Returns
    --------
    var_new : numpy array (2-dim)
        2d variable field with added edge values (Nrows + 2, Ncols + 2)
    
    '''
    
    N1, N2 = var.shape
    
    var_new = np.zeros((N1+2,N2+2))
    
    var_new[1:N1+1,1:N2+1] = var
    
    # upper edge
    dv = var[1,:] - var[0,:]
    var_new[0,1:N2+1] = var_new[1,1:N2+1] - dv

    # lower edge
    dv = var[-2,:] - var[-1,:]
    var_new[-1,1:N2+1] = var_new[-2,1:N2+1] - dv

    # left edge
    dv = var_new[:,2] - var_new[:,1]
    var_new[:,0] = var_new[:,1] - dv

    # right edge
    dv = var_new[:,-3] - var_new[:,-2]
    var_new[:,-1] = var_new[:,-2] - dv

    

    return var_new



######################################################################
######################################################################
######################################################################
#####                          #######################################
##### INTERPOLATION STUFF, ETC #######################################
#####                          #######################################
######################################################################
######################################################################
######################################################################




def get_index(p, lon, lat):
    
    '''
    Given a point p the nearest grid point is returned.


    Parameters
    ----------
    p : numpy array
        point in the 2d grid

    lon : numpy array (2-dim)
        grid values of the 1st dimension (can be longitude, x, etc.)

    lat : numpy array (2-dim)
        grid values of the 2nd dimension (can be latitude, y, etc.)


    Returns
    --------
    ir : numpy array (2-dim)
        row index (the 1st)

    ic : numpy array (2-dim)
        column index (the 2nd)
    '''

    # get shapes of fields 
    Nr, Nc = lon.shape
    
    # start couting the number of rows and columns
    jr = np.arange(Nr)
    jc = np.arange(Nc)
    
    # make a mesh of row and column counts
    ic, ir = np.meshgrid(jc, jr)

    # transfrom the point in a shape used in the interpolation routine
    p = np.asarray([p]).T
    
    # do the interpolation in the row and column mesh
    irp = interpolate_field_to_hor_pos(p, lon, lat, ir, method='nearest')
    icp = interpolate_field_to_hor_pos(p, lon, lat, ic, method='nearest')
    
    return irp, icp

######################################################################
######################################################################


def interpolate_field_to_hor_pos(pos, lon, lat, field, **kwargs):

    '''
    Given a sequence of points (or just one) a field (2d or 3d) given a
    grid is interpolated to the new positions.

        
    Parameters
    ----------
    pos : numpy array (2-dim) 
         sequence of points given as np.array((x,y))
         where x and y can be n-dimensional 

    lon : numpy array (2-dim)
        grid values of the 1st dimension (can be longitude, x, etc.)

    lat : numpy array (2-dim)
        grid values of the 2nd dimension (can be latitude, y, etc.)

    field : numpy array
        2d or 3d field to be interpolated to pos

    **kwargs : dict content
        keyword arguments for the griddata routine of
        the scipy.interpolate package can be passed through  


    Returns
    --------
    f : numpy array 
       field interpolated to the position sequence pos,
       f has the dimension n x Nlev, 
       where n is the dimension of pos, Nlev the vertical dimension

    '''

      
    # set dimensions .................................................
    Nhor = lon.flatten().shape[0]

    if len(field.shape) == 3:
        Nlev = field.shape[2]
    else:
        Nlev = 1


    # set input coordinates ..........................................
    Pin = np.zeros((Nhor,2))
    Pin[:,0] = lon.flatten()
    Pin[:,1] = lat.flatten()
    
    # set input field ................................................
    fin = field.flatten().reshape(Nhor, Nlev)
    
    
    # output coordinates
    NP = pos.flatten().shape[0]/2
    Pout = pos.T
    
    # horizontal interpolation ...........................................
#    print 'Do horizontal interpolation'
    if 'method' not in kwargs:
        kwargs['method'] = 'nearest'
    
    fout = griddata(Pin,fin,Pout, **kwargs)
    
    if NP == 1:
        profile = fout.squeeze()
    else:
        profile = fout
    
    return profile

######################################################################
######################################################################



def make_index_set(nrows, ncols):

    '''
    This simple routine makes an 2d index set for given row and column
    number of a 2d field.


    Parameters
    ----------
    nrows : int
        number of rows

    ncols : int
        number of columns
    
    
    Returns
    -------
    ii : numpy array, dtype = int
        2d field of rows indices

    jj : numpy array, dtype = int
        2d field of column indices
    '''
    
    i = np.arange(0, nrows)
    j = np.arange(0, ncols)
    ii, jj = np.meshgrid(i, j, indexing = 'ij')

    return ii, jj

######################################################################
######################################################################



def create_interpolation_index(lon1, lat1, lon2, lat2, xy = False):

    '''
    Given georeference of two grids in lon and lat, an index is built to map 
    a field on grid 1 to grid 2.

    
    Parameters
    ----------
    lon1 : numpy array (2-dim) 
        longitude of input grid 1

    lat1 : numpy array (2-dim) 
        latitude of input grid 1

    lon2 : numpy array (2-dim) 
        longitude of target grid 2

    lat2 : numpy array (2-dim) 
        latitude of target grid 2

    xy : bool, optional, default = False
        switch if georef field are interpreted as local 
        Cartesian coordinates: Then no internal transformation is used.


    Returns
    --------
    ir : numpy array, 2-dim, dtype = int 
        index field for transforming rows

    ic :  numpy array, 2-dim, dtype = int 
        index field for transforming columns
    '''


    Nr, Nc = lon1.shape
    s2 = lon2.shape

    ir = np.arange(Nr)
    ic = np.arange(Nc)

    ic1, ir1 = np.meshgrid(ic, ir)

    if xy:
        x1, y1 = lon1, lat1
        x2, y2 = lon2, lat2
    else:
        x1, y1 = ll2xy(lon1, lat1)
        x2, y2 = ll2xy(lon2, lat2)

    pos = np.array((x2.flatten(), y2.flatten()))
    
    ir2 = interpolate_field_to_hor_pos(pos, x1, y1, ir1).reshape(s2) 
    ic2 = interpolate_field_to_hor_pos(pos, x1, y1, ic1).reshape(s2)


    return ir2, ic2


######################################################################
######################################################################

def curve_flow_filter(f, numberOfIterations = 5):

    '''
    Smoothing filter depending on isoline curvature. Interface for 
    curvature flow filter from simpleITK toolkit.


    Parameters
    ----------
    f : numpy array (2-dim)
        2d field to be filtered (smoothed)

    numberOfIterations : int, optional, default = 5
        number of iterations, increases smooting effect


    Returns
    --------
    f_sm : numpy array (2-dim)
        smoothed 2d field


    Notes
    ------
    Only works if SimpleITK is installed !!!
    '''

    img = SimpleITK.GetImageFromArray(f)
    img_sm = SimpleITK.CurvatureFlow(img, numberOfIterations = numberOfIterations)
    
    f_sm = SimpleITK.GetArrayFromImage(img_sm)

    return f_sm
    

######################################################################
######################################################################

def spline_smooting(x, y, s = 0.1, fixed_endpoints = True):
    
    '''
    Smoothes a curve with smooting spline.
    
    
    Parameters
    ----------
    x : numpy array (1-dim)
        abscissa values

    y : numpy array (1-dim)
        ordinate values

    s : float, optional, default = 0.1
        smooting parameter 

    fixed_endpoints : nool, optional, default = True
        switch, if endpoinds should be hold fixed
    

    Returns
    --------
    y_smooth : numpy array (1-dim)
        smoothed version of y


    Notes
    ------
    Uses the function ``scipy.interpolate.UnivariateSpline``.
    '''
   
    # check for nan values
    ym = np.ma.masked_invalid( y )

    # mask for clean values
    m = np.logical_not( ym .mask ) 
    
    # set weights for possible endpoint
    w = np.ones_like(x[m])
    
    if fixed_endpoints:
        w[0] = 100.
        w[-1] = 100.
    
    # make interpolating spline
    s = scipy.interpolate.UnivariateSpline(x[m], ym[m], s = s, w = w)


    # write spline-interpolated values on output field (allow for NaNs)
    f = np.nan * np.ones_like( y )
    f = np.ma.masked_invalid( f )
    f[m] = s(x[m])
    
    return f
    
######################################################################
######################################################################

def remap_field(lon, lat, f, dr = 0.5):

    '''
    Do nearest nearbor remapping on a simple equi-distant
    local cartesian coordiante grid.
    
    Uses curvature flow filter for smoothing and pixel edge
    correction.


    Parameters
    ----------
    lon : numpy array, (2-dim)
        longitude

    lat : numpy array, (2-dim)
        latitude
    
    f : numpy array, (2-dim)
        2d field to be interpolated

    dr : float, optional, default = 0.5
        grid spacing in km


    Returns
    --------
    xnew : numpy array, (2-dim)
        regridded, equi-distant x-coordinate in km

    ynew : numpy array, (2-dim)
        regridded, equi-distant y-coordinate in km

    fint : numpy array, (2-dim)
        remapped field
    
    '''

    # convert to local coordiantes
    x, y = ll2xy(lon, lat)


    # make new grid
    xn = np.arange(x.min(), x.max(), dr)
    yn = np.arange(y.min(), y.max(), dr)
    xnew, ynew = np.meshgrid(xn, yn)

    nlon, nlat = xy2ll(xnew, ynew)
    ir, ic = create_interpolation_index(lon, lat, nlon, nlat)

    fint = curve_flow_filter(f[ir, ic], 3)

    return xnew, ynew, fint

######################################################################
######################################################################


def make_vert_cut(p1, p2, lon, lat, vg, **kwargs):


    '''
    Makes a horizontal - vertical cut through a 3d field.


    Parameters
    ----------
    p1 : list or numpy array
       1st end of the cutting line (lon, lat)

    p2 : list or numpy array
       2nd end of the cutting line (lon, lat)

    lon : numpy array (2-dim)
        longitude

    lat : numpy array (2-dim)
        latitude

    vg : numpy array (3-dim)
        3d field to be interpolated

    **kwargs: dict content
         keyword arguments for the griddata routine of
         the scipy.interpolate package can be passed through  


    Returns
    --------
    s : numpy array (1-dim)
        distance of km from one to the other end of the cutting line

    lo : numpy array (1-dim)
        longitude along the cut

    la : numpy array (1-dim)
        latitude along the cut

    v : numpy array (2-dim)
        field interpolated to the cutting line
    
        Vertical levels are kept.

    '''


    # prepare grid and start/end points ------------------------------
    # get grid in x and y
    xg, yg = ll2xy(lon,lat)

    # get points
    x1, y1 = ll2xy(p1[0],p1[1])
    x2, y2 = ll2xy(p2[0],p2[1])


    # get distance between start and end point
    d = abs((x1 - x2) + 1j * (y1 -y2))
    # ================================================================
    

    # prepare subsampling along the cutout ---------------------------

    # set smallest distance between the new grid points
    ds = 2.


    # number of new along-cutout grid points 
    N = int(d / ds)


    # subsample along the cutout
    x = np.linspace(x1,x2,N,endpoint=True)
    y = np.linspace(y1,y2,N,endpoint=True)

    pos = np.array((x,y))

    # get the along-cut distance
    x0 = x[0]
    y0 = y[0]

    lo, la = xy2ll(x, y)
    s = abs((x - x0) + 1j * (y - y0))
    # ================================================================


    
    # do interpolation -----------------------------------------------
    if 'method' not in kwargs:
        kwargs['method'] = 'nearest'
    v = interpolate_field_to_hor_pos(pos, xg, yg, vg, **kwargs)
    # ================================================================
    

    
    return s, lo, la, v

######################################################################
######################################################################


def simple_pixel_area(lon, lat, xy = False, uncertainty = False):

    ''' 
    Approximates the grid box area of a lon - lat grid.
    
    The grid is treated as if it is an equidistant grid in local 
    Cartesian coordinates on Earth's surface (i.e. in x and y). 
    The grid point is placed at each of the for corners of the grid box
    and the corresponding average area ist calculated.


    Parameters
    ----------
    lon : numpy array (2-dim)
        longitude

    lat : numpy array (2-dim)
        latitude
    
    xy : bool, optional, default = False
        switch if ``lon`` & ``lat`` are interpreted as local Cartesian coordinates

    uncertainty : bool, optional, default = False
        switch if std of pixel area is output
    
        Deviation are caused by effects of the not-rectangular grid.
    

    Returns
    --------
    a : numpy array (2-dim)
        grid box area

    da : numpy array (2-dim), optional
        standard deviation of grid box area, if ``uncertainty == True``
    
    '''

    if not np.ndim(lon) == 2 and not np.ndim(lat) == 2:
        print('ERROR: please only use 2d matrices for lon and lat')
        sys.exit(0)

    # convert lon,lat into Cartesian coordinates
    if xy:
        x, y = lon, lat
    else:
        x, y = ll2xy(lon, lat)

    # add some edge via linear extrapolation
    x_ext = make_add_edge(x)
    y_ext = make_add_edge(y)

    # make forward, backward difference assuming that 
    #  *) x is connected to column
    dx1 = x_ext[1:-1, 2:] - x 
    dx2 = x - x_ext[1:-1,:-2]

    #  *) y is connected to row
    dy1 = y - y_ext[2:, 1:-1]  
    dy2 = y_ext[:-2, 1:-1] - y
    
    # calculate the area assuming that the grid point is located at 
    # each of the four corners of the grid box
    ac = []
    ac.append(dx1 * dy1)
    ac.append(dx1 * dy2)
    ac.append(dx2 * dy1)
    ac.append(dx2 * dy2)

    # calculate average area
    am = np.dstack(ac).mean(axis=2)
    

    if uncertainty:
        # calculate deviation of different area 
        # (perhaps a measure of uncertainty that simple method)
        da = np.dstack(ac).std(axis=2)
        return am ,da

    else:
        return np.abs( am )


######################################################################
######################################################################


def i2iset(v,i):

    '''
    Returns an index set for a n-dim numpy array given an index i
    which points to a position within the same but flattened array.
    
    It is assumed that the array was flattened in C-mode.


    Parameters
    ----------
    v : numpy array (n-dim)
        field which determines the rank and shape

    i : int
        input index which points to the position in v.flatten()

    
    Returns
    --------
    iset : tuple
        index set which points to the position in v


    Notes
    ------
    This routine helps you to locate a maximum in a n-dimensional array
    using e.g. the ``np.argmax`` function.
    '''
   

    # get shape of input n-dim array
    s = v.shape
   
    # calculate factors which determine the jumps within 
    # the flattened array
    factors = np.cumprod(s[::-1])[:-1]
 
    # loop over different dimension of the n-dim array
    res = i
    iset = []
    for f in factors[::-1]:
        iset.append(res / f)
        res = i % f

    iset.append(res)

    return tuple(iset)



######################################################################
######################################################################


def ldiff(var, axis = -1):

    '''
    Calculates difference of a field at adjacent levels. Now, wrapper
    for ``numpy.diff`` with default on last axis.


    Parameters
    ----------
    var : numpy array (n-dim where n = [1,2,3])
        field; 1d, 2d or 3d array

    axis : int, optional, default = -1
        axis used for layer avering
    
    
    Returns
    -------
    dvar : numpy array (n-dim where n = [1,2,3])
        level difference
    
    '''

 
    return np.diff(var, axis = axis)
 

######################################################################
######################################################################


def lmean(var, axis = -1):

    '''
    Calculates the layer mean of a field.


    Parameters
    ----------
    var : numpy array (n-dim where n = [1,2,3])
        field; 1d, 2d or 3d array

    axis : int, optional, default = -1
        axis used for layer avering
    
    
    Returns
    -------
    varm : numpy array (n-dim where n = [1,2,3])
        layer mean
    
    '''

    # select the axis for layer-average calculations
    if axis == -1:
        ndim = len(var.shape)
    else:
        ndim = axis + 1

    if ndim == 1:
        dvar = 0.5 * (var[1:] + var[:-1])

    elif ndim == 2:
        dvar = 0.5 * (var[:,1:] + var[:,:-1])

    elif ndim == 3:
        dvar = 0.5 * (var[:,:,1:] + var[:,:,:-1])

    return dvar
 


######################################################################
######################################################################



def cutout_fields(fin, slices, vaxis = 0):

    '''
    Cuts out a field or a list of similar fields given row and column slices.


    Parameters
    ----------
    fin : 2d or 3d numpy array or list of fields
        the fields to be cut out

    slices : tuple or list
        region slice such that the cutout is done as ((row1, row2), (col1, col2))

    vaxis: int, optional, default = 0
        for 3d fields, axis which is still varying
    

    Returns
    --------
    fout : 2d or 3d numpy array or list of fields
        resulting cutout (possibly a list of cutted fields)
    
    '''

    (ir1, ir2), (ic1, ic2) = slices

    fnew = []

    if not type(fin) == type([]):
        flist = [fin]
    else:
        flist = fin

    for f in flist:
        if np.ndim(f) == 2:
            fcut = f[ir1:ir2, ic1:ic2]
        elif np.ndim(f) == 3:
            if vaxis == 0:
                fcut = f[:, ir1:ir2, ic1:ic2]

            elif vaxis == 2:
                fcut = f[ir1:ir2, ic1:ic2, :]


        fnew.append(fcut)


    if not type(fin) == type([]):
        return fnew[0]
    else:
        return fnew


######################################################################
######################################################################

def cutout_field4box(f, ind, bsize, **kwargs):

    '''
    Cuts out a field based on center pix index and box size.


    Parameters
    ----------
    f : 2d or 3d numpy array
        the field to be cut out

    ind : tuple or list
        center index of horizontal cutout

    bsize : int
        size of the quadratic cutout box


    Returns
    --------
    fcut : 2d or 3d numpy array
        resulting quadratic cutout 

    '''

    # extract the vaxis keyword if there
    vaxis = kwargs.get('vaxis', 0)

    # get the fill value
    fill_value = kwargs.pop('fill_value', 0)


    # get field dimensions and prepare stacking if needed for adding edges

    if np.ndim(f) == 2:
        nrows, ncols = f.shape

    elif np.ndim(f) == 3:
        if vaxis == 0:
            ntime, nrows, ncols = f.shape

        if vaxis == 2:
            nrows, ncols, nlev = f.shape

    # get index values
    irow, icol = ind


    # set box edges
    ir1 = irow -  bsize // 2
    ir2 = irow +  bsize // 2  +  1

    ic1 = icol -  bsize // 2
    ic2 = icol +  bsize // 2  +  1

    # check for box which extend the field dimensions
    dc_left, dc_right = (0, 0)
    dr_left, dr_right = (0, 0)

    if ic1 < 0:
        dc_left =  - ic1
        ic1 = 0

    if ic2 > ncols:
        dc_right = ic2 - ncols
        ic2 = ncols  

    if ir1 < 0:
        dr_left =  - ir1
        ir1 = 0

    if ir2 > nrows:
        dr_right = ir2 - nrows
        ir2 = nrows 
   
    #set region 
    reg = ((ir1, ir2), (ic1, ic2))

    # do the cut
    fcut =  cutout_fields(f, reg, **kwargs)

    if np.ndim(f) == 2:
        fcut = np.ma.reshape(fcut, fcut.shape + (1,))

    elif np.ndim(f) == 3 and vaxis == 0:
        fcut = fcut.transpose(1,2,0)


    # and now add edges if needed
    s = fcut.shape
    if dc_left != 0:
        z = np.ma.zeros( (s[0], dc_left, s[2] ) )
        fcut = np.column_stack([z, fcut])

    if dc_right != 0:
        z = np.ma.zeros( (s[0], dc_right, s[2]) )
        fcut = np.column_stack([fcut, z])

    s = fcut.shape
    if dr_left != 0:
        z = np.ma.zeros( (dr_left,) + s[1:] )
        fcut = np.row_stack([z, fcut])

    if dr_right != 0:
        z = np.ma.zeros( (dr_right,) + s[1:] )
        fcut = np.row_stack([fcut, z])

    if np.ndim(f) == 3 and vaxis == 0:
        fcut = fcut.transpose(2,0,1)


    return fcut.squeeze().astype(f.dtype)

    
####################################################################
####################################################################


def tube_cutout4box(v3d, i0, boxsize):

    '''
    Performs tube cutout with given index set. The center index might not
    be constant.


    Parameters
    ----------
    v3d : numpy array (3-dim)
        3d field, 1st dimension is time (or variable)

    i0 : list or numpy array (2-dim)
        center indices
        (irvec, icvec) index tuple, (row index vector, co index vector)

    boxsize : int
        size of the quadratic cutout box


    Returns
    --------
    vtube : numpy array (3-dim)
        3d tube, quadratic cutout in the last two dimensions
    

    Notes
    ------
    This function can be used top make vertical cutouts, but a transpose command
    has to be applied before and after.
    '''

    # get number of time steps
    ntime = v3d.shape[0]

    # extract index vector
    irvec, icvec = i0

    # initialize tube field
    vtube = np.ma.zeros((ntime, boxsize, boxsize))
    
    # and loop over time for cutout
    for n in range(ntime):
 
        # row and column index
        ir = irvec.flatten()[n]
        ic = icvec.flatten()[n]

        # cutout box
        vtube[n] =  cutout_field4box(v3d[n], (ir, ic), boxsize)

    return vtube


######################################################################
######################################################################


def cutout_cluster(c, nc, 
                   nedge = 30):
    
    '''
    Makes a cutout of a categorial field c for an object of class / number nc.


    Parameters
    ----------
    c : numpy array, 2-dim, dtype = int
        2d categorial field

    nc : int
        category number / class which is chosen

    nedge : int, optional, default = 30
        number of edge pixels added / included 


    Returns
    --------
    ccut : numpy array, 2-dim, dtype = int 
        cutout of categorial field c
    '''


    # get index field ------------------------------------------------
    nrow, ncol = c.shape
    irow, icol = make_index_set(nrow, ncol)
    # ================================================================


    # make mask and mask the indices ---------------------------------
    mask = (c== nc)

    ir = irow[mask]
    ic = icol[mask]
    # ================================================================

    
    # determine cutout region ----------------------------------------
    ir1 = ir.min() - nedge
    ir2 = ir.max() + nedge
    
    ic1 = ic.min() - nedge
    ic2 = ic.max() + nedge

    # check bounds
    if ir1 < 0:
        ir1 = 0

    if ir2 >= nrow:
        ir2 = nrow 

    # check bounds
    if ic1 < 0:
        ic1 = 0

    if ic2 >= ncol:
        ic2 = ncol 

    region = ((ir1, ir2), (ic1, ic2))
    # ================================================================


    # do the cutout --------------------------------------------------
    ccut = cutout_fields(c, region)
    # ================================================================

    return ccut

######################################################################
######################################################################

def region2slice(lon, lat, region):

    '''
    Given a region of ((lon1,lon2),(lat1, lat2)) the routine provides 
    a slicing ((row1, row2), (col1, col2)) than accomodates the region.


    Parameters
    ----------
    lon : numpy array (2-dim)
        longitude

    lat : numpy array (2-dim)
        latitude

    region : tuple of two tubles
        ((lon1, lon2), (lat1, lat2)) 
    

    Returns
    --------
    slice : tuple of two tubles
        region slice for cuuting a field ((row1, row2), (col1, col2))
    
    '''


    # get corner points
    (lon1, lon2), (lat1, lat2) = region

    # get the correspoding index of the corners 
    ir1, ic1 = get_index((lon1, lat2), lon, lat)
    ir2, ic2 = get_index((lon2, lat1), lon, lat)

    return ((ir1, ir2), (ic1, ic2))


######################################################################
######################################################################



def low2hres_region(region):

    '''
    Based on low-res. region of MSG channels an correspoding HRV region is 
    calculated.


    Parameters
    ----------
    region : tuple of two tubles
        region slice for cuuting a field ((row1, row2), (col1, col2))

    
    Returns
    --------
    hrv_region : corresponding HRV region as ((row1, row2), (col1, col2))
    

    Notes
    ------

    Calculation
    
    It is assumed that the hrv region is a dependent variable which 
    depends on the low res region attribute
    That means: Given the low res region cutout hrv region is determined 
    as corresponding cutout which 

      (i) exactly fits with the low res region and 
      (ii) refers to the artifical hrv full disk of 11136 x 11136
    
    low res pixel with the index (I,J) = (0,0)
    has a high res pixel with index (i_m, j_m) = (2,2) in the middle
    and (i_u, j_u) = (1,1) in the upper corner
    see doc "MSG Level 1.5 Image Data Format Description", Fig.8
    which then leads to the relation between the low res pixel (I,J)
    and its corresponding upper corner high res pixel (i_u, j_u):
    (i_u, j_u) = 3 * (I, J ) + 1
    '''
        
        
    (r1, r2), (c1, c2) = region 
    
    hrv_region = ((3*r1 + 1, 3*r2 + 1), (3*c1 + 1, 3*c2 + 1))
    
    return hrv_region
    

######################################################################
######################################################################


if __name__ == '__main__':
    # test cutout_field4box
 
    test = np.arange(25).reshape(( 5, 5)).astype(np.int)

    c = cutout_field4box(test, (2,2), 15)
    print(c)
    print(c.shape)


    test = np.arange(25*2).reshape((5, 5, 2)).astype(np.int)

    c = cutout_field4box(test, (2,2), 7, vaxis = 2)
    print(c[:,:,0])
    print()
    print(c[:,:,1])
    print(c.shape)
