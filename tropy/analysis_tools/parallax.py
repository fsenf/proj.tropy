#!/usr/bin/env python

import sys, os, glob, copy
import numpy as np

import grid_and_interpolation as gi
import datetime, time
import scipy.ndimage

import tropy.io_tools.radolan as radolan

######################################################################
######################################################################

def parallax_correct_on_radolan(lon, lat, cth, f, 
                                do_masking =True,
                                satlon = 9.5, 
                                smoothing = True,
                                gauss_sig = 3,
                                do_sequential_hole_filling = False): 

    '''
    Do parallax-correction of field f based on cth (m). Corrected f 
    is remapped onto Radolan grid.

    INPUT
    =====
    lon: longitude
    lat: latitude
    cth: cloud-top height (m)
    f: field that should be corrected
    do_masking: OPTIONAL, switch which determines if clear-sky and invalid
                are also placed onto the output field, if True only
                cloudy values are remapped.
    do_sequential_hole_filling: OPTIONAL, if True median values of an 
        sequentially increasing neighborhood are taken to fill holes.

    OUTPUT
    ======
    fpar: parallax-corrected field
    '''


    # do smoothing ---------------------------------------------------
    if smoothing or gauss_sig != 0:
        print '..use cth smoothing with sigma', gauss_sig
        cth_sm =  scipy.ndimage.percentile_filter(cth, 90, size = 4 * gauss_sig)
        cth_sm =  scipy.ndimage.gaussian_filter(cth_sm, gauss_sig)
    else:
        cth_sm = cth
    # ================================================================


    # make masking ---------------------------------------------------
    if do_masking:
        f = np.ma.masked_invalid(f)
        mnot = f.mask
        mf = np.logical_not(mnot)
    
        m = (mf) & ( cth_sm != 0)
    else:
        m = np.ones_like(cth_sm).astype(np.bool)
    # ================================================================

        

    # parallax correction --------------------------------------------
    plon, plat = parallax_correction(lon[m], lat[m], cth_sm[m], satlon = satlon)

    # and remapping
    xp, yp = radolan.rado_ll2xy(plon, plat)
    # ================================================================


    # create indices - not sure if round is right ... ----------------
    R = radolan.Radolan()
    xnew, ynew = R.xg, R.yg
    fac = 1.
    dr = 1.

    ic = np.round((xp - xnew.min()) / (fac*dr)).astype(np.int)
    ir = np.round((yp - ynew.min()) / (fac*dr)).astype(np.int)

    # check bounds
    Nr, Nc = xnew.shape
    ic = np.ma.masked_outside(ic, 0, Nc - 1)
    ir = np.ma.masked_outside(ir, 0, Nr - 1)

    # save additional mask
    madd = np.logical_not(ic.mask) & np.logical_not(ir.mask)
    ic = ic[madd]
    ir = ir[madd]
    # ================================================================


    # successive index-based interpolation ---------------------------
    Nsteps = 10
    h = np.linspace(0, cth_sm.max(), Nsteps + 1)  

    # initial field
    fpar = -9999 * np.ones_like( f )
    fpar[~m] = f[~m]


    # iteration
    for n in range(Nsteps):
        h_low = h[n]
        h_upp = h[n + 1]

        # make mask
        ma = np.logical_and(cth_sm[m][madd] >= h_low, cth_sm[m][madd] < h_upp)

        fpar[ir[ma],ic[ma]] = f[m][madd][ma]

    # ================================================================
    

    # do sequential hole filling -------------------------------------
    if do_sequential_hole_filling:
        fpar = sequential_hole_filling(fpar, filtertype='max')
    # ================================================================



    return np.ma.masked_where(fpar == -9999, fpar)

######################################################################
######################################################################

def parallax_correct_field_remap(lon, lat, cth, f, 
                                    dr = 0.3, fac = 3):

    '''
    Do parallax-correction of field f based on cth (m). f and cth are
    remapped onto an equi-distant (x,y) grid before.

    INPUT
    =====
    lon: longitude
    lat: latitude
    cth: cloud-top height (m)
    f: field that should be corrected
    dr: OPTIONAL grid distance of intermediate equi-distant (x,y) grid
    fac: factor such that fac*dr is grid distance of output grid

    OUTPUT
    ======
    xout: output x-coordinate
    yout: output y-coorfinate
    fpar: parallax-corrected field
    '''

    # first do remaping ----------------------------------------------
    xnew, ynew, fnew = gi.remap_field(lon, lat, f, dr = dr)
    xnew, ynew, cth_new = gi.remap_field(lon, lat, cth, dr = dr)

    xout, yout, fout = gi.remap_field(lon, lat, f, dr = fac*dr)
    xout, yout, cout = gi.remap_field(lon, lat, cth, dr = fac*dr)
    # ================================================================


    # parallax correction --------------------------------------------
    nlon, nlat = gi.xy2ll(xnew, ynew)
    plon, plat = parallax_correction(nlon, nlat, cth_new)
    xp, yp = gi.ll2xy(plon, plat)
    # ================================================================


    # create indices - not sure if round is right ... ----------------
    ic = np.round((xp - xnew.min()) / (fac*dr)).astype(np.int)
    ir = np.round((yp - ynew.min()) / (fac*dr)).astype(np.int)

    # check bounds
    Nr, Nc = xnew.shape
    ic = np.ma.masked_outside(ic, 0, Nc - 1)
    ir = np.ma.masked_outside(ir, 0, Nr - 1)

    # save additional mask
    madd = np.logical_not(ic.mask) & np.logical_not(ir.mask)
    ic = ic[madd]
    ir = ir[madd]
    # ================================================================


    # successive index-based interpolation ---------------------------
    Nsteps = 10
    h = np.linspace(0, cth.max(), Nsteps + 1)  

    # initial field
    fpar = np.where(cout == 0, fout, -9999)

    # iteration
    for n in range(Nsteps):
        h_low = h[n]
        h_upp = h[n + 1]

        # make mask
        ma = np.logical_and(cth_new[madd] > h_low, cth_new[madd] <= h_upp)

        fpar[ir[ma],ic[ma]] = fnew[madd][ma]

    # ================================================================

    
    # do sequential hole filling -------------------------------------
    fpar = sequential_hole_filling(fpar)
    # ================================================================

    return xout, yout, gi.curve_flow_filter(fpar,3) #scipy.ndimage.median_filter(fpar, (5,3))


######################################################################
######################################################################


def sequential_hole_filling(f, NaN = -9999, Nmax = 20, 
                            filtertype = 'median'):

    '''
    Fills holes that appear in parallax correction with
    median values of sequetially increasing neigborhoods.

    INPUT
    =====
    f: field with holes (to be filled)
    NaN: OPTIONAL, not-a-number value that is masked out
    Nmax: OPTIONAL, numberof iterations

    OUTPUT
    ======
    fout: output field with holes filled
    '''


    fout = copy.copy(f)

    print '... start with %d Nan values' % len(f[f == NaN])
    fmin = NaN
    nint = 1
    while fmin == NaN and nint < Nmax:

        print nint

        if filtertype == 'median':
            # calculate local median
            fmed = scipy.ndimage.median_filter(fout, 2 * nint + 1)

        elif filtertype == 'max':
            # calculate local median
            fmed = scipy.ndimage.maximum_filter(fout, 2 * nint + 1)


        # create mask
        ma = (fout == NaN)

        # update field with median values
        fout[ma] = fmed[ma]

        # calculate new minimum
        fmin = fout.min()

        nint += 1

    print '... end with %d Nan values' % len(fout[fout == NaN])
    
    return fout 



######################################################################
######################################################################


def extremely_simple_cth(bt108):

    return 12. / ( 210. - 300. ) * (bt108 - 300.)



######################################################################
######################################################################


def parallax_correction(lon, lat, height, satlon = 9.5):

    '''
    Calculates parallax-corrected longitude and latitude values. Adapted 
    from NWCSAF.

    INPUT
    =====
    lon: longitude
    lat: latitude
    height: cloud-top height (m)
    datlon: OPTIONAL, longitude position of satellite 

    OUTPUT
    ======
    plon: corrected longitude
    plat: corrected latitude
    '''


    # earth ellipsoid
    POLE_EARTH_RADIUS = 6356.58
    EQUATOR_EARTH_RADIUS = 6378.17
    aorbit=42165390

    aplr = POLE_EARTH_RADIUS * 1000
    aeqr = EQUATOR_EARTH_RADIUS * 1000
    aratio = aeqr / aplr

    # Convert angles from degrees to radians and height from ft to m
    # --------------------------------------------------------------
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)

    # Get x-y-z coordinates for satellite */
    # slat1 is the geometric latitude for geodetic lat satlat
    # ------------------------------------------------------- */
    satlon = np.deg2rad(satlon)
    satlat = 0
    slat1 = np.arctan( np.tan(satlat)*(aratio*aratio))
    
    xs = aorbit * np.cos(slat1) * np.sin(satlon)
    ys = aorbit * np.sin(slat1)
    zs = aorbit * np.cos(slat1) * np.cos(satlon)

    # Get x-y-z coordinates for surface point */
    # --------------------------------------- */
    lat1 = np.arctan( np.tan(lat) * (aratio*aratio) )
    
    ri = aeqr/np.sqrt((np.cos(lat1) * np.cos(lat1))\
                          +(aratio*aratio) * (np.sin(lat1)*np.sin(lat1)))

    xi = ri * np.cos(lat1) * np.sin(lon)
    yi = ri * np.sin(lat1)
    zi = ri * np.cos(lat1) * np.cos(lon)

    # b is the new aratio */
    # ------------------- */
    b = (((aeqr + height) / (aplr + height)) * ((aeqr + height) / (aplr + height)))

    xsmxi = xs - xi
    ysmyi = ys - yi
    zsmzi = zs - zi

    e = (xsmxi*xsmxi) + b*(ysmyi*ysmyi) + (zsmzi*zsmzi)
    ef = 2*(xi*xsmxi+b*yi*ysmyi+zi*zsmzi)
    eg = (xi*xi)+(zi*zi)+b*(yi*yi)-((aeqr+height)*(aeqr+height))
    a = (np.sqrt((ef*ef)-4*e*eg) - ef)/2/e

    # Corrected xyz */
    # ------------- */
    xc = xi + a * xsmxi
    yc = yi + a * ysmyi
    zc = zi + a * zsmzi

    # Convert back to lat/lon */
    # ----------------------- */
    aux = np.arctan(yc / np.sqrt((xc*xc)+(zc*zc)))

    lat_new  = np.arctan( np.tan(aux) / (aratio*aratio))
    lon_new = np.arctan2(xc,zc)

    return np.rad2deg(lon_new), np.rad2deg(lat_new)


######################################################################
######################################################################

if __name__ == '__main__':

    from l15_msevi.msevi import MSevi
    import MSGtools
    import plotting_tools.bmaps
    import pylab as pl

    t = datetime.datetime(2012,5,23,16,0)
 #   reg = ((200, 500), (1700, 2000))
    s = MSevi(time = t)
    s.lonlat()
    s.load('IR_108')
    s.rad2bt()
    # pl.imshow(s.bt['IR_108'])
    R = radolan.Radolan()
    R.read(t)
    cth = MSGtools.get_nwcsaf_prod('CTTH_HEIGHT', t, calibrate = True)
    cth [cth<0] = 0

    btp = parallax_correct_on_radolan(s.lon, s.lat, cth, s.bt['IR_108'],
                do_masking =False, do_sequential_hole_filling=True) 
    pl.figure()
    btp = np.where(np.isfinite(btp), btp,0)
    pl.pcolormesh(R.xg, R.yg,
            btp, vmin = 210, vmax = 300)

    pl.contour(R.xg, R.yg, R.mask(), linewidths = 3, colors = 'k')

    pl.figure()
    (x,y), mp = plotting_tools.bmaps.make_map(s.lon, s.lat, region='germ')
    mp.pcolormesh(x,y,s.bt['IR_108'],vmin = 210, vmax = 300)

    xr, yr = mp(R.lon, R.lat)
    mp.contour(xr, yr, R.mask(), linewidths = 3, colors = 'k')
