#!/usr/bin/env python

# load libraries -----------------------------------------------------
import sys, os

import numpy as np
import datetime
import cv2
import scipy.ndimage

import analysis_tools.grid_and_interpolation as gi


######################################################################
######################################################################
#
#  === HERE THE OPTICAL FLOW ROUTINES START ===
#
######################################################################
######################################################################

def displacement_from_opt_flow(f1, f2, 
                               vmin = 0, 
                               vmax = 1.,
                               nlayers = 5,
                               boxsize = 15):

    '''
    Derives displacement vector between to fields f1 and f2 
    using optical flow method.

    INPUT
    =====
    f1, f2: two 2d fields
    vmin: optional, lower bound of fields for rescaling
    vmax: optional, upper bound of fields for rescaling
    nlayers: optional, number of layers for pyramidal coarse-graining
    boxsize: optional, size of avering box from which flow is estimated
    
    

    OUTPUT
    ======
    flow: displacement vector, 1st component is related to u, i.e. 
          displacement along row, 2nd to v, i.e. along column

    '''

    
    # field rescaling ------------------------------------------------
    fr1 = (f1 - vmin) / (vmax - vmin)
    fr1 = np.clip( fr1, 0., 1.)

    fr2 = (f2 - vmin) / (vmax - vmin)
    fr2 = np.clip( fr2, 0., 1.)
    # ================================================================



    # convert fields to images and applied optical flow --------------
    im1 = (fr1*255).astype(np.uint8)
    im2 = (fr2*255).astype(np.uint8)

    flow = cv2.calcOpticalFlowFarneback(im1,im2,  0.5, nlayers, boxsize, 3, 5, 1.2, 0)
    # ================================================================

    return flow


######################################################################
######################################################################


def morph_trans_opt_flow(f, flow, method = 'forward'):
    
    '''
    Applies morphological transformation of field f given a displacement 
    field.
    
    INPUT
    =====
    f: 2d field that is transformed (source)
    flow: displacement vector between source and target stage
    

    OUTPUT
    ======
    ftrans: transformed field
    
    '''


    # get shape of field and corresponding index set -----------------
    nrows, ncols = f.shape
    irow, icol = gi.make_index_set(nrows, ncols)
    # ================================================================
    

    # get index shift from real-values flow field --------------------

    ishift = np.round(flow).astype(np.int)

    # care that row and column are transposed here!!!!
    ishift_row = ishift[:, :, 1]
    ishift_col = ishift[:, :, 0]
    
    if method == 'forward':
        ir = irow - ishift_row
        ic = icol - ishift_col
    elif method == 'backward':
        ir = irow + ishift_row
        ic = icol + ishift_col
        

    ir = np.clip(ir, 0, nrows - 1)
    ic = np.clip(ic, 0, ncols - 1)
    # ================================================================

    return f[ir, ic]

######################################################################
######################################################################


def flow_velocity(lon, lat, f3d, 
                  dt = 5.,
                  flow_input = None,
                  vmin = 0, 
                  vmax = 1.):

    '''
    Calculates flow field for Lagrangian displacements of tracers
    in the field.

    INPUT
    =====
    lon: longitude
    lat: latidute
    f3d: 3d fields (time on 1st axis)
    dt: optional, time interval in minutes
    vmin: optional, lower bound of fields for rescaling
    vmax: optional, upper bound of fields for rescaling


    OUTPUT
    ======
    u: zonal velocity, 
    v: meridional velocity 


    COMMENTS
    ========
    Velocities are derived from displacements between subsequent time slots
    and hence, time dimension is ntimes - 1.

    '''


    # # field dimensions -----------------------------------------------
    ntimes, nrows, ncols = f3d.shape
    # # ================================================================

    
    # use local Cartesian coordinates on sphere ----------------------
    x, y = gi.ll2xy(lon, lat)

    x3d = np.expand_dims(x, 0).repeat(ntimes, axis = 0)
    y3d = np.expand_dims(y, 0).repeat(ntimes, axis = 0)
    # ================================================================


    dx = Lagrangian_change(x3d, 
                           flow_input = flow_input,
                           tracer_field = f3d, 
                           vmin = vmin, 
                           vmax = vmax)

    dy = Lagrangian_change(y3d, 
                           flow_input = flow_input,
                           tracer_field = f3d, 
                           vmin = vmin, 
                           vmax = vmax)

    return dx * 1000. / (dt * 60.), dy  * 1000. / (dt * 60.)


######################################################################
######################################################################


def displacement_vector(tracer_field,
                      vmin = 0, 
                      vmax = 1.):

    '''
    Calculates displacement vector for several times.

    INPUT
    =====
    tracer_field: 3d fields (time on 1st axis)


    OUTPUT
    ======
    flow: displacement vector

    '''



    # field dimensions -----------------------------------------------
    ntimes, nrows, ncols = tracer_field.shape

    flow = np.zeros((ntimes - 1, nrows, ncols, 2))
    # ================================================================

    


    # loop over time -------------------------------------------------
    for i in range(ntimes - 1):

        f1 = tracer_field[i]       # first 
        f2 = tracer_field[i + 1]   # next

        flow[i] = displacement_from_opt_flow(f1, f2, 
                                          vmin = vmin, 
                                          vmax = vmax)

    # ================================================================


    return flow

######################################################################
######################################################################


def Lagrangian_change(f3d, 
                      tracer_field = None, 
                      gauss_sigma = 1.,
                      flow_input = None,
                      vmin = 0, 
                      vmax = 1.):

    '''
    Calculates Lagrangian change of a field using possibly another
    field to generate optical flow.

    INPUT
    =====
    f3d: 3d fields (time on 1st axis)


    OUTPUT
    ======
    df: Lagrangian change of field f


    COMMENTS
    ========
    Lagrangian changes are derived from displacements between subsequent time slots
    and hence, time dimension is ntimes - 1.

    '''



    # field dimensions -----------------------------------------------
    ntimes, nrows, ncols = f3d.shape

    df = np.zeros((ntimes - 1, nrows, ncols))
    # ================================================================

    
    # set tracer field -----------------------------------------------
    if tracer_field == None:
        tracer_field = f3d
    # ================================================================


    # loop over time -------------------------------------------------
    for i in range(ntimes - 1):

        f1 = f3d[i]       # first 
        f2 = f3d[i + 1]   # next

        if flow_input == None:
            flow = displacement_from_opt_flow(tracer_field[i], 
                                              tracer_field[i + 1], 
                                              vmin = vmin, 
                                              vmax = vmax)
        else:
            flow = flow_input[i]


        # the first is shifted by advection
        ft = morph_trans_opt_flow(f1, flow) 

        # the shifted field is compared to next
        df[i] = (f2 - ft) 
    # ================================================================


    if gauss_sigma > 0:
        df = scipy.ndimage.gaussian_filter(df, (0, gauss_sigma, gauss_sigma) )

    return df


######################################################################
######################################################################
