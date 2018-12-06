#!/usr/bin/env python

# load libraries -----------------------------------------------------
import sys, os

import numpy as np
import datetime
import cv2
import scipy.ndimage

import grid_and_interpolation as gi
import statistics as stats


######################################################################
######################################################################
#
# Set default OFLOW parameters
#
######################################################################
######################################################################

_farneback_default_parameters = {'flow': None,
                                 'pyr_scale':0.5, 
                                 'levels': 5,
                                 'win_size':15,
                                 'iterations':3,
                                 'poly_n':5,
                                 'poly_sigma':1.2,
                                 'flags':0 }


_tvl1_default_parameters = {'flow': None,
                            'epsilon':0.01,
                            'lambda':0.2,
                            'outer_iterations':20,#40,
                            'inner_iterations':5,#7,
                            'gamma':0.4,
                            'scales_number':3,#5,
                            'tau':0.25,
                            'theta':0.8,
                            'warpings_number':3,#5,
                            'scale_step':0.5,
                            'median_filtering':1,
                            'use_initial_flow':0}

######################################################################
######################################################################
#
#  === HERE THE OPTICAL FLOW ROUTINES START ===
#
######################################################################
######################################################################

def displacement_from_opt_flow(f1, f2, 
                               method = 'farneback',
                               **kwargs):

    '''
    Derives displacement vector between to fields f1 and f2 
    using optical flow method (either farnebaeck, or tvl1).


    Parameters
    ----------

    f1 : numpy array, 2dim
         field for past time step

    f2 : numpy array, 2dim
         field for actual time step

    method : str, optional, default = 'farneback'
         method selected for optical flow calculations
         possible options: 'farneback' and 'tvl1'


    kwargs : dict
        parameters for opencv optical flow algorithm 
        if not given, default is taken from  _farnebaeck_default_parameters
    
    

    Returns
    --------
    flow : numpy array, 3dim
        displacement vector as index shift
        1st component is related to u, i.e. displacement along row, 
        2nd to v, i.e. along column
    '''


    if method == 'farneback':
        return  displacement_from_opt_flow_farneback(f1, f2, 
                                                     **kwargs)
    elif method == 'tvl1':
        return  displacement_from_opt_flow_tvl1(f1, f2, 
                                                **kwargs)
    else:
        raise ValueError('Unknown method,  possible options: "farneback" and "tvl1"')
    


######################################################################
######################################################################


def displacement_from_opt_flow_farneback(f1, f2, 
                                         vmin = None, 
                                         vmax = None,
                                         **kwargs):

    '''
    Derives displacement vector between to fields f1 and f2 
    using optical flow method after Farneback (2003).


    Parameters
    ----------

    f1 : numpy array, 2dim
         field for past time step

    f2 : numpy array, 2dim
         field for actual time step

    vmin : float, optional, default = None
        lower bound of fields for rescaling

    vmax : float, optional, default = None
        upper bound of fields for rescaling

    kwargs : dict
        parameters for opencv optical flow algorithm 
        if not given, default is taken from  _farnebaeck_default_parameters
    

    Returns
    --------
    flow : numpy array, 3dim
        displacement vector as index shift
        1st component is related to u, i.e. displacement along row, 
        2nd to v, i.e. along column
    '''


    # get parameters from keywords -----------------------------------

    # get default
    flow_parameters =  _farneback_default_parameters.copy()
    
    # and overwrite with input paras
    flow_parameters.update( kwargs )
    # ================================================================


    
    # field rescaling ------------------------------------------------
    fr1 = stats.normalize_field( f1, vmin = vmin, vmax = vmax)
    fr2 = stats.normalize_field( f2, vmin = vmin, vmax = vmax)
    # ================================================================



    # convert fields to images and applied optical flow --------------
    im1 = (fr1*255).astype(np.uint8)
    im2 = (fr2*255).astype(np.uint8)

    flow = cv2.calcOpticalFlowFarneback(im1, im2,  
                                        flow_parameters['flow'],
                                        flow_parameters['pyr_scale'],
                                        flow_parameters['levels'],
                                        flow_parameters['win_size'],
                                        flow_parameters['iterations'],
                                        flow_parameters['poly_n'],
                                        flow_parameters['poly_sigma'],
                                        flow_parameters['flags'])

    # ================================================================

    return flow


######################################################################
######################################################################

def displacement_from_opt_flow_tvl1(f1, f2, 
                                    vmin = None, 
                                    vmax = None,
                                    **kwargs):

    '''
    Derives displacement vector between to fields f1 and f2 
    using optical flow method after  Zach et al (2007).


    Parameters
    ----------

    f1 : numpy array, 2dim
         field for past time step

    f2 : numpy array, 2dim
         field for actual time step

    vmin : float, optional, default = None
        lower bound of fields for rescaling

    vmax : float, optional, default = None
        upper bound of fields for rescaling

    kwargs : dict
        parameters for opencv optical flow algorithm 
        if not given, default is taken from  _farnebaeck_default_parameters
    

    Returns
    --------
    flow : numpy array, 3dim
        displacement vector as index shift
        1st component is related to u, i.e. displacement along row, 
        2nd to v, i.e. along column
    '''


    # get parameters from keywords -----------------------------------

    # get default
    flow_parameters =  _tvl1_default_parameters.copy()
    
    # and overwrite with input paras
    flow_parameters.update( kwargs )
    # ================================================================


    
    # field rescaling ------------------------------------------------
    fr1 = stats.normalize_field( f1, vmin = vmin, vmax = vmax)
    fr2 = stats.normalize_field( f2, vmin = vmin, vmax = vmax)
    # ================================================================



    # convert fields to images and applied optical flow --------------
    im1 = (fr1*255).astype(np.uint8)
    im2 = (fr2*255).astype(np.uint8)


    optflow=cv2.createOptFlow_DualTVL1()
    
    optflow.setEpsilon( flow_parameters['epsilon'] )
    optflow.setLambda( flow_parameters['lambda'] )
    optflow.setOuterIterations( flow_parameters['outer_iterations'] )
    optflow.setInnerIterations( flow_parameters['inner_iterations'] )
    optflow.setGamma( flow_parameters['gamma'] )
    optflow.setScalesNumber( flow_parameters['scales_number'] )
    optflow.setTau( flow_parameters['tau'] )
    optflow.setTheta( flow_parameters['theta'] )
    optflow.setWarpingsNumber( flow_parameters['warpings_number'] )
    optflow.setScaleStep( flow_parameters['scale_step'] )
    optflow.setMedianFiltering( flow_parameters['median_filtering'] )
    optflow.setUseInitialFlow( flow_parameters['use_initial_flow'] )

    # calculate flow
    flow = optflow.calc( im1, im2, 
                         flow_parameters['flow'] )
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
