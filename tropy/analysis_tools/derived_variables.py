#!/usr/bin/python


'''
This module contains / will contain variables
derived from some basic input fields.
The name "basic" refers to the typical set of
variables stored in climate or weather model
output.
'''


import sys, os

#import tropy.ECMWFtools as ECMWFtools
import numpy as np
from .grid_and_interpolation import ldiff, lmean

######################################################################
######################################################################


def calc_lwp(p, T, qv,  qc, g = 9.8, axis = -1):

    '''
    Calculates the liquid water path.


    Parameters
    ----------
    p : numpy array
        atmospheric pressure [hPa]

    T : numpy array
        temperature [K]

    qv : numpy array
        water vapor mixing ratio [kg / kg]

    qc : numpy array
        liquid water mixing ratio [kg / kg]

    
    Returns
    --------
    lwp : numpy array
        liquid water path [g / m**2]


    Notes
    ------
    * hydrostatic approximation is used to convert vertical integration from
      pressure to height levels
    * impact of condesate mass on air density is ignored.
    '''
    
    # pressure units conversion (hPa --> Pa) !!!!

    if p.max() < 1200.:
        p = p * 1e2

    # get ideal gas constants
    RL = 287.    # gas constant of dry air
    RV = 461.5   # gas constant of water vapor
    
    R = RL + (RV - RL) * qv  # gas constant of humid air

    # approximated gas density
    rho = p / (R * T)
    
    # liquid cloud water content
    lwc = rho * qc

    # pressure level distance
    dp = ldiff(p, axis = axis)

    # height level distance 
    dz = dp / (g * lmean(rho, axis = axis))

#    print dz.min(),dz.max()

    # vertical integral of LWC
    lwp = (lmean(lwc, axis = axis) * dz).sum(axis = axis)

    return lwp * 1e3


######################################################################
######################################################################
