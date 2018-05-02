#!/usr/bin/python

import sys, os

import ECMWFtools
import numpy as np
from analysis_tools.grid_and_interpolation import ldiff, lmean

######################################################################
######################################################################


def calc_lwp(p, T, qv,  qc, g = 9.8, axis = -1):

    '''
    Calculates the liquid water path.


    USAGE
    =====
    lwp = calc_lwp(p, T, qv,  qc, g = 9.8, axis = -1)


    INPUT
    =====
    p: atmospheric pressure [hPa]
    T: temperature [K]
    qv: water vapor mixing ratio [kg / kg]
    qc: liquid water mixing ratio [kg / kg]

    
    OUTPUT
    ======
    lwp: liquid water path [g / m**2]

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
