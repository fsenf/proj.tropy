#!/usr/bin/env python



######################################################################
######################################################################
# 
# Description
# ===========
#
# Package for the calculation of several thermodynamical properties
# of moist air (including condensate).


# ToDo:
# ====
#
# * Some final tests of the routines should be done. 
# * There was an issue with mixing ratio vs. specific humidity, 
#   but might be resolved ... 
#

######################################################################
######################################################################



import sys
import numpy as np


######################################################################
######################################################################

def thermodynamic_constants():
    
    '''
    Description:
    ============
    Sets a lot of thermodynamic constants.
    
    '''

    # gas constants from Davies-Jones, 2009, MWR, 137, p. 3137

    const = {}
    const['R_d'] = 287.04 # specific gas constant of dry air [J kg-1 K-1]
    const['R_v'] = 461.5 # specific gas constant of water vapour [J kg-1 K-1]
    const['eps'] = const['R_d'] / const['R_v']
    const['c_pd'] = 1005.7 # Specific heat at constant pressure of dry air [J kg-1 K-1]
    const['c_pv'] = 1875.  # Specific heat at constant pressure of water vapor [J kg-1 K-1]
    const['kappa_d'] = const['R_d'] / const['c_pd']
    
    return const
    


######################################################################
######################################################################

def moist_gas_constant(r):
    
    '''
    Description:
    ===========
    Calculates gas constant of moist air.


    Arguments:
    =========
    r: mixing ratio i.e. vapor mass per dry air mass in [kg / kg]

    
    Return:
    =======
    R_m: gas constant of moist air [J / kg / K]

    '''

    # get constants
    c = thermodynamic_constants()
    R_d, eps = c['R_d'], c['eps']

    R_m = R_d * ( 1 + r / eps) / (1 + r)

    return R_m

######################################################################
######################################################################

def moist_specific_heat(r):
    
    '''
    Description:
    ===========
    Calculates specific heat at constant pressure of moist air.


    Arguments:
    =========
    r: mixing ratio i.e. vapor mass per dry air mass in [kg / kg]

    
    Return:
    =======
    c_pm: specific heat at constant pressure of moist air

    '''

    # get constants
    c = thermodynamic_constants()
    c_pd, c_pv = c['c_pd'], c['c_pv']

    c_pm = ( c_pd + r * c_pv) / (1 + r)

    return c_pm 

######################################################################
######################################################################

def dry_air_density(r, p, T):

    '''
    Description:
    ===========
    Calculates dry air density.


    Arguments:
    =========
    r: mixing ratio i.e. vapor mass per dry air mass in [kg / kg]
    p: total gas pressure [Pa]
    T: temperature [K]

    
    Return:
    =======
    rho_d: dry air density in [kg / m**3]

    '''

    # get constants
    c = thermodynamic_constants()
    R_d, R_v = c['R_d'], c['R_v']


    # calculate dry air density
    rho_d = p / ( R_d * T * ( 1 + r * R_v / R_d) )

    return rho_d

######################################################################
######################################################################

def total_density(r, r_w, p, T):

    '''
    Description:
    ===========
    Calculates total density of moist air with hydrometeors.


    Arguments:
    =========
    r: mixing ratio i.e. vapor mass per dry air mass in [kg / kg]
    r: mixing ratio of condensed water i.e. condensed water mass per dry air mass in [kg / kg]
    p: total gas pressure [Pa]
    T: temperature [K]

    
    Return:
    =======
    rho_d: dry air density in [kg / m**3]

    '''

    # get dry air density
    rho_d =  dry_air_density(r, p, T)

    rho = rho_d * ( 1. + r + r_w)

    return rho

######################################################################
######################################################################

def absolute_humidity(r, p, T):

    '''
    Description:
    ===========
    Calculates absolute humidity.


    Arguments:
    =========
    r: mixing ratio i.e. vapor mass per dry air mass in [kg / kg]
    p: total gas pressure [Pa]
    T: temperature [K]

    
    Return:
    =======
    rho_v: absolute humidity / m**3]

    '''


    # get dry air density
    rho_d = dry_air_density(r, p, T)

    rho_v = r * rho_d

    return rho_v

######################################################################
######################################################################


def specific_humidity(r):

    '''
    Description:
    ===========
    Calculates specific humidity.


    Arguments:
    =========
    r: mixing ratio i.e. vapor mass per dry air mass in [kg / kg]

    
    Return:
    =======
    q: specific humidity i.e. vapor mass per moist air mass [kg / kg] 

    '''


    return r / ( 1. + r )

######################################################################
######################################################################



def water_vapor_pressure(r, p, T):

    '''
    Description:
    ===========
    Calculates water vapor pressure.


    Arguments:
    =========
    r: mixing ratio i.e. vapor mass per dry air mass in [kg / kg]
    p: total gas pressure [Pa]
    T: temperature [K]

    
    Return:
    =======
    e : water vapor pressure in [Pa]

    '''

    # get constants
    c = thermodynamic_constants()
    R_v = c['R_v']


    rho_d = dry_air_density(r, p, T)

    p_v = r * rho_d * R_v * T


    return p_v
    

######################################################################
######################################################################

def saturation_pressure(T):

    '''
    Description:
    ===========
    Calculates saturation water vapor pressure 
    after Bolton (1980) MWR 108, p.1046, eq. (10)


    Argument:
    =========
    T: temperature [K]

    
    Return:
    =======
    es : water vapor pressure in [Pa]

    '''
    
    # in degree celcius
    T = T - 273.15

    # Markowski book p 13 eq (2.16)
    es = 6.112 * np.exp(17.67 * T / (T + 243.5) )
    
    return es * 1e2 # unit conversion in Pa

######################################################################
######################################################################


def saturation_over_ice(T, method = 'PK'):
    
    '''
    Description:
    ===========
    Calculates saturation water vapor pressure over ice.


    Argument:
    =========
    T: temperature [K]

    
    Return:
    =======
    es : water vapor pressure in [Pa]

    '''

    if method == 'PK':
        # constants
        e0 = 610.64
        T0 = 273.15
        Ai = 21.88
        Bi = 7.65
        
        es = e0 * np.exp( Ai * (T - T0) / ( T - Bi))

    elif method == 'Murphy':
        es = np.exp(9.550426 - 5723.265 / T + 3.53068 * np.log(T) - 0.00728332 * T)

    return es


######################################################################
######################################################################


def relative_humidity(r, p, T):

    '''
    Description:
    ===========
    Calculates relative humidity.


    Arguments:
    =========
    r: mixing ratio i.e. vapor mass per dry air mass in [kg / kg]
    p: total gas pressure [Pa]
    T: temperature [K]

    
    Return:
    =======
    H: relative humidity in [%]

    '''



    e = water_vapor_pressure(r, p, T)
    es = saturation_pressure(T)

    H = e / es

    return H * 100.
    

######################################################################
######################################################################


def dew_point(r, p, T):
    '''
    Description:
    ===========
    Calculates dew point temperature
    after Markowski book p 13 eq (2.25)


    Arguments:
    ==========
    r: mixing ratio i.e. vapor mass per dry air mass in [kg / kg]
    p: total gas pressure [Pa]
    T: temperature [K]


    
    Return:
    =======
    Td: dew point temperature in [K]

    '''
    
    # get vapor pressure
    e = water_vapor_pressure(r, p, T) * 1e-2 # unit conversion to hPa


    # Markowski book p 13 eq (2.25)
    Td = 243.5 / (17.67 / np.log(e / 6.112) - 1.)
    
    return Td + 273.15

######################################################################
######################################################################


def H2r(hrel, p, T):
    '''
    Description:
    ===========
    Converts relative humidity H into mixing ratio r.


    Argument:
    =========
    hrel: relative humidity in [%]
    p: total gas pressure [Pa]
    T: temperature [K]

    
    Return:
    =======
    es : water vapor pressure in [Pa]

    '''

    # get constants
    c = thermodynamic_constants()
    eps = c['eps']

    # unit conversion
    r = hrel / 100.

    # saturation pressure
    es = saturation_pressure(T)

    # vapor pressure
    e = r * es

    # mixing ratio
    r = eps * e / (p  -  e)

    return r

######################################################################
######################################################################


def watermass_fraction2rw(qw, r):

    '''
    Description:
    ===========
    Converts water mass fraction to condensed water mixing ratio.


    Argument:
    =========
    qw: condensed water mass fraction [kg / kg]
    r: mixing ration of water vapor [kg / kg]
    p: total gas pressure [Pa]
    T: temperature [K]

    
    Return:
    =======
    rw: condensed water mixing ratio  [kg / kg]

    '''


    rw = qw / (1. - qw) * (1 + r) 

    return rw

######################################################################
######################################################################

def both_mass_fractions2mixing_ratios(qv, qw):

    '''
    Description:
    ===========
    Converts water and the vapor mass fraction to 
    condensed water and vapor mixing ratio.


    Argument:
    =========
    qv: water vapor mass fraction [kg / kg]
    qw: condensed water mass fraction [kg / kg]

    
    Return:
    =======
    r:  water vapor mixing ratio  [kg / kg]
    rw: condensed water mixing ratio  [kg / kg]

    '''
    
    r = qv / ( 1 - qv - qw)

    rw = r * qw / qv

    return r, rw

######################################################################
######################################################################



def lifting_condensation_level_temperature(r, p, T, eq = 15):

    '''
    Description:
    ===========
    Calculates lifting condensation level temperature.
    after Bolton (1980) MWR 108, p.1046, eq. (15), (21) or (22)


    Arguments:
    ==========
    r: mixing ratio i.e. vapor mass per dry air mass in [kg / kg]
    p: total gas pressure [Pa]
    T: temperature [K]


    
    Return:
    =======
    TL: lifting condensation level temperature in [K]

    '''
    
    
    if eq == 15:
        
        TD = dew_point(r, p, T)

        TL = 1./ ( 1 / ( TD - 56.) + np.log(T / TD) / 800.) + 56.

    elif eq == 21:
        
        e = water_vapor_pressure(r, p, T) * 1e-2 # unit conversion to hPa
        
        TL = 2840. / ( 3.5 * np.log(T) - np.log(e) - 4.805) + 55.

    elif eq == 22:
        hrel = relative_humidity(r, p, T)
        
        TL = 1./ ( 1 / ( T - 55.) - np.log( hrel / 100.) / 2840.) + 55.
    else:
        print 'wrong equation number'

    return TL

######################################################################
######################################################################
    
def dry_air_potential_temperature(r, p, T):

    '''
    Description:
    ===========
    Calculates dry air potential temperature.


    Arguments:
    ==========
    r: mixing ratio i.e. vapor mass per dry air mass in [kg / kg]
    p: total gas pressure [Pa]
    T: temperature [K]


    
    Return:
    =======
    theta_d: dry air potential temperature (K)
    
    '''

    # get constants
    c = thermodynamic_constants()
    kappa_d = c['kappa_d']

    e = water_vapor_pressure(r, p, T)
    
    theta_d = T * (1000e2 / (p - e)) ** kappa_d
     
    
    return theta_d

######################################################################
######################################################################
    
def moist_potential_temperature(r, p, T):

    '''
    Description:
    ===========
    Calculates moist air potential temperature.


    Arguments:
    ==========
    r: mixing ratio i.e. vapor mass per dry air mass in [kg / kg]
    p: total gas pressure [Pa]
    T: temperature [K]


    
    Return:
    =======
    theta_m: moist air potential temperature (K)
    
    '''

    # get exponent
    R_m = moist_gas_constant(r)
    c_pm = moist_specific_heat(r)
    kappa_m = R_m / c_pm

    
    theta_m = T * (1000e2 / p) ** kappa_m
     
    
    return theta_m

######################################################################
######################################################################


def equivalent_potential_temperature(r, p, T, eq = 15):

    '''
    Description:
    ===========
    Calculates equivalent potential temperature.
    after Bolton (1980) MWR 108, p.1046, eq. (43) 

    Davies-Jones (2009) MWR137, eq. (6.3)

    also using either one of (15), (21) or (22)


    Arguments:
    ==========
    r: mixing ratio i.e. vapor mass per dry air mass in [kg / kg]
    p: total gas pressure [Pa]
    T: temperature [K]


    
    Return:
    =======
    theta_E: 

    '''

    # get constants
    c = thermodynamic_constants()
    c_pd = c['c_pd']

    
    TL = lifting_condensation_level_temperature(r, p, T, eq = eq)
    theta = moist_potential_temperature(r,p,T)


    # Theta_E = (
    #     T * (1000e2 / p) ** (0.2854* ( 1. - 0.28e-3 * r)) * 
    #     np.exp((3.376 / TL - 0.00254) * r * (1. + 0.81e-3 * r))
    #     )

    
    L0 = 2.711e6
    L1 = 1109
    C = 273.15

    kern = ( L0 - L1 * (TL - C) ) * r / (  c_pd * TL ) 

    theta_E = theta * np.exp(kern)

    return theta_E

######################################################################
######################################################################
    


if __name__ == '__main__':
    
    p = 1000e2
    T = 30 + 273.15
    r = hrel2qv(100., p, T)

    print dry_air_density(r, p, T)
    print absolute_humidity(r, p, T)
    print water_vapor_pressure(r, p, T)
    print saturation_pressure(T)
    print dew_point(r, p, T)
    hrel = relative_humidity(r, p, T)
    print hrel2qv(hrel, p, T)
    
    print lifting_condensation_level_temperature(r, p, T, eq = 15)
    print lifting_condensation_level_temperature(r, p, T, eq = 21)
    print lifting_condensation_level_temperature(r, p, T, eq = 22)
    print equivalent_potential_temperature(r, p, T)

