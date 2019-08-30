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
    Sets some thermodynamic constants.


    Parameters
    ----------
    None

    
    Returns
    --------
    const : dict
        set of thermodynamic constans saved in dictionary
    
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
    Calculates gas constant of moist air.


    Parameters
    ----------
    r :  float or numpy array
        mixing ratio i.e. vapor mass per dry air mass in [kg / kg]

    
    Returns
    --------
    R_m : float or numpy array
        gas constant of moist air [J / kg / K]
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
    Calculates specific heat at constant pressure of moist air.


    Parameters
    ----------
    r : numpy array
        mixing ratio i.e. vapor mass per dry air mass in [kg / kg]

    
    Returns
    --------
    c_pm : numpy array
        specific heat at constant pressure of moist air

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
    Calculates dry air density.


    Parameters
    ----------
    r : numpy array
        mixing ratio i.e. vapor mass per dry air mass in [kg / kg]

    p : numpy array
        total gas pressure [Pa]

    T : numpy array
        temperature [K]

    
    Returns
    --------
    rho_d : numpy array
        dry air density in [kg / m**3]

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
    Calculates total density of moist air with hydrometeors.

    Parameters
    ----------
    r : numpy array
        mixing ratio i.e. vapor mass per dry air mass in [kg / kg]

    r_w : numpy array
        mixing ratio of condensed water i.e. condensed water mass per dry air mass in [kg / kg]

    p : numpy array
        total gas pressure [Pa]

    T : numpy array
        temperature [K]

    
    Returns
    --------
    rho_d : numpy array
        dry air density in [kg / m**3]

    '''

    # get dry air density
    rho_d =  dry_air_density(r, p, T)

    rho = rho_d * ( 1. + r + r_w)

    return rho

######################################################################
######################################################################

def absolute_humidity(r, p, T):

    '''
    Calculates absolute humidity.


    Parameters
    ----------
    r : numpy array
        mixing ratio i.e. vapor mass per dry air mass in [kg / kg]

    p : numpy array
        total gas pressure [Pa]

    T : numpy array
        temperature [K]

    Returns
    --------
    rho_v : numpy array
        absolute humidity / m**3]

    '''


    # get dry air density
    rho_d = dry_air_density(r, p, T)

    rho_v = r * rho_d

    return rho_v

######################################################################
######################################################################


def specific_humidity(r):

    '''
    Calculates specific humidity.


    Parameters
    ----------
    r : numpy array
        mixing ratio i.e. vapor mass per dry air mass in [kg / kg]

    
    Returns
    --------
    q : numpy array
        specific humidity i.e. vapor mass per moist air mass [kg / kg] 

    '''


    return r / ( 1. + r )

######################################################################
######################################################################



def water_vapor_pressure(r, p, T):

    '''
    Calculates water vapor pressure.


    Parameters
    ----------
    r : numpy array
        mixing ratio i.e. vapor mass per dry air mass in [kg / kg]

    p : numpy array
        total gas pressure [Pa]

    T : numpy array
        temperature [K]

    
    Returns
    --------
    e : numpy array
        water vapor pressure in [Pa]

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
    Calculates saturation water vapor pressure 
    after Bolton (1980) MWR 108, p.1046, eq. (10)


    Parameters
    ----------
    T : numpy array
        temperature [K]

    
    Returns
    --------
    es : numpy array
         water vapor pressure in [Pa]

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
    Calculates saturation water vapor pressure over ice.


    Parameters
    ----------
    T : numpy array
        temperature [K]

    method : str, optional, default = 'PK'
        use an equation which approximates the vapor pressure over ice
        method in ['PK', 'Murphy']
    
    Returns
    --------
    es : numpy array
        water vapor pressure in [Pa]


    Notes
    ------
    'PK' refers to Pruppacher and Klett (1997)
    'Murphy' refers to Murphy and Koop (2005) eq. (7)


    References
    ----------
    Murphy, D. M., and T. Koop (2005), Review of the vapour pressures of ice and supercooled water for atmospheric applications, Quart. J. Roy. Meteor. Soc., 131(608), 1539-1565, doi:10.1256/qj.04.94.
    

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
    Calculates relative humidity.

    
    Parameters
    ----------
    r : numpy array
        mixing ratio i.e. vapor mass per dry air mass in [kg / kg]

    p : numpy array
        total gas pressure [Pa]

    T : numpy array
        temperature [K]

    
    Returns
    --------
    H : numpy array
        relative humidity in [%]

    '''



    e = water_vapor_pressure(r, p, T)
    es = saturation_pressure(T)

    H = e / es

    return H * 100.
    

######################################################################
######################################################################


def dew_point(r, p, T):

    '''
    Calculates dew point temperature
    after Markowski book p 13 eq (2.25)


    Parameters
    ----------
    r : numpy array
        mixing ratio i.e. vapor mass per dry air mass in [kg / kg]

    p : numpy array
        total gas pressure [Pa]

    T : numpy array
        temperature [K]


    Returns
    --------
    Td : numpy array
        dew point temperature in [K]

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
    Converts relative humidity H into mixing ratio r.


    Parameters
    ----------
    hrel : numpy array
        relative humidity in [%]

    p : numpy array
        total gas pressure [Pa]

    T : numpy array
        temperature [K]

    
    Returns
    --------
    es : numpy array
        water vapor pressure in [Pa]

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
    Converts water mass fraction to condensed water mixing ratio.


    Parameters
    ----------
    qw : numpy array
        condensed water mass fraction [kg / kg]

    r : numpy array
        mixing ration of water vapor [kg / kg]

    p : numpy array
        total gas pressure [Pa]

    T : numpy array
        temperature [K]

    
    Returns
    --------
    rw : numpy array
        condensed water mixing ratio  [kg / kg]

    '''


    rw = qw / (1. - qw) * (1 + r) 

    return rw

######################################################################
######################################################################

def both_mass_fractions2mixing_ratios(qv, qw):

    '''
    Converts water and the vapor mass fraction to 
    condensed water and vapor mixing ratio.

    
    Parameters
    ----------
    qv : numpy array
        water vapor mass fraction [kg / kg]

    qw : numpy array
        condensed water mass fraction [kg / kg]

    
    Returns
    --------
    r :  numpy array
        water vapor mixing ratio  [kg / kg]

    rw : numpy array
        condensed water mixing ratio  [kg / kg]

    '''
    
    r = qv / ( 1 - qv - qw)

    rw = r * qw / qv

    return r, rw

######################################################################
######################################################################



def lifting_condensation_level_temperature(r, p, T, eq = 15):

    '''
    Calculates lifting condensation level temperature.
    after Bolton (1980) MWR 108, p.1046, eq. (15), (21) or (22)

    
    Parameters
    ----------
    r : numpy array
        mixing ratio i.e. vapor mass per dry air mass in [kg / kg]

    p : numpy array
        total gas pressure [Pa]

    T : numpy array
        temperature [K]

    eq : int, optional, default = 15
        which eq. of Bolton (1980) to be chosen 
        eq in [15, 21, 22]

    
    Returns
    --------
    TL : numpy array
        lifting condensation level temperature in [K]

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
        print('wrong equation number')

    return TL

######################################################################
######################################################################
    
def dry_air_potential_temperature(r, p, T):

    '''
    Calculates dry air potential temperature.


    Parameters
    ----------
    r : numpy array
        mixing ratio i.e. vapor mass per dry air mass in [kg / kg]

    p : numpy array
        total gas pressure [Pa]

    T : numpy array
        temperature [K]

    
    Returns
    --------
    theta_d : numpy array
        dry air potential temperature (K)
    
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
    Calculates moist air potential temperature.


    Parameters
    ----------
    r : numpy array
        mixing ratio i.e. vapor mass per dry air mass in [kg / kg]

    p : numpy array
        total gas pressure [Pa]

    T : numpy array
        temperature [K]

    
    Returns
    --------
    theta_m : numpy array
        moist air potential temperature (K)
    
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
    Calculates equivalent potential temperature.
    after Bolton (1980) MWR 108, p.1046, eq. (43) 

    Davies-Jones (2009) MWR137, eq. (6.3)

    also using either one of (15), (21) or (22)


    Parameters
    ----------
    r : numpy array
        mixing ratio i.e. vapor mass per dry air mass in [kg / kg]

    p : numpy array
        total gas pressure [Pa]

    T : numpy array
        temperature [K]

    eq : int, optional, default = 15
        which eq. of Bolton (1980) to be chosen for lifting condensation level temperature
        eq in [15, 21, 22]

    
    Returns
    --------
    theta_E : numpy array
        equivalent potential temperature

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

    print(dry_air_density(r, p, T))
    print(absolute_humidity(r, p, T))
    print(water_vapor_pressure(r, p, T))
    print(saturation_pressure(T))
    print(dew_point(r, p, T))
    hrel = relative_humidity(r, p, T)
    print(hrel2qv(hrel, p, T))
    
    print(lifting_condensation_level_temperature(r, p, T, eq = 15))
    print(lifting_condensation_level_temperature(r, p, T, eq = 21))
    print(lifting_condensation_level_temperature(r, p, T, eq = 22))
    print(equivalent_potential_temperature(r, p, T))

