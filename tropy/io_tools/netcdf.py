#!/usr/bin/env python

# ####################################################################
#
# description
# ===========
# collection of routines to read netcdf files 
#
# to do
# =====
# * mainly developed for icon files: generalization needed
#
# ####################################################################


import os, sys
import numpy as np
import netCDF4
import datetime

import hdf as hio

from standard_config import *

######################################################################
######################################################################
# 
# input forecast data 
# ====================
#
######################################################################
######################################################################

def read_icon_georef(fname):

    dset = {}

    # input dimensions
    f = netCDF4.Dataset(fname, 'r')
 
    lon = f.variables['lon'][:]
    lat = f.variables['lat'][:]

    if np.ndim(lon) == 1 and np.ndim(lat) == 1:
        dset['clon'], dset['clat'] = np.meshgrid(lon, lat, indexing = 'xy')

        dset['lon'] = dset['clon'] # link also the more common variable names
        dset['lat'] = dset['clat'] # link also the more common variable names
        
    else:
        dset['lon'], dset['lat'] = lon, lat
        
    f.close()


    return dset

######################################################################
######################################################################


def save_icon_georef(fname, geopath = None):

    geo = read_icon_georef(fname)

    for k in geo.keys():
        geo[k] = reshape_icon_fields( geo[k] )

    if geopath == None:
        geopath = os.environ['LOCAL_DATA_PATH'] + '/icon/aux'

    gname = '%s/icon_coarse_georef.h5' % geopath
    print '... save georef to %s' % gname
    hio.save_dict2hdf(gname, geo)

    return

######################################################################
######################################################################

  
def read_icon_2d_data(fname, var_list, itime = 0):


    dset = {}
    
    # input dimensions
    f = netCDF4.Dataset(fname, 'r')

 
    # input data
    for vname in var_list:
        dset[vname] = np.array( f.variables[vname][itime] )
    
                
    f.close()

    return dset


######################################################################
######################################################################
  
def roundTime(dt=None, roundTo=60):
    """Round a datetime object to any time laps in seconds
       dt : datetime.datetime object, default now.
      roundTo : Closest number of seconds to round to, default 1 minute.
      Author: Thierry Husson 2012 - Use it as you want but don't blame me.
    """
    if dt == None : dt = datetime.datetime.now()
    
    seconds = (dt - dt.min).seconds
            
    # // is a floor division, not a comment on following line:
    rounding = (seconds+roundTo/2) // roundTo * roundTo
                
    return dt + datetime.timedelta(0,rounding-seconds,-dt.microsecond)


######################################################################
######################################################################
 

def read_icon_time(fname, itime = 0):


    dset = {}
    
    # input dimensions
    f = netCDF4.Dataset(fname, 'r')

    t = f.variables['time'][itime]
    t_int = np.int(t)

    date_str = str( t_int )
    day_fraction = t - t_int
            
    tobj = datetime.datetime.strptime(date_str, '%Y%m%d') 
    tobj += datetime.timedelta(days = day_fraction)

    tobj = roundTime(tobj)
    f.close()

    return tobj


######################################################################
######################################################################
    
def read_icon_4d_data(fname, var_list, itime = 0, itime2 = None):

    '''
    Read netcdf data.

    USAGE:
    ======
    dset = read_icon_4d_data(fname, var_list, itime = 0)


    INPUT:
    ======
    fname: netcdf filename
    var_list: variable name or list of variables
    itime (OPTIONAL): index of 1st dimension (assumed to be time) to be read only
                      if itime = None, the full field is read

    OUTPUT:
    =======
    dset: dictionary containing the fields
    '''



    # check if var_list is one variable or list
    if type(var_list) == type(''):
        var_list = [var_list,]

    dset = {}
    
    # input dimensions
    f = netCDF4.Dataset(fname, 'r')

    # input data
    for vname in var_list:
        v = f.variables[vname]
        
        if itime == None:
            dset[vname] =  np.array( v[:] )
        elif itime2 == None:
            dset[vname] =  np.array( v[itime] )
        else:
            dset[vname] =  np.array( v[itime:itime2] )
                
    f.close()

    return dset


######################################################################
######################################################################

def read_icon_dimension(fname, dim_name):

    f = netCDF4.Dataset(fname, 'r')
    ndim = len( f.dimensions[dim_name] )
    f.close()

    return ndim



######################################################################
######################################################################

def save_icon_time_reference(fname, outfile = None):


    tstamp = os.path.splitext(fname)[0].split('_')[-1]

    print tstamp

    times = {}
    ntimes = read_icon_dimension(fname, 'time')
    for n in range(ntimes):
        tstamp_n = '%s_%s' % (tstamp, str(n).zfill(2))
        t = read_icon_time(fname, n)
        tstr = t.strftime('%Y%m%d_%H%M')
        times[tstamp_n] = tstr

    if outfile == None:
        outfile = '%s/icon/aux/icon_coarse_timeref.h5' % os.environ['LOCAL_DATA_PATH']

    print '... save timeref to %s' % outfile
    hio.update_dict_in_hdf(outfile, times)


    return

######################################################################
######################################################################


if __name__ == '__main__':
  
    print 'start input test'

    try:
        fname = sys.argv[1]
    except:
        print 'ERROR: no input file given!'

    # save_icon_georef(fname)
 
    save_icon_time_reference(fname)
