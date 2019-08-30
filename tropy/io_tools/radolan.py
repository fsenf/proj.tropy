#!/usr/bin/python


import sys, os, glob, copy
import datetime
import tempfile

import numpy as np
import pylab as pl

import tropy.io_tools.netcdf as ncio

######################################################################
######################################################################

'''
Description
===========
The routine collection can 
   * generate filenames for Radolan RX and RW products based on time objects, 
   * read Radolan RX and RW data (also from .gz files) into a Radolan class  
   * provide georeference


User adjustments:
================
* arch in  radoname_from_time might also be changed
'''


#
#
######################################################################
######################################################################

from tropy.standard_config import *

######################################################################


def parse_DWD_quant_composite_header(header):
    """Parses the ASCII header of a DWD quantitative composite file

    Parameters
    ----------
    header : string (ASCII header)

    Returns
    -------
    output : dictionary of metadata retreived from file header

    """
    # empty container
    out = {}
    # RADOLAN product type def
    out["producttype"] = header[0:2]
    # file time stamp as Python datetime object
    out["datetime"] = datetime.datetime.strptime(header[2:8]+header[13:17]+"00",
                                           "%d%H%M%m%y%S")
    # radar location ID (always 10000 for composites)
    out["radarid"] = header[8:13]
    pos_VS = header.find("VS")
    pos_SW = header.find("SW")
    pos_PR = header.find("PR")
    pos_INT = header.find("INT")
    pos_GP = header.find("GP")
    pos_MS = header.find("MS")
    if pos_VS > -1:
        out["maxrange"] = {0:"100 km and 128 km (mixed)",
                           1: "100 km",
                           2:"128 km",
                           3:"150 km" }[int(header[(pos_VS+2):pos_VS+4])]
    else:
        out["maxrange"] = "100 km"
    out["radolanversion"] = header[(pos_SW+2):pos_SW+11]
    out["precision"] = 10**int(header[pos_PR+4:pos_PR+7])
    out["intervalseconds"] = int(header[(pos_INT+3):pos_INT+7])*60
    dimstrings = header[(pos_GP+2):pos_GP+11].strip().split("x")
    out["nrow"] = int(dimstrings[0])
    out["ncol"] = int(dimstrings[1])
    locationstring = header[(pos_MS+2):].strip().split("<")[1].strip().strip(">")
    out["radarlocations"] = locationstring.split(",")
    return out

######################################################################
######################################################################


def read_RADOLAN_composite(fname, missing=-9999, scale = 0.5, offset = -32.5 ):
    """Read quantitative radar composite format of the German Weather Service

    The quantitative composite format of the DWD (German Weather Service) was
    established in the course of the `RADOLAN project <http://www.dwd.de/radolan>`
    and includes several file types, e.g. RX, RO, RK, RZ, RP, RT, RC, RI, RG and
    many, many more (see format description on the project homepage, [DWD2009]).

    At the moment, the national RADOLAN composite is a 900 x 900 grid with 1 km
    resolution and in polar-stereographic projection.

    **Beware**: This function already evaluates and applies the so-called PR factor which is
    specified in the header section of the RADOLAN files. The raw values in an RY file
    are in the unit 0.01 mm/5min, while read_RADOLAN_composite returns values
    in mm/5min (i. e. factor 100 higher). The factor is also returned as part of
    attrs dictionary under keyword "precision".

    Parameters
    ----------
    fname : path to the composite file

    missing : value assigned to no-data cells

    Returns
    -------
    output : tuple of two items (data, attrs)
        - data : numpy array of shape (number of rows, number of columns)
        - attrs : dictionary of metadata information from the file header

    References
    ----------

    .. [DWD2009] Germany Weather Service (DWD), 2009: RADLOAN/RADVO-OP -
        Beschreibung des Kompositformats, Version 2.2.1. Offenbach, Germany,
        URL: http://dwd.de/radolan (in German)

    """
    mask = 4095 # max value integer
    NODATA = missing
    header = '' # header string for later processing
    # open file handle
    f = open(fname, 'rb')
    # read header
    while True :
        mychar = f.read(1)
        if mychar == chr(3) :
            break
        header = header + mychar
    attrs = parse_DWD_quant_composite_header(header)
    attrs["nodataflag"] = NODATA
    if not attrs["radarid"]=="10000":
        warnings.warn("WARNING: You are using function e" +
                      "wradlib.io.read_RADOLAN_composit for a non " +
                      "composite file.\n " +
                      "This might work...but please check the validity " +
                      "of the results")
    if attrs["producttype"] == "RX":
        # read the actual data
        indat = f.read(attrs["nrow"]*attrs["ncol"])
        # convert from 8-bit integers
        # and upgrade to 32-bit ints, so that nodata values may be inserted
        arr = np.frombuffer(indat, np.uint8).astype(np.int)
        arr = np.where(arr==250,NODATA,arr)
        clutter = np.where(arr==249)[0]
    else:
        # read the actual data
        indat = f.read(attrs["nrow"]*attrs["ncol"]*2)
        # convert to 16-bit integers
        arr = np.frombuffer(indat, np.uint16).astype(np.int)
        # evaluate bits 14, 15 and 16
        nodata   = np.where(arr & int("10000000000000",2))
        negative = np.where(arr & int("100000000000000",2))
        clutter  = np.where(arr & int("1000000000000000",2))
        # mask out the last 4 bits
        arr = arr & mask


        # consider negative flag if product is RD (differences from adjustment)
        if attrs["producttype"]=="RD":
            # NOT TESTED, YET
            arr[negative] = -arr[negative]

        # apply precision factor
        arr = arr.astype(np.float) * attrs["precision"]

        # set nodata value
        arr[nodata] = NODATA

    # bring it into shape
    arr = arr.reshape( (attrs["nrow"], attrs["ncol"]) )

    # append clutter mask
    attrs['cluttermask'] = clutter

    # close the file
    f.close()

    
    
    # do the rescaling
    if attrs["producttype"] == "RX":
        attrs['scale'] = scale
        attrs['offset'] = offset
        arr = np.where(arr != missing, scale * arr  +  offset, missing)
    
    return np.ma.masked_equal(arr,missing), attrs

######################################################################
######################################################################


class Radolan(object):

    '''
    Simple class for reading and holding Radolan data and georeference.
    '''

    def __init__(self, do_lonlat = True):

        # calculate georeference
        if do_lonlat:
            self.lonlat()
        

    # ...............................................................


    def read(self, name_or_time, rproduct = 'rx_hdcp2'):

        '''
        Reads Radolan data either 
         (i) from DWD raw data format.
         (ii) or from HDCP2 netcdf format
        '''

        if rproduct in ['rx', 'rw', 'sf']:
            self.rawread(name_or_time, rproduct = rproduct)

        elif rproduct == 'rx_hdcp2':
            self.read_nc(name_or_time)
    # ...............................................................


    def read_nc(self, time):

        '''
        Reads Radolan data from HDCP2 project.
        '''
        
        
        # get name of radolan file
        fname = radoname_from_time(time, 
                                   rproduct = 'rx_hdcp2')

        # determine itime from time object
        # TODO: SIMPLE METHOD TO BE IMPROVED
        itime = 12 * time.hour + time.minute / 5
    
        vname = 'dbz'
        d = ncio.read_icon_4d_data(fname, [vname], itime = itime)[vname]
        d = np.ma.masked_invalid( d )
 
        self.data = np.ma.masked_greater_equal( d, 92)

        # for testing our geo-reference -> HAPPY, HAPPY, JOY, JOY -> it works
        # geo = ncio.read_icon_4d_data(fname, ['lon', 'lat'], itime = None)
        # self.rlon = geo['lon']
        # self.rlat = geo['lat']
        
        
    # ...............................................................


    def rawread(self, name_or_time, rproduct = 'rx'):

        '''
        Reads Radolan data from DWD raw data format.
        '''
        
        # check if input is filename or time object
        if type(name_or_time) == type(''):
            fname = name_or_time

        elif type(name_or_time) == type(datetime.datetime(2012,12,1,12,0)):

            t = name_or_time
            fname = radoname_from_time(t, rproduct = rproduct)

        
        # test filename existence somehow ...
        flist = glob.glob(fname + '*')
        print(flist) 
            
        if len(flist ) == 0:
            print('ERROR: file %s does not exist' % fname)
            sys.exit()

        if len(flist) == 1:

            # filename manipulations
            dir_and_base, ext = os.path.splitext(flist[0])
            base = os.path.basename(dir_and_base)


            if ext == '.gz':
                
                # create gizp output name
                tmpdir = tempfile.mkdtemp()
                newname = '%s/%s' % (tmpdir, base) 

                # uncompress data and redirect it 
                os.system('gunzip -c %s > %s' % (flist[0], newname) )
                
            else:
                newname = fname

        # read the data ..............................................
        self.data, self.meta = read_RADOLAN_composite(newname)


        # cleanup again .............................................. 
        # os.system('gzip %s &' % fname)
        if os.path.isdir(tmpdir):
            os.system('rm -rf %s' %  tmpdir)

        return

    # ................................................................

    def dbz2rr(self):
        
        '''
        Simple Z to R conversion.
        '''

        dbz = self.data

        a = 256.
        b = 1.42
        
        z = dbz / (10*b)
        f = 1. / (a**(1./b))
        
        self.rr =  f * 10**z

        return self.rr

    # ................................................................

    def lonlat(self):

        '''
        Calculates longitude and latitude.
        '''

        x,y = radolan_xy()
        self.xg, self.yg = np.meshgrid(x,y)

        
        self.lon, self.lat = rado_xy2ll(self.xg, self.yg)

        return
    # ................................................................


    def mask(self, thresh=35):
        '''
        Performs masking with defined Z threshold.
        '''
        
        return np.ma.masked_less(self.data, thresh)

######################################################################
######################################################################

def radolan_xy(x0 = -523.4622, y0 = -4658.645, shift = 0, edge = 0):
 
    '''
    Local coordinates of Radolan projection are calculated as vectors.
    

    Input
    =====
    x0: offset for x-direction in km (lower left corner) OPTIONAL
    y0: offset for y-direction in km (lower left corner) OPTIONAL
    shift: subpixel shift (center vs. edge values) OPTIONAL


    Output
    ======
    x: vector of local Radolan x-coordinate (in km)
    y: vector of local Radolan y-coordinate (in km)

    '''

    # grid indices
    i = np.arange(-edge, 900 + edge)
    j = np.arange(-edge, 900 + edge)

    # grid distance
    d = 1.

    x = x0  + d * i + shift
    y = y0  + d * j + shift

    return x, y


######################################################################
######################################################################

def rado_xy2ll(x,y):
  
    '''
    Transforms local Radolan coordinates into lon/lat.
    

    Input
    =====
    x: 2d grid of local Radolan x-coordinate (in km)
    y: 2d grid of local Radolan y-coordinate (in km)


    Output
    =====
    lon: 2d grid of longitude (deg)
    lat: 2d grid of latitude (deg)

    '''


    # earth radius
    R = 6370.04

    # projection parameters
    lon0 = np.deg2rad(10.)
    lat0 = np.deg2rad(60.)


    Dq = x**2 + y**2
    Rfac = R**2. * ( 1 + np.sin(lat0))**2.

    lon = np.arctan(-x/y)  + lon0
    lat = np.arcsin( (Rfac - Dq) / (Rfac + Dq))


    return np.rad2deg(lon), np.rad2deg(lat)




######################################################################
######################################################################


def rado_ll2xy(lon,lat):
       
    '''
    Transforms lon/lat into local Radolan coordinates.
    

    Input
    =====
    lon: 2d grid of longitude (deg)
    lat: 2d grid of latitude (deg)


    Output
    =====
    x: 2d grid of local Radolan x-coordinate (in km)
    y: 2d grid of local Radolan y-coordinate (in km)


    '''

    fac = np.pi / 180. 
    
# transform: deg to rad
    lat_0 = 60. * fac
    lon_0 = 10. * fac
    
    lon = lon * fac
    lat = lat * fac


# stereographic scaling factor
    M = (1. + np.sin(lat_0)) / (1. + np.sin(lat))


# earth radius
    R = 6370.04

# transformation
    x =  R * M * np.cos(lat) * np.sin(lon - lon_0)
    y = -R * M * np.cos(lat) * np.cos(lon - lon_0)


    return x,y


######################################################################
######################################################################    



def rado_xy2ij(x, y, x0 = -523.4622, y0 = -4658.645):

    
    '''
    Transforms local coordinates given in Radolan projection into an index
    field that gives the corresponding nearest neighbor index.

    Assumes grid size of 1 x 1 km.

    Input
    =====
    x: 2d grid of local Radolan x-coordinate (in km)
    y: 2d grid of local Radolan y-coordinate (in km)
    x0: offset for x-direction in km (lower left corner) OPTIONAL
    y0: offset for y-direction in km (lower left corner) OPTIONAL
 

    Output
    =====
    ir: row index
    ic: column index

    '''

    x_rel = x - x0
    y_rel = y - y0

    irow = np.round(y_rel).astype(np.int)
    icol = np.round(x_rel).astype(np.int)

    return irow, icol


######################################################################
######################################################################


def radoname_from_time(tin, 
                       fmt = '%Y%m%d_%H%M', 
                       rproduct = 'rx',
                       arch = '/%s/radolan' % global_data_path):

    '''
    Generates Radolan filename from time object.


    Input
    =====
    tin: time as datetime object or as time string with format fmt
    fmt: format of time string (OPTIONAL)
    rproduct: Radolan product type (OPTIONAL), e.g. 'rx' or 'rw'
    arch: directory where data are stored (OPTIONAL)


    Output
    =====
    fname: Radolan file name


    '''

    # convert input to datetime object -------------------------------
    try:
        t = datetime.datetime.strptime(tin,fmt)

    except:
        t = tin

    if type(t) != datetime.datetime:
        raise ValueError
    
    # ================================================================


    if rproduct in ['rx', 'rw', 'sf']:
        date_str = t.strftime('%Y/%m/%d')
        time_str = t.strftime('%y%m%d%H%M')

        fmt = {'arch':arch, 'date':date_str, 'time':time_str, 'rprod':rproduct}
        fname = '{arch}/{rprod}/{date}/raa01-{rprod}_10000-{time}-dwd---bin'.format(**fmt) 
        
    elif rproduct == 'rx_hdcp2':
        date_str = t.strftime('%Y')
        time_str = t.strftime('%Y%m%d')
        fmt = {'arch':arch, 'date':date_str, 'time':time_str, 'rprod':rproduct}

        fname = '{arch}/{rprod}/{date}/hdfd_miub_drnet00_l3_dbz_v??_{time}000000.nc'.format(**fmt) 
        fname = glob.glob(fname)[0]
        

    return fname

######################################################################
######################################################################

class RadoSequence(Radolan):
    '''
    Reads a sequence of Radolan slots.
    '''

    def __init__(self, **kwargs):

        super(RadoSequence, self).__init__(**kwargs)

        return 
    # ================================================================
    
    def dt_from_rprod(self, **kwargs):

        rproduct = kwargs.get('rproduct', 'rx')

        if rproduct == 'rx':
            dt_default = 5.

        elif rproduct == 'rw':
            dt_default = 60.

        dt = kwargs.get('dt', dt_default)

        return dt

    ##################################################################

    def sload(self, t1, t2, **kwargs):

        '''
        Load a time sequence of Radolan data.

        INPUT
        =====
        t1: start time (time object)
        t2: end time, exclusive/not included (time object)
        dt: optional, time step in minutes
        rproduct: Radolan product either 'rx' or 'rw', set default dt to 5 or 60 min
        '''


        # determine time sequence ------------------------------------
        rproduct = kwargs.get('rproduct', 'rx')
        dt = self.dt_from_rprod(**kwargs)
        self.tlist = self.time_sequence(t1, t2, dt = dt)
        tlist = self.tlist
        ntimes = len(tlist)
        # ============================================================


        # read data from list ----------------------------------------
        for n, t in enumerate(tlist):

            self.read(t, rproduct = rproduct)

            d = self.data
            
            try:
                data_seq[n]  = d
            except:
                nrows, ncols = d.shape
                data_seq = np.ma.zeros((ntimes, nrows, ncols))
                data_seq[n] = d
            

        self.data = data_seq
        # ============================================================
            
        return 


    ##################################################################
    
    def supdate(self, t1, t2, **kwargs):

        '''
        Updates an Radolan sequence class object.
        
        It keeps old data if contained in time range, and reads new data.
        '''

        # determine time sequence ------------------------------------
        rproduct = kwargs.get('rproduct', 'rx')
        dt = self.dt_from_rprod(**kwargs)
        newlist = self.time_sequence(t1, t2, dt = dt)
        ntimes =  len(newlist)
        # ============================================================

        
        # read data from list ----------------------------------------
        for n, t in enumerate(newlist):


            # try to get data from existing array
            try:
                d_index = self.tlist.index( t ) 
                d = self.data[d_index]

            except:
                
                self.read(t, rproduct = rproduct)
                d = self.data
            
            try:
                data_seq[n]  = d
            except:
                nrows, ncols = d.shape
                data_seq = np.ma.zeros((ntimes, nrows, ncols))
                data_seq[n] = d
            

        self.data = data_seq
        self.tlist = newlist
        # ============================================================
   
        return

    ##################################################################
  
    def time_sequence(self, t1, t2, dt = 5.):
        
        # test for time format ---------------------------------------
        if type(t1) == type(''):
            t1_obj = datetime.datetime.strptime(t1, '%Y%m%d_%H%M')
        else:
            t1_obj = t1

        if type(t2) == type(''):
            t2_obj = datetime.datetime.strptime(t2, '%Y%m%d_%H%M')
        else:
            t2_obj = t2
        # ============================================================
        

        # loop over times and make time list -------------------------
        tlist = []        
        t = copy.copy(t1_obj)

        while t < t2_obj:
            tlist.append( t )
            t += datetime.timedelta( minutes = dt)
        # ============================================================

        return tlist


######################################################################
######################################################################

if __name__=='__main__':

    t = datetime.datetime(2014, 6, 8, 0, 0)
    fname = radoname_from_time(t, 
                       rproduct = 'rx_hdcp2')

    r = Radolan()
    r.read(t, rproduct = 'rx_hdcp2')


if False:
    t1 = datetime.datetime(2014, 6, 9, 0, 0)
    t2 = t1 + datetime.timedelta(minutes = 15)

    t3 = t1 + datetime.timedelta(minutes = 5)
    t4 = t3 + datetime.timedelta(minutes = 20)


    r = RadoSequence()
    r.sload(t1, t2)


    r.supdate(t3, t4)

if False:
    r = Radolan()
    r.read(radoname_from_time(t, rproduct='rw'))

    r = Radolan()
    r.read(radoname_from_time(t))

    r = Radolan()
    r.read(t)
