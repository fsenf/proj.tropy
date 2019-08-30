#!/usr/bin/env python
######################################################################
# 
# DESCRIPTION
# ===========
#
# Library for reading and manipulating MSG data from hdf files in the 
# TROPOS Archive.
#
# ####################################################################


# load libraries -----------------------------------------------------
import sys, os, glob
import numpy as np
import h5py
import datetime
import subprocess
# ====================================================================

import tropy.io_tools.netcdf as ncio


######################################################################
######################################################################

def get_highres_cma(t, region = 'de', 
                       scan_type = 'rss', 
                       arch = "none"):


    '''
    Reads product of High-Res Cloud Mask product (Bley et al., 2013) for a given date.


    Parameters
    ----------
    t : datetime object 
        day and time of MSG time slot

    region : str, optional, default = 'de'
        string that defines regional cutout

    scan_type : str, optional, default = 'rss'

    arch : str, optional, default = "none"
        archive directory where NWC-SAF files are saved (OPTIONAL)


    Returns
    --------
    var : numpy array
        Cma product as field
    '''


    # set date and time strings --------------------------------------
    date_string = t.strftime('%Y/%m/%d')
    tstr = t.strftime('%Y%m%dt%H%Mz')
    time_string = tstr
    # ================================================================


   
    # check scan_type ................................................
    scan =  scan_type.lower()

    if not scan in ('rss', 'pzs'):
        print('ERROR: scan_type %s not available')
        print('       use either "rss" or "pzs"')
        return None
    # ================================================================



    # set archive directory ------------------------------------------
    base_dir = "/vols/altair/datasets/eumcst/"

    if arch == 'none':
        arch_dir =  "%s/msevi_%s/l2_hrv_products/%s/%s" % (base_dir, 
                                                           scan, 
                                                           region, 
                                                           date_string)


    else: 
        arch_dir = arch
    # ================================================================

     
    # set file name --------------------------------------------------
    prod_file = arch_dir + "/msg?-sevi-%s-l2-hrvcma-%s-%s.c2.nc" % (time_string, 
                                                                      scan, 
                                                                      region)

    

    # get filename while MSG number is unspecified ...................
    pfile_list =  glob.glob(prod_file)

    if pfile_list:
        prod_file = pfile_list[0]
    else:
        raise Exception('ERROR: ', prod_file,' does not exist!')
    # ================================================================


    # read variable --------------------------------------------------
    # TODO function call should updated !!!
    prod = 'clm'
    var = ncio.read_icon_4d_data(prod_file, prod, itime = None)
    # ================================================================

    return var[prod]

######################################################################
######################################################################


def get_berendes_prod(prod, t, 
                      region = 'de', 
                      scan_type = 'rss', 
                      arch = "none"):


    '''
    Reads product of Berendes Cloud type product for a given date.



    Parameters
    ----------
    prod : str
        name of Berendes product, either ccm or texture

    t : datetime object 
        day and time of MSG time slot

    region : str, optional, default = 'de'
        string that defines regional cutout

    scan_type : str, optional, default = 'rss'

    arch : str, optional, default = "none"
        archive directory where NWC-SAF files are saved (OPTIONAL)


    Returns
    --------
    var : numpy array
        Berendes product as field
    '''


    # set date and time strings --------------------------------------
    date_string = t.strftime('%Y/%m/%d')
    tstr = t.strftime('%Y%m%dt%H%Mz')
    time_string = tstr
    # ================================================================


   
    # check scan_type ................................................
    scan =  scan_type.lower()

    if not scan in ('rss', 'pzs'):
        print('ERROR: scan_type %s not available')
        print('       use either "rss" or "pzs"')
        return None
    # ================================================================



    # set archive directory ------------------------------------------
    base_dir = "/vols/altair/datasets/eumcst/"

    if arch == 'none':
        arch_dir =  "%s/msevi_%s/l2_hrv_products/%s/%s" % (base_dir, 
                                                           scan, 
                                                           region, 
                                                           date_string)


    else: 
        arch_dir = arch
    # ================================================================

     
    # set file name --------------------------------------------------
    prod_file = arch_dir + "/msg?-sevi-%s-l2-berendes_ccm-%s-%s.c2.nc" % (time_string, 
                                                                      scan, 
                                                                      region)

    

    # get filename while MSG number is unspecified ...................
    pfile_list =  glob.glob(prod_file)

    if pfile_list:
        prod_file = pfile_list[0]
    else:
        raise Exception('ERROR: ', prod_file,' does not exist!')
    # ================================================================


    # read variable --------------------------------------------------
    # TODO function call should updated !!!
    var = ncio.read_icon_4d_data(prod_file, prod, itime = None)
    # ================================================================

    return var[prod]

######################################################################
######################################################################

def prod_reference(prod): # LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

    '''
    Given a NWC-SAF product, the routine returns the correspondig
    substring to address the product file.


    Parameters
    ----------
    prod :  str
        name of NWC-SAF product, e.g. CMa, CT, CRR


    Returns
    --------
    name : str
        substring in the NWC-SAF product file


    Notes
    -----
    Further products with its corresponding substrings should be 
    included in future.

    '''

    

# set connection between product name and product file
    prod_list={}

# CMa
    prod_list['CMa'] = 'CMa__'
    prod_list['CMa_DUST'] = 'CMa__'


# CT
    prod_list['CT'] = 'CT___'
    prod_list['CT_PHASE'] = 'CT___'


# CTTH
    prod_list['CTTH_PRESS'] = 'CTTH_'
    prod_list['CTTH_HEIGHT'] = 'CTTH_'
    prod_list['CTTH_TEMPER'] = 'CTTH_'
    prod_list['CTTH_EFFECT'] = 'CTTH_'
    prod_list['CTTH_QUALITY'] = 'CTTH_'

# CRR
    prod_list['CRR'] = 'CRR__'
    prod_list['CRR_ACCUM'] = 'CRR__'
    prod_list['CRR_INTENSITY'] = 'CRR__'


    try:
        name = prod_list[prod]
    except:
        print('no file associated with product')
        return

    return name # TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT



######################################################################
######################################################################


def get_nwcsaf_prod(prod, day, 
                    scan_type = 'rss', 
                    arch = 'none', 
                    region = 'eu', 
                    calibrate = False,
                    is_region_shifted = False):


    '''
    Reads product of NWC-SAF for a given date.


    Parameters
    ----------
    prod : str
        name of NWCSAF product
        prod in ['CMa', 'CMa_DUST', 'CT', 'CTTH_EFFECT', 'CTTH_HEIGHT', 'CTTH_PRESS', 'CTTH_QUALITY', 'CTTH_TEMPER', 'CT_PHASE']


    day : datetime object 
        day and time of MSG time slot

    region : str, optional, default = 'eu'
        string that defines regional cutout

    scan_type : str, optional, default = 'rss'
        one of the two Meteosat scan types

    arch : str, optional, default = "none"
        archive directory where NWC-SAF files are saved 

    calibrate : bool, optional, default = False
        switch for slope/offset calibration 

    is_region_shifted : bool, optional, default = False
        switch if bugfix for (mistakingly) shifted regions has to be applied


    Returns
    --------
    product : numpy array, optional
        calibrated NWCSAF product data (if calibrate == True)

    counts : numpy array, optional
        NWCSAF product count data (if calibrate == False)
    
    scal : float, optional
        scale factor for NWCSAF products (if calibrate == False)

    offset : float, optional
        offset for NWCSAF products (if calibrate == False)


    Notes
    ------
    Data are saved as integer via: DATA = scal * count +  offset


    The following structures can be returned:

    counts, scal, offset = get_nwcsaf_prod(prod, day, calibrate = False) 

     OR

    product = get_nwcsaf_prod(prod, day, calibrate = True)

    '''

    
    # get strings of filename creation -------------------------------
    time_string = day.strftime('%Y%m%d%H%M')
    date_string = day.strftime('%Y/%m/%d')

    # check scan_type ................................................
    scan =  scan_type.lower()

    if not scan in ('rss', 'pzs'):
        print('ERROR: scan_type %s not available')
        print('       use either "rss" or "pzs"')
        return None
    
    # check product name .............................................
    prodname = prod_reference(prod)

    if not prodname:
        print('ERROR: key %s is not available!' % prod)
        return
    # ================================================================

    
    
    # create hdf filename --------------------------------------------

    # set archive directory ..........................................
    base_dir = "/vols/altair/datasets/eumcst/"

    if arch == 'none':
        arch_dir =  "%s/msevi_%s/l2_nwcsaf/eu/%s" % (base_dir, scan, date_string)
    else: 
        arch_dir = arch


    s = '________'
    l = 8 - len(region)
    reg_str = region + s[:l]
    prod_file = arch_dir + "/SAFNWC_MSG?_%s%s_%s-%s.c?.h5" % (prodname, time_string, scan,reg_str)


    # get filename while MSG number is unspecified ...................
    pfile_list =  glob.glob(prod_file)

    if pfile_list:
        prod_file = pfile_list[0]
    else:
        print('ERROR: ', prod_file,' does not exist!')
        return None
    # ================================================================


# I/O ----------------------------------------------------------------
    if not h5py.is_hdf5(prod_file):
        print('ERROR: searched for ',prod_file)
        print('No hdf-file available!')
        return

# open files .........................................................
    fprod = h5py.File(prod_file,"r")

# counts
    counts = fprod[prod].value


# BUGFIX
# ======
# Unfortunately, there was a bug in the region defintion of the 
# nwcsaf 2012 configuration which led to a pixel displacement :-(
    if is_region_shifted:
        newcounts = np.zeros_like(counts)
        newcounts[:-1, :-1] = counts[1:, 1:]

        counts = newcounts
#--- END BUGFIX


# scaling factor
    scal =  fprod[prod].attrs['SCALING_FACTOR']

# offset
    offset =  fprod[prod].attrs['OFFSET']


    fprod.close()
# ====================================================================


    if calibrate:
        return np.where(counts==0, 0, scal * counts + offset)
    else:
        return counts, scal, offset # TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT




######################################################################
######################################################################


def get_cmsaf_prod(prod, day, 
                   scan_type = 'rss', 
                   region = 'eu', 
                   arch = 'none', 
                   calibrate = False,
                   switch2new = True): # LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

    '''
    Reads product of CM-SAF (KNMI retrieval) for a given date.



    Parameters
    ----------    
    prod : str
        name of CMSAF product

    day : datetime object 
        day and time of MSG time slot

    region : str, optional, default = 'eu'
        string that defines regional cutout

    scan_type : str, optional, default = 'rss'
        one of the two Meteosat scan types

    arch : str, optional, default = "none"
        archive directory where NWC-SAF files are saved 

    calibrate : bool, optional, default = False
        switch for slope/offset calibration 

    switch2new : bool, optional, default = True
        if the "new" set of CMSAF products should be selected 
        "new" means collection two: "c2"


    Returns
    --------
    product : numpy array, optional
        calibrated CMSAF product data (if calibrate == True)

    counts : numpy array, optional
        CMSAF product count data (if calibrate == False)
    
    scal : float, optional
        scale factor for CMSAF products (if calibrate == False)

    offset : float, optional
        offset for CMSAF products (if calibrate == False)


    Notes
    ------
    Data are saved as integer via: DATA = scal * count +  offset


    The following structures can be returned:

    counts, scal, offset = get_cmsaf_prod(prod, day, calibrate = False) 

     OR

    product = get_cmsaf_prod(prod, day, calibrate = True)

    '''

    

    time_string = day.strftime('%Y%m%d%H%M')
    if switch2new:
        time_string = day.strftime('%Y%m%dt%H%Mz')

    date_string = day.strftime('%Y/%m/%d')

    
    # check scan_type ................................................
    scan =  scan_type.lower()

    if not scan in ('rss', 'pzs'):
        print('ERROR: scan_type %s not available')
        print('       use either "rss" or "pzs"')
        return None
    # ================================================================


    # create hdf filename --------------------------------------------

    # set archive directory ..........................................
    base_dir = "/vols/altair/datasets/eumcst/"

    if arch == 'none':
        arch_dir =  "%s/msevi_%s/l2_cpp/%s" % (base_dir, scan, date_string)

        if switch2new:
            arch_dir =  "%s/msevi_%s/l2_cpp_new/%s/%s" % (base_dir, 
                                                          scan, 
                                                          region, 
                                                          date_string)


    else: 
        arch_dir = arch


    s = '________'
    l = 8 - len(region)
    reg_str = region + s[:l]
    prod_file = arch_dir + "/msg?-sevi-"+time_string+"-cpp-%s-eu.h5" % scan
        
    if switch2new:
        prod_file = arch_dir + "/msg?-sevi-%s-l2cpp-%s-eu.c2.h5" % (time_string, scan)



    # get filename while MSG number is unspecified ...................
    pfile_list =  glob.glob(prod_file)

    if pfile_list:
        prod_file = pfile_list[0]
    else:
        print('ERROR: ', prod_file,' does not exist!')
        return None
    # ================================================================


# I/O
# ====
    if not h5py.is_hdf5(prod_file):
        print('ERROR: searched for ', prod_file)
        print('No hdf-file available!')
        return

# open files .........................................................
    fprod = h5py.File(prod_file,"r")


    if prod in fprod:
        # counts
        counts = fprod[prod].value

    # scaling factor
        slope =  fprod[prod].attrs['gain']

    # offset
        offset = fprod[prod].attrs['intercept']

    else:
        counts = 0
        slope = 0
        offset = 0


    fprod.close()

    if calibrate:
        var = np.where(counts!=-1, slope * counts + offset,0)
        return var

    else:
        return counts, slope, offset # TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT



######################################################################
######################################################################


def list_cmsaf_prod(day, arch='none'): # LLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

    '''
    Lists products contained in a specific CM-SAF (KNMI retrieval) file
    at a given date.


    Parameters
    ----------
    day : datetime object 
        day and time of MSG time slot

    arch : str, optional, default = "none"
        archive directory where NWC-SAF files are saved 


    Returns
    --------
    L : list 
        list of cmsaf products

    '''


    time_string = day.strftime('%Y%m%d%H%M')
    date_string = day.strftime('%Y/%m/%d')


    # hdf5 file ..........................................................
    if arch == 'none':
        arch_dir="/mnt/zephyr/u1/EUMETCAST_ARCHIVE/msg-sevi/rss/l2_cmsaf/" + date_string
    else: 
        arch_dir = arch

    prod_file = arch_dir + "/msg1-sevi-"+time_string+"-cpp-rss-eu.h5"



# I/O
# ====
    if not h5py.is_hdf5(prod_file):
        print('ERROR: searched for ', prod_file)
        print('No hdf-file available!')
        return

# open files .........................................................
    fprod = h5py.File(prod_file,"r")

    L = list(fprod.keys())
    
    fprod.close()

    return L


######################################################################
######################################################################


def get_seviri_chan(ch_name, day, 
                    scan_type = 'rss', 
                    arch = 'none', 
                    calibrate=True): # LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

    '''
    Reads MSG - SEVIRI radiance of a given channel for a given date.


    Parameters
    ----------
    ch_name : str
        name of MSG SEVIRI channel, e.g. ir_108 or IR_108 (case does not matter)

    day : datetime object 
        day and time of MSG time slot

    scan_type : str, optional, default = 'rss'
        one of the two Meteosat scan types

    arch : str, optional, default = "none"
        archive directory where NWC-SAF files are saved 

    calibrate : bool, optional, default = False
        switch for slope/offset calibration 
        if True, radiances are output

    
    Returns
    --------
    out_field : numpy array
        radiance/counts (depending on calibrate option) of channel <ch_name> at date <day>


    info : dict
        meta info list which includes calibration coefficients
    '''

    # get strings of filename creation -------------------------------
    # change channel name to lower case, e.g. VIS008 --> vis008
    ch_name = ch_name.lower()


    # set time and date strings for file selection
    time_string = day.strftime('%Y%m%dt%H%Mz')
    date_string = day.strftime('%Y/%m/%d')

    # check scan_type ................................................
    scan =  scan_type.lower()

    if not scan in ('rss', 'pzs'):
        print('ERROR: scan_type %s not available')
        print('       use either "rss" or "pzs"')
        return None
    # ================================================================


    # create hdf filename --------------------------------------------
    
    # set archive directory ..........................................
    base_dir = "/vols/altair/datasets/eumcst/"

    if arch == 'none':
        arch_dir =  "%s/msevi_%s/l15_hdf/eu/%s" % (base_dir, scan, date_string)
    else: 
        arch_dir = arch


    # set hdf file name
    hfile = arch_dir  +  "/msg?-sevi-%s-l15hdf-%s-eu.c2.h5" % (time_string, scan)


    # get filename while MSG number is unspecified ...................
    hfile_list =  glob.glob(hfile)

    if hfile_list:
        hfile = hfile_list[0]
    else:
        print('ERROR: ', hfile,' does not exist!')
        return None
    # ================================================================


    # I/O ------------------------------------------------------------
    # check if hdf file exists 
    if not h5py.is_hdf5(hfile):
        print('No hdf-file available!')
        return None
    
    
    # open hdf file
    f = h5py.File(hfile,"r")
    

    # extract meta info for all channels
    ch_info_all=f['meta/channel_info'].value


    # get info of requested channel
    for ch_info_part in ch_info_all:
        if ch_name in list(ch_info_part):
            info=ch_info_part.copy()


    # read offset and slope
    offset = info['cal_offset']
    slope = info['cal_slope']
    
    # Read counts of sat channel
    image_name='l15_images/image_' + ch_name
    counts = f[image_name].value
    

    # calculate radiances
    radiance = np.where(counts!=0,slope * counts + offset,0)
    
    # close hdf file
    f.close()
    # ================================================================


    # do calibration if needed ---------------------------------------
    if calibrate:
        out_field = radiance
    else:
        out_field = counts
    # ================================================================

    return out_field, info # TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT



######################################################################
######################################################################


    

def get_cos_zen(day, 
                scan_type = 'rss',
                arch = 'none'): 

    '''
    Reads sun zenith angle from the partical satellite hdf file
    and calculates the cosine of the sun zenith.


    Parameters
    ----------
    day : datetime object 
        day and time of MSG time slot

    scan_type : str, optional, default = 'rss'
        one of the two Meteosat scan types


    Returns
    --------
    coszen : numpy array, optional
        cosine of sun zenith angle

    '''


    # get strings of filename creation -------------------------------
    # set time and date strings for file selection
    time_string = day.strftime('%Y%m%dt%H%Mz')
    date_string = day.strftime('%Y/%m/%d')

    # check scan_type ................................................
    scan =  scan_type.lower()

    if not scan in ('rss', 'pzs'):
        print('ERROR: scan_type %s not available')
        print('       use either "rss" or "pzs"')
        return None
    # ================================================================


    # create hdf filename --------------------------------------------
    
    # set archive directory ..........................................
    base_dir = "/vols/altair/datasets/eumcst/"

    if arch == 'none':
        arch_dir =  "%s/msevi_%s/l15_hdf/eu/%s" % (base_dir, scan, date_string)
    else: 
        arch_dir = arch


    # set hdf file name
    hfile = arch_dir  +  "/msg?-sevi-%s-l15hdf-%s-eu.c2.h5" % (time_string, scan)


    # get filename while MSG number is unspecified ...................
    hfile_list =  glob.glob(hfile)

    if hfile_list:
        hfile = hfile_list[0]
    else:
        print('ERROR: ', hfile,' does not exist!')
        return None
    # ================================================================



    
    
    # open hdf file
    frss = h5py.File(hfile,"r")
    
    z = frss['geometry/sun_zenith']
    
    # satellite zenith
    zen = np.array(z)
    scale = z.attrs['scale_factor'][0]
    off = z.attrs['add_offset'][0]

    zen = zen * scale + off

    zen = np.deg2rad(zen)

    return np.cos(zen)


######################################################################
######################################################################


    

def get_msg_sat_zen(day, sat_type = 'rss'): 

    '''
    Reads satellite zenith angle from the partical satellite hdf file.


    Parameters
    ----------
    day : datetime object 
        day and time of MSG time slot

    sat_type : str, optional, default = 'rss'
        one of the two Meteosat scan types


    Returns
    --------
    satzen : numpy array
        satellite zenith angle

    '''

    # get the right satellite specification
    if sat_type in ('msg1','meteosat8','rss'):
        scan_type = 'rss'
        sat_name = 'msg1'


    elif sat_type in ('msg2','meteosat9','hrs'):
        scan_type = 'hrs'
        sat_name = 'msg2'

    else:
        'ERROR: unknown sat_type'
        sys.exit()
        

    # set time and date strings for file selection
    time_string = day.strftime('%Y%m%d%H%M')
    date_string = day.strftime('%Y/%m/%d')


    # set hdf file name
    arch_dir="/mnt/zephyr/u1/EUMETCAST_ARCHIVE/msg-sevi/%s/l15/hdf/eu/"% scan_type
    rssfile =arch_dir + date_string + "/%s-sevi-%s-l15-%s-eu.h5" % (sat_name, time_string, scan_type)


    # check if hdf file exists 
    if not h5py.is_hdf5(rssfile):
        print('No hdf-file available!')
        return
    
    
    # open hdf file
    frss = h5py.File(rssfile,"r")
    
    z = frss['geometry/satellite_zenith']
    
    # satellite zenith
    zen = np.array(z)
    scale = z.attrs['scale_factor'][0]
    off = z.attrs['add_offset'][0]

    satzen = zen * scale + off

    return satzen


######################################################################
######################################################################


def rad2refl(rad, info, day, coszen):

    '''
    Conversion of satellite radiances in reflectances.


    Parameters
    ----------
    rad : numpy array
        radiance of a visible channel

    info : dict
        dictionary which contains information about calibration coefficients

    day : datetime object 
        day and time of MSG time slot

    coszen : numpy array
        cosine of sun zenith angle


    Returns
    --------
    refl : numpy array
        reflectance

    '''

    # get the Julian day
    jday = int(day.strftime('%j'))

    # approximative formular for earth sun distance
    ESD = 1. - 0.0167 * np.cos( 2*np.pi*(jday-3.) / 365 )

    # top of atmosphere radiance
    TOARAD = info['f0'] / ESD**2 / np.pi

    return rad / TOARAD / coszen

######################################################################
######################################################################

def rad2bt(rad, info): # LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    
    '''
    Calculates brightness temperature of infrared channels given the 
    radiance and info list containing calibration coefficients.


    Parameters
    ----------
    rad : numpy array
        radiance of an infrared channel

    info : dict
        dictionary which contains information about calibration coefficients. 


    Returns
    --------
    bt : numpy array
        brightness temperature

    '''


    # get coefficients from meta info
    A = info['alpha']
    B = info['beta']
    nu_c = info['nu_c']

    # radiation constants
    C1 = 1.19104e-5
    C2 = 1.43877

    bt = (C2 * nu_c / np.log(C1*nu_c**3./rad + 1) - B) / A

    return bt # TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT



######################################################################
######################################################################


def get_msg_lsm(region,  scan_type = 'rss', arch_dir = 'none'): 

    # LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

    '''
    Reads land sea mask of Meteosat8 pixels saved in an
    auxiliary file in the MSG archive.


    Parameters
    ----------
    region :  str
        predefined cutout of MSG fulldisk, e.g. 'eu'

    scan_type : str
        which type of msg is used, 
        * ) allowed keys are: 'pzs',rss'
    
    arch_dir : str, optional, default = "none"
        archive directory where NWC-SAF files are saved 


    Returns
    --------
    lsm : numpy array
        land sea mask

    '''



# get the right file suffix ..........................................
    if scan_type == 'rss':
        file_suff = 'rss'
        eu_col = 1356

    elif scan_type == 'hrs':
        file_suff = 'hrs'
        eu_col = 1556

    elif scan_type == 'pzs':
        file_suff = 'pzs'
        eu_col = 1556
    else:
        'ERROR: unknown scan_type'
        sys.exit()
        


# hdf5 file ..........................................................
    if arch_dir == 'none':
        arch_dir = '%s/SEVIRI/' % os.environ['LOCAL_DATA_PATH']

    fname = arch_dir + "/auxdata/msevi-geolocation-" + file_suff + ".h5"



# open file ..........................................................
    f=h5py.File(fname,"r")


    lsm = f['land_sea_mask'].value


    if region == 'eu':
# dimensions of europe
        r1=156
        r2=r1+600
        c1=eu_col
        c2=c1+800

        lsm = lsm[r1:r2,c1:c2]

    elif region == 'full':
        pass

    else:
        print('UNKNOWN REGION!')


    return lsm # TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT


 
######################################################################
######################################################################



def get_msg_lon_lat(region, scan_type='rss', arch_dir = 'none'): # LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

    '''
    Reads longitudes and latitudes of Meteosat8 pixels saved in an
    auxiliary file in the MSG archive.


    Parameters
    ----------
    region : str
        predefined cutout of MSG fulldisk, e.g. 'eu'

    scan_type : str, optional, default = 'rss'
        one of the two Meteosat scan types

    arch_dir : str, optional, default = "none"
        archive directory where NWC-SAF files are saved 

    
    Returns
    --------
    lon : numpy array
        longitudes of pixels

    lat : numpy array
        latitudes of pixels

    '''
# get the right file suffix ..........................................
    if scan_type == 'rss':
        file_suff = 'rss'
        eu_col = 1356

    elif scan_type == 'hrs':
        file_suff = 'hrs'
        eu_col = 1556

    elif scan_type == 'pzs':
        file_suff = 'pzs'
        eu_col = 1556
    else:
        'ERROR: unknown scan_type'
        sys.exit()
        


# hdf5 file ..........................................................
    if arch_dir == 'none':
        arch_dir = '%s/SEVIRI/' % os.environ['LOCAL_DATA_PATH']

    fname = arch_dir + "/auxdata/msevi-geolocation-" + file_suff + ".h5"


# open file ..........................................................
    f=h5py.File(fname,"r")
    
# list content .......................................................
    lon = f['longitude'].value
    lat = f['latitude'].value


    if region == 'eu':
# dimensions of europe
        r1=156
        r2=r1+600
        c1=eu_col
        c2=c1+800
        
        lon = lon[r1:r2,c1:c2]
        lat = lat[r1:r2,c1:c2]

    elif region == 'full':
        pass
    else:
        print('UNKNOWN REGION!')

    return lon, lat # TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT



######################################################################
######################################################################


def channel_segment_sets(set_name):

    '''
    Outputs a predefined set of channels and segments'


    Parameters
    ----------
    set_name : str
        name of a predefined set
    

    Returns
    --------
    chan_seg : dict
        Dictionary of channels with a list of segments. 

    '''


    # standard segment range
    eu_segs = list(range(7,9))
    full_segs = list(range(1,9))

    sets = {}

    # European Area
    sets['ir-eu'] = {'IR_039':eu_segs, 'WV_062':eu_segs, 'WV_073':eu_segs, 'IR_087':eu_segs, \
                  '  IR_097':eu_segs, 'IR_108':eu_segs, 'IR_120':eu_segs, 'IR_134':eu_segs}
    
    sets['eu'] = {'VIS006':eu_segs, 'VIS008':eu_segs, 'IR_016':eu_segs,\
                    'IR_039':eu_segs, 'WV_062':eu_segs, 'WV_073':eu_segs, 'IR_087':eu_segs, \
                  '  IR_097':eu_segs, 'IR_108':eu_segs, 'IR_120':eu_segs, 'IR_134':eu_segs,
                    'HRV':list(range(20,24))}

    sets['nc-eu'] = {'VIS006':eu_segs, 'VIS008':eu_segs, 'IR_016':eu_segs}

    
    
    # full disk service

    sets['ir-full'] = {'IR_039':full_segs, 'WV_062':full_segs, 'WV_073':full_segs, \
                           'IR_087':full_segs, 'IR_097':full_segs, 'IR_108':full_segs,\
                           'IR_120':full_segs, 'IR_134':full_segs}
    
    sets['full'] = {'VIS006':full_segs, 'VIS008':full_segs, 'IR_016':full_segs,\
                    'IR_039':full_segs, 'WV_062':full_segs, 'WV_073':full_segs, 'IR_087':full_segs, \
                  '  IR_097':full_segs, 'IR_108':full_segs, 'IR_120':full_segs, 'IR_134':full_segs,
                    'HRV':list(range(1,24))}

    
    sets['nc-full'] = {'VIS006':full_segs, 'VIS008':full_segs, 'IR_016':full_segs}

    sets['ir108-full'] = { 'IR_108':full_segs }

    
    if set_name in sets:
        chan_seg = sets[set_name]
    else:
        print('ERROR: unknown channel-segment set!')
        print('   available keys:')
        print(list(sets.keys()))
        chan_seg = None
        
    
    return chan_seg


######################################################################
######################################################################



def get_HRIT_from_arch(day, chan_set = 'nc-full', scan_type = 'pzs' ,
                       arch_dir = None, out_path = None):
    
    '''
    Gets and decompresses the original HRIT files from archive. A
    typical set of channels and segments can be chosen.


    Parameters
    ----------
    day : datetime object 
        day and time of MSG time slot

    chan_set : str, optional, default = 'nc-full'
        name of a predefined channels set

    scan_type : str, optional, default = 'pzs'
        one of the two Meteosat scan types

    arch_dir : str, optional, default = None
        archive directory where NWC-SAF files are saved 

    out_path : str, optional, default = None
        directory where the HRIT files are saved.

    
    Returns
    --------
    file_list : list 
        list of extracted HRIT files

    '''

    

    # get the right satellite specification
    if scan_type not in ('rss', 'hrs'):
#        raise ValueError
        print('ERROR: ',scan_type,' is not a valid scan_type option') 
        return None

    # set default archive directory
    if not arch_dir:
        arch_dir = '/mnt/zephyr/u1/EUMETCAST_ARCHIVE/msg-sevi/%s/l15_hrit/' % scan_type

    # set output directory
    if not out_path:
        out_path = '/tmp/hrit' + str(np.random.randint(1e10))

    if not os.access(out_path, os.F_OK):
        os.mkdir(out_path)

    # set time and date strings for file selection
    time_string = day.strftime('%Y%m%d%H%M')
    date_string = day.strftime('%Y/%m/%d')


    # create filename template (NOTE that number of MSG have been left out!)
    filename_template = arch_dir + date_string + "/msg?-sevi-%s-l15-%s.tar" % (time_string, scan_type)

    # try to match template with existing file
    arch_file_list = glob.glob(filename_template) 
         

    if arch_file_list:
        arch_file = arch_file_list[0]
    else:
        print('ERROR: ', filename_template,' does not exist!')
        return None


    # build tar command of reading content
    tar_com = 'tar -tf '+ arch_file

    # get content list of tar file
    cont_list = subprocess.getoutput(tar_com).split()
    print(cont_list)
    
    # get channel-segment list for requested set
    chan_seg = channel_segment_sets(chan_set)

    # extract tar files
    extract_list =[]
    for HRIT_file in cont_list:
        
        # get prolog file
        if 'PRO' in HRIT_file:
            extract_list.append(HRIT_file)
            
        # get epilog file
        if 'EPI' in HRIT_file:
            extract_list.append(HRIT_file)
            
        # get channel segments
        for ch in chan_seg:
            for seg in chan_seg[ch]:
                if str(seg).zfill(6) in HRIT_file and ch in HRIT_file:
                    extract_list.append(HRIT_file)

    old_path = os.getcwd()

    os.chdir(out_path)

    hrit_list = []
    for hrit_file in extract_list:
        os.system('tar -xf '+ arch_file+' '+ hrit_file)

        # decompress the segments
        if '-C_' in hrit_file:
            os.system('xRITDecompress ' + hrit_file)
            os.system('rm -f ' + hrit_file)
        
            hrit_list.append(out_path + '/' + hrit_file.replace('C_','__'))

        # or do nothing for epi- and prolog
        else:
            hrit_list.append(out_path + '/' + hrit_file)
            

#    print hrit_list

    os.chdir(old_path)

    return hrit_list


######################################################################
######################################################################


def read_rad_from_hdf(hfile): # LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

    '''
    Reads radiances from synthetic satellite calculation which are 
    intermediately saved in an hdf file.


    Parameters
    ----------
    hfile : str
        file name of the hdf file where synthetic radiances are stored.

    
    Returns
    --------
    RAD : dict
        dictionary of radiances
 
        the channel name, e.g. IR_108, is used as key to access the
           data array (e.g. saved on COSMO grid) 
    '''


    # test if file exists
    if not h5py.is_hdf5(hfile):
        print('No hdf-file available!')
        return

    hf = h5py.File(hfile,'r')


    # read radiance field for each channel
    RAD = {}
    for chan in list(hf):
        RAD[chan] =  np.array(hf[chan])
    
    return RAD # TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT



######################################################################
######################################################################

