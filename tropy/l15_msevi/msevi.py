#!/usr/bin/env python

######################################################################
######################################################################
#
# DESCRIPTION
# ===========
# This file contains the MSevi Container Class.
# 
# It allows reading MSG-SEVIRI either from HRIT or from hdf files
# for a given time slot, region and scan_type.
#
#
######################################################################
######################################################################
   

# load libraries -----------------------------------------------------
import sys, os, glob

import numpy as np
import datetime as dt
import pylab as pl
# import Image
from .msevi_config import _predefined_regions, _narrow_channels,\
    _channel_config, _calibration_constants, _hdf_cutout_regions,\
    _arch_config
# use pyorbital from pytroll group for sun zenith angle calculations
import pyorbital.astronomy

import tropy.MSGtools as MSGtools
import tropy.io_tools.HRIT as hio
import tropy.io_tools.hdf as iohdf
from tropy.io_tools.find_latest_slot import find_latest_slot
from tropy.analysis_tools.grid_and_interpolation import make_hrv_upscaling


# ====================================================================
   


class MSevi(object):
    ''' 
     Container for MSG SEVIRI data. 
    '''


    def __init__(self, \
        time = dt.datetime(2012, 12, 20, 12, 0), \
        region = 'eu',\
        scan_type = 'rss', \
        chan_list = [],\
        **kwargs):

        '''
        Description:
        ============

        Initialization of the MSevi Data Container Class.


        Input keywords arguments:
        =========================

        time: either a predefined string ('now', newest', 'latest') or 
              a datetime object containing the time slot

        region: either a prefined string ('eu', 'full', 'de') or
                a tuple integers ((nr1, nr2), (nc1, nc2)) defining the region which 
                should be retrieved

                ((nr1, nr2), (nc1, nc2)) define the rows and columns of a cutout 
                of a 2d field f as f[nr1:nr2, nc1:nc2]

        scan_type: 'rss' or 'hrs' are used as short cuts from the rapid and the 
                   operational scan, respectively

        chan_list: list of channel names which should be retrieved, 
                   e.g. ['IR_108', 'IR_120']

                   Names are defined in the dict _narrow_chanels.

        Comments:
        ========
        Channels are loaded per default by initialization.

        '''

        #--- [1] Init variables -------------------------------------- 


        # requested scan type, i.e. 'hrs' or 'rss'
        self.scan_type  = scan_type.lower()
        self.sat_type = 'Unknown'


        # requested time slot as datetime object
        if time in ('latest', 'newest', 'now'):
            self.time = find_latest_slot(scan_type=self.scan_type)
        else:
            self.time = time


        # init channel-segment-list
        self.chan_list = []
        self.chan_seg = {}


        # init data container ........................................
        #  rad = radiance
        #  ref = reflectance
        #  bt  = brightness temperature

        self.rad = {}
        self.ref = {}
        self.bt = {}     

        # meta data container ........................................
        self.meta = {}
        #=============================================================


        #---  [2] get region -----------------------------------------
   
        # check for user-defined region 
        if type(region) == type(()):
            self.region = region

        # check for predefined regions
        elif type(region) == type(''):
            self.reg_name = region
            reg_type = region + '-' + self.scan_type
 
            if reg_type in _predefined_regions:
                self.region = _predefined_regions[reg_type]
  
            else:
                print('unknown region');return # sys.exit()

        # calculate HRV region
        # ====================
        # It is assumed that the hrv region is a dependent variable which 
        # depends on the low res region attribute
        # That means: Given the low res region cutout hrv region is determined 
        # as corresponding cutout which 
        # (i) exactly fits with the low res region and 
        # (ii) refers to the artifical hrv full disk of 11136 x 11136
        #
        # low res pixel with the index (I,J) = (0,0)
        # has a high res pixel with index (i_m, j_m) = (2,2) in the middle
        # and (i_u, j_u) = (1,1) in the upper corner
        # see doc "MSG Level 1.5 Image Data Format Description", Fig.8
        # 
        #  which then leads to the relation between the low res pixel (I,J)
        # and its corresponding upper corner high res pixel (i_u, j_u):
        #  (i_u, j_u) = 3 * (I, J ) + 1


        self.hrv_region = low2hres(self.region)
        #=============================================================


        #---  [3] initially load some channels -----------------------
        if chan_list:
            self.load(chan_list, **kwargs)
        #=============================================================

        return

######################################################################

    def __repr__(self):
        return '''
    MSevi Data Container
    ====================
    Time Slot:       %s
    Region:          %s
    Scan Type:       %s
    Satellite Type:  %s
    Loaded Channels: %s
    ''' % (str(self.time), str(self.region),\
          str(self.scan_type), str(self.sat_type),\
          str(self.chan_list))

######################################################################

    def check_channels(self, chan_list):

        # get channel which have been already loaded
        old_chan_list = self.chan_list 

        # 1st check if list is used for channel argument
        if type(chan_list) == type(''):
            chan_list = [chan_list]

        # 2nd check if channel naming is correct
        new_chan_list = []
        for ch in chan_list:

            if ch in _channel_config:
                new_chan_list.append(ch)

            else:
                print('ERROR: Channel name is wrong!')
                print('       Use one of these:')
                print(sorted(_channel_config.keys()))
                return

        # create for additional channels to be loaded
        list_for_load = [c for c in new_chan_list if c not in old_chan_list]

        return list_for_load

######################################################################

    def load(self, chan_list, **kwargs):
        '''
        The load method checks if channels are loaded and then loads
        the ones which are not already there.


        The load method is an interface to an hdf and hrit load method.
        If the region matches, it tries 1st to load data from the hdf 
        file. In the case of no success, the hrit loader is started.

        Input arguments:
        ================
        chan_list: list of channel names

        
        '''

        #--- [1] selection of channels for loading -------------------
        list_for_load = self.check_channels(chan_list)
        #=============================================================


        #--- [2] load channels ---------------------------------------
        if list_for_load:
            try:
                self.load_hdf(list_for_load, **kwargs)

            except IOError:
                self.load_hrit(list_for_load, **kwargs)

        else:
            print(chan_list, 'is already loaded!')
        #=============================================================


        #--- [3] book-keeping of loaded channels ---------------------
        self.chan_list = list(self.rad.keys())
        self.chan_list.sort()
        #=============================================================

        return


######################################################################

    def load_hdf(self, list_for_load, **kwargs):
        
        # testb ++++++++++++++++++++++++++++++
        # work around for stupid hdf bug !!!! 
        # if self.scan_type == 'hrs':
        #    raise IOError
        # teste ++++++++++++++++++++++++++++++ 

        #--- [1] check if requested region is within the hdf region range
        
        (nr1, nr2), (nc1, nc2) = self.region

        (hr1, hr2), (hc1, hc2) = _hdf_cutout_regions[self.scan_type]

        row_cond = hr1 <= nr1 and hr2 >= nr2
        col_cond = hc1 <= nc1 and hc2 >= nc2

        within_hdf_region = row_cond and col_cond


        if not within_hdf_region:
            raise IOError
        else:
            print('Region suggests use of hdf file') 
        #=============================================================


        #--- [2] load channels ---------------------------------------
        self.fname = self.msevi_fname(data_type = 'hdf')

        rad = iohdf.get_seviri_chan(list_for_load, 
                                    self.time, 
                                    fname = self.fname,
                                    scan_type = self.scan_type, 
                                    add_meta=True)
        if rad != None:
            self.rad.update(rad)
        else:
            raise IOError
        #=============================================================


        #--- [3] write meta data in corresponding attribute ----------
        self.meta  = self.rad['meta']
            
        del self.rad['meta']

        # get sat_type from prolog file ..............................
        if self.sat_type == 'Unknown':
            
            ch_name = list(self.rad.keys())[0]
            
            # now I try to get the name of the attribute which 
            # specifies the satellite id, this is a bit more 
            # complicated than necessary because Hartwig incooperated
            # some spelling error in older hdf files
            attr_keys = self.meta['attributes'][ch_name]

            for k in attr_keys:
                if 'sa' in k and 'id' in k:
                    sat_id_attr = k

            
            # retrieve sat id from attributes list of 1st channel
            try:
                sat_id = self.meta['attributes']['general']['satellite_id']
            except:
                sat_id = self.meta['attributes'][ch_name][sat_id_attr]


            # convert sat id into MSG name
            self.sat_type = 'MSG' + str(sat_id[0] - 320)
        #=============================================================
        
       
        #--- [4] cutout region ---------------------------------------
        for ch in list_for_load:
            self.cutout_region_hdf(ch)
        #=============================================================

        return


######################################################################

    def load_hrit(self, list_for_load, **kwargs):

        

        #--- [1] create channel-segment sets to be loaded ------------
        chan_seg_for_load = {}
        for c in list_for_load:
             if c == 'HRV':
                 chan_seg_for_load[c] = hio.segments_for_region(self.hrv_region,\
                     **_channel_config[c])          
             else:
                 chan_seg_for_load[c] = hio.segments_for_region(self.region,\
                     **_channel_config[c])          
          
        self.chan_seg.update(chan_seg_for_load)
        #=============================================================


        #--- [2] read HRIT data --------------------------------------

        self.fname = self.msevi_fname(data_type = 'hrit')

        out_path = '/tmp/hrit' + str(np.random.randint(1e10))
        
        try:
            rad = hio.read_HRIT_data(self.time, chan_seg_for_load,
                                     scan_type = self.scan_type, 
                                     add_meta = True, 
                                     tarname =  self.fname,
                                     out_path = out_path,
                                     **kwargs)
        except:
            os.system('rm -rf %s' % out_path)
            raise IOError

        
        if rad != None:
            self.rad.update(rad)
        else:
            raise IOError

        #=============================================================
 

        #--- [3] write meta data in corresponding attribute ----------
        meta_new  = self.rad['meta']
            

        # rewrite meta data ..........................................
        
        # TODO:
        # * combine meta data of different loads
 
        self.meta = meta_new
        
        del self.rad['meta']

        # get sat_type from prolog file ..............................
        if self.sat_type == 'Unknown':
            pfile = self.meta['prolog_file']

        # position of MSG string in filename
            imsg = pfile.find('MSG')
 
        # extract MSG? name
            self.sat_type = pfile[imsg:imsg+4]
        #=============================================================

       
        #--- [4] cutout region ---------------------------------------
        for ch in list_for_load:
            self.cutout_region_hrit(ch)
        #=============================================================

        return

   #####################################################################

    def msevi_fname(self, data_type = 'hdf'):

        d = {}
        d['scan_type'] = self.scan_type
        d['time_string'] = self.time.strftime(_arch_config['time_string'])
        d['date_string'] = self.time.strftime(_arch_config['date_string'])
        try:
            d['region'] = self.reg_name
        except:
            d['region'] = 'eu'
    
        arch_dir = _arch_config['arch_dir'].format(**d)

        if data_type == 'hdf':
            subpath = _arch_config['hdf_subpath'].format(**d)
            fname = _arch_config['hdf_fname'].format(**d)

        elif data_type == 'hrit':
            subpath = _arch_config['hrit_subpath'].format(**d)
            fname = _arch_config['hrit_tarname'].format(**d)
    
         

        return '%s/%s/%s' % (arch_dir, subpath, fname)


#####################################################################

    def cutout_region_hrit(self, ch):
 
 
        # select HRV or non-HRV region ------------------------------- 
        if ch == 'HRV': 
            region = self.hrv_region
        else:
            region = self.region

        (rowmin,rowmax),(colmin,colmax) = region
        # ============================================================

     
        # calculate the row shift due to missing top segments --------

        # get the most upper segment used
        iseg_max = np.array(self.chan_seg[ch]).max()

        # get the number of segments at all
        Nseg = _channel_config[ch]['Nseg']
        
        # and segment size
        seg_size = _channel_config[ch]['seg_size']

        # calculate row shift
        row_shift = (Nseg - iseg_max) * seg_size

        # and apply it
        rowmin -= row_shift
        rowmax -= row_shift
        # ============================================================
        
        # cut out region --------------------------------------------- 
        self.rad[ch] = self.rad[ch][rowmin:rowmax,colmin:colmax]
        # ============================================================

        return



#####################################################################

    def cutout_region_hdf(self, ch):

 
        # select HRV or non-HRV region ------------------------------- 
        if ch == 'HRV': 
            region = self.hrv_region
            hdf_cutout = low2hres(_hdf_cutout_regions[self.scan_type])
        else:
            region = self.region
            hdf_cutout = _hdf_cutout_regions[self.scan_type]
        # ============================================================

        
        # get definition of region 
        (rowmin,rowmax),(colmin,colmax) =region


        # get hdf region defintion
        (hr1, hr2), (hc1, hc2) = hdf_cutout


        # convert absolute region definition to relative defintion
        r1 = rowmin - hr1
        r2 = rowmax - hr1

        c1 = colmin - hc1
        c2 = colmax - hc1


        # cut out region --------------------------------------------- 
        self.rad[ch] = self.rad[ch][r1:r2,c1:c2]
        # ============================================================

        return

######################################################################


    def lonlat(self):

        # retrieve geo reference from aux file
        lon, lat = MSGtools.get_msg_lon_lat('full', scan_type = self.scan_type)

        # get definition of region 
        (rowmin,rowmax),(colmin,colmax) = self.region

        # cut out region --------------------------------------------- 
        self.lon = lon[rowmin:rowmax,colmin:colmax]
        self.lat = lat[rowmin:rowmax,colmin:colmax]
        # ============================================================

        self.hlon = make_hrv_upscaling(self.lon)
        self.hlat = make_hrv_upscaling(self.lat)

        return

#####################################################################

    def landsea(self):

        # retrieve geo reference from aux file
        lsm = MSGtools.get_msg_lsm('full', scan_type = self.scan_type)

        # get definition of region 
        (rowmin,rowmax),(colmin,colmax) = self.region

        # cut out region --------------------------------------------- 
        self.lsm = lsm[rowmin:rowmax,colmin:colmax]
        # ============================================================

        return

#####################################################################

    def sunzen(self):
        
        d =  self.__dict__

        if 'lon' in d and 'lat' in d:
            self.szen = pyorbital.astronomy.sun_zenith_angle(self.time, self.lon, self.lat)
        
        if 'hlon' in d and 'hlat' in d:
            self.hszen = pyorbital.astronomy.sun_zenith_angle(self.time, self.hlon, self.hlat)


        return

#####################################################################

    def rad2refl(self, chan = None):

        '''
        Calculates reflectance of solar channels given the 
        radiance and info list containing calibration coefficients.
    
        '''
     
        # get the Julian day
        jday = int(self.time.strftime('%j'))
 
        # approximative formular for earth sun distance
        esd = 1. - 0.0167 * np.cos( 2*np.pi*(jday-3.) / 365 )
 
 
        f0 = _calibration_constants[self.sat_type]['f0']
 
        # top of atmosphere radiance
        # toa_rad = f0 / esd**2 
 
 
        # check for which list reflectance is calculated --------------
        if chan:
           if type(chan) == type(''):
               work_list = [chan]
           elif type(chan) == type([]):
               work_list = chan
        else:
           work_list = list(f0.keys())
        # =============================================================
 
 
        # do calculation ----------------------------------------------
        for ch in work_list:
            if ch in self.chan_list:
                self.ref[ch] = self.rad[ch] / (f0[ch] / esd**2)
        # =============================================================
   
        return
######################################################################
 
    def rad2bt(self, chan = None):
        
        '''
        Calculates brightness temperature of infrared channels given the 
        radiance and info list containing calibration coefficients.
    
        '''
    
    
        # get coefficients from meta info
        A = _calibration_constants[self.sat_type]['A']
        B = _calibration_constants[self.sat_type]['B']

        nu_c = _calibration_constants[self.sat_type]['nu_c']


        # radiation constants
        C1 = 1.19104e-5
        C2 = 1.43877
 
        # check for which list reflectance is calculated --------------
        if chan:
           if type(chan) == type(''):
               work_list = [chan]
           elif type(chan) == type([]):
               work_list = chan
        else:
           work_list = list(A.keys())
        # =============================================================
 
 
        # do calculation ----------------------------------------------
        for ch in work_list:
            if ch in self.chan_list:
                self.bt[ch] = np.where(self.rad[ch] == 0, 0,
                                       (C2 * nu_c[ch] / \
                                        np.log(C1 * nu_c[ch]**3./self.rad[ch] + 1) - B[ch]) / A[ch])
        # =============================================================
     
    
        return
    
    

       
######################################################################
######################################################################

def low2hres(region):

        # calculate HRV region
        # ====================
        # It is assumed that the hrv region is a dependent variable which 
        # depends on the low res region attribute
        # That means: Given the low res region cutout hrv region is determined 
        # as corresponding cutout which 
        # (i) exactly fits with the low res region and 
        # (ii) refers to the artifical hrv full disk of 11136 x 11136
        #
        # low res pixel with the index (I,J) = (0,0)
        # has a high res pixel with index (i_m, j_m) = (2,2) in the middle
        # and (i_u, j_u) = (1,1) in the upper corner
        # see doc "MSG Level 1.5 Image Data Format Description", Fig.8
        # 
        #  which then leads to the relation between the low res pixel (I,J)
        # and its corresponding upper corner high res pixel (i_u, j_u):
        #  (i_u, j_u) = 3 * (I, J ) + 1


        (r1, r2), (c1, c2) = region 
        
        hrv_region = ((3*r1 + 1, 3*r2 + 1), (3*c1 + 1, 3*c2 + 1))
    
        return hrv_region

######################################################################
######################################################################


if __name__ == '__main__':

        cin = { \
         'time' : dt.datetime(2012,5,23,12,0), \
        'region': 'eu',\
         'scan_type' : 'pzs',\
        'chan_list' : ['VIS008'],\
    #['VIS006', 'VIS008','IR_016']
             }
    
        cin.update( {     'time' : dt.datetime(2013,6,20,12,0), \
         'scan_type' : 'rss',\
                    })

        scene = MSevi(**cin)

       

