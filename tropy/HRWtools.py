#!/usr/bin/python

# Libraries
# ==========

import sys
import h5py
import numpy as np
import datetime



class HRW:
    ''' Class for High Resolution Winds (HRW) objects. 
    '''
    def __init__(self,day,color):
        self.day = day
        self.lon = np.array([])
        self.lat = np.array([])
        self.u =  np.array([])
        self.v =  np.array([])
        self.p =  np.array([])
        self.T = np.array([])
        self.col = color
        self.des = str()
        self.cloud_type =  np.array([])
        self.qi = np.array([])
        self.confidence = np.array([])



    def append(self,lon,lat,u,v,p,T,ct,qi,conf):
        self.lon = np.append(self.lon,lon)
        self.lat = np.append(self.lat,lat)
        self.u   = np.append(self.u,u)
        self.v   = np.append(self.v,v)
        self.p   = np.append(self.p,p)
        self.T   = np.append(self.T,T)
        self.cloud_type   = np.append(self.cloud_type,ct)
        self.qi   = np.append(self.qi,qi)
        self.confidence   = np.append(self.confidence,conf)


    def mask(self, m):

        # get key set
        ks = self.__dict__.keys()

        # get length
        l = len(self.lon)
        tp = type(self.lon)
        
        
        for k in ks:
            var = self.__dict__[k]

            if type(var) == tp:
                if len(var) == l:
                    self.__dict__[k] = var[m==True]

       
    def __getitem__(self,attr):
        return getattr(self,attr) 
######################################################################
######################################################################


def plot_HRW_to_map(hrw, seq, m):
    '''
    Plots a sorted HRW object onto a map using colored wind barbs.


    USAGE:
    =====
    plot_HRW_to_map(hrw, seq, m)


    INPUT:
    ======
    hrw: level-sorted dictionary of HRW object
    seq: sequence of levels to be plotted
    m: Basemap map object
 

    OUTPUT:
    =======
    None

    '''

    # rows 
    r = [0.075, 0.045, 0.015]
  
    # and columns
    c = [0.3, 0.5, 0.7]
    

    text_pos={150:(c[0],r[0]),250:(c[1],r[0]),350:(c[2],r[0]),\
              450:(c[0],r[1]),550:(c[1],r[1]),650:(c[2],r[1]),\
              750:(c[0],r[2]),850:(c[1],r[2]),950:(c[2],r[2])}
    
    for l in seq:
         
        try:
            h = hrw[l]

            x,y = m(h.lon,h.lat)

            # plot barbs
            m.barbs(x,y,h.u,h.v,length=4.5, pivot='middle',barbcolor=h.col,lw=0.5)

            # print shadow of text labels
            d = 0.001
            figtext(text_pos[l][0]+d,text_pos[l][1]-d,h.desc,color='black',size=12,ha='center')

            # print colored text labels
            figtext(text_pos[l][0],text_pos[l][1],h.desc,color=h.col,size=12,ha='center')
        except:
            pass
    
    
    return
#
######################################################################
######################################################################

def read_HRW(day, chan_list = 'all',\
                 arch_dir='/u1/home/fabian/proj/2012-06_nwcsaf_exp/hrw-tests/export/PGE09/', \
                 sat_type='MSG1',\
                 lon_min=0, lon_max=20, lat_min=40, lat_max=60):

    '''
    Reads an HRW object from the NWCSAF output hdf5 file.


    USAGE:
    =====
    hrw = read_HRW(day, chan_list = 'all', 
             arch_dir='/u1/home/fabian/proj/2012-06_nwcsaf_exp/hrw-tests/export/PGE09/', 
             lon_min=0, lon_max=20, lat_min=40, lat_max=60)


    INPUT:
    ======
    day = time slot as datetime object
    chan_list = list of HRW file keys (look into the hdf5 for key names, OPTIONAL
    arch_dir = directory where NWCSAF-HRW file is saved, OPTIONAL
    sat_type = type of satellite used, i.e. either MSG1 or MSG2, OPTIONAL
    lon_min = minimal longitude saved in HRW object, OPTIONAL
    lon_max = maximal longitude saved in HRW object, OPTIONAL
    lat_min = minimal latitude  saved in HRW object, OPTIONAL
    lat_max = maximal latitude  saved in HRW object, OPTIONAL
    
     

    OUTPUT:
    =======
    hrw: HRW object

    '''
    
    time_string = day.strftime('%Y%m%d%H%M')
    date_string = day.strftime('%Y/%m/%d')

    if sat_type == 'MSG1':
        region_string = '_rss-eu______'
    elif sat_type == 'MSG2':
        region_string = '_hrs-eu______'
    else:
        print 'ERROR: unknown sat_type in read_HRW!'
        sys.exit()
        

        

    prod_file=arch_dir + 'SAFNWC_'+ sat_type +'_HRW__'+time_string + region_string + '.buf.h5'


    # check if file exists
    if not h5py.is_hdf5(prod_file):
        print 'No hdf-file ',prod_file ,' available!'
        return


    # open file
    fprod=h5py.File(prod_file,"r")

    
    hrw = HRW(day, 'black')

    if chan_list == 'all':
        chan_list = fprod.keys()


    for ch in chan_list:
        # get values
        val = fprod[ch].value
        

        # get longitude and latitude
        lon = val['lon']
        lat = val['lat']

        # get speed and direction
        U, al =  val['wind_speed'], val['wind_direction']

        # pressure
        p = val['pressure']

        # temperature
        T = val["temperature"]

        # cloud type
        ct = val["cloud_type"]
        
        
        # quality index
        qi = val["applied_QI"]

        # quality index
        conf = val["confidence"]


    
        # get region mask
        lon_mask = np.logical_and(lon_min < lon, lon < lon_max)
        lat_mask = np.logical_and(lat_min < lat, lat < lat_max)
        mask = np.logical_and(lon_mask,lat_mask)

        # mask values        
        lon = lon[mask]
        lat = lat[mask]

        p = p[mask]

        al = al[mask]
        U = U[mask]

        # angle: deg to rad
        al = np.pi * al / 180.

        T = T[mask]
        ct = ct[mask]
        qi = qi[mask]
        conf = conf[mask]

        # zonal and meridional wind
        u = -U * np.sin(al)
        v = -U * np.cos(al)

        hrw.append(lon,lat,u,v,p,T,ct,qi,conf)               

    fprod.close()
    return hrw
#

######################################################################
######################################################################



def sort_HRW(hrw_all):

    '''
    Sorts the HRW object regarding to different pressure levels.


    USAGE:
    =====
    hrw = sort_HRW(hrw_all)


    INPUT:
    ======
    hrw_all : HRW object containing all HRW points 
    
     

    OUTPUT:
    =======
    hrw: dictionary of HRW objects sorting by pressure levels 

    '''
    

    clev = [ 150,  250,  350,  450,  550,  650,  750,  850,  950]
    ccol = ['magenta','red','orange','yellow','YellowGreen','DarkCyan','cyan','blue','navy']


    # get values ........................
        

    # get longitude and latitude
    day = hrw_all.day

    lon, lat = hrw_all.lon, hrw_all.lat

    # winds
    u, v = hrw_all.u, hrw_all.v

    # pressure
    p = hrw_all.p 

    T = hrw_all.T
    ct = hrw_all.cloud_type
    qi = hrw_all.qi
    conf = hrw_all.confidence
        
    # get corresponding level index, so that p [hPa] is in a 100hPa around level l
    ilev = np.array(np.floor((p/100-100)/100),dtype='int8')
        

    hrw={}
    for i,l in enumerate(clev):
        hrw[l] = HRW(day, ccol[i])
        hrw[l].desc = str(l-50) + '-' + str(l+49) + ' hPa'


    for n, i in enumerate(ilev):
        lev = clev[i]
        hrw[lev].append(lon[n],lat[n],u[n],v[n],p[n],T[n],ct[n],qi[n],conf[n])
        

    return hrw
#
