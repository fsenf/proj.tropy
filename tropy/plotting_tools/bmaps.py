#!/usr/bin/env python
# Collection of Basemap maps

import sys
from mpl_toolkits.basemap import Basemap
import numpy as np
import pylab as pl
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

from standard_config import *

#####################################################################
#####################################################################


def germ_map(lon, lat, mp = None, xy = None):

    if not mp:
        mp = Basemap(projection='stere', lon_0=10, lat_0=90, lat_ts=60, \
                    llcrnrlon = 3.5889, llcrnrlat = 46.9526, \
                    urcrnrlon = 15.7208, urcrnrlat = 54.7405, \
                    resolution='i')
    
    
    # Convert the lat/lon values to x/y projections
    if not xy:
        xy = mp(lon,lat)
    
    mp.drawcoastlines(linewidth=1.5,color='black')
    mp.drawcountries(linewidth=1.5,color='black')
    mp.drawparallels(np.arange(40.,60.,1.),labels=[1,0,0,0])
    mp.drawmeridians(np.arange(-10.,20.,2.),labels=[0,0,0,1])
    
    return xy, mp
   
#####################################################################
#####################################################################


def narval_map(lon, lat, mp = None, xy = None):

    if not mp:
        mp = Basemap(projection='mill',
                     llcrnrlat = -12,
                     urcrnrlat = 22,
                     llcrnrlon = -70,
                     urcrnrlon = 17,
                     resolution='i')
            
   
    
    # Convert the lat/lon values to x/y projections
    if not xy:
        xy = mp(lon,lat)
    
    mp.drawcoastlines(linewidth = 1, color='black')

    mp.drawparallels( np.arange(-10., 20., 5.), linewidth = 0, labels=[1,0,0,0])
    mp.drawmeridians( np.arange(-60., 30., 15.),linewidth = 0, labels=[0,0,0,1])
    
    return xy, mp
   
#####################################################################
#####################################################################


def make_map(lon, lat, mp = None, xy = None, region = 'germ', flex_kws = None, 
             states = False, kreise = False, merid = True, paral = True,
             leipzig = False, color = 'black', lw = 1.5):


    if region == 'germ':
        if not mp:
            mp = Basemap(projection='stere', lon_0=10, lat_0=90, lat_ts=60, \
                             llcrnrlon = 3.5889, llcrnrlat = 46.9526, \
                             urcrnrlon = 15.7208, urcrnrlat = 54.7405, \
                             rsphere = (6378169., 6356583.8),\
                             resolution='i')
        dlat = 1.
        dlon = 2.

    elif region == 'meu':
        if not mp:
            mp = Basemap(width=1500e3, height=1500e3, resolution='i',
                         projection='stere', lat_0=52, lon_0=10, lat_ts=52,
                         rsphere = (6378169., 6356583.8),\
                             )
    
        dlat = 3.
        dlon = 3.


    elif region == 'europe':
        if not mp:
            mp = Basemap(width=4000e3, height=4000e3, resolution='i',
                         projection='stere', lat_0=52, lon_0=10, lat_ts=52,
                         rsphere = (6378169., 6356583.8),\
                             )
    
        dlat = 10.
        dlon = 15.

    elif region == 'cosmo-de':
        if not mp:
            mp = Basemap(projection='stere', lon_0=10, lat_0=90, lat_ts=60, \
                             llcrnrlon = 2, llcrnrlat = 44, \
                             urcrnrlon = 21.5, urcrnrlat = 56.5, \
                             rsphere = (6378169., 6356583.8),\
                             resolution='i')
        dlat = 3.
        dlon = 3.


    elif region == 'sax':
        if not mp:
            mp = Basemap(projection='stere', lon_0=10, lat_0=90, lat_ts=60, \
                             llcrnrlon = 10, llcrnrlat = 50, \
                             urcrnrlon = 16., urcrnrlat = 52.5, \
                             rsphere = (6378169., 6356583.8),\
                             resolution='h')
        merid = False
        paral = False
        states = True
        leipzig = True

    elif region == 'flex':
            
        if not mp:
            mp = Basemap(projection='stere',  lat_ts = 60, \
                         rsphere = (6378169., 6356583.8),\
                         resolution='h', **flex_kws)

        dlat = 1.
        dlon = 2.


    elif region == 'narval':
        
        if not mp:
            mp = Basemap(projection='cyl',
                         llcrnrlat = -12,
                         urcrnrlat = 22,
                         llcrnrlon = -70,
                         urcrnrlon = 17,
                         resolution='i')
            
        dlat = 5.
        dlon = 15.




    else:
        print 'ERROR: unknown map region'
        sys.exit()


    # Convert the lat/lon values to x/y projections
    if not xy:
        xy = mp(lon,lat)
    
    mp.drawcoastlines(linewidth=lw, color=color)
    mp.drawcountries(linewidth=lw, color=color)
    if states:
        shfile = '%s/gshhs/vg2500_geo84/vg2500_bld' % local_data_path
        mp.readshapefile(shfile, 'states', linewidth = lw, color=color)

    if kreise:
        shfile = '%s/gshhs/vg2500_geo84/vg2500_krs' % local_data_path
        mp.readshapefile(shfile, 'kreise', linewidth = 0.5, color='darkblue')
        


    if leipzig:

        shfile = '%s/gshhs/vg2500_geo84/vg2500_krs' % local_data_path
        mp.readshapefile(shfile, 'kreise', linewidth = 0.5, drawbounds = False)

        patches   = []

        for info, shape in zip(mp.kreise_info, mp.kreise):
            shnum =  info['SHAPENUM']
            if shnum == 363:
                patches.append( Polygon(np.array(shape), True) )

            if  shnum >= 353 and shnum <= 365:
                x, y = zip(*shape) 
                mp.plot(x, y, marker = None, color = 'k', lw = 0.5)
                    

        ax = pl.gca()
        ax.add_collection(PatchCollection(patches, facecolor= 'darkred', 
                                          edgecolor='k', linewidths=1., zorder=2, alpha = 0.3))

        

    if paral:
        mp.drawparallels(np.arange(-90.,90.,dlat),labels=[1,0,0,0])
    if merid:
        mp.drawmeridians(np.arange(-180.,120.,dlon),labels=[0,0,0,1])
    
    return xy, mp

#####################################################################
#####################################################################



if __name__ == '__main__':
    
    import pylab as pl
    import numpy as np  

    pl.clf()

    ilon = np.linspace(0, 17, 200)
    jlat = np.linspace(45, 55, 250)

    lat, lon = np.meshgrid(jlat, ilon) 

#    (x,y),mp = germ_map(lon, lat)
    (x,y),mp = make_map(lon, lat, region = 'sax')

    pl.show()
