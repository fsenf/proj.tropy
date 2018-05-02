#!/usr/bin/python

import numpy as np
import sys
sys.path.insert(0,'/u1/home/fabian/lib/tropy')
import MSGtools
from pyresample import geometry, image, plot
import datetime as dt
import pylab as pl
import Image
from msevi_rgb import MSeviRGB


rgb_list = ['pytroll_nc', 'am']

time = dt.datetime(2011,6,6,12,30)
cin = { \
    'time' : time,\
#	'scan_type' : 'rss',\
	'scan_type' : 'hrs',\
        'rgb_type':rgb_list,\
        'region':'eu'\
        }

scene = MSeviRGB(**cin)

scene.lonlat()
lon, lat = scene.lon, scene.lat

grid_def = geometry.GridDefinition(lons=lon, lats=lat)

nx, ny = 1200, 1200
area_def = geometry.AreaDefinition('areaD', 'Europe (3km, HRV, VTC)', 'areaD',\
                                       {'a': '6378144.0', 'b': '6356759.0',\
                                            'lat_0': '50.00', 'lat_ts': '50.00',\
                                            'lon_0': '8.00', 'proj': 'stere'},\
                                       nx, ny,\
                                [-1370912.72, -909968.64,  1029087.28, 1490031.36])



bmap = plot.area_def2basemap(area_def, resolution = 'i')

print 'Start Resampling'
for rgb_str in scene.images.keys():
    img = scene.images[rgb_str]

    rgb = np.array(img)
    rgb_new = np.zeros((nx,ny,3))

    for i in range(3):
#        print i
        data = rgb[:,:,i]
        msg_con_nn = image.ImageContainerNearest(data, grid_def, radius_of_influence=50000)
        area_con_nn = msg_con_nn.resample(area_def)
    
        rgb_new[:,:,i] = area_con_nn.image_data

    pl.figure()

    bmap.imshow(rgb_new[::-1]/255.)
    bmap.drawcoastlines()

    pname = 'test_%s.png' % rgb_str 

    pl.title(scene.time.strftime('%Y-%m-%d %H:%M UTC'))
    pl.savefig(pname, dpi = 200)
    pl.close()
