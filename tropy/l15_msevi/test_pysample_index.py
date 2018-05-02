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
import io_tools.hdf as iohdf



scan_type = 'rss'
scan_type = 'hrs'


rgb_list = ['pytroll_nc', 'am']

time = dt.datetime(2011,6,6,12,30)
cin = { \
    'time' : time,\
	'scan_type' : scan_type,\
        'rgb_type':rgb_list,\
        'region':'eu'\
        }

scene = MSeviRGB(**cin)

scene.lonlat()
lon, lat = scene.lon, scene.lat

grid_def = geometry.GridDefinition(lons=lon, lats=lat)


projection = {
    'a': '6378144.',
    'b': '6356795.',
    'lat_0': '90.',
    'lon_0': '5.',
    'lat_ts': '50.',
    'proj': 'stere',
    'ellps': 'bessel'
    }

area_extent = (-155100.436345, -4441495.37946, 868899.563655, -3417495.37946)

nx, ny = 1024, 1024

area_def = geometry.AreaDefinition('germ', 'Germany Stereographic Projection', 'germ_stere',\
                                       projection, nx, ny, area_extent)


bmap = plot.area_def2basemap(area_def, resolution = 'i')


ifile='proj_index_germ_%s.h5' % scan_type
ind = iohdf.read_dict_from_hdf(ifile)

ic, ir = ind['ic'], ind['ir']

print 'Start Resampling'
for rgb_str in scene.images.keys():
    img = scene.images[rgb_str]

    rgb = np.array(img)
    rgb_new = rgb[ir, ic, :]

    pl.figure()

    bmap.imshow(rgb_new[::-1]/255.)
    bmap.drawcoastlines()

    pname = 'test_index_%s.png' % rgb_str 

    pl.title(scene.time.strftime('%Y-%m-%d %H:%M UTC'))
    pl.savefig(pname, dpi = 200)
    pl.close()
