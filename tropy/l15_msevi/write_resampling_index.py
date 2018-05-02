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


# scan_type = 'rss'
scan_type = 'hrs'

lon, lat = MSGtools.get_msg_lon_lat('eu', scan_type)

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

nr, nc = lon.shape
ir = np.arange(nr)
ic = np.arange(nc)

ind = {}
ind['ic'], ind['ir'] = np.meshgrid(ic, ir)

print 'Start Resampling'

ind_new = {}
for k in ind:

    data = ind[k]
    msg_con_nn = image.ImageContainerNearest(data, grid_def, radius_of_influence=50000)
    area_con_nn = msg_con_nn.resample(area_def)
    
    ind_new[k] =  area_con_nn.image_data


ifile='proj_index_germ_%s.h5' % scan_type
iohdf.save_dict2hdf(ifile, ind_new)
