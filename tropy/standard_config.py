#!/usr/bin/env python

import os

try:
    local_home_path  =  os.environ['LHOME']
    proj_path = '%s/proj' % local_home_path
    pics_path = '%s/pics' % local_home_path
except:
    pass

try:
    global_data_path =  os.environ['GLOBAL_DATA_PATH']
    cosmo_path = '%s/cosmo/de' % global_data_path
except:
    pass


try:
    local_data_path  =  os.environ['LOCAL_DATA_PATH']
except:
    pass

try:
    meteosat_archive_path =  os.environ['METEOSAT_ARCHIVE_PATH']
except:
    pass



