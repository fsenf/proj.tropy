#!/usr/bin/env python

import os

local_home_path  =  os.environ['LHOME']

global_data_path =  os.environ['GLOBAL_DATA_PATH']
local_data_path  =  os.environ['LOCAL_DATA_PATH']
meteosat_archive_path =  os.environ['METEOSAT_ARCHIVE_PATH']


cosmo_path = '%s/cosmo/de' % global_data_path
proj_path = '%s/proj' % local_home_path
pics_path = '%s/pics' % local_home_path

