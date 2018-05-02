#!/usr/bin/python

import sys

local_libs = ['numpy', 'scipy']

for lib in local_libs:
    try:
        exec('print %s.__version__' % lib) 
        exec('del %s' % lib)
    except:
        pass

    for k in sys.modules.keys():
        if lib in k:
            del sys.modules[k]

local_path = '/u1/home/fabian/.local/lib/python2.6/site-packages'

if local_path in sys.path:
    sys.path.remove(local_path)

sys.path.insert(1,local_path)
