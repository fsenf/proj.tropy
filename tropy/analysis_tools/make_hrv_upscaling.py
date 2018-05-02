#/usr/bin/env python

import numpy as np
import grid_and_interpolation as gi

######################################################################
######################################################################
# this routine is obsolete and only left for backwards compatibility #
######################################################################
######################################################################


def make_hrv_upscaling(v):
    
    '''
    A 2d field is upscaled to 3-fold higher resolution using linear 
    interpolation. Edge values are determined by linear extrapolation.
   

    USAGE:
    =====
    vhigh = make_hrv_upscaling(v)

    
    INPUT
    =====
    v: 2d variable field (Nrow, Ncols)

    
    OUTPUT
    ======
    vhigh: 2d variable field with 3-times higher resolution (3*Nrows, 3*Ncols)
    
    '''

    return gi.make_hrv_upscaling(v)


######################################################################
######################################################################


def make_add_edge(var):

    '''
    Adds egde region to 2d data field. Values in added left and right column, 
    and lower and upper row are linearily extrapolated from their neighbors.
   


    USAGE:
    =====
    var_new = make_add_edge(var)

    
    INPUT
    =====
    var: 2d variable field (Nrow, Ncols)

    
    OUTPUT
    ======
    var_new: 2d variable field with added edge values (Nrows + 2, Ncols + 2)
    
    '''
    

    return gi.make_add_edge(var)


######################################################################
######################################################################


if __name__ == "__main__":


    import sys
    import MSGtools



    lon,lat = MSGtools.get_msg_lon_lat('eu')
    
    lat_new = make_hrv_upscaling(lat)
    lon_new = make_hrv_upscaling(lon)
