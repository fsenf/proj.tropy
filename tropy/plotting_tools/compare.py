import sys
from shaded import shaded
import numpy as np
import pylab as py


######################################################################
######################################################################
def compare_shaded(x, y, *args, **kwargs):
    '''
    The function compare can be used to fastly compare several 2d fields.

    The coordinates of the x and y axis are needed. An arbitrary number
    of fields is plotted using the function shaded. Corresponding option
    keywords can be used.

    '''

    # get the number of fields
    n = len(args)

    # assume a square to optimally place figures
    nrow = int(np.sqrt(n))
    ncol = int(np.ceil(float(n) / nrow))

    # aspect ratio of the figure
    aspect = float(nrow) / float(ncol)


    if aspect < 0.5: 
        aspect = 0.5

    if aspect < 2./3.:
        nx = 9.
        ny = aspect * nx
    else:
        ny = 6.
        nx = ny / aspect

    fontsize = 13  -  int(ncol*nrow)
    if fontsize < 8 :
        fontsize = 8

    py.figure(figsize=(nx,ny))

    for i, v in enumerate(args):
        py.subplot(nrow, ncol, i + 1)
        shaded(x,y,v,**kwargs)
        py.xticks(fontsize=fontsize)
        py.yticks(fontsize=fontsize)
        py.title(str(i+1))

        cb = py.colorbar()
        py.axes(cb.ax)
        py.yticks(fontsize=fontsize)
        
    py.show()

    return


######################################################################
######################################################################

######################################################################
######################################################################

def compare(*args, **kwargs):
    '''
    The function compare can be used to fastly compare several 2d fields.

    The coordinates of the x and y axis are needed. An arbitrary number
    of fields is plotted using the function shaded. Corresponding option
    keywords can be used.

    '''

    if kwargs.has_key('georef'):
        georef = True
        del kwargs['georef']
    else:
        georef = False

    # get the number of fields
    n = len(args)

    # assume a square to optimally place figures
    nrow = int(np.sqrt(n))
    ncol = int(np.ceil(float(n) / nrow))

    # aspect ratio of the figure
    aspect = float(nrow) / float(ncol)


    if aspect < 0.5: 
        aspect = 0.5

    if aspect < 2./3.:
        nx = 9.
        ny = aspect * nx
    else:
        ny = 6.
        nx = ny / aspect

    fontsize = 13  -  int(ncol*nrow)
    if fontsize < 8 :
        fontsize = 8

    py.figure(figsize=(nx,ny))

    if georef:
        x = args[0]
        y = args[1]
        args.pop(0)
        args.pop(0)


    for i, v in enumerate(args):
        py.subplot(nrow, ncol, i + 1)

        if georef:
            py.pcolormesh(x,y,v,**kwargs)
        else:
            py.pcolormesh(v,**kwargs)

        py.xticks(fontsize=fontsize)
        py.yticks(fontsize=fontsize)
        py.title(str(i+1))

        cb = py.colorbar()
        py.axes(cb.ax)
        py.yticks(fontsize=fontsize)
        
    py.show()

    return
