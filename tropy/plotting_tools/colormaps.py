#!/usr/bin/env python


import numpy as np
import matplotlib.colors as mcol
import pylab as pl
import copy

######################################################################
######################################################################

def colorname_based_cmap(colname, 
                         start_col = None, 
                         final_col = None,
                         reverse = False):
    '''
    A matplotlib colormap is constructed based on one color'

    USAGE
    =====
    cmap = colorname_based_cmap(colname, 
                         start_col = None, 
                         final_col = None,
                         reverse = False)

    
    INPUT
    =====
    colname : supported colorname as basis for colormap 
    start_col : color at one end of the colormap
    final_col : color at the other end of the colormap
    reverse: option if order is reversed

    OUTPUT
    ======
    cmap: resulting colormap

    '''

    cname = colname.lower()
    
    col_dict = pl.matplotlib.colors.cnames

    # check if color is known
    if cname not in col_dict:
        print('the color name is not available')

        raise NameError

    chex1 = pl.matplotlib.colors.cnames[cname]
    rgb1 = pl.matplotlib.colors.hex2color(chex1)
    h,s,v = pl.matplotlib.colors.rgb_to_hsv(np.array(rgb1).reshape(1,1,3))[0,0]

    # determine starting color
    if start_col == None:
        s0 = s - 0.3
        if s0 <= 0.1:
            s0 = 0.1

        v0 = v - 0.3
        if v0 <= 0.1:
            v0 = 0.1

        rgb0 = pl.matplotlib.colors.hsv_to_rgb(np.array([h,s0,v0]).reshape(1,1,3)).squeeze()

    else:
        chex0 = pl.matplotlib.colors.cnames[start_col.lower()]
        rgb0 = pl.matplotlib.colors.hex2color(chex0)
        

    # determine final color
    if final_col == None:
        s2 = s - 0.3
        if s2 <= 0.1:
            s2 = 0.1

        v2 = v + 0.3
        if v2 >= 0.9:
            v2 = 0.9
        rgb2 = pl.matplotlib.colors.hsv_to_rgb(np.array([h,s2,v2]).reshape(1,1,3)).squeeze()
    else:
        chex2 = pl.matplotlib.colors.cnames[final_col.lower()]
        rgb2 = pl.matplotlib.colors.hex2color(chex2)


    if reverse:
        rgb = copy.deepcopy(rgb0)

        rgb0 = copy.deepcopy(rgb2)
        rgb2 = copy.deepcopy(rgb)


    # set resulting cmap
    cdict = {}
    for i, col in enumerate(['red', 'green', 'blue']): 
        cdict[col] = [(0., 0., rgb0[i]), (0.5, rgb1[i], rgb1[i]), (1., rgb2[i], 1.)]

        cmap = pl.matplotlib.colors.LinearSegmentedColormap(cname + 'based',cdict,256, 
                                                            gamma = 0.5)
    return cmap

######################################################################
######################################################################



def nice_cmaps(cmap_name):

    if cmap_name == 'red2blue_disc':
       cmap_new =  pl.matplotlib.colors.LinearSegmentedColormap.from_list('new', [ '#d0220a', '#dc2506', '#e16618', '#edb526', '#fdf738', '#96c41b', '#3c9d10', '#007a30', '#57c3f7', '#48afe8', '#1c89d3', '#0e5dae', '#0d4780'], gamma = 1.2, N = 12)

    if cmap_name == 'red2blue':
       cmap_new =  pl.matplotlib.colors.LinearSegmentedColormap.from_list('new', [ '#d0220a', '#dc2506', '#e16618', '#edb526', '#fdf738', '#96c41b', '#3c9d10', '#007a30', '#57c3f7', '#48afe8', '#1c89d3', '#0e5dae', '#0d4780'], gamma = 1.2)

    elif cmap_name =='white-green-orange':
        cmap_new =  pl.matplotlib.colors.LinearSegmentedColormap.from_list('new', ['white', 'yellowgreen', 'forestgreen',  'gold','orangered', 'brown', 'black'], gamma = 1.4)

    elif cmap_name == 'white-blue-green-orange':
        cmap_new =  pl.matplotlib.colors.LinearSegmentedColormap.from_list('new', 
                                        ['white', '#1ca8e9','#097cc4', '#095389', 
                                         'forestgreen','yellowgreen',  'gold','orangered', 
                                         'brown', 'black'], gamma = 1.)
    elif cmap_name =='white-purple-orange':
         cmap_new =  pl.matplotlib.colors.LinearSegmentedColormap.from_list('new', ['white', 'powderblue', 'blue', 'purple', 'orangered', 'orange', 'gold', 'brown', 'black'], gamma = 1.4)


    elif cmap_name =='ocean-clouds':
         cmap_new =  pl.matplotlib.colors.LinearSegmentedColormap.from_list('ocean', ['teal', 'turquoise',  'paleturquoise', 'white', 'lightgray', '#15212c', '#0d141b' ], gamma = 1.3)


    elif cmap_name =='land-clouds':
         cmap_new =  pl.matplotlib.colors.LinearSegmentedColormap.from_list('ocean', ['teal', 'turquoise', 'paleturquoise', 'white', 'lightgray', '#81995f', '#5b6440'], gamma = 1.3)



    return cmap_new


######################################################################
######################################################################


def enhanced_colormap(vmin = 200., vmed = 240., vmax = 300.):

    
    nfull = 256

    ngray = int( nfull * (vmax - vmed) / (vmax - vmin) )
    ncol = nfull - ngray

    colors1 = pl.cm.gray_r(np.linspace(0., 1., ngray))
    colors2 = pl.cm.jet_r(np.linspace(0, 1., ncol))

    # combine them and build a new colormap
    colors = np.vstack((colors2, colors1))
    mymap = mcol.LinearSegmentedColormap.from_list('enhanced_colormap', colors)


    return mymap


######################################################################
######################################################################


def enhanced_wv62_cmap(vmin = 200., vmed1 = 230., vmed2 = 240., vmax = 260.):

    
    nfull = 256

    ncopp = int( nfull * (vmax - vmed2) / (vmax - vmin) )
    ngray = int( nfull * (vmed2 - vmed1) / (vmax - vmin) )
    ncol = nfull - (ncopp + ngray)

#    colors1 = pl.cm.copper(np.linspace(0., 1., ncopp))
    colors1 = pl.cm.afmhot(np.linspace(0, 1., ncopp))
    colors2 = pl.cm.gray_r(np.linspace(0., 1., ngray))
    colors3 = pl.cm.jet_r(np.linspace(0, 1., ncol))

    # combine them and build a new colormap
    colors = np.vstack((colors3, colors2, colors1))
    mymap = mcol.LinearSegmentedColormap.from_list('enhanced_cmap', colors)


    return mymap


######################################################################
######################################################################
    
def dwd_sfmap():

    rrlevs = [0, 0.1, 1, 2, 5, 10, 15, 20, 30, 50, 80, 100, 150, 200, 500]
    
    colvals = np.array([[ 0.98823529,  1.        ,  0.75686275],
                        [ 0.98431373,  1.        ,  0.36078431],
                        [ 0.8745098 ,  0.98823529,  0.14901961],
                        [ 0.62745098,  0.83921569,  0.14901961],
                        [ 0.27058824,  0.76470588,  0.4745098 ],
                        [ 0.        ,  0.83921569,  0.84705882],
                        [ 0.06666667,  0.63137255,  0.83921569],
                        [ 0.02745098,  0.00784314,  0.98823529],
                        [ 0.57254902,  0.19607843,  0.71764706],
                        [ 0.85490196,  0.15686275,  0.77647059],
                        [ 0.90588235,  0.05098039,  0.04705882],
                        [ 0.53333333,  0.05490196,  0.05098039],
                        [ 0.30980392,  0.05490196,  0.05098039]])
    
    cmap = pl.matplotlib.colors.ListedColormap(colvals)

    return rrlevs, cmap

######################################################################
######################################################################

if __name__ == '__main__':
    
    # cm_test = colorname_based_cmap('DarkGreen', start_col ='maroon', final_col = 'yellow')

    # x = np.linspace(0,1,100)
    # y = x

    # xx,yy = np.meshgrid(x,y)

    # f = xx**2 + yy**2

    
    # pl.imshow(f, cmap = cm_test)
    # pl.show()


    r = 240 + 100 * np.random.randn(100,100)
    import scipy.ndimage
    r = scipy.ndimage.gaussian_filter(r, 3)

    pl.figure(figsize=(16,6))
    pl.subplot(121)
    pl.imshow(r, vmin = 240, vmax = 300, cmap = pl.cm.gray_r)
    pl.imshow(np.ma.masked_greater(r, 240), vmax = 240, vmin = 210, cmap = pl.cm.jet_r)
    pl.colorbar()

    pl.subplot(122)
    cmap =  enhanced_colormap(vmin = 210., vmed = 240., vmax = 300.)
    pl.imshow(r, vmin = 210, vmax = 300., cmap = cmap)
    pl.colorbar()

