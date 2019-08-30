#!/usr/bin/env python

# load libraries -----------------------------------------------------
import sys, os, glob

from .msevi import MSevi
import numpy as np
import datetime as dt
import scipy.ndimage

import pylab as pl
import matplotlib.colors as col
import tropy.plotting_tools.colormaps as colormaps

from PIL import Image, ImageDraw, ImageFont

# ====================================================================
   


######################################################################
######################################################################


class MSeviRGB(MSevi):
    ''' 
     Container for MSG SEVIRI RGB images. 
    '''


    def __init__(self, rgb_type = 'pytroll_nc', msevi = None, tstamp = True, **kwargs):


        # init variables ---------------------------------------------

        self.rgb_type = rgb_type
 
        # images .....................................................
        self.images = {}
        # ============================================================


        # initialized MSG Data Container object part -----------------
        if msevi:
            # write content of MSevi instance on MSeviRBG 
            for att in msevi.__dict__:
               self.__setattr__(att, msevi.__dict__[att])

        else:
            # if no MSevi instance is imported then initialize MSeviRGB
            super(MSeviRGB, self).__init__(**kwargs)
        #=============================================================

        # make rgb images --------------------------------------------
        self.create_rgb(self.rgb_type, tstamp = tstamp, **kwargs)
        #=============================================================

        return

######################################################################

    def __repr__(self):
       return super(MSeviRGB, self).__repr__() + \
             'RGB-Images:      %s'  % str(list(self.images.keys())) 
 
#####################################################################

    def create_rgb(self, rgb_type, **kwargs):

        
        tstamp = kwargs.get('tstamp', True)


        # check if several composites should be made -----------------
        if type(rgb_type) == type([]):

            rgb_list = rgb_type

            for rgb in rgb_list:
                self.create_rgb(rgb, **kwargs)

            return
        else:
            rgb_job = rgb_type.lower()
        #=============================================================


        
        # natural color composite ------------------------------------
        if rgb_job in ('nc'):  #, 'all'):  
             
            rgb = 'nc'
    
            # load required channels ................................. 
            required_channels =  ['IR_016', 'VIS008', 'VIS006']

            self.load(required_channels)

            # calculate reflectance ..................................
            self.rad2refl(required_channels)

           
            # set red, green, blue ...................................    
            r = self.ref['IR_016']
            g = self.ref['VIS008']
            b = self.ref['VIS006']
 
            imgr = Image.fromarray(to255(r))
            imgg = Image.fromarray(to255(g))
            imgb = Image.fromarray(to255(b))
    
            self.images[rgb] =  Image.merge('RGB',(imgr,imgg,imgb))
        # ============================================================



        # pytrolls natural color composite ---------------------------
        if rgb_job in ('pytroll_nc', 'all'):  
             
            rgb = 'pytroll_nc'
    
            # load required channels ................................. 
            required_channels =  ['IR_016', 'VIS008', 'VIS006']

            self.load(required_channels)

            # if self.check_channels(required_channels):
            #    print 'Not all here!'

            # calculate reflectance ..................................
            self.rad2refl(required_channels)

           
            # set red, green, blue ...................................    

            r = self.ref['IR_016'] / 0.9
            g = self.ref['VIS008'] / 0.9
            b = self.ref['VIS006'] / 0.9
 
            imgr = Image.fromarray(to255(r, gamma = 1.8))
            imgg = Image.fromarray(to255(g, gamma = 1.8))
            imgb = Image.fromarray(to255(b, gamma = 1.8))
    
            self.images[rgb] =  Image.merge('RGB',(imgr,imgg,imgb))
        # ============================================================


        # high res natural color composite ---------------------------
        if rgb_job in ('nc_hrv', 'all'):  
             
            rgb = 'nc_hrv'
    
            # load required channels ................................. 
            required_channels =  ['IR_016', 'VIS008', 'VIS006', 'HRV']

            self.load(required_channels)

            # if self.check_channels(required_channels):
            #    print 'Not all here!'

            # calculate reflectance ..................................
            self.rad2refl(required_channels)

            
            # get hres channel reflectance ...........................
            hfac = 2.  # unphysical enhancement factor
            hgamma = 0.8
            w0 = 0.3
            w1 = 0.3

            hrv = self.ref['HRV']
            h = hrv ** (1. / hgamma)
            hrv_low = scipy.ndimage.convolve(hrv, np.ones((3, 3)) / 9.)
            hrv_high = hrv - hrv_low

            lchans = ['VIS006', 'VIS008', 'IR_016']
            c_high = {}
            for lc in lchans:
                c = self.ref[lc].repeat(3, axis = 0).repeat(3, axis = 1)
                c_low =  scipy.ndimage.convolve(c, np.ones((3,3)) / 9.)
                c_high[lc] = c_low + hfac * hrv_high


            # weights
            w = w0 + w1 * h 

           
            # set red, green, blue ...................................    
            r = c_high['IR_016'] / 0.9 * (1 - w)  + h * w
            g = c_high['VIS008'] / 0.9 * (1 - w)  + h * w
            b = c_high['VIS006'] / 0.9 * (1 - w)  + h * w
 
            imgr = Image.fromarray(to255(r, gamma = 1.8))
            imgg = Image.fromarray(to255(g, gamma = 1.8))
            imgb = Image.fromarray(to255(b, gamma = 1.8))
    
            self.images[rgb] =  Image.merge('RGB',(imgr,imgg,imgb))
        # ============================================================


        # composite between infrad and hres visible channels ---------
        if rgb_job in ('ir_hrv', 'all'):  
             
            rgb = 'ir_hrv'
    

            # load required images ...................................
            required_images = ['nc_hrv', 'col_108']
            
            for im in required_images:
                if im not in self.images:
                    self.create_rgb(im, **kwargs)

            img_h = self.images['nc_hrv']
            img_ir = self.images['col_108'].resize(img_h.size)

            btmin = 210.
            btmax = 250.
            m1 = (self.bt['IR_108'] - btmin) / (btmax - btmin)
            m1 = m1.repeat(3, axis = 0).repeat(3, axis = 1)

            m = m1

            mask = Image.fromarray(to255(m, gamma = 1.), mode='L')
    
            self.images[rgb] =  Image.composite(img_h, img_ir, mask)
        # ============================================================


        # composite between infrad and hres visible channels ---------
        if rgb_job in ('nc_merge', 'all'):  
             
            rgb = 'nc_merge'
    

            # load required images ...................................
            required_images = ['pytroll_nc', 'ir_natcol']
            
            for im in required_images:
                if im not in self.images:
                    self.create_rgb(im, **kwargs)

            img_vis = self.images['pytroll_nc']
            img_ir = self.images['ir_natcol']


            xm = 80
            dx = 5

            self.sunzen()
            m =  0.5*(1 + np.tanh(( self.szen -xm ) / dx))
            mask = Image.fromarray(to255(m, gamma = 1.), mode='L')
    
            self.images[rgb] =  Image.composite(img_ir, img_vis, mask)
        # ============================================================



        
        # air mass composite -----------------------------------------
        if rgb_job in ('am', 'all'):  
             
            rgb = 'am'
    
            # load required channels ................................. 
            required_channels =  ['WV_062', 'WV_073', 'IR_097', 'IR_108']

            self.load(required_channels)


            # calculate brightness temperature .......................
            self.rad2bt(required_channels)


            # minimum and maximum values .............................
            (rmin, rmax) = (-25, 0)     
            (gmin, gmax) = (-40, 5)     
            (bmin, bmax) = (-243, -208)     
             
 
            # set red, green, blue ...................................    
            r = (self.bt['WV_062'] - self.bt['WV_073'] - rmin) / (rmax - rmin)
            g = (self.bt['IR_097'] - self.bt['IR_108'] - gmin) / (gmax - gmin)
            b = (                 -  self.bt['WV_062'] - bmin) / (bmax - bmin)
 
            imgr = Image.fromarray(to255(r))
            imgg = Image.fromarray(to255(g))
            imgb = Image.fromarray(to255(b))
    
            self.images[rgb] =  Image.merge('RGB',(imgr,imgg,imgb))
        # ============================================================

        
        # convective storms composite --------------------------------
        if rgb_job in ('severe_storms', 'all'):  
             
            rgb = 'severe_storms'
    
            # load required channels ................................. 
            required_channels =  ['WV_062', 'WV_073', 'IR_039', 'IR_016', 'IR_108', 'VIS006']

            self.load(required_channels)


            # calculate brightness temperature .......................
            self.rad2bt( ['WV_062', 'WV_073', 'IR_039',  'IR_108'])
            self.rad2refl( [ 'IR_016',  'VIS006'])



            # minimum and maximum values .............................
            (rmin, rmax) = (-35, 5)     
            (gmin, gmax) = (-5, 60)     
            (bmin, bmax) = (-0.7, 0.25)     
             
 
            # set red, green, blue ...................................    
            r = (self.bt['WV_062'] - self.bt['WV_073'] - rmin) / (rmax - rmin)
            g = (self.bt['IR_039'] - self.bt['IR_108'] - gmin) / (gmax - gmin)
            b = (self.ref['IR_016']-  self.ref['VIS006'] - bmin) / (bmax - bmin)
 
            imgr = Image.fromarray(to255(r))
            imgg = Image.fromarray(to255(g, gamma=0.5))
            imgb = Image.fromarray(to255(b))
    
            self.images[rgb] =  Image.merge('RGB',(imgr,imgg,imgb))
        # ============================================================



        # ir108 color enhanced ---------------------------------------
        if rgb_job in ('ir_natcol', 'all'):  
             
            rgb = 'ir_natcol'
    
            # load required channels ................................. 
            required_channels =  ['IR_108',]

            self.load(required_channels)


            # calculate brightness temperature .......................
            self.rad2bt(required_channels)

            if 'lsm' not in self.__dict__:
                self.landsea()

            if 'lon' not in self.__dict__:
                self.lonlat()

            # prepare masking ........................................
            # get brightness temperature 
            tsurf = kwargs.get( 'ir_natcol_tsurf', 270)
            lowcloud_delta = kwargs.get( 'ir_natcol_lowcloud_delta', 15)
            lowcloud_factor  = kwargs.get( 'ir_natcol_lowcloud_factor', 0.8)

            bmax = tsurf + 35 * np.cos(np.deg2rad(self.lat))
            bmed = bmax - lowcloud_delta
            bmin = 210
            bt108 = self.bt['IR_108']


            # try two ranges

            v1 = (bt108 - bmin) / (bmed - bmin)
            v2 = (bt108 - bmed) / (bmax - bmed)

            fac = lowcloud_factor
            v = np.where(bt108 < bmed, v1 * fac, fac + v2 * (1 - fac))
            v = np.clip(v, 0., 1.)


            lsm = self.lsm
            m = ( lsm == 2 )

            bsea = np.ma.masked_where(m, v)
            bland = np.ma.masked_where(m == False, v)


            # create mask for 3 rgb + 1 alpha = 4 dimensions
            msk = np.dstack([m,]*4)
    

            # get color schemes ......................................
            cmap_sea = colormaps.nice_cmaps('ocean-clouds')
            cmap_land = colormaps.nice_cmaps('land-clouds')


            
            # do plotting ............................................
            # make the bw plot with the bone colormap
            pcm_sea = pl.pcolormesh(bsea,
                                    vmin = 0, 
                                    vmax = 1., 
                                    cmap = cmap_sea)

            pcm_land = pl.pcolormesh(bland, 
                                    vmin = 0, 
                                    vmax = 1., 
                                    cmap = cmap_land)



            # convert plots to rgba values ...........................
            rgb_sea = pcm_sea.to_rgba(bsea, bytes = True)
            rgb_land = pcm_land.to_rgba(bland, bytes = True)
            pl.close()

            self.images[rgb] =  Image.fromarray(np.where(msk, rgb_land, rgb_sea),mode='RGBA')


        # ir108 color enhanced ---------------------------------------
        if rgb_job in ('col_108', 'all'):  
             
            rgb = 'col_108'
    
            # load required channels ................................. 
            required_channels =  ['IR_108']

            self.load(required_channels)


            # calculate brightness temperature .......................
            self.rad2bt(required_channels)

            
            # prepare masking ........................................
            # get brightness temperature 
            bt108 = self.bt['IR_108']

            # create masked array from the colored part of the image
            tmin = 210
            tmax = 240
            btm = np.ma.masked_greater(bt108, tmax)
            btm_ext = np.ma.masked_greater(bt108, 210)


            # create mask for 3 rgb + 1 alpha = 4 dimensions
            msk = np.dstack([btm.mask, btm.mask, btm.mask, btm.mask])
            msk_ext = np.dstack([btm_ext.mask, btm_ext.mask, btm_ext.mask, btm_ext.mask])


            
            # do plotting ............................................
            # make the bw plot with the bone colormap
            pcm_bw = pl.pcolormesh(bt108, vmin = tmax, vmax = 300, cmap = pl.cm.gray_r)

            # make the color plot with jet
            pcm_col = pl.pcolormesh(btm, vmin = tmin, vmax = tmax, cmap = pl.cm.jet_r)

            # extended color scale
            cmap_new = col.LinearSegmentedColormap.from_list('new', 
                                                             ['darkred', 'black',  'pink', 'black'][::-1])    
            
            pcm_ext = pl.pcolormesh(btm_ext,
                                    cmap = cmap_new, vmax = 210, vmin = 180)


            # convert plots to rgba values ...........................
            rgb_bw = pcm_bw.to_rgba(bt108, bytes = True)
            rgb_col = pcm_col.to_rgba(btm, bytes = True)
            rgb_ext = pcm_ext.to_rgba(btm_ext, bytes = True)

            img1 =  Image.fromarray(np.where(msk, rgb_bw, rgb_col),mode='RGBA')
            

            self.images[rgb] = Image.fromarray(np.where(msk_ext, img1, rgb_ext), mode='RGBA')


        # dust composite ---------------------------------------------
        if rgb_job in ('dust', 'all'):  
             
            rgb = 'dust'
    
            # load required channels ................................. 
            required_channels =  ['IR_087', 'IR_108', 'IR_120']

            self.load(required_channels)


            # calculate brightness temperature .......................
            self.rad2bt(required_channels)


            # minimum and maximum values .............................
            (rmin, rmax) = (-4, 2)     
            (gmin, gmax) = (0, 15)     
            (bmin, bmax) = (261, 289)     
             
 
            # set red, green, blue ...................................    
            r = (self.bt['IR_120'] - self.bt['IR_108'] - rmin) / (rmax - rmin)
            g = (self.bt['IR_108'] - self.bt['IR_087'] - gmin) / (gmax - gmin)
            b = (self.bt['IR_108'] - bmin) / (bmax - bmin)
 
            imgr = Image.fromarray(to255(r))
            imgg = Image.fromarray(to255(g, gamma = 2.5))
            imgb = Image.fromarray(to255(b))
    
            self.images[rgb] =  Image.merge('RGB',(imgr,imgg,imgb))



        # high res clouds composite ----------------------------------
        if rgb_job in ('hrv_clouds', 'all'):  
             
            rgb = 'hrv_clouds'
    
            # load required channels ................................. 
            required_channels =  ['IR_108', 'HRV']

            self.load(required_channels)


            # calculate brightness temperature .......................
            self.rad2bt('IR_108')
            self.rad2refl('HRV')


            # minimum and maximum values .............................
            # (bmin, bmax) = (self.bt['IR_108'].min(), self.bt['IR_108'].max())     
            (bmin, bmax) = (203, 323)
 
            # set red, green, blue ...................................    
            rg = (self.ref['HRV'])
            b = (bmax - self.bt['IR_108']) / (bmax - bmin)
 
            imgrg = Image.fromarray(to255(rg))
            imgb = Image.fromarray(to255(b)).resize(imgrg.size)
    
            self.images[rgb] =  Image.merge('RGB',(imgrg,imgrg,imgb))
            

        # high res fog composite -------------------------------------
        if rgb_job in ('hrv_fog', 'all'):  
             
            rgb = 'hrv_fog'
    
            # load required channels ................................. 
            required_channels =  ['IR_016', 'HRV']

            self.load(required_channels)


            # calculate brightness temperature .......................
            self.rad2refl(required_channels)


 
            # set red, green, blue ...................................    
            r = (self.ref['IR_016']) / 0.7
            gb = (self.ref['HRV'])
 
            imggb = Image.fromarray(to255(gb, gamma = 1.5))
            imgr = Image.fromarray(to255(r)).resize(imggb.size)
    
            self.images[rgb] =  Image.merge('RGB',(imgr,imggb,imggb))


        # high res severe storms composite ---------------------------
        if rgb_job in ('hrv_severe_storms', 'all'):  
             
            rgb = 'hrv_severe_storms'
    
            # load required channels ................................. 
            required_channels =  ['IR_108', 'IR_039', 'HRV']

            self.load(required_channels)


            # calculate brightness temperature .......................
            self.rad2refl('HRV')
            self.rad2bt(['IR_108', 'IR_039'])


 
            # minimum and maximum values .............................
            (bmin, bmax) = (-60, -40)
 
            # set red, green, blue ...................................    
            rg = (self.ref['HRV'])
            b = (self.bt['IR_108'] - self.bt['IR_039'] - bmin) / (bmax - bmin)
 
            imgrg = Image.fromarray(to255(rg))
            imgb = Image.fromarray(to255(b,gamma=2)).resize(imgrg.size)
    
            self.images[rgb] =  Image.merge('RGB',(imgrg,imgrg,imgb))
            


        # high res severe storms composite ---------------------------
        if rgb_job in ('hrv_vis_ros', 'all'):  
             
            rgb = 'hrv_vis_ros'
    
            # load required channels ................................. 
            required_channels =  ['IR_108', 'VIS008', 'HRV']

            self.load(required_channels)


            # calculate brightness temperature .......................
            self.rad2refl(['HRV', 'VIS008'])
            self.rad2bt('IR_108')


 
            # minimum and maximum values .............................
            (bmin, bmax) = (220, 300)
 
            # set red, green, blue ...................................    
            r = (self.ref['HRV'])
            g = (self.ref['VIS008'])
            b = (bmax - self.bt['IR_108']) / (bmax - bmin)
 
            imgr = Image.fromarray(to255(r))
            imgg = Image.fromarray(to255(g)).resize(imgr.size)
            imgb = Image.fromarray(to255(b)).resize(imgr.size)

            self.images[rgb] =  Image.merge('RGB',(imgr, imgg, imgb))
            

        # time stamps ------------------------------------------------
        if tstamp:
            if rgb_job == 'all':
                for rgb in list(self.images.keys()):
                    self.time_stamp(rgb)
    
            else:
                self.time_stamp(rgb_job)        
        # ============================================================
        
        return

######################################################################

    def show(self, rgb_str,**kwargs):
        
        if 'command' not in kwargs:
            kwargs['command'] = 'display'

        overlay = kwargs.pop('overlay', False)


        if not type(rgb_str) == type(''):
            print('Use only one rgb-type string as argument')
            return

        if rgb_str not in self.images:
            self.create_rgb(rgb_str, **kwargs)

        if overlay:
            self.overlay()

        if rgb_str == 'all':
            for rgb in list(self.images.keys()):
                self.images[rgb].show(**kwargs)
        
        else:
            rgb = rgb_str
            self.images[rgb].show(**kwargs)
            
        return

    
######################################################################
    
    def overlay(self):

        if self.reg_name != 'eu':
            return
        else:
            # hard coded ATTENTION! and only for EU
            ofile = '/vols/talos/home/fabian/pics/borders_hrv.png'

            o = Image.open(ofile)

            for rgb_str in list(self.images.keys()):
                i = self.images[rgb_str].resize((2400,1800)).convert(mode='RGBA')
                self.images[rgb_str] = Image.composite(o,i,o)

            return

######################################################################



    def oshow(self, rgb_str, **kwargs):
        self.show(rgb_str, overlay = True)
    
######################################################################

    def save(self, rgb_str,**kwargs):


        # check if images exist --------------------------------------
        if rgb_str not in self.images:
            self.create_rgb(rgb_str, **kwargs)
        # ============================================================

 

        # do the save for all rgb images -----------------------------
        if rgb_str == 'all':

            for rgb in list(self.images.keys()):

                self.save(rgb, **kwargs)

            return
        # ============================================================
        

        # generate default filename ----------------------------------
        time_str = self.time.strftime('%Y-%m-%d_%H%M')
        default_name ='%s_%s_%s.png' % (self.sat_type.lower(), rgb_str, time_str)
        # ============================================================


        # get picture filename ---------------------------------------
        if 'filename' not in kwargs:
            fname = default_name 
        else:
            fname = kwargs['filename']
            del kwargs['filename']
        # ============================================================

        # get picture directory --------------------------------------
        if 'pic_dir' not in kwargs:
            pic_dir ='.'
        else:
            pic_dir = kwargs['pic_dir']    

        if not os.access(pic_dir, os.F_OK):
            os.makedirs(pic_dir)
        # ============================================================

        
        fname = pic_dir + '/' + fname

        print('image saved in ',fname)
        self.images[rgb_str].save(fname, **kwargs)
            
        return

        

######################################################################

    def time_stamp(self, rgb_str, offset = 38, fontsize = 16):
        
        img = self.images[rgb_str]
        nx, ny = img.size

        # create an ImageDraw object from image
        dimg = ImageDraw.Draw(img)

    # load font
        fontfile='/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf'
        sans = ImageFont.truetype(fontfile, fontsize)


    # include time stamp 
        offset = offset
        dimg.text((offset, ny-offset ),self.time.strftime('%Y-%m-%d  %H:%M UTC'), \
                      font=sans, fill='black')
        
        white_offset = offset + 2
        dimg.text((offset, ny - white_offset),self.time.strftime('%Y-%m-%d  %H:%M UTC'), \
                      font=sans, fill='white')
    
        return


######################################################################
######################################################################


def to255(d, gamma = 1.0):

    dat = d.copy()
    # remove outliners for dat>1. or dat<0.
    dat[dat>1.] = 1.
    dat[dat<0.] = 0.

    if not gamma == 1.0:
       dat = dat ** (1./gamma) 
   
    # transform (0,1)-data to (0,255)-image-scale 
    rgb = 255 * dat
   
    return rgb.astype('uint8')

######################################################################
######################################################################



if __name__ == '__main__':

    cin = { \
     'time' : 'now',\
#        'time': dt.datetime(2015,9,1,8,0),\
     'scan_type' : 'rss',\
        'rgb_type':'pytroll_nc',\
        'region':'eu'\
#        'region': 'germ'\
         }

    scene = MSeviRGB(**cin)

    scene.oshow('ir_hrv')
#    scene.oshow('nc_hrv')
#    scene.oshow('col_108')

#    scene.oshow('pytroll_nc')
