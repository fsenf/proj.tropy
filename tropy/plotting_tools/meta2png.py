#!/usr/bin/env python

import sys, os
from PIL import Image, PngImagePlugin
import pylab as pl

######################################################################
######################################################################


def meta2png(fname, meta):

    # add meta data to image
    im = Image.open(fname)
    im.info.update(meta)

    # these can be automatically added to Image.info dict                                                                              
    # they are not user-added metadata
    reserved = ('interlace', 'gamma', 'dpi', 'transparency', 'aspect')

    # undocumented class
    meta = PngImagePlugin.PngInfo()

    # copy metadata into new object
    for k,v in im.info.items():
        if k in reserved: continue
        meta.add_text(k, v, 0)

    # and save
    im.save(fname, "PNG", pnginfo=meta)

######################################################################
######################################################################



class pngsave(object):
    '''
    That class is designed to save pylab figures in png and add meta data.
    '''
 
    def __init__(self, *args, **kwargs):

        self.__call__(*args, **kwargs)

        return

    def __call__(self, fname, **kwargs):

        fig = kwargs.pop('fig', pl.gcf())

        print('... save image to ', fname)

        fig.savefig(fname, **kwargs)

        meta2png(fname, self.meta())
        return


    def meta(self):
        
        meta = {'Author': 'Fabian Senf'}
        meta['Source'] =  '{cwd}/{sname}'.format(cwd = os.getcwd(), sname =  sys.argv[0])

        return meta

######################################################################
######################################################################


if __name__ == '__main__':

    meta = {'Author': 'Fabian Senf'}
    meta['Source'] =  '{cwd}/{sname}'.format(cwd = os.getcwd(), sname =  sys.argv[0])

    fname = 'msg_today.png'

    meta2png(fname, meta)
