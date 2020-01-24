#!/usr/bin/env python

'''
Module designed to write an Author Name and the Source Filename 
into the Meta-Data of an PNG file saved with matplotlib.
'''


import sys, os
from PIL import Image, PngImagePlugin
import pylab as pl
from IPython.core.display import Javascript
from IPython.display import display

######################################################################
######################################################################

def get_notebook_name():
    """
    Returns the name of the current notebook as a string
    
    From From https://mail.scipy.org/pipermail/ipython-dev/2014-June/014096.html
    """

    display(Javascript('IPython.notebook.kernel.execute("theNotebook = " + \
    "\'"+IPython.notebook.notebook_name+"\'");'))
    #
    nb_full_path = os.path.join(os.getcwd(), theNotebook)
    
    return os.path.join(os.getcwd(), theNotebook)

######################################################################
######################################################################


def meta2png(fname, meta):

    '''
    Saves meta data into image file.


    Parameters
    ----------
    fname : str
        name of png file

    meta : dict
        collection of extra meta data to be stored in image file


    
    '''

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


    Parameters
    ----------
    *args : list
        other positional arguments passed to `plt.savefig`

    **kwargs : dict
        other optional arguments passed to `plt.savefig``

        'author' : Place yeur name into author keyword 
        
        'source' : Specify your source filename if needed,
                   if not set, it tries to automatically find the filename
        
    '''
 
    def __init__(self, *args, **kwargs):

        # getting the author attribute
        self.author = kwargs.pop('author', 'Fabian Senf')
        self.source = kwargs.pop('source', 'Auto')


        self.__call__(*args, **kwargs)

        return


    def __call__(self, fname, **kwargs):

        fig = kwargs.pop('fig', pl.gcf())

        print('... save image to ', fname)

        fig.savefig(fname, **kwargs)

        meta2png(fname, self.meta())
        return


    def meta(self):

        '''
        Set the meta data, esp. the Author name and Source file information.
        '''

        meta = {'Author': self.author }

        if self.source == 'Auto':
            meta['Source'] =  '{cwd}/{sname}'.format(cwd = os.getcwd(), sname =  sys.argv[0])
        else:
            meta['Source'] = self.source

        return meta

######################################################################
######################################################################


if __name__ == '__main__':

    meta = {'Author': 'Fabian Senf'}
    meta['Source'] =  '{cwd}/{sname}'.format(cwd = os.getcwd(), sname =  sys.argv[0])

    fname = 'msg_today.png'

    meta2png(fname, meta)
