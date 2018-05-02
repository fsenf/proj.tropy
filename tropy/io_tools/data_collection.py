#!/usr/bin/env python

import numpy as np
import h5py
import file_status
import time
from data_settings import settings

######################################################################
######################################################################


class Dataset(object):
    '''
    Description:
    ============
    The Dataset class is a building block of the DataCollection class
    which will be used to easy input and output dataset sets.
    '''
    
    # ----------------------------------------------------------------

    def __init__(self, name, data = None, setting = settings('bt')):

        self.name = name

        if type(setting) == type({}):

            for k in setting.keys():
                self.__dict__[k] = setting[k]

        if data == None:
            self.data = np.array([])
        else:
            self.add(data)

    # ----------------------------------------------------------------
   

    # def description(self):
    #     print
    #     print ' Data Container '
    #     print ' ============== '

    #     c = self.__dict__
    #     for a in sorted(self.attrs()):
    #         print '%10s :' % a, c[a] 
            
    #     return ''
    # # ----------------------------------------------------------------

    # def __repr__(self):
    #     return self.description()
    # ----------------------------------------------------------------


    def add(self, array):
        
        a = array
        a0 = self.add_offset
        s  = self.scale_factor
        nan = self._FillValue
        
        counts = np.where (a != nan,  ( a - a0 ) / s, nan)

        self.data = np.array(counts, dtype = self.dtype)
        return
    # ----------------------------------------------------------------
    
    def array(self):
        
        counts = self.data
        a0 = self.add_offset
        s  = self.scale_factor
        nan = self._FillValue
        
        array = np.where (counts != nan,  s* counts + a0, nan)

        return array

    # ----------------------------------------------------------------

    def __getitem__(self, attr):
        
        if attr == self.name:
            return self.array()

        else:
            print 'Dataset has name: ', self.name
            return

    # ----------------------------------------------------------------


    def save(self, fname, subpath = None, mode = 'overwrite', time_out=300.):
        

        # check if file is open
        t_start = time.time()
        t = 0
        # while file_status.test(fname) and t < time_out:
        #     t = time.time() - t_start
        #     time.sleep( 0.1 )

        # if file_status.test(fname):
        #     print 'sorry, it is timeout, no data saved ...'
        #     return


        hfile = h5py.File(fname, 'a')
        
        # start with root group
        g = hfile


        # create subdirectory structure if needed
        if subpath != None:
            g = g.require_group(subpath)


        # create dataset
        if mode == 'overwrite':
            try:
                del g[self.name]
 #               print 'overwrite dataset'
            except:
                pass
#                print 'dataset does not exist'
        else:
            try:
                g[self.name]
                hfile.close()
                return
            except:
                pass

        d = g.create_dataset(self.name, data = self.data , compression="gzip")
        

        # create attributes
        for a in self.attrs():
            val = self.__dict__[a]

            d.attrs.create(a, val)
        

        hfile.close()

        return

    # ----------------------------------------------------------------

    def load(self, fname, subpath = None):
        
        hfile = h5py.File(fname, 'r')
        
        if subpath != None:
            vname = '%s/%s' % (subpath, self.name)
        else:
            vname = self.name

        try:
            d = hfile[vname]

            # get dataset
            self.data = np.array(d)
            
            # get attributes
            for a in d.attrs.keys():
                self.__dict__[a] = d.attrs[a]
        except:
            print 'dataset does not exist'
 

        hfile.close()

        return

    # ----------------------------------------------------------------


    def attrs(self):

        no_attrs = [ 'data']

        keys = self.__dict__.keys()

        for na in no_attrs:
            keys.remove(na)

        return keys


######################################################################
######################################################################
        

class DataCollection(dict):

    # ----------------------------------------------------------------

    def __init__(self):
        pass
    # ----------------------------------------------------------------

    def add(self, vname, array, setting = 'bt', subpath=None):
        
        dset = Dataset(vname, data = array, setting = setting)
        
        if subpath != None:
            if not self.has_key(subpath):
                self[subpath] = {}


            self[subpath][vname] = dset
        else:
            self[vname] = dset

        return

    # ----------------------------------------------------------------

    def save(self, fname):
        
        for vname, dset, subpath in self.list():
            dset.save(fname, subpath = subpath)
        
        return

    # ----------------------------------------------------------------

    def list(self):
        
        l = []
        for k in self.keys():
            
            if type(self[k]) == type({}):
                for vname in self[k].keys():
                    l.append([vname, self[k][vname], k])

            else:
                l.append([k, self[k], None])
                
        return l
######################################################################
######################################################################


if __name__ == '__main__':


    z = np.random.randn(100,100,100)
    z1 = np.random.randn(100,100,10)
    z2 = np.random.randn(100,100)

    d = DataCollection()

    d.add('z', z, subpath='adsk')
    d.add('z1', z1, subpath='adsk')
    d.add('z2', z2, subpath = 'sdfj/dsf/kg')
    d.add('z2', z2)


    d.list()
    d.save('test.h5')

