#!/usr/bin/env python

# ####################################################################
# ####################################################################
# Library for the interaction with hdf files #########################
# ####################################################################
# ####################################################################

import h5py
import numpy as np
import os, glob
import datetime
import copy

######################################################################
######################################################################

def list_hdf_groups(fname):

    '''
    Makes a list of hdf groups.

    Only the 1st level is implemented, TDB: recursive listing.
    


    USAGE
    =====
    glist = list_hdf_groups(fname)

    
    INPUT
    =====
    fname: filename of hdf file


    
    OUTPUT
    ======
    glist: list of group names

    '''

    group_type = h5py._hl.group.Group


    f = h5py.File(fname, 'r')

    glist = []
    for kname in sorted(f.keys()):
        if type(f[kname]) == group_type:
            glist += [ kname, ]
    f.close()

    return glist


######################################################################
######################################################################


def save_dict2hdf(fname, d, mode = 'w'):

    '''
    The content of nested dictionaries of arbritrary depth is saved 
    into an hdf file.

    The key, subkeys, etc. are mapped into the hdf-groups  directory
    structure.

    USAGE
    =====
    save_dict2hdf(fname, d)

    
    INPUT
    =====
    fname: filename of output hdf file
    d: dictionary which contains the data

    
    OUTPUT
    ======
    None

    '''

    # check if directory exists ......................................
    fdir = os.path.dirname(fname)
    if fdir and not os.access(fdir, os.F_OK):
        os.makedirs(fdir)

    # open file and write ............................................
    f = h5py.File(fname, mode = mode)
    
    save_dict_cont(f, d)
    
    f.close()

    return

######################################################################
######################################################################

def save_dict_cont(f, d):

   """ Recursively saves nested dictionaries."""

   for key, value in d.items():

       if isinstance(value, dict):
           
           g = f.create_group(str(key))

           save_dict_cont(g, value)

       else:

           # try to use compression for data saving
           try:
               f.create_dataset(str(key), data = value, compression="gzip")
           except:
               f.create_dataset(str(key), data = value)
           
   return


######################################################################
######################################################################

def read_var_from_hdf(fname, vname, subpath = None):

    '''
    A specific variable is read from an hdf file. Group and subgroups can be 
    specified.

    USAGE
    =====
    v = read_var_from_hdf(fname, vname, subpath = None)

    
    INPUT
    =====
    fname: filename of input hdf file
    vname: variable name as given in hdf file
    
    optional:
    --------
    subpath: group or subgroup path with /

    
    OUTPUT
    ======
    v: variable as numpy array


    '''

    
    print('.. open ', fname)
    f = h5py.File(fname, 'r')

    # handle groups and subgroups
    g = copy.copy(f)
    if subpath != None:
       grouplist = subpath.split('/')
       for gname in grouplist:
           if len(gname) != 0:
               g = g[gname]

    v = np.array( g[vname] )

    f.close()

    return v


######################################################################
######################################################################



def read_dict_from_hdf(fname):

    '''
    The content of an hdf file with arbitrary depth of subgroups is saved in
    nested dictionaries.

    USAGE
    =====
    d = read_dict_from_hdf(fname)

    
    INPUT
    =====
    fname: filename of output hdf file

    
    OUTPUT
    ======
    d: nested dictionary which contains the data


    '''

    f = h5py.File(fname, 'r')
    
    d = {}
    read_dict_cont(f, d)
    
    f.close()

    return d


######################################################################
######################################################################

def read_dict_cont(f, d):

   """ Recursively read nested dictionaries."""

   # set simple types of hdf content
   data_type = h5py._hl.dataset.Dataset
   group_type = h5py._hl.group.Group


   for key, value in f.items():

       if isinstance(value, group_type):
           
           d[key] = {} 

           read_dict_cont(value, d[key])

       else:
          try:
             d[key] = np.array(value)
          except:
             d[key] = value.value
           
   return


######################################################################
######################################################################

def dict_merge(a, b):
    '''
    Recursively merges dict's. not just simple a['key'] = b['key'], if
    both a and bhave a key who's value is a dict then dict_merge is called
    on both values and the result stored in the returned dictionary.

    from https://www.xormedia.com/recursively-merge-dictionaries-in-python/ 
    '''


    if not isinstance(b, dict):
        return b
    result = copy.deepcopy(a)
    for k, v in b.items():
        if k in result and isinstance(result[k], dict):
                result[k] = dict_merge(result[k], v)
        else:
            result[k] = copy.deepcopy(v)
    return result

######################################################################
######################################################################


def update_dict_in_hdf(fname, din):

    '''
    The content of nested dictionaries of arbritrary depth is updated 
    in an hdf file.

    The key, subkeys, etc. are mapped into the hdf-groups  directory
    structure.

    USAGE
    =====
    update_dict_in_hdf(fname, d)

    
    INPUT
    =====
    fname: filename of output hdf file
    d: dictionary which contains the data

    
    OUTPUT
    ======
    None

    '''


    # if hdf file does not exists than just write data to file
    if not os.access(fname, os.F_OK) or not h5py.is_hdf5(fname):
        save_dict2hdf(fname, din)

    else:

        tmpname = '/tmp/' + os.path.basename(fname)
        os.system('mv %s %s' % (fname, tmpname))

        # read original data
        dfile = read_dict_from_hdf(tmpname)

    # merge old and new data
        dout = dict_merge(dfile, din)

    # save the merge data set again
        save_dict2hdf(fname, dout, mode = 'w')

        os.system('rm %s' %tmpname)
    return

######################################################################
######################################################################

def get_seviri_chan(chan_list, day, 
                    scan_type = 'rss', 
                    fname = None,
                    calibrate = True, 
                    add_meta = False,
                    add_geo = False, 
                    arch_dir = None):

    '''
    Reads MSG - SEVIRI radiance of a given channel for a given date.


    USAGE:
    =====
    var = get_seviri_chan(chan_list, day, scan_type = 'rss', calibrate=True, add_meta = False)


    INPUT:
    ======
    chan_list: name of channel, e.g. ir_108 or IR_108 (case does not matter)
               or list of channel names

    day: datetime object which include day and time of MSG time slot

    scan_type: sets of scanning modus, i.e. 'rss' or 'pzs' (DEFAULT: 'rss')

    calibrate : optional, decides if output is radiance (True) or counts (False)

    add_meta: meta data of the MSG-SEVIRI hdf file


    OUTPUT:
    =======
    var: dictionary including radiance/counts of channel <ch_name> at date <day>
         and if chosen meta data

    '''


    # check input ----------------------------------------------------
    if type(chan_list) == type(''):
        chan_list= [chan_list]


    # create filename or take given one ------------------------------
    if fname != None:
        hfile = fname

    else:
        # check scan_type ................................................
        scan =  scan_type.lower()

        if not scan in ('rss', 'pzs'):
            print('ERROR: scan_type %s not available')
            print('       use either "rss" or "pzs"')
            return None



    # set time and date strings for file selection ...................
        time_string = day.strftime('%Y%m%d%H%M')
        date_string = day.strftime('%Y/%m/%d')



    # create hdf filename ............................................
    # set archive directory 
        if not arch_dir:
            base_dir = "/mnt/zephyr/u1/EUMETCAST_ARCHIVE/msg-sevi/"
    
            arch_dir = base_dir + "%s/l15_hdf/eu" % scan
    

    # set hdf file name
        hfile = "%s/%s/msg?-sevi-%s-l15-%s-eu.h5" % (arch_dir, date_string, 
                                                     time_string, scan)
    # ================================================================


    # get filename while MSG number is unspecified ...................
    hfile_list =  glob.glob(hfile)

    if hfile_list:
        hfile = hfile_list[0]
    else:
        print('ERROR: ', hfile,' does not exist!')
        raise IOError
    # ================================================================



    # I/O ------------------------------------------------------------
    # check if hdf file exists 
    if not h5py.is_hdf5(hfile):
        print('No hdf-file available!')
        raise IOError
    
    
    # open hdf file
    f = h5py.File(hfile,"r")
    

    # Read counts of sat channel .....................................
    counts = {}
    for ch in chan_list:
        fimg = f['l15_images']
        image_name = 'image_' + ch.lower()
        if image_name not in list(fimg.keys()):
            print(image_name,' not in ',list(fimg.keys()))
            raise KeyError 
        counts[ch] = fimg[image_name].value


    # Do calibration .................................................
    if calibrate:
        # extract meta info for all channels
        ch_info_all = f['meta/channel_info'].value

        info = {}
        for info_line in ch_info_all:
            ch_name = str( info_line['name'], 'utf-8' ).lower()
            info[ch_name] = info_line
            
        rad = {}
        for ch in list(counts.keys()):
            ch_name = ch.lower()
            offset = info[ch_name]['cal_offset']
            slope = info[ch_name]['cal_slope']
            c = counts[ch]
            rad[ch] = np.where(c != 0, slope * c + offset, 0)
    
        var = rad
    else:
        var = counts

    #pin

    if add_meta:
        # get meta data from the meta data group in file
        m = f['meta']
        meta = {}
        read_dict_cont(m, meta)

        attr = {}
        # get general file attributes
        attr['general'] = {}
        for k in list(f.attrs.keys()):
            attr['general'][k] = f.attrs[k]

        # and get additional meta data from images attributes
        for ch in chan_list:
            image_name = 'l15_images/image_' + ch.lower()
            im = f[image_name]
            
            attr[ch] = {}
            
            for k in list(im.attrs.keys()):
                attr[ch][k] = im.attrs[k]
            
        meta['attributes'] = attr
        
        var['meta'] = meta

    if add_geo:
        # get meta data from the meta data group in file
        g = f['geometry']
        geo = {}
        read_dict_cont(g, geo)
        
        var['geometry'] = geo


    # close hdf file
    f.close()
    # ================================================================


    return var



######################################################################
######################################################################


if __name__ == '__main__':

    t = datetime.datetime(2012,1,1,12,0)
    ch = ['IR_108', 'IR_120']

    sev = get_seviri_chan(ch,t, add_meta=True)
