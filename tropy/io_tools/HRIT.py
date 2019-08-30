#!/usr/bin/env python
######################################################################
# 
# DESCRIPTION
# ===========
#
# The module is part of the INPUT / OUTPUT tool box and concentrates
# on the HRIT file format, and issues concerning file handling within
# the archive structure.
#
# ####################################################################


# load libraries -----------------------------------------------------
import sys, os, glob
import numpy as np
import h5py
import datetime
import subprocess
import array
import struct
from . import bit_conversion as bc
import time
# ====================================================================

######################################################################
######################################################################



def SEVIRI_channel_list():

    scl = {1:'VIS006', 2:'VIS008', 3: 'IR_016', 4: 'IR_039',\
               5: 'WV_062', 6: 'WV_073', 7: 'IR_087', 8: 'IR_097',\
               9: 'IR_108', 10: 'IR_120', 11: 'IR_134', 12: 'HRV'}
    return scl

######################################################################
######################################################################
           


def bit_conversion(input_set, inbit=8, outbit=10):
    '''
    An integer array which is based on a n-bit representation is
    converted to an integer array based on a m-bit representation.

    

    USAGE:
    =====
    output_set  = bit_conversion(input_set, inbit=8, outbit=10)


    INPUT:
    =====
    input_set : numpy integer array in a n-bit representation
    inbit     : number of bits for input_set (n)
    outbit    : number of bits for output_set (m)


    OUPUT:
    =====
    output_set : numpy integer array in a m-bit representation

    
    COMMENT:
    ========
    The size of input_set should be a multiple of inbit / outbit!


    '''
    
    # convert input integer array to its bit representation ----------

    # get dimension of input array
    in_dim = input_set.shape[0]

    # base of the input array
    base_in = pow(2,np.arange(inbit,dtype='uint16'))

    # initialize bit martix 
    bits_in = np.zeros((inbit,in_dim)).astype( np.uint16 )

    # start the residual part 
    res = input_set

    # loop over the base 
    for n,b in enumerate(base_in[::-1]):
    
        # get the bit belonging to the base b
        bits_in[n] = res // b      

        # calculate the residual 
        res = res % b
   
    # ================================================================

    

    # derive output array from bit representation --------------------
  
    # dimension of the output array
    out_dim = in_dim * inbit / outbit
        
    # get the bits of the output by reshaping the input bits
    bits_out = bits_in.reshape((outbit,out_dim),order='F')
    
    # base of the output array
    base_out = pow(2,np.arange(outbit,dtype=base_in.dtype))

    # initialize the output array
    output_set = np.zeros((out_dim),dtype=input_set.dtype)

    # project the output bits on the output base
    for n,b in enumerate(base_out[::-1]):
        output_set += bits_out[n] * b 

    # ================================================================

    return output_set

######################################################################
######################################################################


def read_HRIT_seg_data(hrit_file): # LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    '''
    Reads directly the counts data for an HRIT segment.
    

    USAGE:
    =====
    seg_data = read_HRIT_seg_data(hrit_file)


    INPUT:
    =====
    hrit_file : HRIT file of the considered segment


    OUPUT:
    =====
    seg_data: integer segment data field of radiance counts 

    '''


    # open HRIT as binary file .......................................
    hf = open(hrit_file,'rb')

    # HRIT primary header record -------------------------------------
    h1_type = array.array('B')
    h1_type.fromfile(hf, 1)

    h1_length = array.array('B')
    h1_length.fromfile(hf, 2)
    n1_length = struct.unpack_from('>H', h1_length)[0]

    h1_file_type_code = array.array('B')
    h1_file_type_code.fromfile(hf, 1)

    h1_total_hlength = array.array('B')
    h1_total_hlength.fromfile(hf, 4)
    n1_total_hlength = struct.unpack_from('>L', h1_total_hlength)[0]

    h1_data_length = array.array('B')
    h1_data_length.fromfile(hf, 8)
    n1_data_length = struct.unpack_from('>Q', h1_data_length)[0]
    # ================================================================


    # HRIT secondary header record -----------------------------------
    h2_type = array.array('B')
    h2_type.fromfile(hf, 1)

    h2_length = array.array('B')
    h2_length.fromfile(hf, 2)
    n2_length = struct.unpack_from('>H', h2_length)[0]

    h2_nb = array.array('B')
    h2_nb.fromfile(hf, 1)
    NB = h2_nb[0]

    h2_nc = array.array('B')
    h2_nc.fromfile(hf, 2)
    NC = struct.unpack_from('>H', h2_nc)[0]

    h2_nl = array.array('B')
    h2_nl.fromfile(hf, 2)
    NL = struct.unpack_from('>H', h2_nl)[0]

    h2_comp_flag = array.array('B')
    h2_comp_flag.fromfile(hf, 1)

    if h2_comp_flag[0] != 0:
        print('ERROR: only uncompressed data can be retrieved')
        return
    # ================================================================


    # check dimensions -----------------------------------------------
    if NB*NC*NL != n1_data_length :
        print('dimension error!')
        print(NB*NC*NL)
        print(n1_data_length)
        return
    else:
        Ndat = n1_data_length
    # ================================================================


    # read binary values ---------------------------------------------
    hf.seek(n1_total_hlength)
   
    bv = np.fromfile(hf, dtype=np.uint8)

        
    # convert the byte (8-bit) array to 10-bit integer array
    datvalues = bc.bit_conversion(bv)

    # get the field
    seg_data = datvalues.reshape(NL,NC,order='C')[::-1,::-1]

    # ================================================================


    hf.close()


    return seg_data # TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT


######################################################################
######################################################################


def write_HRIT_seg_data(seg_data, hrit_file): # LLLLLLLLLLLLLLLLLLLLLL
    '''
    Reads directly the counts data for an HRIT segment.
    

    USAGE:
    =====
    write_HRIT_seg_data(seg_data, hrit_file)


    INPUT:
    =====
    hrit_file : HRIT file of the considered segment
    seg_data: integer segment data field of radiance counts 

    '''

    # open HRIT as binary file .......................................
    hf = open(hrit_file,'r+b')

    # HRIT primary header record -------------------------------------
    h1_type = array.array('B')
    h1_type.fromfile(hf, 1)

    h1_length = array.array('B')
    h1_length.fromfile(hf, 2)
    n1_length = struct.unpack_from('>H', h1_length)[0]

    h1_file_type_code = array.array('B')
    h1_file_type_code.fromfile(hf, 1)

    h1_total_hlength = array.array('B')
    h1_total_hlength.fromfile(hf, 4)
    n1_total_hlength = struct.unpack_from('>L', h1_total_hlength)[0]

    h1_data_length = array.array('B')
    h1_data_length.fromfile(hf, 8)
    n1_data_length = struct.unpack_from('>Q', h1_data_length)[0]
    # ================================================================



    # HRIT secondary header record -----------------------------------
    h2_type = array.array('B')
    h2_type.fromfile(hf, 1)

    h2_length = array.array('B')
    h2_length.fromfile(hf, 2)
    n2_length = struct.unpack_from('>H', h2_length)[0]

    h2_nb = array.array('B')
    h2_nb.fromfile(hf, 1)
    NB = h2_nb[0]

    h2_nc = array.array('B')
    h2_nc.fromfile(hf, 2)
    NC = struct.unpack_from('>H', h2_nc)[0]

    h2_nl = array.array('B')
    h2_nl.fromfile(hf, 2)
    NL = struct.unpack_from('>H', h2_nl)[0]

    h2_comp_flag = array.array('B')
    h2_comp_flag.fromfile(hf, 1)

    if h2_comp_flag[0] != 0:
        print('ERROR: only uncompressed data can be retrieved')
        return
    # ================================================================


    # check dimensions -----------------------------------------------
    if NB*NC*NL != n1_data_length :
        print('dimension error!')
        print(NB*NC*NL)
        print(n1_data_length)
        return
    else:
        Ndat = n1_data_length
    # ================================================================

        
    # transform 10-bit segment data to 8-bit -------------------------
    s_10bit = seg_data[::-1,::-1].reshape(NL*NC,order='C')
    s_8bit  = bit_conversion(s_10bit,inbit=10,outbit=8) 
    # ================================================================


    # read binary values ---------------------------------------------
    hf.seek(n1_total_hlength)
    hf.write(np.array(s_8bit,dtype='uint8'))
    # ================================================================


    hf.close()


    return seg_data # TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

######################################################################
######################################################################


def read_slope_offset_from_prolog(pro_file, NJUMP = 386993+72): # LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    '''
    Reads directly slope and offset from HRIT PRO file.
    

    USAGE:
    =====
    slope, offset = read_slope_offset_from_prolog(pro_file)


    INPUT:
    =====
    pro_file : PROLOG file of the considered time slot


    OUPUT:
    =====
    slope: slope of the data
    offset: offset of the data

    
    COMMENT:
    ========
    Following transformation is applied:

    DATA = SLOPE * COUNTS  +  OFFSET
    
    if COUNTS are valid data points!

    
    !! ATTENTION !! Jump in binary file is hard-coded!

    '''


    # read SEVIRI CHANNEL LIST
    sev_chan = SEVIRI_channel_list()

    # open HRIT as binary file .......................................
    hf = open(pro_file,'rb')

    # HRIT primary header record -------------------------------------
    h1_type = array.array('B')
    h1_type.fromfile(hf, 1)

    h1_length = array.array('B')
    h1_length.fromfile(hf, 2)
    n1_length = struct.unpack_from('>H', h1_length)[0]

    h1_file_type_code = array.array('B')
    h1_file_type_code.fromfile(hf, 1)

    h1_total_hlength = array.array('B')
    h1_total_hlength.fromfile(hf, 4)
    n1_total_hlength = struct.unpack_from('>L', h1_total_hlength)[0]

    h1_data_length = array.array('B')
    h1_data_length.fromfile(hf, 8)
    n1_data_length = struct.unpack_from('>Q', h1_data_length)[0]
    # ================================================================



    # read binary values ---------------------------------------------
    # jump to position
    hf.seek(n1_total_hlength + NJUMP)

    # read slope
    slope_offset_le = array.array('d')
    # slope_offset_le.read(hf,24)
    slope_offset_le.fromfile(hf,24)
    slope_offset = np.array(struct.unpack_from('>24d',slope_offset_le),'float').reshape(12,2)

    slope={}
    for i, sl in enumerate(slope_offset[:,0]):
        slope[sev_chan[i + 1]] = sl

    offset={}
    for i, os in enumerate(slope_offset[:,1]):
        offset[sev_chan[i + 1]] = os
    # ================================================================

    hf.close()


    return slope, offset # TTTTTTT

######################################################################
######################################################################


def combine_segments(seg_dict, missing = -1):
    ''' 
    Combines individual segments of SEVIRI scans. Segments will be sorted
    by key.
    

    USAGE:
    =====
    combined = combine_two_segments(seg_list)


    INPUT:
    =====
    seg_list: dictionary of segment, sorted by key from bottom to top
    
    OUTPUT:
    ======
    combined: combination of both segments
    '''

    
    # get the ste of keys
    ks = list(seg_dict.keys())

    # and sort them
    ks.sort()
    ks.reverse()

    combined = []

    for k in ks:
        
        # get segment data
        seg = seg_dict[k]

        # and combine them
        combined.append(seg) 
        
    
    return np.ma.masked_equal(np.concatenate(combined), missing)

######################################################################
######################################################################



def combine_two_segments(seg1,seg2):
    ''' 
    Combines two segments where seg1 is the upper and seg2 the 
    lower segment.
    

    USAGE:
    =====
    seg = combine_two_segments(seg1,seg2)


    INPUT:
    =====
    seg1: upper segment
    seg2: lower segment
    
    OUTPUT:
    ======
    seg: combination of both segments
    '''

    
    (NL,NC) = seg1.shape
    
    # combine both segments
    seg = append(seg1,seg2).reshape(NL*2,NC)

    
    return seg

######################################################################
######################################################################


def divide_two_segments(seg):
    ''' 
    Divides one segment into two with half the line size.
    

    USAGE:
    =====
    seg1, seg2 = divide_two_segments(seg)


    INPUT:
    ======
    seg: combination of two segments


    OUTPUT:
    =====
    seg1: upper segment
    seg2: lower segment
    
     '''

    # get shape of the combination
    (NL2,NC) = seg.shape
    
    # line number of only one segment
    NL2 = NL/2
    
    # upper segment
    seg1 = seg[:NL,:]
    
    # lower segment
    seg2 = seg[NL:,:]

    
    return seg1, seg2

######################################################################
######################################################################


def channel_segment_sets(set_name):

    '''
    Outputs a predefined set of channels and segments'

    USAGE:
    =====
    chan_seg = channel_segment_sets(set_name)


    INPUT:
    ======
    set_name: name of a predefined set
    

    OUTPUT:
    =======
    chan_seg: Dictionary of channels with a list of segments. 

    '''


    # standard segment range
    eu_segs = list(range(7,9))
    full_segs = list(range(1,9))

    sets = {}

    # European Area
    sets['ir-eu'] = {'IR_039':eu_segs, 'WV_062':eu_segs, 'WV_073':eu_segs, 'IR_087':eu_segs, \
                  '  IR_097':eu_segs, 'IR_108':eu_segs, 'IR_120':eu_segs, 'IR_134':eu_segs}
    
    sets['eu'] = {'VIS006':eu_segs, 'VIS008':eu_segs, 'IR_016':eu_segs,\
                    'IR_039':eu_segs, 'WV_062':eu_segs, 'WV_073':eu_segs, 'IR_087':eu_segs, \
                  '  IR_097':eu_segs, 'IR_108':eu_segs, 'IR_120':eu_segs, 'IR_134':eu_segs,
                    'HRV':list(range(20,24))}

    sets['nc-eu'] = {'VIS006':eu_segs, 'VIS008':eu_segs, 'IR_016':eu_segs}

    sets['ir108-eu'] = { 'IR_108': eu_segs }
    
    
    # full disk service

    sets['ir-full'] = {'IR_039':full_segs, 'WV_062':full_segs, 'WV_073':full_segs, \
                           'IR_087':full_segs, 'IR_097':full_segs, 'IR_108':full_segs,\
                           'IR_120':full_segs, 'IR_134':full_segs}
    
    sets['full'] = {'VIS006':full_segs, 'VIS008':full_segs, 'IR_016':full_segs,\
                    'IR_039':full_segs, 'WV_062':full_segs, 'WV_073':full_segs, 'IR_087':full_segs, \
                  '  IR_097':full_segs, 'IR_108':full_segs, 'IR_120':full_segs, 'IR_134':full_segs,
                    'HRV':list(range(1,24))}

    
    sets['nc-full'] = {'VIS006':full_segs, 'VIS008':full_segs, 'IR_016':full_segs}

    sets['ir108-full'] = { 'IR_108':full_segs }

    
    if set_name in sets:
        chan_seg = sets[set_name]
    else:
        print('ERROR: unknown channel-segment set!')
        print('   available keys:')
        print(list(sets.keys()))
        chan_seg = None
        
    
    return chan_seg


######################################################################
######################################################################


def segments_for_region(region, Nrow=3712, Nseg=8, seg_size = 464):

    # extract row and column range from region definition
    (rowmin,rowmax),(colmin,colmax) = region

    # get range of the different segment
    #   indices are counted from top to bottom,
    #   but segments are counted the other-way-around
    segmin = Nrow - (np.arange(Nseg) + 1) * seg_size
    segmax = Nrow - np.arange(Nseg) * seg_size


    # get array of segments
    segs = np.arange(1,Nseg + 1)[np.logical_and(segmin < rowmax, segmax > rowmin)]

    return segs

######################################################################
######################################################################



def get_HRIT_from_arch(day, chan_seg =  channel_segment_sets('nc-full'), 
                       scan_type = 'pzs',
                       tarname = None,
                       arch_dir = None,
                       out_path = None):
    
    '''
    Gets and decompresses the original HRIT files from archive. A
    typical set of channels and segments can be chosen.


    USAGE:
    =====
    file_list =  get_HRIT_from_arch(day, chan_set = 'nc-full', scan_type = 'pzs' ,\
                           arch_dir = None, out_path = None)


    INPUT:
    ======
    day: datetime object which include day and time of MSG time slot
    chan_set = name for a channel-segment set
    scan_type: which type of msg is used, 
                    * ) allowed keys are: 'pzs' and 'rss'

    arch_dir = archive directory
    out_path = directory where the HRIT files are saved.



    OUTPUT:
    =======
    file_list : list of extracted HRIT files

    '''

 
    # create filename or take given one ------------------------------
    if tarname != None:
         filename_template = tarname
    else:
    # get the right satellite specification
        if scan_type not in ('rss', 'pzs', 'hrs'):
            #        raise ValueError
            print('ERROR: ',scan_type,' is not a valid scan_type option') 
            return None

    # set default archive directory
        if not arch_dir:
            arch_dir = '/mnt/zephyr/u1/EUMETCAST_ARCHIVE/msg-sevi/%s/l15_hrit/' % scan_type


    # set time and date strings for file selection
        time_string = day.strftime('%Y%m%d%H%M')
        date_string = day.strftime('/%Y/%m/%d/')


    # create filename template (NOTE that number of MSG have been left out!)
        filename_template = arch_dir + date_string + "msg?-sevi-%s-l15-%s.tar" % (time_string, scan_type)
    # ================================================================


    # check if tar file exists ---------------------------------------
    # try to match template with existing file 
    arch_file_list = glob.glob(filename_template)          


    if arch_file_list:
        arch_file = arch_file_list[0]
    else:
        print('ERROR: ', filename_template,' does not exist!')
        return None
    # ================================================================




    # get content list of tar file -----------------------------------

    # build tar command of reading content
    tar_com = 'tar -tf '+ arch_file

    cont_list = subprocess.getoutput(tar_com).split()
#    print cont_list
    

    # extract tar files
    extract_list =[]
    for HRIT_file in cont_list:
        
        # get prolog file
        if 'PRO' in HRIT_file:
            extract_list.append(HRIT_file)
            
        # get epilog file
        if 'EPI' in HRIT_file:
            extract_list.append(HRIT_file)
            
        # get channel segments
        for ch in chan_seg:
            for seg in chan_seg[ch]:
                if str(seg).zfill(6) in HRIT_file and ch in HRIT_file:
                    extract_list.append(HRIT_file)
    # ================================================================


    # set output directory -------------------------------------------
    old_path = os.getcwd()

    if not out_path:
        out_path = '/tmp/hrit' + str(np.random.randint(1e10))

    if not os.access(out_path, os.F_OK):
        os.mkdir(out_path)

    os.chdir(out_path)
    # ================================================================


    hrit_list = []
    for hrit_file in extract_list:
        os.system('tar -xf '+ arch_file+' '+ hrit_file)

        # decompress the segments
        if '-C_' in hrit_file:
            os.system('xRITDecompress ' + hrit_file)
            os.system('rm -f ' + hrit_file)
        
            hrit_list.append(out_path + '/' + hrit_file.replace('C_','__'))

        # or do nothing for epi- and prolog
        else:
            hrit_list.append(out_path + '/' + hrit_file)
            

#    print hrit_list

    os.chdir(old_path)

    return hrit_list

######################################################################
######################################################################


def read_HRIT_data(day, chan_seg, calibrate=True, scan_type = 'pzs', \
                   arch_dir = None, standard_seg = (464, 3712), add_meta = False, **kwargs):

    '''
    Read HRIT data given a time slot, satellite type and channel-segment set


    USAGE:
    =====
    combined =   read_HRIT_data(day, chan_seg, calibrate=True, scan_type = 'pzs', \
                   arch_dir = None, add_meta = False, **kwargs)


    INPUT:
    ======
    day =  datetime object which include day and time of MSG time slot
    chan_seg = a channel-segment set, dictionary of channel names with lists of segments
    calibrate = flag, if slope / offset calibration should be done  
    scan_type =  which type of msg is used, 
                    * ) allowed keys are: 'pzs' and 'rss'

    arch_dir = archive directory
    add_meta = flag, if meta data should be added
    out_path = directory where the HRIT files are saved.



    OUTPUT:
    =======
    combined = dictionary of retrieved channels with combined segments 
    '''

    kwargs['scan_type'] = scan_type
    kwargs['arch_dir'] = arch_dir


    # get HRIT data from archive
    hrit_list = get_HRIT_from_arch(day, chan_seg = chan_seg, **kwargs)


    # get channel segments
    chan_data = {}

    for ch in chan_seg:
        chan_data[ch] = {}

        for seg in chan_seg[ch]:

           # look for hrit file in file list, empty for missing segments
           hchan_seg = [hr for hr in hrit_list if ch in hr and str(seg).zfill(6) in hr] 

           if hchan_seg:
               hrit_file = hchan_seg[0]
               print('... reading ', hrit_file)
               chan_data[ch][seg] = read_HRIT_seg_data(hrit_file)
           else:
               print('... missing segment %s of channel %s' % (seg, ch))
               chan_data[ch][seg] = -np.ones(standard_seg)

    print() 
    print('Combine segments')
    # combine all segments to a "full" images 
    combined = {}
    for ch in chan_seg:
        combined[ch] = combine_segments(chan_data[ch])
        
    # get name of PRO file
    for hrit_file in hrit_list: 
        if 'PRO' in  hrit_file: 
            pro_file = hrit_file
 
    # calibrate satellite counts with slope and offset


    if calibrate:
        print()                        
        print('Do calibration')  
        slope, offset = read_slope_offset_from_prolog(pro_file)

        for ch in chan_seg:

            s = slope[ch]
            o = offset[ch]
            c = combined[ch]

            combined[ch] = np.where(c == 0, 0, s*c + o)
        
 
    # add meta data
    if add_meta:
        combined['meta'] = {}
        combined['meta']['time']   = day
        combined['meta']['hrit_list'] = hrit_list
        combined['meta']['prolog_file'] = pro_file

        if calibrate:
            combined['meta']['offset'] = offset
            combined['meta']['slope']  = slope

   
       

    # clean up
    for hrit_file in hrit_list:
        os.system('rm -f ' + hrit_file)

    out_dir = os.path.dirname(hrit_file)
    if not os.listdir(out_dir):
        os.rmdir(out_dir)

    return combined
    

######################################################################
######################################################################


if __name__ == '__main__':
     
    hfile='/tmp/hrit5481570323/H-000-MSG2__-MSG2________-IR_108___-000007___-201206081200-__'
   
    for i in range(24):
        seg = read_HRIT_seg_data(hfile)
