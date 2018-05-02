#!/usr/bin/python

# load libraries -----------------------------------------------------
import os
import datetime
from l15_msevi.msevi_config import _arch_config
# ====================================================================
   

def find_latest_file(directory, partial_file_name):
    '''
    Finds the most recent file in a directory given a piece of filename


    USAGE:
    =====
    lastest_file =  find_latest_file(directory, partial_file_name)


    INPUT:
    ======
    directory: directory name
    partial_file_name: a piece of the search filename


    OUTPUT:
    =======
    latest_file:  name of the latest file
                  -> None if directory of file does not exist


    COMMENT:
    ========
    No recursive search is done.

    '''


    # appends a slash if necessary
    if not directory[-1] == '/': 
        directory = directory + '/'

    # checks if directory exists and then list its content
    if os.access(directory, os.F_OK):
        # list all the files in the directory
        files = os.listdir(directory)
    else:
        return None
    
 
    # remove all file names that don't match partial_file_name string
    files = filter(lambda x: x.find(partial_file_name) > -1, files)

    if not files:
        return None

    # create a dict that contains list of files and their modification timestamps
    name_n_timestamp = dict([(x, os.stat(directory+x).st_mtime) for x in files])
 

    # get the file with the latest timestamp
    latest_file = max(name_n_timestamp, key=lambda k: name_n_timestamp.get(k))

    return latest_file


######################################################################
######################################################################



def find_latest_slot(scan_type='rss', arch_dir = None, time_string = None, arch_type='new', aparam = None):

    if arch_type == 'old':
        return find_latest_slot_oldarch(scan_type = scan_type, arch_dir = arch_dir, time_string = time_string)
    else:
        
        # new config
        c = _arch_config.copy()
        c.update({'scan_type':scan_type})

        
        return find_latest_slot_oldarch(scan_type = scan_type, 
                                        arch_dir = '{arch_dir}/msevi_{scan_type}/l15_hrit/'.format(**c), 
                                        time_string = _arch_config['time_string'])
	    


######################################################################
######################################################################



def find_latest_slot_oldarch(scan_type='hrs', arch_dir = None, time_string = None):

    # search for the latest MSG tar file -----------------------------

    to_day = datetime.datetime.today() + datetime.timedelta(days=1)
    date = to_day

    # main archive directory
    if arch_dir == None:
        main_arch_path = '/mnt/zephyr/u1/EUMETCAST_ARCHIVE/msg-sevi/%s/l15_hrit/' % scan_type
    else:
        main_arch_path = arch_dir


    # set default time string argument 
    if time_string == None:
        time_string = '%Y%m%d%H%M'
  

    latest_file = None

    while not latest_file:

        # set search path
        date_string = date.strftime('%Y/%m/%d/')
        search_path = main_arch_path +  date_string

        # search for latest tar file
        latest_file = find_latest_file(search_path,'tar')


        # if nothing is found then take the day before
        if not latest_file:
            date -=datetime.timedelta(days=1)
    # ================================================================

    print 'latest file:', search_path + latest_file

    # properties of arch file ----------------------------------------

    # CAUTION THIS DEPEND ON FILE NOMENCLATURE !!!!
    arch_time_string = latest_file.split('-')[2]

    date = datetime.datetime.strptime(arch_time_string, time_string)
    # ================================================================

    return date

######################################################################
######################################################################


if __name__ == '__main__':
    
    print find_latest_slot(scan_type = 'hrs')
    print find_latest_slot(scan_type = 'rss')

