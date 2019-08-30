#!/usr/bin/env python

import os

######################################################################
######################################################################

def test(fname):

    ''' 
    Returns status of a file:

       True: if open
       False: if not open

    '''
   
    user = os.environ['USER'] 
    ret = os.system('lsof -u %s | grep %s > /dev/null' % (user, fname))

    if ret == 0:
        return True
    else:
        return False

######################################################################
######################################################################


if __name__ == '__main__':
    
    fname = 'test.h5'

    print(test(fname))

    
    
