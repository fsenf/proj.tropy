#!/usr/bin/env python

import os

# first try platform independent solution
try:
    from _f90_bit_conversion import bit_conversion
except:

    try:
        hostname = os.uname()[1]
    
        if 'satyr' in hostname:
            from _f90_bit_conversion_satyr import bit_conversion
    
        elif 'lenovo-senf' in hostname:
            from _f90_bit_conversion_lenovosenf import bit_conversion
    
        else:
            from _f90_bit_conversion_altair  import bit_conversion
    except:
        print('bit_conversion error')
