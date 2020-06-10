#!/usr/bin/env python3

"""@

Main Program:  CVersion.py 

Functions:     get_Version
               get_now
               get_lastUpdate

Author:        Wayne Dawson
Creation date: 190817
Last Update:   200218 python3 compatible (from RVersion)
Version:       0.0

Purpose:

To unify the version numbers in a systematic way. 

"""

import sys
from datetime import datetime



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
__version_info__  = ("1.0", "1")
__lastUpdate_info__ = ('2020','02','18')
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
__lastUpdate__ = '-'.join(__lastUpdate_info__)        
__version__    = "version %s, release %s" % (__version_info__[0], __version_info__[1]) 

def get_Version():
    a = "(%s, last update %s)" % (__version__, __lastUpdate__)
    return a
#

def get_lastUpdate():
    return __lastUpdate__
#


def get_now():
    #v = datetime.now().strftime("%Y-%m-%d %H:%M:%S ")    
    v = datetime.now().strftime("%a %b %d %H:%M:%S %Y")    
    return v
#

def main(cl):
    print (get_now())
    print (get_lastUpdate())
    print (get_Version())
#



if __name__ == '__main__':
    main(sys.argv)
#
