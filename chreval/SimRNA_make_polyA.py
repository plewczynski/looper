#!/usr/bin/env python

"""@@@

Main Module:   SimRNA_make_polyA.py 

Classes:       

Functions:     make_polyA_seq

Author:        Wayne Dawson
creation date: parts 2016, made into a separate object 170314
last update:   170314
version:       0

Purpose:

After obtaining the meta-structure of chromatin, we need to build a 3D
structure. This program generates a polyA sequence of a specified
length that can be used to begin a SimRNA simulation using a zeroed
out statistical potential

"""

import sys
SHOWMAIN = False
PROGRAM = "SimRNA_make_polyA.py"

def usage():
    print "USAGE: %s integer" % PROGRAM
#

def make_polyA_seq(cl):
    #
    n = -1
    
    if SHOWMAIN:
        print "cl: ", cl
        print "number of args: %d" % (len(cl))
    #
    
    if len(cl) < 2:
        emsg = "ERROR: too few arguments"
        usage()
        sys.exit(1)
    #
    
    for clk in cl:
        if clk == "-h" or clk == "--help":
            usage()
            sys.exit(0)
        #
    #
    try:
        n = int(sys.argv[1])
    except ValueError:
        print "ERROR: command requires an integer argument"
        print "       entered argument '%s'" % sys.argv[1]
        usage()
        sys.exit(1)
    #
    
    print 'a'*n

#






if __name__ == '__main__':
    make_polyA_seq(sys.argv)
#
