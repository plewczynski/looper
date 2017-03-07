#!/usr/bin/env python

# builds a polyA sequence of a specified length

import sys
SHOWMAIN = False
PROGRAM = "SimRNA_make_polyA.py"

def usage():
    print "USAGE: %s integer" % PROGRAM
#

def main(cl):
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








if __name__ == '__main__':
    main(sys.argv)
