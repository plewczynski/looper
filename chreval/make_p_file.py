#!/usr/bin/env python3

from math import exp, log
import sys
import os
import string
import argparse


# main tool objects
from Functions import KahanSumExp
from FileTools import getHeadExt

# for entropy calculations
from Constants import kB # [kcal/molK] (Boltzmann constant)

# A hardwired infinity. Somewhat hate this, but it seems there is no
# way around this darn thing.
from Constants import INFINITY

# for the KahanSumExp
from Functions import flag_KahanSumExp_USE_SORT

kBT = kB * 310.15


def calc_Z(dGdata):
    debug_Boltzmann = True
    
    # Result: $Z = 1 + \sum_k {\exp(-dG_k / kBT) }$
    
    # The method for calculating the partition function may seem a
    # bit odd at first brush. After all, the mathematical formula
    # is really quite straight forward.
    
    # The problem is that floats have a limit in size
    # that can be easily exceeded. So we have first construct the
    # exponents and then analyze them through a subsidiary
    # function KahanSumExp() to add the exponentials together.
    
    # math.exp
    explist = []
    for dG_k in dGdata:
        explist += [ -dG_k / kBT ]
    #
    
    Z, shift = KahanSumExp(explist)
    re_wt = float(2**shift)
    
    
    if shift > -1: # self.debug:
        if debug_Boltzmann:
            print ("Z = ", Z)
            print ("shift: %4d,  new exp wt: %.1f" % (shift, re_wt))
        #
        
        # print out the first ten arguments of on the partition
        # function exponential term
        
        #markings                     #        #
        print ("ndx    arg       dG/kBT          dG        p(ndx)")
        for k in range(0, len(explist)):
            if k > 10: # could be a huge list
                break
            #
            
            print ("%2d  %8.2f   %8.2f  %12.2f    %8.3g" \
                   % (k,
                      explist[k],
                      (-dGdata[k]/kBT),
                      dGdata[k],
                      calc_p(dGdata[k], Z, shift)))
        #
        
        #sys.exit(0)
    
    #

    return Z, shift

#

def calc_p(dG, Z, shift):
    
    # math.exp
    if shift > 0:
        # print ((float(self.shift)*log(2)))
        exponent = - (dG / kBT)
        if exponent - shift * log(2) < -708.396:
            p = 0.0 # surely, this is neglectable
        else:
            p = exp( exponent - float(shift)*log(2) ) / Z
        #
        
    else:
        exponent = - (dG / kBT)
        if exponent - shift * log(2) < -708.396:
            p = 0.0 # surely, this is neglectable
        else:
            p = exp( exponent ) / Z
        #
        
    #
    
    return p
#

class Make_pFile(object):
    def __init__(self):
        self.addAB = False
    #
    
    
    def make_p_file(self, flhd):
        
        fp = open(flhd + "_dG.dat", 'r')
        lfp_dG = fp.readlines()
        fp.close()
        
        fp = open(flhd + ".dat", 'r')
        lfp = fp.readlines()
        fp.close()
        print (len(lfp), len(lfp_dG))
        
        pdataset = []
        #kk = 0
        w = 100
        header = ''
        pflnm = flhd + "_pe.dat"
        
        for j in range(0, len(lfp_dG)):
            st = lfp_dG[j].strip()
            # print (st)
            # find the cutoff at length
            """
            print ((' '*w) + 'v')
            print (st[:w])
            sys.exit(0)
            """
            if j < 23:
                header += st + '\n'
            #
            
            if j == 23:
                if self.addAB:
                    print (header)
                    tt = st[:w] + "    A           B           p_k  ----------->>>>"
                    print (tt)
                    header += tt + '\n'
                    fp = open(pflnm, 'w')
                    fp.write(header)
                    fp.close()
                    #sys.exit(0)
                else:
                    
                    print (header)
                    tt = st[:w] + "    p_k  ----------->>>>"
                    print (tt)
                    header += tt + '\n'
                    fp = open(pflnm, 'w')
                    fp.write(header)
                    fp.close()
                    #sys.exit(0)
                #
                
            #
            
            if len(st) == 0:
                continue
            #
            
            if st[:1] == '#':
                continue
            #
            
            srv = lfp[j].strip().split()
            #print (j)
            #print (srv[25], srv[26])
            #sys.exit(0)
            
            #print (header)
            #sys.exit(0)
            # extract this dG dataset
            stv = st.split()
            print (stv[10:])
            
            leader = st[:w]
            if self.addAB:
                leader += "  %8.4f  " % (float(srv[25].strip()))
                leader += "  %8.4f  " % (float(srv[26].strip()))
            #
            
            dGdata = []
            for sdG_k in stv[10:]:
                dGdata += [ float(sdG_k) ]
            #|endfor
            
            Z, shift = calc_Z(dGdata)
            
            pdata = []
            for dG_k in dGdata:
                pdata += [ calc_p(dG_k, Z, shift) ]
            #|endfor
            
            
            pdataline = leader
            
            for p_k in pdata:
                pdataline += "%10.4g  " % p_k
            #|endfor
            
            pdataline += '\n'
            fp = open(pflnm, 'a')
            fp.write(pdataline)
            fp.close()
            
            
        #|endfor
        
    #
    
#


def main(cl):
    
    """@
    
    NOTE: you must specify the root file, not the *_dG.dat, *_p.dat,
    etc. files.
    
    originally, this was built to fix a problem in the middle of a
    calculation where the *_p.dat file was the same as *_dG.dat
    file. However, I also added the A and B columns, so then we had
    the problem that the data was being produced this way. I probably
    should have kept the original format, but I didn't. Now what
    happens is a new file *_pe.dat is made.
    
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', default="", type=str,
                        dest='inflnm',
                        help='read in main analyze_loops.py output file \
                        and rebuild the *_p.dat file with A and B columns. \
                        All files must have the extension \"dat\".')

    parser.add_argument('-addAB', action='store_true', default=False,
                        dest='addAB',
                        help='include the A and B info in file.')
    
    args = parser.parse_args()
    
    flnm  = args.inflnm
    optAB = args.addAB

    if not os.path.exists(flnm):
        print ("ERROR: %s does not exist")
        sys.exit(1)
    #
    
    flhd, ext = getHeadExt(flnm)
    if not ext == "dat":
        print ("ERROR: all extensions for analyze_loops files should have")
        print ("       the extension \"dat\" != \"%s\" <- input file." % ext)
        sys.exit(1)
    #

    a = Make_pFile()
    a.addAB = optAB
    
    a.make_p_file(flhd)
    
    
#
    

if __name__ == '__main__':
    # running the program
    main(sys.argv)
#    
