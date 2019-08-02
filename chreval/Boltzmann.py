#!/usr/bin/env python

"""Main Program:  Boltzmann.py 

Classes:       Boltzmann

Author:        Wayne Dawson
creation date: mostly 2016 and up to March 2017 (separated from chreval 180709)
last update:   190522
version:       0

Purpose:

calculates the the effective Boltzmann distribution


Comments

None presently. It is one of the few modules that I have not had to
change drastically over time.

"""

from math import exp, log
import sys
import os
import string

# main tool objects
from Functions import KahanSumExp


# LThread object representation
from LThread import DispLThread
from LThread import LThread


# ################################################################
# ######################  Global constants  ######################
# ################################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# for entropy calculations
from Constants import kB # [kcal/molK] (Boltzmann constant)

# A hardwired infinity. Somewhat hate this, but it seems there is no
# way around this darn thing.
from Constants import INFINITY

# for the KahanSumExp
from Functions import flag_KahanSumExp_USE_SORT

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ################################################################
# ################################################################


# ################################################################
# ###########  global functions/objects and constants  ###########
# ################################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# 1. control parameters for the program

PROGRAM      = "Boltzmann.py"  # name of the program

# tests
TEST5        = True  # test the Kahan algorithm

# debugging options 
SHOWMAIN     = False # Debugging in main()


# 2. special local global function
def usage():
    print "USAGE: %s" % PROGRAM
#
# At present, there is not much reason to use this at all. The only
# thing that can be done with this object independent of the programs
# that use it is to run a test, and that test is not so interesting
# presently. Nevertheless, maybe it will be useful at some point.

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ################################################################
# ################################################################

#####################################################################
################  CALCULATE BOLTZMANN DISTRIBUTION  #################
#####################################################################

class Boltzmann:
    def __init__(self, calc):
        self.debug_Boltzmann = False
        self.Z    = 1.0
        self.explist = [] # to ensure that the exp doesn't over extend
        self.shift   = 0  #
        self.rewt    = 1.0
        self.kB   = kB
        self.calc = calc
        self.fe   = calc.fe
        self.N    = calc.N
        self.T    = calc.T
        self.flag_set_Z = False
        
        self.debug = False
        
        # to ensure that Z is calculated before doing other operations.
    #
    
    def calc_Z(self, lt, T):
        kBT = self.kB * T
        
        # Result: $Z = 1 + \sum_k {\exp(-dG_k / kBT) }$
        
        # The method for calculating the partition function may seem a
        # bit odd at first brush. After all, the mathematical formula
        # is really quite straight forward.
        
        # The problem is that floats have a limit in size
        # that can be easily exceeded. So we have first construct the
        # exponents and then analyze them through a subsidiary
        # function KahanSumExp() to add the exponentials together.
        
        # math.exp
        self.explist = []
        for ltk in lt:
            self.explist += [ -ltk.dG / kBT ]
        #
        
        self.Z, self.shift = KahanSumExp(self.explist)
        re_wt = float(2**self.shift)
        self.flag_set_Z = True # calculation of Z is now done!
        
        # 
        if self.shift > -1: # self.debug:
            if self.debug_Boltzmann:
                print "Z = ", self.Z
                print "shift: %4d,  new exp wt: %.1f" % (self.shift, re_wt)
            #

            # print out the first ten arguments of on the partition
            # function exponential term
            #                             #        #
            print "ndx    arg       dG/kBT          dG        p(ndx)"
            for k in range(0, len(self.explist)):
                if k > 10: # could be a huge list
                    break
                print "%2d  %8.2f   %8.2f  %12.2f    %8.3g" \
                    % (k, self.explist[k], (-lt[k].dG/kBT), lt[k].dG, self.calc_p(lt[k], T))
            #
            # sys.exit(0)
        #
    #
    
    
    def set_LThread_TdS(self, lt, T):
        if not self.flag_set_Z:
            print "ERROR(Boltzmann.calc_p()): Z is undefined"
            sys.exit(1)
        #
        
        for m in range(0, len(lt)):
            lt[m].TdS = 0.0
            for tr in lt[m].thread:
                v = tr.ij_ndx
                i = v[0]
                j = v[1]
                # print "ij= ", i, j
                lt[m].TdS += self.fe.TdS(i,j,T)
            #
        #
        if self.debug_Boltzmann:
            dt = DispLThread(self.calc.N)
            for m in range(0, len(lt)):
                ss_m = dt.makeLThreadDotBracket_1b(lt[m], True)
                print m, lt[m].TdS, ss_m, lt[m].p
            sys.exit(0)
        #
        return lt
    #
    
    def set_LThread_p(self, lt, T):
        if not self.flag_set_Z:
            print "ERROR(Boltzmann.calc_p()): Z is undefined"
            sys.exit(1)
        # Note:
        # $p = exp( -dG_k / kBT ) / Z$ 
        # $p = exp( -dG_k / kBT ) / \[ 1 + \sum_k {\exp(-dG_k / kBT)} \]$ 
        
        kBT = self.kB * T
        # math.exp
        if self.shift > 0:
            for k in range(0, len(lt)):
                exponent = - (lt[k].dG / kBT)
                if exponent - self.shift * log(2) < -708.396:
                    lt[k].p = 0.0
                else:
                    lt[k].p = exp( exponent - float(self.shift)*log(2) ) / self.Z
            #
        else:
            for k in range(0, len(lt)):
                exponent = - (lt[k].dG / kBT)
                if exponent - self.shift * log(2) < -708.396:
                    lt[k].p = 0.0 # surely, this is neglectable
                else:
                    lt[k].p = exp( exponent ) / self.Z
                #
            #
        #
        return lt
    #
    
    def calc_p(self, ltk, T):
        if not self.flag_set_Z:
            print "ERROR(Boltzmann.calc_p()): Z is undefined"
            sys.exit(1)
        #
        kBT = self.kB * T
        # math.exp
        if self.shift > 0:
            # print (float(self.shift)*log(2))
            exponent = - (ltk.dG / kBT)
            if exponent - self.shift * log(2) < -708.396:
                p = 0.0 # surely, this is neglectable
            else:
                p = exp( exponent - float(self.shift)*log(2) ) / self.Z
        else:
            exponent = - (ltk.dG / kBT)
            if exponent - self.shift * log(2) < -708.396:
                p = 0.0 # surely, this is neglectable
            else:
                p = exp( exponent ) / self.Z
            #
        #
        return p
    #
#

# this is intended to try various cases where it is a bit important ot
# make adjustments.

def test5():
    global flag_KahanSumExp_USE_SORT
    flag_KahanSumExp_USE_SORT = True
    dG_vr_kBT = [10, 37, 34, 0.1, 0.0004, 34, 37.1, 37.2, 36.9, 709, 710, 711]
    Z, shift = KahanSumExp(dG_vr_kBT)
    print "{0} x 2^{1}".format(Z, shift)
    wt = float(2**shift)
    p = []
    for k in range(0, len(dG_vr_kBT)):
        p += [ exp( dG_vr_kBT[k] - float(shift)*log(2) ) / Z ]
    summation = 0.0
    for pk in p:
        summation += pk
        
    if shift > 1.0: # self.debug:
        print Z
        print "shift: ", shift, "new exp wt: ", wt
        for k in range(0, len(dG_vr_kBT)):
            print "%2d  %8.2f   %12.3g" % (k, dG_vr_kBT[k], p[k])
        print "summation of p: %12.5g" % summation
        # sys.exit(0)
    #
#    



def main(cl):
    # presently doesn't do anything
    print cl
#


if __name__ == '__main__':
    if TEST5: # test the Kahan algorithm
        test5()
    else:
        main(sys.argv)
    #
#
