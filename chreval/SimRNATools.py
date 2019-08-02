#!/usr/bin/env python

"""@@@


Main Module:   SimRNATools.py 

Classes:       SimRes
               SimRNAData

Functions:

Author:        Wayne Dawson
creation date: parts 2016, made into a separate object 170314
last update:   170314
version:       0

Purpose:

Like, ChPairData, I have some reservations about making this a special
object with a data format. Nevertheless, I think for managing input
and outputs, it is still useful to put this in the form of an
object. It would be next to impossible to convert the information from
a SimRNA restraint file into ChPairData, let alone LThread or Motif
data. Nevertheless, for processing, it seems better to group this as a
single object that can "dumb down" higher level data formats to
something useful for building SimRNA outputs.


Comments:

Particularly problematical with SimRNA data outputs is that the only
version of the program insists on only having data inputs with no room
for comments. I always hated that about SimRNA, because 6 months from
now, I have absolutely no information about _why_ I wrote down some
strange thing in the file. Admittedly, most people don't write
comments, and I often neglect to, but I also _do_ write comments. So
SimRNA restraint files are almost next to useless for storing anything
other than the absolutely undeniably irrefutable things that tell
SimRNA what to do, and nothing whatsoever more.

Perhaps later, I might make join the other tools involving SimRNA so
this can at least work in the direction of generating 3D structures
from Chreval outputs.

"""

from FileTools import FileTools
from FileTools import getHeadExt
from ChPair import ChPairData
import sys
import os

PROGRAM = "SimRNATools.py"

def usage():
    print "Usage: %s <file>.chpair" % PROGRAM
#


resdist = { 'N~N' : 9.0,  'P~P' : 18.7 }
res_wt  = { 'N~N' : 0.4,  'P~P' :  0.4 }

class SimRes:
    def __init__(self, i, chainA, atomA, j, chainB, atomB, distAB, wt, xi = 7.0):
        self.i = i
        self.j = j
        self.chainA = chainA
        self.atomA  = atomA
        self.chainB = chainB
        self.atomB  = atomB
        self.ABdist = distAB
        self.wt = wt
        self.xi = xi
    #
    
    def disp_SimRes(self, fmt = "slope"):
        # we must add 1 to the indices because most notation systems
        # start from 1 rather than zero.
        s = ''
        if fmt == "slope" or fmt == "SLOPE":
            s = "SLOPE    %s/%d/%s   %s/%d/%s %6.2f  %6.2f  %8.3f\n" % \
                (self.chainA, self.i + 1, self.atomA, self.chainB, self.j + 1, self.atomB, (self.ABdist - 0.5), (self.ABdist + 0.5), self.wt)
        elif fmt == "cle" or fmt == "CLE":
            s = "CLE     %s/%d/%s   %s/%d/%s %6.2f  %6.2f  %8.3f\n" % \
                (self.chainA, self.i + 1, self.atomA, self.chainB, self.j + 1, self.atomB, self.ABdist, self.xi, self.wt)
        #
        return s
    #
#



class SimRNAData:
    def __init__(self, fmt = 'MJB'):
        # options for employing SimRNA
        self.xi      = 7.0    # [nt] Kuhn length
        self.resType = "slope" # option for "slope" or "CLE"
        self.simres  = []
        self.fmt     = fmt # MJB or WKD (WKD allows comments!!!)
    #
    
    
    def checkBBlist(self, bblist):
        flag_pass = True
        for bb in bblist:
            if not resdist.has_key(bb):
                bbsplit = bb.split('~')
                if len(bbsplit) < 2:
                    print "ERROR: incomplete information on the bead ~ bead information"
                else:
                    print "ERROR: no definition for %s ~ %s interactions" % (bbsplit[0], bbsplit[1])
                #
                print "       input key:               '%s'" % bb
                sbb = ''
                for sr in resdist.keys():
                    sbb += "%s " % sr
                print "       allowed keys for resdist: %s" % sbb
                flag_pass = False
            #
            if not res_wt.has_key(bb):
                bbsplit = bb.split('~')
                if len(bbsplit) < 2:
                    print "ERROR: incomplete information on the bead ~ bead information"
                else:
                    print "ERROR: no definition for %s ~ %s interactions" % (bbsplit[0], bbsplit[1])
                #
                print "       input key:               '%s'" % bb
                sbb = ''
                for sr in res_wt.keys():
                    sbb += "%s " % sr
                print "       allowed keys: %s" % sbb
                flag_pass = False
            #
        #
        return flag_pass 
    #
    
    
    
    def LThread2SimRes(self, lt, bblist = ['N~N']):
        # CLE    A/1/P  A/15/P  18.7  7.0  -0.2
        # CLE    A/1/N  A/15/N   9.0  7.0  -0.2
        if not self.checkBBlist(bblist):
            print "... list cannot be used"
            sys.exit(1)
        #
        for tr in lt.thread:
            v = tr.ij_ndx
            i = v[0]; j = v[1]
            ctp = tr.ctp    #  B, I, M, etc.
            btp = tr.btp    #  s, sp, sa, etc.
            dG  = tr.dGij_B #  the FE of the specific bond at ij
            
            for bb in bblist:
                bbsplit = bb.split('~')
                aA = bbsplit[0]; aB = bbsplit[1]
                # print aA, aB
                self.simres += [SimRes(i, 'A', aA, j, 'A', aB, resdist[bb], res_wt[bb])]
            #
        #
    #
    
    def ChPair2SimRes(self, chpair, bblist = ['N~N']):
        # CLE    A/1/P  A/15/P  18.7  7.0  -0.2
        # CLE    A/1/N  A/15/N   9.0  7.0  -0.2
        if not self.checkBBlist(bblist):
            print "... list cannot be used"
            sys.exit(1)
        #
        for chp in chpair.data:
            i = chp.i; j = chp.j
            for bb in bblist:
                bbsplit = bb.split('~')
                aA = bbsplit[0]; aB = bbsplit[1]
                # print aA, aB
                self.simres += [SimRes(i, 'A', aA, j, 'A', aB, resdist[bb], res_wt[bb])]
            #
        #
    #
    
    
    
    def print_SimRNArestraints(self, flnm, resType = "slope", flag_save = False):
        # CLE    A/1/P  A/15/P  18.7  7.0  -0.2
        # CLE    A/1/N  A/15/N   9.0  7.0  -0.2
        s = ''
        flag_include_comments = False
        if flag_include_comments:
            if self.resType == "CLE":
                s += "# type  res1      res2     r_min     Kuhn    weight\n"   
            else:
                s += "# type  res1      res2      bgn       end    weight\n"
            #
        #
        
        for srk in self.simres:
            s += srk.disp_SimRes(resType)
        #
        
        if flag_save:
            try:
                fpx = open(flnm, 'w')
            except IOError:
                print "ERROR: cannot open %s" % flnm
                sys.exit(1)
            fpx.write(s)
            fpx.close()
        else:
            print s
            
        #
        return s
    #
#
# ###########################
# ###  Main: for testing  ###
# ###########################

def main(cl):
    print cl
    if len(cl) > 1:
        iflnm = cl[1]
        chdt = ChPairData()
        chdt.read_ChPairFile(iflnm)
        print chdt.disp_ChPairData(chdt.data)
        srdt = SimRNAData()
        srdt.ChPair2SimRes(chdt, ['N~N'])
        srdt.print_SimRNArestraints("test.res") 
        print "finished the test:"
        print "Read a chpair file and convert it to SimRNA restraints"
    else:
        usage()
    #
    
#

# Main
if __name__ == '__main__':
    main(sys.argv)
#
