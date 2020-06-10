#!/usr/bin/env python3

"""@@@

Main Module:   ChPair.py

Classes:       ChPair
               ChPairData

Functions:     LThread2ChPair 
               shuffle_ChPairData
               main

Author:        Wayne Dawson
creation date: parts 2016, made into a separate object 170314
last update:   2002010 (some minor format adjustments, upgrade to python3)
version:       0

Purpose:

180709:

   There seem to a variety of data formats that have been created in
   this Chreval package. There is Motif (the most complete), then
   there is LThread/LNode, which is hopelessly less complete and only
   marginally different from Pair. Here I have this ChPair.

   I think this data format was originally developed because Pair or
   LNode were too highly entwined with Calculate to use them for
   writing code to generate 3D structure building information.  The
   remaining utility of this class is for building shuffled data of
   the pairing interactions. The shuffling is used to validate the 3D
   cryo-EM like structure fitting programs. Since this file contains
   only the pairing information, it is fairly easy to to shuffle heat
   maps without changing other information.


"""

from FileTools import FileTools
from FileTools import getHeadExt
import sys
import random
import os

PROGRAM = "ChPair.py"


def LThread2ChPair(lt, title = "undefined"):
    
    CPData = ChPairData()
    CPData.cpdtitle = title
    CPData.sqlen = lt.sqlen
    CPData.set_sqlen = True
    CPData.dG  = lt.dG
    CPData.TdS = lt.TdS
    CPData.p   = lt.p
    for tr in lt.thread:
        v = tr.ij_ndx
        i = v[0]; j = v[1]
        ctp = tr.ctp    #  B, I, M, etc.
        btp = tr.btp    #  s, sp, sa, etc.
        dG  = tr.dGij_B #  the FE of the specific bond at ij
        if not (btp == 'bgn' or btp == 'end'):
            CPData.data += [ChPair(i,j,ctp,btp,dG)]
        #
        
    #|endfor
    
    return CPData
#

def shuffle_ChPairData(CPData):
    debug_shuffle_ChPairData = False
        
    nCPData = ChPairData()
    nCPData.cpdtitle = CPData.cpdtitle
    nCPData.sqlen   = CPData.sqlen
    nCPData.set_sqlen = True
    nCPData.dG  = CPData.dG
    nCPData.TdS = CPData.TdS
    nCPData.p   = CPData.p
    
    newdata = []
    used = {}
    i = -1; j = -1
    for dt in CPData.data:
        ib = dt.i; jb = dt.j; ctp = dt.ctp; btp = dt.btp; dG = dt.dG
        v_ceil   = CPData.sqlen - jb
        v_floor  = -ib
        k = int(random.uniform(v_floor, v_ceil)) # (no  int(x + 0.5))
        i = ib + k; j = jb + k
        flag_match = True
        ll = 0
        while flag_match:
            if not ctp == 'W' and ((i in used) or (j in used)):
                if debug_shuffle_ChPairData:
                    print ("reassign: i,j: ", i, j)
                #
                
                k = int(random.uniform(v_floor, v_ceil)) # (no  int(x + 0.5))
                i = ib + k; j = jb + k
                
            else:
                if debug_shuffle_ChPairData:
                    print ("W: i,j: ", i, j)
                #
                
                flag_match = False
            #
            
        #|endwhile
        
        used.update({ i : True }); used.update({ j : True })
        if not ctp == 'W':
            ctp = 'B'; btp = '-'
        #
        
        newdata +=  [ChPair(i, j, ctp, btp, dG)]
    
    #|endfor
    
    nCPData.data = newdata
    return nCPData
#




class ChPair(object):
    def __init__(self, i, j, ctp, btp, dG):
        self.i   = i
        self.j   = j
        self.ctp = ctp
        self.btp = btp
        self.dG  = dG
    #
    
    def disp_ChPair(self):
        s = "%5d %5d  %2s  %5s  %12.3f\n" % (self.i, self.j, self.ctp, self.btp, self.dG)
        return s
    #
    
    def __str__(self):
        s = "ij(%3d,%3d)[ctp=%2s][btp=%5s][%dG=%12.2f]" \
            % (self.i, self.j, self.ctp, self.btp, self.dG)
        return s
    #

    def __repr__(self):
        return self.__str__()
    #
        
#


class ChPairData(object):
    def __init__(self):
        self.data  = []
        self.sqlen = -1
        self.set_sqlen = False
        self.dG    = 1e10
        self.TdS   = 1e10
        self.p     = 0.0
        self.cpdtitle = ''
    #
    
    def reset_ChPairData(self):
        self.data  = []
        self.sqlen = -1
        self.set_sqlen = False
        self.dG    = 1e10
        self.TdS   = 1e10
        self.p     = 0.0
        self.cpdtitle = ''
    #
    
    def put_sqlen(self, n):
        print ("got to here")
        if isinstance( n, int ):
            self.sqlen = n
            self.set_sqlen = True
        else:
            print ("ERROR: input sequence length type is not a POSITIVE integer")
            print ("       variable: '%s'" % n)
            sys.exit(1)
        #
        
        if n <= 0:
            print ("ERROR: input sequence length type is not a POSITIVE integer")
            print ("       variable: '%s'" % n)
            sys.exit(1)
        #
        
    #
    
    
    
    # ###########
    # read chpair file
    #------------------------------------
    def read_ChPairFile(self, iflnm):
        debug_read_ChPairFile = False
        #self.reset_ChPairData()
        if debug_read_ChPairFile:
            print ("input filename: ", iflnm)
        #
        
        if not os.path.isfile(iflnm):
            print ("ERROR: %s not found" % iflnm)
            sys.exit(1)
        #
        
        fp = open(iflnm, 'r')
        lfp = fp.readlines()
        fp.close()
        
        jmax = -1
        self.data = [] # reset
        for k in range(0, len(lfp)):
            s = lfp[k].strip().split()
            if len(s) > 0:
                if s[0][0] == '#':
                    continue
                #
                
                if len(s) == 5:
                    try:
                        i   = int(s[0])
                        j   = int(s[1])
                        ctp = s[2]
                        btp = s[3]
                        dG  = float(s[4])
                    except ValueError:
                        print ("ERROR: read in data is not recognizable")
                        sys.exit(0)
                    #
                    
                    if debug_read_ChPairFile:
                        print (s)
                    #
                    
                    self.data += [ChPair(i, j, ctp, btp, dG)]
                    
                    if j > jmax:
                        jmax = j
                    #
                    
                else:
                    print ("ERROR: input file '%s' has too many entries. " % iflnm)
                    print ("       I don't understand it." )
                    print ("       line(%d): " % k, s)
                    sys.exit(1)
                #
                
            #
            
        #|endfor
        
        self.info_in_ChPair(lfp, jmax, iflnm)
        
    #
    
    def info_in_ChPair(self, lfp, jmax, iflnm):
        flag_set_N = False
        if self.set_sqlen:
            flag_set_N = True
        #
        
        k = 0
        ll = len(lfp)-1
        while not flag_set_N and k < ll:
            if lfp[k][0] == '#':
                ss = lfp[k].strip().split()
                if ss[1] == 'N':
                    if  not self.set_sqlen:
                        try:
                            self.sqlen = int(ss[2])
                        except ValueError:
                            print ("ERROR: stored sequence length is not an integer")
                            print ("       N = '%s'??" % ss[2])
                            sys.exit(1)
                        except IndexError:
                            print ("ERROR: missing sequence length data in file %s!" % iflnm)
                            sys.exit(1)
                        #
                        
                    #
                    
                    flag_set_N = True
                    
                elif ss[1] == 'dG':
                    try:
                        self.dG = float(ss[2])
                    except ValueError:
                        print ("ERROR: recorded total free energy is not a float")
                        print ("       dG = '%s'??" % ss[2])
                        sys.exit(1)
                    except IndexError:
                        print ("ERROR: missing free energy data in file %s!" % iflnm)
                        sys.exit(1)
                    #
                    
                elif ss[1] == 'TdS':
                    try:
                        self.TdS = float(ss[2])
                    except ValueError:
                        print ("ERROR: recorded total TdS weight is not a float")
                        print ("       TdS = '%s'??" % ss[2])
                        sys.exit(1)
                    except IndexError:
                        print ("ERROR: missing TdS weight data in file %s!" % iflnm)
                        sys.exit(1)
                    #
                    
                elif ss[1] == 'p':
                    try:
                        self.p = float(ss[2])
                    except ValueError:
                        print ("ERROR: recorded Boltzmann probability is not a float")
                        print ("       p = '%s'??" % ss[2])
                        sys.exit(1)
                    except IndexError:
                        print ("ERROR: missing Boltzmann probability in file %s!" % iflnm)
                        sys.exit(1)
                    #
                    
                elif k == 0:
                    if len(ss) > 1:
                        self.cpdtitle = ss[1]
                    else:
                        self.cpdtitle = "None"
                    #
                    
                #
                
            #
            
            k += 1
        #|endwhile
        
        if not flag_set_N:
                
            print ("cannot find the length of the sequence in %s" % iflnm)
            flhd, ext = getHeadExt(iflnm)
            
            print ("flhd = ", flhd)
            altflnm = flhd + ".DBN"
            print (altflnm)
            if os.path.exists(altflnm):
                ffp = open(altflnm, 'r')
                lffp = ffp.readlines()
                ffp.close()
                print (len(lffp[1]))
                self.sqlen = len(lffp[1])
                print ("found length: %d" % N)
            else:
                print ("WARNING: could not find any length information with ")
                print ("         this simple algorithm. Using the maximum ")
                print ("         from the list." )
                self.sqlen = jmax
            #
            
        #
        
        print ("title: ", self.cpdtitle)
        # sys.exit(0)
        return self.sqlen
    #
    
    
    """@
    
    It is a little confusing, but disp_ChPairData() displays a
    specific ChPair data set independent of where it comes from. On
    the other hand, print_ChPairData only prints data that is made
    within this object. There are reasons for that, though I am not
    sure that creating this confusion warrents the excuses for this
    construction. Presently, it is not used enough to concern me.

    """
    def disp_ChPairData(self, chdata):
        s = ''
        for d in chdata:
            s += d.disp_ChPair()
        #
        return s
    #
    
    def print_ChPairData(self, flnm = "Not_to_be_saved_to_a_file"):
        #
        s  = '# %s\n'            % self.cpdtitle
        s += '# N      %6d\n'    % self.sqlen
        s += '# dG   %8.2f\n'    % self.dG
        s += '# TdS  %8.2f\n'    % self.TdS
        s += '# p       %8.5f\n' % self.p
        s += '#  i     j    ctp   btp       dG\n'
        for dtk in self.data:
            s += dtk.disp_ChPair()
        #|endfor
        
        if not flnm == "Not_to_be_saved_to_a_file":
            try:
                fp = open(flnm, 'w')
            except IOError:
                print ("ERROR: cannot open %s" % flnm)
                sys.exit(1)
            fp.write(s)
            fp.close()
        else:
            print (s)
            
        #
        
        return s
    #
#


# ###########################
# ###  Main: for testing  ###
# ###########################

def main(cl):
    print (cl)
    if len(cl) > 1:
        iflnm = cl[1]
        dt = ChPairData()
        dt.read_ChPairFile(iflnm)
        print (dt.disp_ChPairData(dt.data))
    else:
        print ("nothing to do")
    #
    
#

# Main
if __name__ == '__main__':
    main(sys.argv)
#
