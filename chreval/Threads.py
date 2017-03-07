#!/usr/bin/env python

# package name: Threads
# creater:      Wayne K Dawson
# creation date: 170106
# last update:   170126


# This module is used primarily by "chreval.py". However, as with many
# packages that are building up, the package has become complicated
# enough, that the objects in this module are used by more than one program.

# Therefore, since this is part of a common set of objects, it is
# better to share the code among all the components that need
# it. Therefore, I have separated this package from chreval to allow
# independent testing and so that can be used by other programs for
# various things (for example, building heatmaps from input data).

# Moreover, I would rather work with it as an independent module so
# that I can test various things without having to run the whole
# program.

# Threads is a pair type representation intermediate between the Motif
# and Vienna representations for pairing; where Vienna is the most
# compact, but the least informative and Motif is the most detailed
# but not particularly compact.


import sys
import string
from MatrixTools import MatrixTools

# contains free energy parameters
import FreeEnergy

# for Vienna object representation
from Vienna import Vstruct

# for Motif object representation
from Motif import Motif # a core object of these building programs
from Motif import Link
from Motif import LGroup
from Motif import Map # main map for all the FE and structures


# ################################################################
# ############  Local-global functions and constants  ############
# ################################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# 1. control parameters tests etc in the program

# Test   function
#  0     a _very basic_ test of this module
#  1     parser for simple structures
#  2     insert into a list
#  3     convert a sequence to Motif

TEST = 3


# 2. special local global function

PROGRAM      = "Threads.py"  # name of the program

# This is not used so much now that I have introduced GetOpts.py,
# but it may still be useful at the very beginning of the program and
# possibly in other parts.

def usage():
    print "USAGE: %s" % PROGRAM
#

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ################################################################
# ################################################################



# ################################################################
# ################################################################
# #################################################################
# ###############  General configuration CONSTANTS  ###############
# ###############     settings used in Threads     ################
# #################################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# A hardwired infinity. Somewhat hate this, but it seems there is no
# way around this darn thing.
from Constants import INFINITY

# Pseudoknot and parallel stem 1D notation operators labels.  This
# notation at least works with VARNA
from Constants import num2lpr
from Constants import num2rpr
from Constants import lpr2num
from Constants import rpr2num

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ################################################################
# ################################################################




class LSegment:
    # this is basically a structure
    def __init__(self, ij_ndx, l_ndx, dGij_B, ctp, btp):
        self.ij_ndx = ij_ndx # (i, j)
        self.l_ndx  = l_ndx  # 0, 1, 2 etc.
        self.dGij_B = dGij_B # the FE of the specific bond at ij
        self.ctp    = ctp    # connection type
        self.btp    = btp    # bond type
        #
    #
    
    def get_ij(self):
        return self.ij_ndx[0], self.ij_ndx[1]
    #
    
    def disp_lseg(self):
        s = " (%3d, %3d)[%2d][%s-%5s] dGij_B = %8.3f" \
            % (self.ij_ndx[0], self.ij_ndx[1], self.l_ndx, self.ctp, string.ljust(self.btp, 5), self.dGij_B)
        return s
    #
#

class LThread:
    def __init__(self, N):
        self.sqlen  = N
        self.thread = []
        self.dG     = INFINITY
        self.TdS    =   0.0
        self.p      = -99.99
    #
    
    def add_lsegment(self, ij_ndx, l_ndx, Vij_B, ctp = 'B', btp = 's'):
        self.thread +=  [LSegment(ij_ndx, l_ndx, Vij_B, ctp, btp)]
        self.compute_dG()
    #
    
    def compute_dG(self):
        dG = 0.0
        for tk in self.thread:
            dG += tk.dGij_B
        self.dG = dG
    #
    
#


class DispThreads:
    
    def __init__(self, N):
        self.N      = N
        # calc.fe.mtools would also work, but I want this module to be
        # a little more independent and it doesn't cost so much to
        # initialize a separate object of class MatrixTools.
        self.mtools = MatrixTools()
        
        # options for employing SimRNA
        self.add_P    = False  # include P to P in restraint file
        self.add_N    = True   # include N to N in restraint file
        self.xi       = 7.0    # [nt] Kuhn length
        self.NNdist   = 9.0    # [A] bp N to N distance
        self.PPdist   = 18.7   # [A] bp P to P distance
        self.Nweight  = -0.2   # the N weight for SimRNA restraint
        self.Pweight  = -0.2   # the P weight for SimRNA restraint
        self.restType = "slope" # option for "slope" or "CLE"
    #
    
    # display a given thread
    def disp_LThread(self, lt):
        lt.compute_dG()
        print "total free energy: %8.3f" % lt.dG
        for tr in lt.thread:
            print tr.disp_lseg()
        #
    #
    
    def set_strand_direction(self, ctp, btp, counter):
        btpp = list(btp)
        # print btpp
        i_lab = '('; j_lab = ')'
        if ctp == 'S':
            if btpp[0] == 's':
                if btpp[1] == 'a':
                    i_lab = '('
                    j_lab = ')'
                else:
                    i_lab = num2lpr[counter]
                    j_lab = num2rpr[counter]
                    counter += 1
                #
            else:
                if btp == 'c':
                    # 
                    i_lab = '('
                    j_lab = ')'
                elif btp == 't' or btp == 'r' or btp == 'l':
                    # 't', 'r', 'l'
                    i_lab = num2lpr[counter]
                    j_lab = num2rpr[counter]
                    counter += 1
                else:
                    print "problems in makeLThreadDotBracket_1b():"
                    sys.exit(1)
                #
            #
        else:
            # probably a PK
            i_lab = num2lpr[counter]
            j_lab = num2rpr[counter]
            counter += 1
            
        #
        return i_lab, j_lab, counter
    #
    def makeBlank(self):
        bbb = []
        for i in range(0, self.N):
            bbb += ['.']
        return bbb

    
    def makeLThreadDotBracket_1b(self, lt, structure_layout = 0):
        flag_debug = False # True
        seqv  = []
        sssv  = []
        ctcfv = []
        pkv   = []
        npks  = -1
        k_pointer = 0
        nws   = -1
        w_pointer = 0
        counter = 3
        
        for i in range(0, self.N):
            seqv += ['c']
        #
        
        for i in range(0, self.N):
            sssv   += ['.']
            ctcfv  += ['.']
        #
        
        for tr in lt.thread:
            v = tr.ij_ndx
            ctp = tr.ctp
            btp = tr.btp
            if flag_debug:
                print v, ctp, btp, tr.dGij_B        
            #
            
            if not counter < len(num2lpr): 
                # this is _not_ good, and it can lead to major issues
                # for long sequences
                counter = 3
            #
            
            # remove references to PKs
            if ctp == 'K' and btp == 'bgn':
                npks += 1
                pkv += [self.makeBlank()] # add another linkage channel
                #print "begin PK(%d)" % npks
                continue
            #
            if ctp == 'K' and (btp == 'end' or (npks > k_pointer and btp == 'l')):
                #print "end of PK(%d)" % npks
                #print "k_pointer = %d, npks = %d, btp = %s" % (k_pointer, npks, btp)
                k_pointer += 1
                continue
            #
            
            if ctp == 'W' and (btp == 'bgn' or btp == 'end'):
                if btp == 'end':
                    ctcfv[v[0]] = '{'
                    ctcfv[v[1]] = '}'
                    sssv[v[0]] = '{'
                    sssv[v[1]] = '}'
                #
                continue
            #
            
            if ctp == 'S' and (btp == 'bgn' or btp == 'end'):
                continue
            #
            
            
            if ctp == 'l':
                #print v
                if btp == 'sp':
                    i_lab, j_lab, counter = self.set_strand_direction(ctp, btp, counter)
                    sssv[v[0]] = i_lab
                    sssv[v[1]] = j_lab
                else:
                    sssv[v[0]] = '['
                    sssv[v[1]] = ']'
                #
                #print npks, k_pointer
                if npks > k_pointer:
                    pkv[k_pointer][v[0]] = '('
                    pkv[k_pointer][v[1]] = ')'
                    #print pkv[k_pointer]
                else:
                    pkv[npks][v[0]] = '('
                    pkv[npks][v[1]] = ')'
                    #print pkv[npks]
                #
            elif ctp == 'W': # island (wyspa)
                # print v
                sssv[v[0]] = '|'
                sssv[v[1]] = '|'
                ctcfv[v[0]] = '|'
                ctcfv[v[1]] = '|'
            elif ctp == 'S':
                if btp == 'sp':
                    sssv[v[0]] = num2lpr[counter]
                    sssv[v[1]] = num2rpr[counter]
                    counter += 1
                else:
                    sssv[v[0]] = '('
                    sssv[v[1]] = ')'
                #
            else:
                # print v
                sssv[v[0]] = '('
                sssv[v[1]] = ')'
            #
            
            if btp == 's' or btp == 'sa' or btp == 'sp':
                seqv[v[0]] = 'x'
                seqv[v[1]] = 'y'
            elif btp == 'c':
                seqv[v[0]] = 'W'
                seqv[v[1]] = 'Z'
                sssv[v[0]] = '{'
                sssv[v[1]] = '}'
            elif btp == 't' or btp == 'r' or btp == 'l':
                seqv[v[0]] = 'w'
                seqv[v[1]] = 'z'
                sssv[v[0]] = '{'
                sssv[v[1]] = '}'
            elif btp == 'wyspa':
                seqv[v[0]] = 'I'
                seqv[v[1]] = 'I'
            elif btp == '-':
                if ctp == 'P' or ctp == 'J':
                    # 'P' and 'J' are real possibilities now.
                    seqv[v[0]] = 'c'
                    seqv[v[1]] = 'c'
                else:
                    # It currently has to go here!!!!
                    print "ERROR: undefined btp = %s" % btp
                    print "       connection type = %s" % ctp
                    print "       ij = %d,%d        " % (v[0], v[1])
                    sys.exit(1)
                #
            else:
                if not (btp == 'bgn' or btp == 'end'):
                    print "ERROR: undefined btp = %s" % btp
                    print "       connection type = %s" % ctp
                    print "       ij = %d,%d        " % (v[0], v[1])
                    sys.exit(1)
                #
            #
            
                    
        #
        
        
        seq = string.join(seqv, '') 
        sss = string.join(sssv, '') 
        # also: python> ''.join(sssv); python> ''.join(seqv) 
        
        s = ''
        if structure_layout == 0:
            s = sss # this will only return the ss string
        elif structure_layout == 1:
            s = seq + '\n' + sss + '\n'
        #
        else:
            s  = seq + '\n'
            s += sss + '\n'
            # print pkv
            if len(pkv) > 0:
                for pkvk in pkv:
                    s += ''.join(pkvk) + '\n'
            s += ''.join(ctcfv) + '\n'

            if flag_debug:
                print "target exit"
                print s
                sys.exit(0)
            #
        #
        return s
    #
    
    def makeLThreadDotBracket_VARNA(self, lt, structure_layout = 0):
        seqv  = []
        sssv  = []
        counter = 3
        
        for i in range(0, self.N):
            seqv += ['c']
            sssv += ['.']
        #
        
        
        # print s
        for tr in lt.thread:
            v = tr.ij_ndx
            ctp = tr.ctp
            btp = tr.btp
            
            # remove references to PKs
            if ctp == 'K' and (btp == 'bgn' or btp == 'end'):
                continue
            #
            
            if ctp == 'W' and (btp == 'bgn' or btp == 'end'):
                if btp == 'end':
                    sssv[v[0]] = '{'
                    sssv[v[1]] = '}'
                #
                continue
            #
            
            if ctp == 'S' and (btp == 'bgn' or btp == 'end'):
                continue
            #

            if not counter < len(num2lpr):
                # this is _not_ good, and it can lead to major issues
                # for long sequences
                counter = 3
            
            
            if ctp == 'l':
                # print v
                if btp == 'sp':
                    i_lab, j_lab, counter = self.set_strand_direction(ctp, btp, counter)
                    sssv[v[0]] = i_lab
                    sssv[v[1]] = j_lab
                else:
                    sssv[v[0]] = '['
                    sssv[v[1]] = ']'
                #
            elif ctp == 'W': # island (wyspa)
                # print v
                sssv[v[0]] = '|'
                sssv[v[1]] = '|'
            elif ctp == 'S':
                # print v
                if btp == 'sp':
                    sssv[v[0]] = num2lpr[counter]
                    sssv[v[1]] = num2rpr[counter]
                    counter += 1
                else:
                    sssv[v[0]] = '('
                    sssv[v[1]] = ')'
                #
            else:
                # print v
                sssv[v[0]] = '('
                sssv[v[1]] = ')'
            #
            
            if btp == 's' or btp == 'sa' or btp == 'sp':
                seqv[v[0]] = 'x'
                seqv[v[1]] = 'y'
            elif btp == 'c':
                sssv[v[0]] = '{'
                sssv[v[1]] = '}'
                seqv[v[0]] = 'W'
                seqv[v[1]] = 'Z'
            elif btp == 't' or btp == 'r' or btp == 'l':
                sssv[v[0]] = '{'
                sssv[v[1]] = '}'
                seqv[v[0]] = 'w'
                seqv[v[1]] = 'z'
            elif btp == 'wyspa':
                seqv[v[0]] = 'I'
                seqv[v[1]] = 'I'
            elif btp == '-':
                if ctp == 'P' or ctp == 'J':
                    # 'P' and 'J' are real possibilities now.
                    seqv[v[0]] = 'c'
                    seqv[v[1]] = 'c'
                else:
                    print "ERROR: undefined btp = %s" % btp
                    print "       connection type = %s" % ctp
                    print "       ij = %d,%d        " % (v[0], v[1])
                    sys.exit(1)
                #
            else:
                if not (btp == 'bgn' or btp == 'end'):
                    print "ERROR: undefined btp = %s" % btp
                    print "       connection type = %s" % ctp
                    print "       ij = %d,%d        " % (v[0], v[1])
                    sys.exit(1)
                #
            #
            
                    
        #
        
        
        seq = string.join(seqv, '') 
        sss = string.join(sssv, '') 
        # also: python> ''.join(sssv); python> ''.join(seqv) 
        
        s = ''
        if structure_layout == 0:
            s = sss # this will only return the ss string
        elif structure_layout == 1:
            s = seq + '\n' + sss + '\n'
        #
        else:
            s  = seq + '\n'
            s += sss + '\n'
        return s
    #
    
    def makeLThreadHeatMap(self, lt, flag_no_header = False):
        
        hmap = []
        hmap = self.mtools.initialize_matrix(hmap, self.N, 0)
        
        
        # print s
        for tr in lt.thread:
            v    = tr.ij_ndx
            ctp = tr.ctp
            btp = tr.btp
            # remove references to PKs
            if ctp == 'K' and (btp == 'bgn' or btp == 'end'):
                continue
            #
            if ctp == 'W' and (btp == 'bgn' or btp == 'end'):
                continue
            #
            
            if ctp == 'S' and (btp == 'bgn' or btp == 'end'):
                continue
            #
            i = v[0]; j = v[1]
            hmap[i][j] = 1
            hmap[j][i] = 1
        #
        return self.mtools.make_heatmap(hmap, flag_no_header)
    #
    
    def printLThreadHeatMap(self, flnm, lt, flag_no_header = False):
        s = self.makeLThreadHeatMap(lt, flag_no_header)
        try:
            fp = open(flnm, 'w')
        except IOError:
            print "ERROR: cannot open %s" % flnm
            sys.exit(1)
        fp.write(s)
        fp.close()
    #
    
    def print_ijPairList(self, flnm, lt, flag_save = False):
        #
        ff = flnm.split('.')
        ext = ff[len(ff)-1]
        title = flnm[0:len(flnm) - len(ext) - 1]
        s  = '# %s\n'            % title
        s += '# dG   %8.2f\n'    % lt.dG
        s += '# TdS  %8.2f\n'    % lt.TdS
        s += '# p       %8.5f\n' % lt.p
        s += '#  i     j    ctp   btp       dG\n'
        for tr in lt.thread:
            v = tr.ij_ndx
            i = v[0]; j = v[1]
            ctp = tr.ctp    #  B, I, M, etc.
            btp = tr.btp    #  s, sp, sa, etc.
            ndx = tr.l_ndx  #  0, 1, 2 etc.
            dG  = tr.dGij_B #  the FE of the specific bond at ij
            if not (btp == 'bgn' or btp == 'end'): 
                s += '%4d  %4d  %5s %5s %8.2f\n' % (i, j, ctp, btp, dG) 
        #
        if flag_save:
            try:
                fp = open(flnm, 'w')
            except IOError:
                print "ERROR: cannot open %s" % flnm
                sys.exit(1)
            fp.write(s)
            fp.close()
        else:
            print s
            
        #
    #
    
    def print_SimRNArestraints(self, flnm, lt, restType = "slope", flag_save = False):
        # CLE    A/1/P  A/15/P  18.7  7.0  -0.2
        # CLE    A/1/N  A/15/N   9.0  7.0  -0.2
        s = ''
        flag_include_comments = False
        if flag_include_comments:
            if self.restType == "CLE":
                s += "# type  res1      res2     r_min     Kuhn    weight\n"   
            else:
                s += "# type  res1      res2      bgn       end    weight\n"
            #
        #
        
        for tr in lt.thread:
            v = tr.ij_ndx
            i = v[0]; j = v[1]
            ctp = tr.ctp    #  B, I, M, etc.
            btp = tr.btp    #  s, sp, sa, etc.
            ndx = tr.l_ndx  #  0, 1, 2 etc.
            dG  = tr.dGij_B #  the FE of the specific bond at ij
            if self.restType == "CLE":
                if self.add_N:
                    s += "CLE     %s/%d/N   %s/%d/N %6.2f  %6.2f  %8.3f\n" % \
                        ('A', i + 1, 'A', j + 1, self.NNdist, self.xi, self.Nweight)
                #
                if self.add_P:
                    s += "CLE     %s/%d/P   %s/%d/P %6.2f  %6.2f  %8.3f\n" % \
                        ('A', i + 1, 'A', j + 1, self.PPdist, self.xi, self.Pweight)
                #
            elif self.restType == "slope":
                if self.add_N:
                    s += "SLOPE    %s/%d/N   %s/%d/N %6.2f  %6.2f  %8.3f\n" % \
                        ('A', i + 1, 'A', j + 1, \
                         (self.NNdist - 0.5), (self.NNdist + 0.5), -self.Nweight)
                #
                if self.add_P:
                    s += "SLOPE    %s/%d/P   %s/%d/P %6.2f  %6.2f  %8.3f\n" % \
                        ('A', i + 1, 'A', j + 1, \
                         (self.PPdist - 0.5), (self.PPdist + 0.5), -self.Pweight)
                #
            #
        #
        if flag_save:
            try:
                fp = open(flnm, 'w')
            except IOError:
                print "ERROR: cannot open %s" % flnm
                sys.exit(1)
            fp.write(s)
            fp.close()
        else:
            print s
            
        #
    #
    
    def printLThreadDotBracket_VARNA(self, flnm, lt):
        ff = flnm.split('.')
        if not lt.p == -99.99:
            s = "> %s    %8.3f   %10.8f\n" % (ff[0], lt.dG, lt.p)
        else:
            s = "> %s    %8.3f\n" % (ff[0], lt.dG)
        s += self.makeLThreadDotBracket_VARNA(lt, 1)
        try:
            fp = open(flnm, 'w')
        except IOError:
            print "ERROR: cannot open %s" % flnm
            sys.exit(1)
        fp.write(s)
        fp.close()
    #
#



class Vienna2LThread:
    def __init__(self, vs):
        self.N        = vs.N        # sequence length
        self.vstr     = vs.vstr     # the 1D structure
        self.BPlist   = vs.BPlist   # the decomposed 1D ss structure
        self.PKlist   = vs.PKlist   # the decomposed 1D apPK structure
        self.CTCFlist = vs.CTCFlist # the decomposed 1D CTCF structure
        self.vs       = vs          # for employing functions
        self.lt = LThread(self.N)
        self.fe = FreeEnergy.FreeEnergy(FreeEnergy.FESetUp())
        self.initialized_apss = False
        self.initialized_apPK = False
        self.initialized_CTCF = False
        self.apPK_r = {}
        self.apPK_l = {}
        self.ppPK_r = {}
        self.ppPK_l = {}
        self.initialize_v2t   = True
    #

    def vienna2thread(self):
        if not self.initialize_v2t:
            print "ERROR: Vienna2LThread() is not initialized!"
            sys.exit(1)
        #
        debug = False
        self.scan_for_ss(0, self.N, 0)
        if debug:
            self.disp_lt()
        #
        self.scan_for_PK(0, self.N, 0)
        if debug:
            self.disp_lt()
        #
        self.scan_for_CTCFs()
        if debug:
            self.disp_lt()
        #
    #
    
    def initialize(self, Xlist):
        debug = False
        for k in range(0, len(Xlist)):
            i, j = Xlist[k].get_Pair()
            # presently, the weight is 5.0 because I don't have a way
            # to enter this info
            dG = self.fe.calc_dG(i, j, 5.0, self.fe.T)
            self.lt.thread += [LSegment((i,j), 0, dG, '-', 's')] 
        #
        if debug:
            print "xxx"
            print self.vstr
            for xk in Xlist:
                print " %4d %4d " % (xk.i, xk.j)
            #
        #
    #
    
    def scan_for_ss(self, ib, jb, ndx):
        debug = False
        if debug:
            print "ndx: ", ndx, len(self.BPlist)
        #
        ndx_pntr = ndx
        if ndx == 0:
            # initialize LThread with items from BPlist
            self.initialize(self.BPlist)
            self.initialized_apss = True
        #
        
        if ndx >= len(self.BPlist):
            # Regardless of what list we have created, when the
            # condition ndx >= len(self.BPlist)!!! occurs, we have run
            # out of items in the list. Therefore search should
            # (indeed must) terminate here.
            return (-1, -1), ndx_pntr-1
        #
        i, j = self.BPlist[ndx].get_Pair()
        if i >= jb:
            # we required an upper boundary of jb, and we have an
            # ordered list up the heirarchy, therefore, i >= jb means
            # the search terminates
            return (-1, -1), ndx_pntr-1
        else:
            # at least one more potential item is present that we must
            # look at.
            v, ndx_pntr = self.scan_for_ss(i,j, ndx+1)
            if v[0] == -1:
                # Recall that this is an ordered list derived from the
                # secondary structure as indicated in a Vienna
                # formatted 1D string. In this case, our level of
                # recursion has reached a point where we can search no
                # further along the trajectory of a given antiparallel
                # stem. Based on the way Vienna derives these lists,
                # we deduce that we have reached the end of a stem
                # where a hairpin loop closes the stem (or we just
                # have a hairpin).
                self.lt.thread[ndx].ctp = 'B'
                if debug:
                    print ndx, "B"
                #
            else:
                # Here, we know we have a branch at (i,j). We must now
                # explore that branch until we either discover a
                # branch adjacent to it (in an M-loop), or we find
                # that we have a single branch inside of this first
                # branch, or we find a branch directly adjacent to
                # this one at (i+1,j-1) in the next step.
                
                # First find out if this is an M-loop.
                branches = [(i, j)]
                k = j+1
                while k < jb:
                    
                    if self.vstr[k] == '(':
                        v = self.find_this_ij(k, self.BPlist)
                        if v[0] == -1:
                            # shouldn't happen because we observe a '(', which
                            # means there should have been a ')'!
                            print "ERROR!"
                            sys.exit(1)
                        #
                        branches += [v]
                        k = v[1] + 1
                    else:
                        k += 1
                    #
                #
                if debug:
                    print branches
                #
                
                # now expand the additional branches if found in the
                # same way as the current one.
                for m in range(1, len(branches)):
                    
                    vv = branches[m]
                    if debug:
                        print m, vv
                    #
                    v, ndx_pntr = self.scan_for_ss(vv[0],vv[1], ndx_pntr+1)
                #
                
                # assuming we found more than one branch, we need to
                # figure out what connection type (ctp) should be
                # assigned to the base of this branch; whether the
                # branch we found should be designated as a stem or a
                # I-loop or a hairpin
                if len(branches) > 1:
                    if self.is_ap_stem(ndx, self.BPlist):
                        # This is only satisfied if the next item on
                        # the list references (i+1,j-1).
                        self.lt.thread[ndx].ctp = 'S'
                        self.lt.thread[ndx].btp = 'sa'
                        if debug:
                            print ndx, "S"
                        #
                    else:
                        # !!!hmm, I think there is a problem here!!!
                        # !!!('I' or 'B'). However, it seems to come
                        # !!!out ok anyway. Maybe the earlier step
                        # !!!overrides this.
                        self.lt.thread[ndx].ctp = 'I'
                        if debug:
                            print ndx, "I"
                        #
                    #
                    if ndx > 0:
                        # for LThread, we can only assign a 'P' to the
                        # very first structure on the list, assuming
                        # that we found branches there.
                        self.lt.thread[ndx-1].ctp = 'M'
                        ndx_pntr = ndx-1
                        if debug:
                            print "%4d %s, %4d  %s" % ((ndx-1), "M", ndx, "I")
                        #
                    else:
                        if debug:
                            # the output will look a bit strange
                            # because ndx-1 = -1!!!
                            print "%4d %s, %4d  %s" % ((ndx-1), "P", ndx, "I")
                        #
                    #
                elif len(branches) == 1:
                    # In this case, we found that we only have one
                    # branch that turned up.
                    if self.is_ap_stem(ndx, self.BPlist):
                        # This is only satisfied if the next item on
                        # the list references (i+1,j-1).
                        self.lt.thread[ndx].ctp = 'S'
                        self.lt.thread[ndx].btp = 'sa'
                        if debug:
                            print ndx, "S"
                        #
                    else:
                        # !!!hmm, I think there is a problem here!!!
                        # !!!('I' or 'B'). However, it seems to come
                        # !!!out ok anyway. Maybe the earlier step
                        # !!!overrides this.
                        self.lt.thread[ndx].ctp = 'I'
                        if debug:
                            print ndx, "I"
                        #
                    #
                #
            #
        #
        return (i,j), ndx_pntr
    #
    
    
    
    def scan_for_PK(self, ib, jb, ndx):
        debug = False
        
        apBPtails = self.find_apStem_tails(self.BPlist)
        if debug:
            print "apBPtails:  ", apBPtails
        #
        apPKtails = self.find_apStem_tails(self.PKlist)
        if debug:
            print "apPKtails:  ", apPKtails
        #
        apPKjoints = self.compress_PKlist(apPKtails)
        if debug:
            print "apPKjoints: ", apPKjoints
        #
        self.apPK_r = {}
        for pkj in apPKjoints:
            for bpt in apBPtails:
                i = bpt[0]; j = bpt[1]
                #i, j = bpt.get_Pair()
                if (pkj[0] < i and i < pkj[1] and pkj[1] < j) \
                or (i < pkj[0] and pkj[0] < j and j < pkj[1]):
                    #  ...[...(...]....)....
                    if self.apPK_r.has_key(pkj):
                        self.apPK_r[pkj] += [(i,j)]
                    else:
                        self.apPK_r.update({pkj : [(i,j)]})
                    #
                #
            #
        #
        
        self.apPK_l = {}
        for pkj in apPKjoints:
            for pkt in self.PKlist:
                if pkt.v == 'p':
                    continue
                else:
                    i_pk, j_pk = pkt.get_Pair()
                    if pkj[0] <= i_pk and j_pk <= pkj[1]:
                        if self.apPK_l.has_key(pkj):
                            self.apPK_l[pkj] += [(i_pk,j_pk)]
                        else:
                            self.apPK_l.update({pkj : [(i_pk,j_pk)]})
                        #
                    #
                #
            #
        #
        if debug:
            print "apPK_r: ", self.apPK_r
            print "apPK_l: ", self.apPK_l
        #
        
        ppPKjoints = self.find_ppStem_tails(self.PKlist)
        if debug:
            print ppPKjoints
        #
        self.ppPK_r = {}
        for pkj in ppPKjoints:
            for bpt in apBPtails:
                i = bpt[0]; j = bpt[1]
                #i, j = bpt.get_Pair()
                if (pkj[0] < i and i < pkj[1] and pkj[1] < j) \
                or (i < pkj[0] and pkj[0] < j and j < pkj[1]):
                    #  ...[...(...]....)....
                    if self.ppPK_r.has_key(pkj):
                        self.ppPK_r[pkj] += [(i,j)]
                    else:
                        self.ppPK_r.update({pkj : [(i,j)]})
                    #
                #
            #
        #
        
        self.ppPK_l = {}
        npkl = len(self.PKlist)
        debug4 = False
        for pkj in ppPKjoints:
            i_pkr = pkj[0]; j_pkr = pkj[1]
            if debug4:
                print "next ij_pkr: ", pkj
            #
            for k in range(0, npkl):
                if self.PKlist[k].v == 'a':
                    continue
                else:
                    i_kp0, j_kp0 = self.PKlist[k].get_Pair()
                    if debug4:
                        print "ppfit (%2d,%2d) <- (%2d,%2d)" % (i_pkr, j_pkr, i_kp0, j_kp0)
                    #
                    if i_pkr == i_kp0 and j_pkr == j_kp0:
                        if debug4:
                            print "assign ij_kp0: ", i_kp0, j_kp0
                        #
                        self.ppPK_l.update({pkj : [(i_kp0,j_kp0)]})
                        if k + 1 >= npkl:
                            continue
                        i_kp1, j_kp1 = self.PKlist[k+1].get_Pair()
                        if i_kp0 + 1 == i_kp1 and j_kp0 + 1 == j_kp1:
                            self.ppPK_l[pkj] += [(i_kp1,j_kp1)]
                    #
                    elif i_pkr < i_kp0 and i_kp0 < j_kp0 and j_pkr < j_kp0:
                        if k + 1 >= npkl:
                            continue
                        i_kp1, j_kp1 = self.PKlist[k+1].get_Pair()
                        if i_kp0 + 1 == i_kp1 and j_kp0 + 1 == j_kp1:
                            self.ppPK_l[pkj] += [(i_kp1,j_kp1)]
                        else:
                            break
                        #
                    #
                #
            #
        #
        if debug:
            print "ppPK_r: ", self.ppPK_r
            print "ppPK_l: ", self.ppPK_l
            print "Finished scanning the dataset"
        #
        
        # first we build the ap PKs.
        if debug:
            print "build the antiparallel PKs"
            print "apPK_r: ", self.apPK_r
            print "apPK_l: ", self.apPK_l
        #
        rkeys = self.apPK_r.keys()
        rkeys.sort(key=self.getKey_0) # set the order
        lkeys = self.apPK_l.keys()
        lkeys.sort(key=self.getKey_0) # set the order
        
        i_r = -1; j_r = -1
        i_l = -1; j_l = -1
        i_K = -1; j_K = -1
        for rk in rkeys:
            if not self.apPK_l.has_key(rk):
                print "error: no linkage for ", rk
                sys.exit(1)
            #
            i_l = rk[0]; j_l = rk[1]
            n = len(self.apPK_r[rk])
            if n > 1:
                i_r = self.apPK_r[rk][0][0]; j_r = self.apPK_r[rk][n-1][1]
                i_K = i_r;              j_K = j_r
                if i_l < i_r or j_l > j_r:
                    print "doesn't make sense!, n = ", n 
                    print self.apPK_r[rk]
                    print self.apPK_l[rk]
                    sys.exit(1)
                #
            else:
                i_r = self.apPK_r[rk][0][0]; j_r = self.apPK_r[rk][0][1]
                if i_l < i_r and j_l < j_r:
                    i_K = i_l;  j_K = j_r
                elif i_r < i_l and j_r < j_l:
                    i_K = i_r;  j_K = j_l
                else:
                    print "doesn't make sense!, n = ", n
                    print self.apPK_r[rk]
                    print self.apPK_l[rk]
                    sys.exit(1)
                #
            #
            if debug:
                print "ij_K: ", i_K, j_K
            #
            start = self.locate_start_insert(i_K, j_K, debug)
            self.lt.thread.insert(start, LSegment((i_K,j_K), 0, 0.0, 'K', 'bgn'))
            start = self.locate_end_insert(i_K, j_K, start, debug)
            for pkl in self.apPK_l[rk]:
                i_l = pkl[0]; j_l = pkl[1]
                # presently, the weight is 5.0 because I don't have a way
                # to enter this info into the set of data
                dG = self.fe.calc_dG(i_l, j_l, 5.0, self.fe.T)
                self.lt.thread.insert(start, LSegment((i_l,j_l), 0, dG, 'l', 'sa'))
                start += 1
            #
            self.lt.thread.insert(start, LSegment((i_K,j_K), 0, 0.0, 'K', 'end'))
        #
        # self.disp_lt()
        
        
        # next we build the pp PKs.
        if debug:
            print "build the parallel (pp) PKs"
            print "ppPK_r: ", self.ppPK_r
            print "ppPK_l: ", self.ppPK_l
        #
        rkeys = self.ppPK_r.keys()
        rkeys.sort(key=self.getKey_0)
        lkeys = self.ppPK_l.keys()
        lkeys.sort(key=self.getKey_0)
        
        # parallel stems
        for lk in lkeys:
            if not self.ppPK_r.has_key(lk):
                if debug:
                    print "has not ", lk
                #
                i_sp0 = lk[0]; j_sp0 = lk[1]
                start = self.locate_start_insert(i_sp0, j_sp0, debug)
                n = len(self.lt.thread)
                if debug:
                    print "start, n: ", start, n
                    if start < n:
                        print self.lt.thread[start].disp_lseg()
                    #
                #
                self.lt.thread.insert(start, LSegment((i_sp0,j_sp0), 0, 0.0, 'S', 'bgn'))
                start += 1
                for pstm in self.ppPK_l[lk]:
                    if debug:
                        print pstm
                    #
                    i_sp = pstm[0]; j_sp = pstm[1] 
                    dG = self.fe.calc_dG(i_sp, j_sp, 5.0, self.fe.T)
                    self.lt.thread.insert(start, LSegment((i_sp,j_sp), 0, dG, 'S', 'sp'))
                    start += 1
                #
                self.lt.thread.insert(start, LSegment((i_sp0,j_sp0), 0, 0.0, 'S', 'end'))
            #
        #
        
        
        # parallel chain interactions mixed with secondary structure.
        for rk in rkeys:
            if not self.ppPK_l.has_key(rk):
                print "error: no linkage for ", rk
                sys.exit(1)
            #
            i_l = rk[0]; j_l = rk[1]
            n = len(self.ppPK_r[rk])
            if n > 1:
                i_r = self.ppPK_r[rk][0][0]; j_r = self.ppPK_r[rk][n-1][1]
                i_K = i_r;              j_K = j_r
                if i_l < i_r or j_l > j_r:
                    print "doesn't make sense!, n = ", n 
                    print self.ppPK_r[rk]
                    print self.ppPK_l[rk]
                    sys.exit(1)
                #
            else:
                i_r = self.ppPK_r[rk][0][0]; j_r = self.ppPK_r[rk][0][1]
                if i_l < i_r and j_l < j_r:
                    i_K = i_l;  j_K = j_r
                elif i_r < i_l and j_r < j_l:
                    i_K = i_r;  j_K = j_l
                else:
                    print "doesn't make sense!, n = ", n
                    print self.ppPK_r[rk]
                    print self.ppPK_l[rk]
                    sys.exit(1)
                #
            #
            
            if debug:
                print i_K, j_K
            #
            
            start = self.locate_start_insert(i_K, j_K, debug)
            if debug:
                print "beginning position of pp pk: ", start
            #
            self.lt.thread.insert(start, LSegment((i_K,j_K), 0, 0.0, 'K', 'bgn'))
            start = self.locate_end_insert(i_K, j_K, start, debug)
            if debug:
                print "end of enclosed region of pp pk: ", start
            #
            
            for pkl in self.ppPK_l[rk]:
                i_l = pkl[0]; j_l = pkl[1]
                # presently, the weight is 5.0 because I don't have a way
                # to enter this info into the set of data
                dG = self.fe.calc_dG(i_l, j_l, 5.0, self.fe.T)
                self.lt.thread.insert(start, LSegment((i_l,j_l), 0, dG, 'l', 'sp'))
                start += 1
            #
            self.lt.thread.insert(start, LSegment((i_K,j_K), 0, 0.0, 'K', 'end'))
            
        #
        self.initialized_apPK = True
        if debug:
            self.disp_lt()
            #sys.exit(0)
        #
    #
    
    def getKey_0(self, item):
        return item[0]
    #
    
    def find_apStem_tails(self, Xlist):
        # test example
        #  0         10        20        30        40        50        60        70        80        90
        #  |         |         |         |         |         |         |         |         |         | 
        # "{((((((.A.AAAAA......))))).((((........a..aaaaa....))))..).|..([)]([)].ABC.abc.}....{.....}"
        tails = []
        k = 0
        flag_cont = True
        n = len(Xlist)
        while flag_cont:
            if k >= n:
                flag_cont = False
            else:
                if Xlist[k].v == 'p':
                    k += 1
                else:
                
                    i, j = Xlist[k].get_Pair()
                    k = self.find_apStem_head(k, Xlist)
                    tails += [(i,j)]
                    k +=1
                #
            #
        #
        return tails
    #
    
    def find_apStem_head(self, k, Xlist):
        flag_cont = True
        if k >= len(Xlist) - 1:
            flag_cont = False
        while flag_cont:
            i_0, j_0 = Xlist[k  ].get_Pair()
            i_1, j_1 = Xlist[k+1].get_Pair()
            if i_0 + 1 == i_1 and j_0 - 1 == j_1:
                k+=1
            else:
                flag_cont = False
            #
        #
        return k
    #
    
    def find_ppStem_tails(self, Xlist):
        # test example
        #  0         10        20        30        40        50        60        70        80        90
        #  |         |         |         |         |         |         |         |         |         | 
        # "{((((((.A.AAAAA......))))).((((........a..aaaaa....))))..).|..([)]([)].ABC.abc.}....{.....}"
        tails = []
        k = 0
        flag_cont = True
        n = len(Xlist)
        while flag_cont:
            if k >= n:
                flag_cont = False
            else:
                if Xlist[k].v == 'a':
                    k += 1
                else:
                
                    i, j = Xlist[k].get_Pair()
                    k = self.find_ppStem_head(k, Xlist)
                    tails += [(i,j)]
                    k +=1
                #
            #
        #
        return tails
    #
    
    def find_ppStem_head(self, k, Xlist):
        flag_cont = True
        if k >= len(Xlist) - 1:
            flag_cont = False
        while flag_cont:
            if k >= len(Xlist) - 1:
                flag_cont = False
            else:
                i_0, j_0 = Xlist[k  ].get_Pair()
                i_1, j_1 = Xlist[k+1].get_Pair()
                if i_0 + 1 == i_1 and j_0 + 1 == j_1:
                    k+=1
                else:
                    flag_cont = False
                #
            #
        #
        return k
    #
    
    def compress_PKlist(self, stems):
        flag_cont = True
        comp_stems = stems
        k = 0
        while flag_cont:
            if k >= len(stems) - 1:
                flag_cont = False
            else:
                i0 = stems[k  ][0]; j0 = stems[k  ][1]
                i1 = stems[k+1][0]; j1 = stems[k+1][1]
                if i0 < i1 and j1 < j0:
                    del comp_stems[k+1]
                else:
                    k += 1
                #
            #
        #
        return comp_stems
    #
    
    
    
    def scan_for_CTCFs(self):
        debug = False
        if not self.initialized_apPK:
            print "ERROR: antiparallel secondary structure or PK data not initialized"
            sys.exit(1)
        #
        
        CTCFisland = []
        count = 0
        for cl in self.CTCFlist:
            #if count == 2:
            #    sys.exit(0)
            i_W = cl.i; j_W = cl.j
            if debug:
                print "search ij_W(%3d,%3d), len(lt.thread) = %d" \
                    % (i_W, j_W, len(self.lt.thread))
            #
            start = self.locate_start_insert(i_W, j_W)
            self.lt.thread.insert(start, LSegment((i_W,j_W), 0, 0.0, 'W', 'bgn'))
            start += 1
            CTCFisland = self.vs.expand_island(cl)
            for cik in CTCFisland:
                if debug:
                    print "insert pair: (%d,%d): " % (cik.i, cik.j)
                #
                dG = -100/float(len(CTCFisland))
                self.lt.thread.insert(start, LSegment((cik.i,cik.j), 0, dG, 'W', 'wyspa'))
                start += 1
            #
            end = self.locate_end_insert(i_W, j_W, start)
            self.lt.thread.insert(end, LSegment((i_W,j_W), 0, 0.0, 'W', 'end'))
            if debug:
                for thk in self.lt.thread:
                    print thk.disp_lseg()
                print "-------"
            #
            count += 1
        #
        
        if debug:
            print "len(self.CTCFlist): ", len(self.CTCFlist)
            for k in range(0, len(self.CTCFlist)):
                print self.CTCFlist[k].disp_Pair()
            #
            print "islands"
            for cl in CTCFisland:
                print cl.disp_Pair()
            #
            print "LSegment"
            for kk in range(0, len(self.lt.thread)):
                print self.lt.thread[kk].disp_lseg()
        #
    #
    
    def locate_start_insert(self, i_X, j_X, debug = False):
        if debug:
            print "locate_start_insert(%3d,%3d):" % (i_X, j_X)
        #
        k = 0
        i_prv = 0;      i_nxt = 0
        j_prv = self.N; j_nxt = 0
        start = -1
        for k in range(0, len(self.lt.thread)):
            i_nxt = self.lt.thread[k].ij_ndx[0]
            j_nxt = self.lt.thread[k].ij_ndx[1]
            if debug:
                print "ij_X(%3d,%3d), ij_prv(%3d,%3d), ij_nxt(%3d,%3d)" \
                    % (i_X, j_X, i_prv, j_prv, i_nxt, j_nxt)
                print "i_prv(%d) <= i_X(%d) <= i_nxt(%d), j_nxt(%d) <= j_X(%d)" \
                    % (i_prv, i_X, i_nxt, j_nxt, j_X)
            #
            if i_prv <= i_X and i_X <= i_nxt and j_nxt <= j_X:
                start = k
                break
            #
            
            i_prv = i_nxt; j_prv = j_nxt
        #
        if start == -1:
            start = len(self.lt.thread)
        #
        
        if debug: 
            print "locate_start_insert: start = ", start
        #
        return start
    #
    
    
    def locate_end_insert(self, i_X, j_X, start, debug = False):
        if start >= len(self.lt.thread):
            start = len(self.lt.thread) - 1
        #
        
        if debug: 
            print "locate_end_insert(%3d,%3d)[%d]:" % (i_X, j_X, start)
        #
        k = start
        i_prv = self.lt.thread[k].ij_ndx[0]; i_nxt = 0
        j_prv = self.lt.thread[k].ij_ndx[0]; j_nxt = 0
        end = len(self.lt.thread)
        for k in range(start+1, len(self.lt.thread)):
            i_nxt = self.lt.thread[k].ij_ndx[0]
            j_nxt = self.lt.thread[k].ij_ndx[1]
            if debug:
                print "ij_X(%3d,%3d), ij_prv(%3d,%3d), ij_nxt(%3d,%3d)" \
                    % (i_X, j_X, i_prv, j_prv, i_nxt, j_nxt)
                print "j_prv(%d) <= j_X(%d) <= j_nxt(%d)" \
                    % (j_prv, j_X, j_nxt)
            #
            if j_X <= i_nxt and j_X <= j_nxt:
                end = k
                break
            #
            i_prv = i_nxt; j_prv = j_nxt
        #
        
        if debug:
            print "locate_end_insert: start = ", start
            print "locate_end_insert: end   = ", end
        #
        return end
    #
    
    def find_this_ij(self, k_target, Xlist):
        i = -1; j = -1
        flag_found = False
        for xk in Xlist:
            i, j = xk.get_Pair()
            # print k_target, i,j
            if i == k_target:
                # print i, j
                flag_found = True
                break
           #
        #
        if not flag_found:
            i = -1; j = -1
        #
        return (i,j)
    #
    
    def is_ap_stem(self, ndx, Xlist):
        flag_is_ap_stem = False
        i, j = Xlist[ndx].get_Pair()
        n = len(Xlist)
        ix = -1; jx = -1
        if ndx < n - 1:
            ix, jx = Xlist[ndx+1].get_Pair()
            if ix == i+1 and jx == j-1:
                flag_is_ap_stem = True
            #
        #
        return flag_is_ap_stem
    #
    
    def is_pp_stem(self, ndx, Xlist):
        flag_is_pp_stem = False
        i, j = Xlist[ndx].get_Pair()
        n = len(Xlist)
        ix = -1; jx = -1
        if ndx < n - 1:
            ix, jx = Xlist[ndx+1].get_Pair()
            if ix == i+1 and jx == j+1:
                flag_is_pp_stem = True
            #
        #
        return flag_is_pp_stem
    #
    
    
    def disp_lt(self):
        for ltk in self.lt.thread:
            print ltk.disp_lseg()
        #
    #
    
#

class LThread2Vienna:
    def __init__(self):
        self.debug = False
        self.vs    = ''
        self.ss_seq = ''
    #
    
    def lthread2vienna(self, lt, flag_disp = False):
        dt = DispThreads(lt.sqlen)
        self.ss_seq = dt.makeLThreadDotBracket_VARNA(lt, 0)
        if self.debug:
            print "\"%s\"" % ss_seq
        #
        self.vs = Vstruct()
        self.vs.convert_CTCFstruct(self.ss_seq, flag_disp)
    #
#


class LThread2Motif:
    def __init__(self, N):
        self.N       = N
        self.smap    = Map(self.N)
        self.initialize_l2m   = True
        print "made LThread2Motif()"
    #
    
    def assign_B(self, lt, tr_k):
        print "assign_B"
        v   = lt.thread[tr_k].ij_ndx
        ctp = lt.thread[tr_k].ctp
        btp = lt.thread[tr_k].btp
        dG  = lt.thread[tr_k].dGij_B
        
        i_B = v[0]; j_B = v[1]
        link_B = Link(i_B, j_B, dG, ctp, btp, [v])
        self.smap.glink[i_B][j_B].add_link(link_B)
        for mm in self.smap.glink[i_B][j_B].lg[0].motif:
            print mm.show_Motif()
    #
    
    def assign_I(self, lt, tr_k):
        print "assign_I"
        i_t, j_t = lt.thread[tr_k].get_ij()
        ctp      = lt.thread[tr_k].ctp
        btp      = lt.thread[tr_k].btp
        dG_t     = lt.thread[tr_k].dGij_B
        
        v_h      = lt.thread[tr_k+1].ij_ndx
        i_h, j_h = lt.thread[tr_k+1].get_ij()
        dG_h = self.smap.glink[i_h][j_h].lg[0].Vij
        dG = dG_t + dG_h
        
        print i_t, j_t, dG_t, dG_h
        link_I = Link(i_t, j_t, dG, ctp, btp, [v_h])
        self.smap.glink[i_t][j_t].add_link(link_I)                
        for mm in self.smap.glink[i_t][j_t].lg[0].motif:
            print mm.show_Motif()
    #
    def assign_apS(self, lt, tr_k):
        print "assign_apS"
        flag_cont = True
        k_h = 1
        ndx = tr_k
        n   = len(lt.thread)
        while flag_cont:
            if ndx >= n-1:
                flag_cont = False
            else:
                i_sk0, j_sk0 = lt.thread[ndx].get_ij()
                i_sk1, j_sk1 = lt.thread[ndx+1].get_ij()
                if i_sk0 + 1 == i_sk1 and j_sk1 == j_sk0 - 1:
                    ndx += 1
                    k_h += 1
                else:
                    flag_cont = False
                #
            #
        #
        stem = []
        i_s0, j_s0 = lt.thread[tr_k].get_ij()
        dG = 0.0
        for k in range(0, k_h):
            stem += [(i_s0 + k, j_s0-k)]
            if k < k_h -1:
                dG += lt.thread[tr_k + k_h].dGij_B
        #
        dG += self.smap.glink[i_s0+k_h-1][j_s0-k_h+1].lg[0].motif[0].Vij
        link_S = Link(i_s0, j_s0, dG, 'S', 'sa', stem)
        self.smap.glink[i_s0][j_s0].add_link(link_S)                
        for mm in self.smap.glink[i_s0][j_s0].lg[0].motif:
            print mm.show_Motif()
        return tr_k
    #
    
    def assign_ppS(self, lt, tr_k):
        print "assign_ppS"
        flag_cont = True
        k_h = 1
        ndx = tr_k
        n   = len(lt.thread)
        while flag_cont:
            i_sk0, j_sk0 = lt.thread[ndx].get_ij()
            ctp = lt.thread[ndx].ctp
            btp = lt.thread[ndx].btp
            if ndx >= n-1 or btp == 'end':
                flag_cont = False
            else:
                i_sk1, j_sk1 = lt.thread[ndx+1].get_ij()
                print ndx, i_sk0, j_sk0, i_sk1, j_sk1
                if i_sk0 + 1 == i_sk1 and j_sk1 == j_sk0 + 1:
                    ndx += 1
                    k_h += 1
                else:
                    flag_cont = False
                #
            #
        #
        stem = []
        i_s0, j_s0 = lt.thread[tr_k].get_ij()
        dG = 0.0
        for k in range(0, k_h):
            stem += [(i_s0 + k, j_s0 + k)]
            if k < k_h -1:
                dG += lt.thread[tr_k + k_h].dGij_B
        #
        dG += self.smap.glink[i_s0+k_h-1][j_s0+k_h-1].lg[0].motif[0].Vij
        link_S = Link(i_s0, j_s0, dG, 'S', 'sa', stem)
        self.smap.glink[i_s0][j_s0].add_link(link_S)                
        for mm in self.smap.glink[i_s0][j_s0].lg[0].motif:
            print mm.show_Motif()
        return tr_k
    #
    
    def scan_lt(self, lt, tr_k):
        print "tr_k: ", tr_k, len(lt.thread), lt.thread[tr_k].ctp, lt.thread[tr_k].btp
        if tr_k < 0:
            return 0
        #
        ctp = lt.thread[tr_k].ctp
        btp = lt.thread[tr_k].btp
        if btp == 'bgn' or btp == 'end':
            tr_k += 1
        #
        if tr_k < len(lt.thread)-1:
            print "evaluate tr_k ", tr_k
            i_k0, j_k0 = lt.thread[tr_k  ].ij_ndx
            i_k1, j_k1 = lt.thread[tr_k+1].ij_ndx
            ctp = lt.thread[tr_k].ctp
            btp = lt.thread[tr_k].btp
            print "ij_k0, ij_k1: ", (i_k0, j_k0), (i_k1, j_k1)
            if i_k0 < i_k1 and j_k1 < j_k0:
                tr_k = self.scan_lt(lt, tr_k + 1)
                ctp = lt.thread[tr_k].ctp
                btp = lt.thread[tr_k].btp
                if ctp == 'I':
                    self.assign_I(lt, tr_k)
                elif ctp == 'S' and btp == 'sa':
                    self.assign_apS(lt, tr_k)
                #
                tr_k -= 1
            elif i_k0 < i_k1 and j_k0 < j_k1:
                print "parallel"
                if ctp == 'I':
                    self.assign_I(lt, tr_k)
                elif ctp == 'S' and btp == 'sp':
                    print "pp"
                    self.assign_ppS(lt, tr_k)
                #
                tr_k -= 1
            else:
                print "inside here"
                ctp = lt.thread[tr_k].ctp
                if ctp == 'B':
                    self.assign_B(lt, tr_k)
                #
                tr_k -= 1
            #
        else:
            print "outside here"
            ctp = lt.thread[tr_k].ctp
            if ctp == 'B':
                self.assign_B(lt, tr_k)
            elif ctp == 'I':
                self.assign_I(lt, tr_k)
            tr_k -= 1
            #
        #
        return tr_k
    #

            
    def thread2motif(self, lt):
        l2v = LThread2Vienna()
        l2v.lthread2vienna(lt)
        # 
        dt = DispThreads(self.N)
        dt.makeLThreadDotBracket_VARNA(lt, 0)
        #
        self.scan_lt(lt, 0)
        print "end"
        
       
            
        # remove references to PKs
        #if ctp == 'K' and (btp == 'bgn' or btp == 'end'):
        #    continue
        ##
            
        #if ctp == 'W' and (btp == 'bgn' or btp == 'end'):
        #    if btp == 'end':
        #        sssv[v[0]] = '{'
        #        sssv[v[1]] = '}'
        #    #
        #    continue
        #
        
        #if ctp == 'S' and (btp == 'bgn' or btp == 'end'):
        #    continue
        ##
            
    #

    #
#

def test0(cl):
    
    N = 100
    dG = -1.0  # this is not important for most tests of this type
    dt = DispThreads(N)
    lt = [] # [LThread()] #
    
    ndx = 0
    lt += [LThread()]
    # add_lsegment(ij_ndx = (i,j), l_ndx = 0, dGij_B = dG, ctp = ctype, btp = btype)
    lt[ndx].add_lsegment((0,99), 0, dG, 'S', 'sa')
    lt[ndx].add_lsegment((1,98), 0, dG, 'B', 'sa')
    for thr in lt[ndx].thread:
        print thr.disp_lseg()
    #
    print dt.makeLThreadDotBracket_VARNA(lt[ndx], 0)
    dt.print_ijPairList("test.chpair", lt[ndx])
#

def test1(cl):
    vs = Vstruct()
    #          0         10        20        30        40        50        60        70        80        90
    #          |         |         |         |         |         |         |         |         |         |
    #ss_seq = "{(((.(((((((.[[...)))))..]].((((...))))..)))))..|..([)]([)]...}.{.(..)..((..)).(..)..}...{.....|...|....|....}"
    #ss_seq = "{((((((.A.AAAAA......))))).BBBB........a..aaaaa....bbbb..).|..([)]([)].ABC.abc.}....{.....}"
    #ss_seq = "{((((((.A.AAAAA......))))).((((........a..aaaaa....))))..).|..([)]([)].ABC.abc.}....{.....}"
    #ss_seq = "{((((((.A.AAAAA......))))).((((........a..aaaaa....))))..).|..([)]([)]...}....{.A.B.C..a.b.c.....}"
    #ss_seq = "{((((((.A.AAAAA......))))).((((........a..aaaaa....))))..).|..([)]([)]...}....{.ABC..DEF..abc..def...}"
    #ss_seq = "{((((((.A.AAAAA.<BC..))))).((((.>bc....a..aaaaa....))))..).|..([)]([)]...}....{.A.B.C..a.b.c.....}"
    #ss_seq = "{((((((.A.AAAAA.<BC..))))).((((.>bc.DE.a..aaaaa..de))))..).|..([)]([)].........}....{.....}"
    ss_seq  = "{((((((.A.AAAAA.<BC..))))).((((.>bc....a..aaaaa....))))..).|..([)]([)]...}....{.ABC..DEF..abc..def...}"
    
    if len(cl) > 1:
        ss_seq = cl[1]
    vs.convert_CTCFstruct(ss_seq, False)
    ds = Vienna2LThread(vs)
    ds.vienna2thread()
    ds.disp_lt()
    dt = DispThreads(len(ss_seq))
    # print dt.makeLThreadDotBracket_VARNA(ds.lt, 0)
    l2v = LThread2Vienna()
    l2v.lthread2vienna(ds.lt, True)
#    

def test2(cl):
    aList = [123, 'xyz', 'zara', 'abc']
    aList.insert( 3, 2009)
    print "Final List : ", aList
#
    
def test3(cl):
    #          0         10        20        30        40        50        60        70        80        90
    #          |         |         |         |         |         |         |         |         |         |
    ss_seq  = ".ABCD.......abcd...."
    #ss_seq = ".((((..(...).)))).."
    vs = Vstruct()
    
    if len(cl) > 1:
        ss_seq = cl[1]
    vs.convert_CTCFstruct(ss_seq, False)
    ds = Vienna2LThread(vs)
    ds.vienna2thread()
    ds.disp_lt()
    l2m = LThread2Motif(ds.lt.sqlen)
    l2m.thread2motif(ds.lt)
#    


def main(cl):
    if TEST == 0:
        test0(cl)
    elif TEST == 1:
        test1(cl)
    #
    elif TEST == 2:
        test2(cl)
    elif TEST == 3:
        test3(cl)
#
   
# Main
if __name__ == '__main__':
    main(sys.argv)


