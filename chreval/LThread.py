#!/usr/bin/env python

"""Main Module:   LThread.py 

Objects:       DispLThread
               LNode
               LThread
               LThread2Vienna


Author:        Wayne Dawson
creation date: 170126
last update:   190705
version:       0.1

Purpose:

This is another data structure like Pair (the class for Vienna
formated pairing information) and ChPair (the class for chromatin
formated pairing information). 

Comments:

190705:

   This originated as a result of refactoring of the original LThread
   developed for the chromatin program chreval in when doing analysis
   of large swaths of chromatin data. Originally, all these Kluges
   were included in the existing module Calculate.

   It would be better that these are all unified under one class, but
   they developed somewhat independently at various stages of
   development of independent packages, so a coherent data
   representation was not part of the scheme at the time. Perhaps it
   was also a matter of convenience at the time to just get something
   working.

"""

import sys
import string

"""historically set up in chromatin"""
from HeatMapTools import HeatMapTools

"""for Vienna Pair object representation"""
from Pair import find_this_ij
from Pair import is_ap_stem
from Pair import is_pp_stem
from Pair import Pair
from Pair import SortPair
from Pair import vsPair2list

from Vienna import Vstruct

"""Other objects and tools"""
from FileTools  import getHeadExt

"""Important constants. """

# Pseudoknot and parallel stem 1D notation operators labels.  This
# notation at least works with VARNA. Useful for data configuration,
# representation, and transformation.

from Constants  import num2lpr
from Constants  import num2rpr
from Constants  import lpr2num
from Constants  import rpr2num

# A hardwired infinity. I still resent this sort of nonsense, but it
# seems there is no way to completely work around this darn thing.
from Constants  import INFINITY

"""basic functions"""
from BasicTools import roundoff
from BasicTools import initialize_matrix

"""Tree structure building/analysis tools"""
from NaryTree   import Node
from NaryTree   import NodeAnalysis
from NaryTree   import TreeBuilder

# ################################################################
# ############  Local-global functions and constants  ############
# ################################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# 1. control parameters tests etc in the program

# Test   function
#  0     a _very basic_ test of this module

TEST = 0


# 2. special local global function

PROGRAM      = "LThread.py"  # name of the program

# This is not used so much now that I have introduced GetOpts.py,
# but it may still be useful at the very beginning of the program and
# possibly in other parts.

def usage():
    print "USAGE: %s" % PROGRAM
#

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ################################################################
# ################################################################



"""
#############################################################
#####################   Class LThread   #####################
#############################################################

LThread is useful in building suboptimal structure because you only
need to copy it. Because I want to minimize the content rather than
complicate it, it is probably better to revise LThread representation
to basically a notation similar to Pair in Vienna. Because it contains
some additional (though minimal) information, I still think LThread
may be more useful than Pair, which carries even less information.

"""

class LNode(object):
    # this is basically a structure
    def __init__(self, ij_ndx, dGij_B, ctp, btp):
        self.ij_ndx = ij_ndx # (i, j) This is very important!
        # probably should have ij_ndx.i, ij_ndx.j calling it a pair.
        
        self.dGij_B = dGij_B # the FE of the specific bond at ij
        self.ctp    = ctp    # connection type
        self.btp    = btp    # bond type

        # 190705: thinking ahead, because it will be necessary
        # eventually anyway. Presently, this is simply extra useless
        # baggage, but it makes it more compatible for converting
        # between LNode data format to Pair data format.
        self.ch_i   = 'A'  # chain i
        self.ch_j   = 'A'  # chain j
        #
    #
    
    def get_ij(self):
        return self.ij_ndx[0], self.ij_ndx[1]
    #
    
    def disp_lnode(self):
        s = " (%3d, %3d)[%s-%5s] dGij_B = %8.2f" \
            % (self.ij_ndx[0], self.ij_ndx[1],
               self.ctp, string.ljust(self.btp, 5), self.dGij_B)
        return s
    #
#

class LThread(object):
    def __init__(self, N):
        self.sqlen  = N
        self.thread = []
        self.dG     = INFINITY
        self.TdS    =   0.0
        self.p      = -99.99
    #
    
    def add_lnode(self, ij_ndx, Vij_B, ctp = 'B', btp = 's'):
        self.thread +=  [LNode(ij_ndx, Vij_B, ctp, btp)]
        self.compute_dG()
    #
    
    def compute_dG(self):
        dG = 0.0
        for tk in self.thread:
            dG += tk.dGij_B
        self.dG = dG
    #
#


class DispLThread(object):
    
    def __init__(self, N):
        self.N      = N

        
        
        # options for employing SimRNA (not used!!!)
        self.add_P    = False  # include P to P in restraint file
        self.add_N    = True   # include N to N in restraint file
        self.xi       = 7.0    # [nt] Kuhn length
        self.NNdist   = 9.0    # [A] bp N to N distance
        self.PPdist   = 18.7   # [A] bp P to P distance
        self.Nweight  = -0.2   # the N weight for SimRNA restraint
        self.Pweight  = -0.2   # the P weight for SimRNA restraint
        self.restType = "slope" # option for "slope" or "CLE"
        self.vSeq = []
    #
    
    # display a given thread
    def disp_LThread(self, lt):
        # lt ==> [LThread()]: list of type LThread()
        #         LThread.thread  ==> list of type LNode()
        #                         
        for ltk in lt:
            ltk.compute_dG()
            print "total free energy: %8.3f" % ltk.dG
            print "number of pairs: ", len(ltk.thread)
            for thk in ltk.thread:
                print thk.disp_lnode() # LThread.thread[k].disp_lnode()
            #
        #
        # sys.exit(0)
    #
    
    
    
    # display a given thread
    def disp_this_LThread(self, lt, n):
        # lt ==> [LThread()]: list of type LThread()
        #         LThread.thread  ==> list of type LNode()
        #                         
        lt[n].compute_dG()
        print "total free energy: %8.3f" % lt[n].dG
        print "number of pairs: ", len(lt[n].thread)
        for thk in lt[n].thread:
            print thk.disp_lnode() # LThread.thread[k].disp_lnode()
        #
        
        # sys.exit(0)
    #
    

    
    
    def set_strand_direction(self, ctp, btp, counter, callprgm = "undefined"):
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
                    print "problems in resolving structure in %s():" % callprgm
                    sys.exit(1)
                #
            #
        else:
            # probably a PK
            i_lab = num2lpr[counter]
            j_lab = num2rpr[counter]
            counter += 1
            
        
        return i_lab, j_lab, counter
    #
    
    def makeBlank(self):
        bbb = []
        for i in range(0, self.N):
            bbb += ['.']
        return bbb
    #
    
    def makeLThreadDotBracket_1b(self, lt, structure_layout = 0, is_chromatin = True):
        flag_debug = False # True # 
        if flag_debug:
            print "Enter makeLThreadDotBracket_1b()", len(lt.thread)
        #
        seqv  = self.vSeq
        sssv  = []
        ctcfv = []
        pkv   = []
        npks  = -1
        k_pointer = 0
        nws   = -1
        w_pointer = 0
        counter = 3

        if is_chromatin:
            for i in range(0, self.N):
                seqv += ['c']
            #
        #
        
        for i in range(0, self.N):
            sssv   += ['.']
            ctcfv  += ['.']
        #
        
        ktr = -1
        while ktr < len(lt.thread) - 1:
            ktr += 1
            tr = lt.thread[ktr]
            v = tr.ij_ndx
            ctp = tr.ctp
            btp = tr.btp
            if flag_debug:
                print v, ctp, btp, tr.dGij_B        
            #
            
            if not counter < len(num2lpr): 
                # this is _not_ good, and it can lead to major issues
                # for long sequences ... here it just arbitrarilly
                # forces it down without consideration of anything, so
                # it is unlikely to be a good situation if it
                # happens. Anyway, presently, there are not any other
                # ways to do a single line sequence with parallel
                # stems using varna.
                counter = 3
            #
            
            # remove references to PKs
            if ctp == 'K' and btp == 'bgn':
                
                npks += 1
                pkv += [self.makeBlank()] # add another linkage channel
                if flag_debug:
                   print "begin PK(%d)" % npks
                #
                continue
            #
            if ctp == 'R' and btp == 'bgn':
                
                npks += 1
                pkv += [self.makeBlank()] # add another linkage channel
                if flag_debug:
                   print "begin PK(%d)" % npks
                #
                continue
            #
            if ctp == 'K' and (btp == 'end' or (npks > k_pointer and btp == 'l')):
                # This seems to be a problematic spot
                
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
                if btp == 'bgn' or btp == 'end':
                    # I think this is probably the point where we can
                    # increment the counter over a succession of of
                    # the same antiparallel stem; end should certainly
                    # kick this to increment up one.

                    # For the case of parallel stems, this is handled
                    # with the step below. So it seems like what needs
                    # to be done is look a little deeper into
                    # things. It might be good to have the label
                    # "bgn-ap" and "bgn-pp". This would allow
                    # incrementing the thing as a unit (ap) or
                    # individually (pp). It would also be possibe to
                    # look ahead to see whether the set is 'ap' or
                    # 'pp'. 
                    continue
                #
                #print v
                if btp == 'sp':
                    i_lab, j_lab, counter \
                        = self.set_strand_direction(ctp,
                                                    btp,
                                                    counter,
                                                    "makeLThreadDotBracket_1b")
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
                if is_chromatin:                
                    seqv[v[0]] = 'x'
                    seqv[v[1]] = 'y'
                #
            elif btp == 'c':
                if is_chromatin:                
                    seqv[v[0]] = 'W'
                    seqv[v[1]] = 'Z'
                #
                sssv[v[0]] = '{'
                sssv[v[1]] = '}'
            elif btp == 't' or btp == 'r' or btp == 'l':
                if is_chromatin:                
                    seqv[v[0]] = 'w'
                    seqv[v[1]] = 'z'
                #
                sssv[v[0]] = '{'
                sssv[v[1]] = '}'
            elif btp == 'wyspa':
                if is_chromatin:                
                    seqv[v[0]] = 'I'
                    seqv[v[1]] = 'I'
                #
            elif btp == '-':
                if ctp == 'P' or ctp == 'J':
                    # 'P' and 'J' are real possibilities now.
                    if is_chromatin:                
                        seqv[v[0]] = 'c'
                        seqv[v[1]] = 'c'
                    #
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
                #sys.exit(0)
            #
        #
        return s
    #
    
    def makeLThreadDotBracket_VARNA(self, lt, structure_layout = 0, is_chromatin = True):
        flag_debug = False # True
        seqv  = self.vSeq
        sssv  = []
        counter = 3
        
        if is_chromatin:
            for i in range(0, self.N):
                seqv += ['c']
                sssv += ['.']
            #
        else:
            for i in range(0, self.N):
                sssv += ['.']
            #
        #
        
        ktr = -1
        #for tr in lt.thread:
        while ktr < len(lt.thread) - 1:
            ktr += 1
            tr = lt.thread[ktr]
            v = tr.ij_ndx
            ctp = tr.ctp
            btp = tr.btp
            
            if flag_debug:
                print ktr, v, ctp, btp, sssv
            #
            
            # remove references to PKs
            if ctp == 'K' and (btp == 'bgn' or btp == 'end'):
                if btp == 'bgn':
                    ctpt = lt.thread[ktr+1].ctp
                    if flag_debug:
                        print "K:  ", ctpt
                    #
                    if ctpt == 'P':
                        ktr += 1
                    #
                    #sys.exit(0)
                continue
            #
            # remove references to PKs
            if ctp == 'R' and (btp == 'bgn' or btp == 'end'):
                if btp == 'bgn':
                    ctpt = lt.thread[ktr+1].ctp
                    if flag_debug:
                        print "R:  ", ctpt
                    #
                    if ctpt == 'P':
                        ktr += 1
                    #
                    #sys.exit(0)
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
                if btp == 'bgn' or btp == 'end':
                    # I think this is probably the point where we can
                    # increment the counter over a succession of of
                    # the same antiparallel stem; end should certainly
                    # kick this to increment up one.

                    # For the case of parallel stems, this is handled
                    # with the step below. So it seems like what needs
                    # to be done is look a little deeper into
                    # things. It might be good to have the label
                    # "bgn-ap" and "bgn-pp". This would allow
                    # incrementing the thing as a unit (ap) or
                    # individually (pp). It would also be possibe to
                    # look ahead to see whether the set is 'ap' or
                    # 'pp'. 
                    continue
                #
                
                # print v
                if btp == 'sp':
                    i_lab, j_lab, counter \
                        = self.set_strand_direction(ctp,
                                                    btp,
                                                    counter,
                                                    "makeLThreadDotBracket_VARNA")
                    sssv[v[0]] = i_lab
                    sssv[v[1]] = j_lab
                else:
                    # I think this will have problems if there is more
                    # than one antiparallel PK. It should also
                    # increment, but the trouble is that we do not so
                    # easily know the beginning and end of this
                    # container; i.e., we work at the micro-level from
                    # the thread list, not from the container level.
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
                if is_chromatin:
                    seqv[v[0]] = 'x'
                    seqv[v[1]] = 'y'
                #
            elif btp == 'c':
                sssv[v[0]] = '{'
                sssv[v[1]] = '}'
                if is_chromatin:
                    seqv[v[0]] = 'W'
                    seqv[v[1]] = 'Z'
                #
            elif btp == 't' or btp == 'r' or btp == 'l':
                sssv[v[0]] = '{'
                sssv[v[1]] = '}'
                if is_chromatin:
                    seqv[v[0]] = 'w'
                    seqv[v[1]] = 'z'
                #
            elif btp == 'wyspa':
                if is_chromatin:
                    seqv[v[0]] = 'I'
                    seqv[v[1]] = 'I'
                #
            elif btp == '-':
                if ctp == 'P' or ctp == 'J':
                    # 'P' and 'J' are real possibilities now.
                    if is_chromatin:
                        seqv[v[0]] = 'c'
                        seqv[v[1]] = 'c'
                    #
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
        
        # calc.fe.mtools would also work, but I want this module to be
        # a little more independent and it doesn't cost so much to
        # initialize a separate object of class HeatMapTools.
        mtools = HeatMapTools() # default setup GenerateHeatMapTools()
        
        hmap = []
        hmap = initialize_matrix(hmap, self.N, 0)
        
        
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
            if ctp == 'l' and (btp == 'bgn' or btp == 'end'):
               continue
            #
            
            i = v[0]; j = v[1]
            hmap[i][j] = 1
            hmap[j][i] = 1
        #
        return mtools.make_heatmap(hmap, flag_no_header)
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
    
    
    def printLThreadDotBracket_VARNA(self, flnm, lt, is_chromatin = True):
        ff = flnm.split('.')
        if not lt.p == -99.99:
            s = "> %s    %8.3f   %10.8f\n" % (ff[0], lt.dG, lt.p)
        else:
            s = "> %s    %8.3f\n" % (ff[0], lt.dG)
        s += self.makeLThreadDotBracket_VARNA(lt, 1, is_chromatin)
        try:
            fp = open(flnm, 'w')
        except IOError:
            print "ERROR: cannot open %s" % flnm
            sys.exit(1)
        fp.write(s)
        fp.close()
    #

    # NOTE(190706): In the original DispThreads, the method
    # getLThread2DotBracket() was originally located here. Because the
    # method depended on TreeNode2DotBracket and Vienna2TreeNode, it
    # would introduce a circular definition were I to call LThread in
    # TreeNode and then call TreeNode in LThread. Therefore, this
    # method was coverted to an independent class LThread2DotBracket
    # to achieve this purpose and can only be accessed this way.
#



class LThread2Vienna(SortPair):
    """@

    This class is designed to convert an object of class LThread to an
    object of class Vienna. Specifically, the input is the list
    "thread" from class LThread (which is a list of objects of class
    LNode). 

    It is maybe important to emphasis that this is a full scale
    reduction of the LThread class (some details are lost) to the
    slightly similar form of various lists objects of type Pair in
    Vienna. Originally, LThread was intended to be simplified
    representation of structures based on information derived from an
    object of class Map (typically smap in these
    programs). Unfortunately, it was only a little more advanced in
    representation compared to Vienna, so it did not really serve that
    purpose so well. I'm still struggling with the question of whether
    to upgrade LThread to provide better details (the most important
    being branching) or merge LThread with the classes Vienna and
    ChPair (another bastard representation that was developed on the
    fly).

    At any rate, the object of class Vienna can be used to generate
    objects from the module Motif, and this, in turn, can be used to
    generate much larger representations such as an object of class
    Map. Therefore, one of the key purposes of this module is
    debugging. With tools like this, we can reduce parts of an actual
    Map to Vienna representation and regenerate that part of the
    original Map as a new object. Most of the time, this should be
    done by writing secondary structure into Vienna and converting it
    to Map, but LThread2Vienna was introduced as proof of principle so
    we could go backwards and forwards between Map and Vienna.

    """
    def __init__(self):
        self.vsBPlist  = [] # ss
        self.vsBProots = []
        self.vsPKlist  = [] # PKs
        self.vsPKroots = []
        self.vsMPlist  = [] # e.g., ctcf etc
        self.N         = -1
        self.vstr      = ''
        self.vseq      = ''
        """
        The most knarly of messes is the islands, so I think the easiest
        thing to do is separate these out first and solve
        them. Afterwards, it is fairly easy to deal with all the rest
        of the structrural types (when converting it to simple Vienna
        Pair representations, of course).
        
        The output should finally be equivalent to self.vsBPlist and
        self.vsPKlist. Presently, I am focusing on self.vsBPlist
        constructions. ... in part because I don't have the PK
        constructions working in Vienna2TreeNode.
        
        Finally, to build a proper conversion to Vienna, we need to
        address the following variables.
        
        Vstruct:   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        self.vstr ='' # structure
        self.vseq ='' # sequence (optional) 
        self.N    = -1  # sequence length
        self.wt   = -1 # used with my_generation to add specified weights to heatmaps
        # standard secondary structure (ss)
        self.vsBPlist  = []
        self.vsBProots = []
        # a both ss and pseudoknot (PK) in same line
        self.vsPKlist  = []
        self.vsPKroots = []
        # CTCF islands
        self.MPlist = []
        
        Comments: 
        
        * self.MPlist is actually what I have already figured out as
          part of the self.vsBPlist. I guess it is not in BPlist but
          separate, but still uses Pair representation.
        
        * We can get self.N from LThreads. 
        
        * From chreval, we can get a "structure" sequence for
          self.vstr and self.vseq can be "generated"
        
        * The main thing I have to consider is how to generate PKroots
          and BProots, but I think that is not so hard.
        
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        """
    #
    
    def reset_LThread2Vienna(self):
        self.vsBPlist  = [] # ss
        self.vsBProots = []
        self.vsPKlist  = [] # PKs
        self.vsPKroots = []
        self.vsMPlist  = [] # e.g., ctcf etc
        self.N         = -1
        self.vstr      = ''
        self.vseq      = ''
    #
    
    def set_vstr(self, vstr):
        self.vstr = vstr
    #
    
    def set_vseq(self, vseq):
        self.vseq = vseq
    #
    
    
    def lt2vs(self, lt):
        debug_lt2vs = False # True # 
        if debug_lt2vs:
            print "Entered lt2vs"
        #
        
        self.reset_LThread2Vienna()
        self.N = lt.sqlen
        
        # first deal with the objects of class MultiPair
        wPairList = []
        for tr in lt.thread:
            ijv = tr.ij_ndx
            ctp = tr.ctp
            btp = tr.btp
            
            if ctp == 'W' and btp == "wyspa":
                #
                i = ijv[0]; j = ijv[1]
                wPairList += [(i,j, ctp, "ctcf")]
            elif ctp == 'w' and btp == "wyspa":
                i = ijv[0]; j = ijv[1]
                wPairList += [(i,j, ctp, "ctcf")]
            #
        #
        
        if debug_lt2vs:
            print "wPairList: "
            for ww in wPairList:
                print ww
            #
        #
        
        if len(wPairList) > 0:
            # now group related W in to "lslands" (wyspa)
            wGrpList = []
            ka = 0
            if debug_lt2vs:
                print "0 wPairList: ", wPairList
            #
            
            while ka < len(wPairList)-1:
                iw1 = wPairList[ka][0]; jw1 = wPairList[ka][1]
                kb = ka + 1
                newgroup = [wPairList[ka]]
                # we can trust that the LThread organization will always
                # place the largest arch first on the list.
                while kb < len(wPairList):
                    iw2 = wPairList[kb][0]; jw2 = wPairList[kb][1]
                    if iw1 <= iw2 and jw2 <= jw1:
                        newgroup += [wPairList[kb]]
                        del wPairList[kb]
                    else:
                        kb += 1
                    #
                #
                wGrpList += [newgroup]
                ka += 1
            #
            ka_max = len(wPairList)-1
            ka_last = len(wGrpList) -1
            if debug_lt2vs:
                print "wGrpList:  ", wGrpList
                print "1 wPairList: ", wPairList
            #
            
            if len(wGrpList) > 0:
                if wGrpList[ka_last][0][1] < wPairList[ka_max][0]:
                    wGrpList += [[wPairList[ka_max]]]
                #
            else:
                wGrpList += [[wPairList[ka_max]]]
            if debug_lt2vs:
                print "wGrpList: "
                for ww in wGrpList:
                    print ww
                #
            #
            
            # sys.exit(0)
            # Now we can build the list of islands in Pair format
            for k in range(0, len(wGrpList)):
                ww = wGrpList[k]
                #print len(ww), ww
                #print ww[0]
                
                iw = ww[0][0]; jw = ww[0][1]
                p = Pair()
                p.put_ssPair(iw, jw, "ctcf", 'W')
                # put_ssPair(self, i, j, nm = 'bp', v = 'a')
                last_k = iw
                for k in range(1, len(ww)):
                    iwk = ww[k][0]; jwk = ww[k][1]
                    if not iwk == last_k:
                        p.contacts += [iwk]
                        last_k = iwk
                    if not jwk == jw:
                        p.contacts += [jwk]
                        last_k = jwk
                    #
                #
                self.vsMPlist += [p]
            #
            if debug_lt2vs:
                print "vsMPlist: "
                for mpk in self.vsMPlist:
                    print mpk.disp_Pair()
                #
                # sys.exit(0)
            #
        #
        
        for tr in lt.thread:
            ijv = tr.ij_ndx
            ctp = tr.ctp
            btp = tr.btp
            
            # need to think about how these are assigned!!!!
            
            # remove references to PKs
            if ctp == 'K' and (btp == 'bgn' or btp == 'end'):
                continue
            #
            if ctp == 'R' and (btp == 'bgn' or btp == 'end'):
                continue
            #
            
            # need to devise something to pick out the linkage stems
            # here with K and R. These would be assigned to vsPKlist.
            
            if ctp == 'W' or ctp == 'w':
                continue
            #
            
            if ctp == 'S' and (btp == 'bgn' or btp == 'end'):
                continue
            #
            if ctp == 'P' or ctp == 'J':
                continue
            #
            if ctp == 'l' and (btp == 'bgn' or btp == 'end'):
                continue
            #
            i = ijv[0]; j = ijv[1]
            p = Pair()
            p.put_ssPair(i, j)  # put_ssPair(self, i, j, nm = 'bp', v = 'a')
            if btp == "sp" or btp == "sa" or btp == 's':
                if debug_lt2vs:
                    print "ij(%2d,%2d)[btp=%s]" % (i, j, btp)
                #
                if btp == "sp":
                    p.v = "p"
                    p.name = "bp"
                else:
                    p.v = "a"
                    p.name = "bp"
                #
                self.vsBPlist += [p]
            elif btp == "lp" or btp == "la":
                if debug_lt2vs:
                    print "ij(%2d,%2d)[btp=%s]" % (i, j, btp)
                #
                if btp == "lp":
                    p.v = "p"
                    p.name = "bp"
                else:
                    p.v = "a"
                    p.name = "bp"
                #
                self.vsPKlist += [p]
            else:
                if debug_lt2vs:
                    print "lt2vs: found an unusual one..."
                    print "ij(%d,%d), ctp = %s, btp = %s" % (i, j, ctp, btp)
                    sys.exit(0)
                    p.v = btp
                    p.name = "mp"
                #
            #
                
        #
        
        self.vsBPlist = self.sortvsList(self.vsBPlist, 'i')
        self.vsPKlist = self.sortvsList(self.vsPKlist, 'i')
        vs = Vstruct()
        self.vsBProots = vs.findroots(self.vsBPlist, False)        
        self.vsPKroots = vs.findroots(self.vsPKlist, False)

        """@
        
        190723: 

           After sorting, check for redundant indicies. Originally,
           BIM were real indicies that were special to the structures;
           however, now that we have RNA where we require dinucleotide
           base pairs, this is not always a valid arrangement. Because
           there are still some compatibility issues with the
           chromatin and RNA programs, I cannot take a unilateral
           approach to the problem without affecting one or the other
           of the current approaches in their current state of
           development.
        
           So the current work around presently is to just remove
           redundant indicies. I don't like this so much, but I cannot
           work around the matter so easily without making the
           definition of Stem exactly (absolutely exactly) the same
           for both programs, which currently they are not.

        """
        
        self.vsBPlist = self.filter_redundant_indices(self.vsBPlist)
        self.vsPKlist = self.filter_redundant_indices(self.vsPKlist)
        
        if debug_lt2vs:
            print "vsBPlist: "
            for vv in self.vsBPlist:
                print vv.disp_Pair()
            #
            print "vsBProots: "
            for vv in self.vsBProots:
                print vv.disp_Pair()
            #
            print "vsPKlist: "
            for vv in self.vsPKlist:
                print vv.disp_Pair()
            #
            print "vsPKroots: "
            for vv in self.vsPKroots:
                print vv.disp_Pair()
            #
            print "vsMPlist: "
            for vv in self.vsMPlist:
                print vv.disp_Pair()
            #
        #

        return self.vsBPlist
    #

    def filter_redundant_indices(self, Xlist):
        k = 0
        while k < len(Xlist)-1:
            i1 = Xlist[k  ].i; j1 = Xlist[k  ].j
            i2 = Xlist[k+1].i; j2 = Xlist[k+1].j
            if i1 == i2 and j1 == j2:
                #print "found redundant index ij1(%2d,%2d) ij2(%2d,%2d)" % (i1, j1, i2, j2)
                #sys.exit(0)
                del Xlist[k]
            else:
                k += 1
            #
            
        #
        return Xlist
    #

#


##
# LThread construction routine


def test0(cl):
    
    N = 100
    dG = -1.0  # this is not important for most tests of this type
    dt = DispLThread(N)
    lt = [] # [LThread()] #
    
    ndx = 0
    lt += [LThread(N)]
    # add_lnode(ij_ndx = (i,j), dGij_B = dG, ctp = ctype, btp = btype)
    lt[ndx].add_lnode((0,99), dG, 'S', 'sa')
    lt[ndx].add_lnode((1,98), dG, 'B', 'sa')
    for thr in lt[ndx].thread:
        print thr.disp_lnode()
    #
    print dt.makeLThreadDotBracket_VARNA(lt[ndx], 0)
    # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    from ChPair     import LThread2ChPair
    chdt = LThread2ChPair(lt[ndx], "test0")
    chdt.p = 1.0
    # the argument "test0" is just used to give this data set a name.
    chdt.print_ChPairData()
    # no argument in print_ChPairData() means that the output will
    # only be displayed.
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
#





def main(cl):
    if TEST == 0:
        test0(cl)
    #
    
#
   
# Main
if __name__ == '__main__':
    main(sys.argv)
#

