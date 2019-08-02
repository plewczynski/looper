#!/usr/bin/env python

"""@@@

Main Module:   Trace.py 

Classes:       Trace

Author:        Wayne Dawson
creation date: mostly 2016 a little bit in 2017 up to March
last update:   190402
version:       0

Purpose:

Trace() does the trace back steps for a particular structure.


Comments:

At this point, Calculate() has been run (there should be a check for
that here!!!). The program runFE() only searches for the structure
with the minimum free energy. However, if we take advantage of the
additional components in the push down stack LGroup(), we can scan a
distribution of structures. Ultimately, there should be a sort
function that organizes the suboptimal structures in some progressive
logical order.

"""

from math import exp
import sys
import os
import string

# Motif object representation
from Motif import Motif # a core object of these building programs
from Motif import Link
from Motif import LGroup
from Motif import Map # main map for all the FE and structures
from Motif import XLoop
from Motif import MBL
from Motif import Stem
from Motif import PseudoKnot

# LThread object representation
from LThread import DispLThread
from LThread import LThread

# Other objects and tools
from FileTools   import getHeadExt

# ################################################################
# ######################  Global constants  ######################
# ################################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# for entropy calculations

# A hardwired infinity. Somewhat hate this, but it seems there is no
# way around this darn thing.
from Constants import INFINITY

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ################################################################
# ################################################################


# ################################################################
# ###########  global functions/objects and constants  ###########
# ################################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# 1. control parameters for the program

PROGRAM      = "Trace.py"  # name of the program

# debugging options

# Debugging in main()
SHOWMAIN             = False #

# debugging Traces()
# the top of the program
DEBUG_get_traces_top = False # True # 
# work all the way into the details
DEBUG_get_traces     = False # True # 

# Generally for checking Trace()
STOP_AT_FIND = False  # stops program when encouters condition

# 2. special local global function

def usage():
    print "USAGE: %s" % PROGRAM
#
# At present, there is not much reason to use this at all. I presently
# don't even have a simple test to run it with. Nevertheless, maybe it
# will be useful at some point.

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ################################################################
# ################################################################


class Trace:
    def __init__(self, calc):
        # calc: the object Calculate()
        self.N        = calc.N
        self.T        = calc.T
        self.hv       = calc.fe.hv
        self.smap     = calc.fe.smap
        self.calc     = calc
        self.dGmin    = INFINITY
        
        
        # the linkage thread variable
        self.lt = [] # [LThread()] # 
        self.debug_traces_top = DEBUG_get_traces_top
        # find hot spot regions for further search
        self.debug_get_traces = DEBUG_get_traces
        # search for suboptimal structures
        self.warning_boundary = 10000
        self.maxlayers = self.calc.fe.maxlayers
        # This is used in traceback routines. See comments in the
        # constructor of Calculate()
        
        # Hot Spots
        self.wt_HS = 0.1
        self.V_HS = INFINITY
        self.hotspot = []
        
        
    #
    
    """@
    
    This should mainly be used for debugging. It displays all the
    glink data for all the data points in the upper triangle, so it
    can print out a lot of data. However, it can be useful for
    detecting problems in the calculated results.
    
    """
    def disp_allFE(self, smap):
        for j in range(0, self.N):
            for i in range(0,j):
                if len(smap.glink[i][j].lg) > 0:
                    for k in range(0, len(smap.glink[i][j].lg)):
                        ctp  = smap.glink[i][j].lg[k].motif[0].get_ctp() # conn type
                        bdt  = smap.glink[i][j].lg[k].motif[0].get_btp() # bond type
                        Vij  = smap.glink[i][j].lg[k].Vij
                        branching = smap.glink[i][j].lg[k].motif[0].get_base()
                        print "(%2d, %2d)[%2d][%s_%s][%8.2f]: " % (i,j, k, ctp, bdt, Vij), branching
                else:
                    print "(%2d, %2d): " % (i,j), "empty"
                #
            #
        #
    #
    
    
    
    def setup_HotSpot(self, dGmin, wt = 0.7):
        self.wt_HS = wt
        self.dGmin = dGmin
        dGmax_ratio = self.wt_HS * self.dGmin
        dGmax_diff  = self.dGmin + self.calc.dGrange
        if dGmax_ratio > dGmax_diff:
            dGratio = dGmax_diff/dGmin
            print "selecting max FE from the difference result +/- %8.2f (ratio: %8.4f)" \
                % (self.calc.dGrange, dGratio)
            self.V_HS = dGmax_diff
            self.wt_HS = dGratio
        else:
            dGrange = dGmax_ratio - dGmin 
            print "selecting max FE from the ratio result %8.4f (+/- %8.2f)" \
                % (self.wt_HS, dGrange)
            self.V_HS = dGmax_ratio
        #
    #
    
    
    # sort the results from HotSpots
    def ins_SortHotSpot(self, s):
        for i in range(1,len(s)):    
            j = i                    
            while j > 0 and s[j][3] < s[j-1][3]: 
                s[j], s[j-1] = s[j-1], s[j] # syntactic sugar: swap the items
                j=j-1 
            #
        #
        return s
    #
    
    """@
    
    The top of the list need not be an 'M', 'K', 'I' or 'B'!!!  So
    this makes sure that the structure that is selected is one of
    these. It is used in I-loop searches; i.e., find_best_I() and
    related structures. So the purpose is to find a structure with a
    single closing point at (k,l), not a pMBL or a J-loop. The main
    question I have now is whether this should also include a
    PK. Currently (160911) by scanNN_MIB(), scan_narrow(),
    find_best_I(), localscan_for_I()
    
    """
    def get_kref(self, k, l, V, ctp):
        #
        flag_debug = False # True # 
        
        kref = 0
        if flag_debug:
            print "get_kref", k, l, V, ctp
        for m in range(0, len(self.smap.glink[k][l].lg)):
            ctpx = self.smap.glink[k][l].lg[m].motif[0].get_ctp()
            btpx = self.smap.glink[k][l].lg[m].motif[0].get_btp()
            Vx = self.smap.glink[k][l].lg[m].Vij
            if flag_debug:
                print "(%2d,%2d)[%d]: %8.3f [%s][%s]" % (k,l, m, Vx, ctpx, btpx)
            if ctpx == ctp and Vx >= V:
                kref = m
                break
                #
        if flag_debug:
            print "kref = ", kref
        #
        return kref
    #
    
    
    def localscan_getBranches(self, p, q, ctp, branching):
        flag_debug = False
        
        if ctp == 'S' or ctp == 'K' or ctp == 'R' or ctp == 'W' or ctp == 'B':
            branching += self.smap.glink[p][q].lg[0].motif[0].get_base()
        #
        elif ctp == 'J' or ctp == 'I':
            for jk in self.smap.glink[p][q].lg[0].motif[0].get_branches():
                branching += [jk]
            #
        #
        elif ctp == 'P' or ctp == 'M':    
            for jk in self.smap.glink[p][q].lg[0].motif[0].get_branches():
                branching += [jk]
            #
        #
        else:
            print "ERROR: encountered an unrecognized symbol (%s)" % (ctp)
            print "in localscan_for_M(%d,%d)" % (p,q)
            sys.exit(1)
            
        #
        if flag_debug:
            print "ctp = %s" % ctp
            print "branching: ", branching
        #
        return branching
    #
    
    """@
    
    This is in a similar style as the M-loop calculation in
    Calculate(). However, here, it saves more of the information and
    it does not care if the search turns up type 'P' or 'J' or any of
    the multitude of other motif types.  Everything that satisfies the
    cutoff is fair game.
    
    An important thing to remember about "list_M" is that "list_M" is
    a place holder (a container) for the particular motif
    object. Therefore, the container motif must point _to_ the motif
    object.
    
    """
    def localscan_for_M(self, i, j, list_M, minFE):
        debug_localscan_for_M1 = False # True # 
        debug_localscan_for_M2 = False # True # 
        best_Mk = ()
        best_dGMk = 0.0
        flag_M  = False
        best_Ik = ()
        best_dGIk = 0.0
        flag_I  = False
        motif = []
        fpp = {(i,j) : 0, (i+1,j) : -1, (i+1,j-1) : -2, (i,j-1) : -3 }
        ipp = {0 : (i,j), -1 : (i+1,j), -2 : (i+1,j-1), -3 : (i,j-1) }
        
        if debug_localscan_for_M1:
            location = "Trace: enter localscan_for_M(%2d,%2d)" % (i, j)
            self.show_list_M(list_M, location)
        #endif
        
        for k in range(1, j-i-1):
            branching = []
            E1 = self.calc.dG[i][i+k]; E2 = self.calc.dG[i+k+1][j]
            V = E1 + E2
            
            ctp1 = 'X'; ctp2 = 'X'
            n1 = len(self.smap.glink[i][i+k].lg)   # n1 > 0 if Lgroup has data
            n2 = len(self.smap.glink[i+k+1][j].lg) # n2 > 0 if Lgroup has data
            if not n1 == 0:
                ctp1 = self.smap.glink[i  ][i+k].lg[0].motif[0].get_ctp()
            #
            if not n2 == 0:
                ctp2 = self.smap.glink[i+k+1][j].lg[0].motif[0].get_ctp()
            #
            if debug_localscan_for_M2:
                print "localscan_for_M: k(%2d), (%2d,%2d)[%8.2f][%s]|(%2d,%2d)[%8.2f][%s] V_HS(%8.2f)" \
                    % (k, i, i+k, E1, ctp1, i+k+1, j, E2, ctp2, self.V_HS)
            #
            
            # The variable "V_HS" is used to select only free energy
            # within a range of the minimum free energy. In effect, it
            # slices off the outside shell of free energy (something
            # like an onion) and records those structures.
            if V < self.V_HS or E1 < self.V_HS or E2 < self.V_HS:
                
                if debug_localscan_for_M2:
                    print "n12: ", n1, n2, (i, j), k
                #
                if not (ctp1 == 'X' or ctp2 == 'X'): # was n1 > 0 and n2 > 0:
                    if not V >= minFE:
                        continue
                    #
                    flag_M = True
                    # found some real M-loop
                    branching = []
                    branching  = self.localscan_getBranches(i,     i+k, ctp1, branching)
                    branching  = self.localscan_getBranches(i+k+1, j,   ctp2, branching)
                    
                    motif_P = Motif(i, j, V, 'P', '-', branching)
                    list_M += [(i,j, [motif_P], V)]
                    if V < best_dGMk: 
                        best_Mk = (i,j, [motif_P], V, k)
                        best_dGMk = V
                        #
                        if debug_localscan_for_M2:
                            n1v = [n1, self.smap.glink[i    ][i+k].lg[0].motif[0].get_base()]
                            n2v = [n2, self.smap.glink[i+k+1][j  ].lg[0].motif[0].get_base()]
                            print "i(%2d), p(%2d), q(%2d), j(%2d), k(%2d), n12(%2d,%2d) E12(%8.2f, %8.2f), " \
                                % (i, i+k, i+k+1, j, k, n1, n2, E1, E2), branching
                        #
                    #
                elif not ctp1 == 'X': # was n1 > 0:
                    if not E1 >= minFE:
                        continue
                    #
                    ctp  = self.smap.glink[i][i+k].lg[0].motif[0].get_ctp()
                    btp  = self.smap.glink[i][i+k].lg[0].motif[0].get_btp()
                    base = self.smap.glink[i][i+k].lg[0].motif[0].get_base()
                    mjn  = self.smap.glink[i][i+k].lg[0].motif[0].get_branches()
                    if debug_localscan_for_M2:
                        print "base: ", base
                    #
                    branching = []
                    branching = self.localscan_getBranches(i, i+k, ctp1, branching)
                    # self.get_kref(i, i+k, E1, ctp)
                    ctpx = 'J'
                    btpx = '-'
                    if len(branching) > 1:
                        flag_M = True
                        ctpx = 'P'
                    else:
                        flag_I = True
                    #
                    
                    motif_J = Motif(i, j, E1, ctpx, btpx, branching)
                    if debug_localscan_for_M2:
                        if ctp == 'S':
                            print "n1: "
                            print "boundaries: ", i, (i + k)
                            print motif_J.show_Motif()
                            #print "stop at 0 in localscan_for_M"; sys.exit(0)
                        #
                    #
                    list_M += [(i,j, [motif_J], E1)]
                    if len(branching) == 1:
                        if E1 < best_dGIk:
                            best_Ik = (i,j, [motif_J], E1, k)
                            best_dGIk = E1
                            if debug_localscan_for_M2:
                                n1x = [n1, ctp, motif_J.get_branches()]
                                print "1: i(%2d), p(%2d), q(%2d), j(%2d), k(%2d), n1(%2d) E1(%8.2f), " \
                                    % (i, i+k, i+k+1, j, k, n1, E1), branching
                            #
                        #
                    else:
                        if E1 < best_dGMk:
                            best_Mk = (i,j, [motif_J], E1, k)
                            best_dGMk = E1
                            if debug_localscan_for_M2:
                                n1x = [n1, ctp, motif_J.get_branches()]
                                print "1: i(%2d), p(%2d), q(%2d), j(%2d), k(%2d), n1(%2d) E1(%8.2f), " \
                                    % (i, i+k, i+k+1, j, k, n1, E1), branching
                            #
                        #
                    #
                    
                elif not ctp2 == 'X': # was n2 > 0:
                    if not E2 >= minFE:
                        continue
                    #
                    ctp  = self.smap.glink[i+k+1][j].lg[0].motif[0].get_ctp()
                    btp  = self.smap.glink[i+k+1][j].lg[0].motif[0].get_btp()
                    mjn  = self.smap.glink[i+k+1][j].lg[0].motif[0].get_branches()
                    base = self.smap.glink[i+k+1][j].lg[0].motif[0].get_base()
                    if debug_localscan_for_M2:
                        print "base: ", base
                    #
                    branching = []
                    branching = self.localscan_getBranches(i+k+1, j,   ctp2, branching)
                    # self.get_kref(i+k+1, j, E2, ctp)
                    ctpx = 'J'
                    btpx = '-'
                    if len(branching) > 1:
                        flag_M = True
                        ctpx = 'P'
                    else:
                        flag_I = True
                    #
                    
                    motif_J = Motif(i, j, E2, ctpx, btpx, branching)
                    if debug_localscan_for_M2:
                        if ctp == 'S':
                            print "n2: "
                            print "boundaries: ", (i+k+1), j
                            print motif_J.show_Motif()
                            #print "stop at 1 in localscan_for_M"; sys.exit(0)
                        #
                    #
                    list_M += [(i,j, [motif_J], E2)]
                    if len(branching) == 1:
                        if E2 < best_dGIk:
                            best_Ik = (i,j, [motif_J], E2, k)
                            best_dGIk = E2
                            if debug_localscan_for_M2:
                                n2x = [n2, ctp, motif_J.get_branches()]
                                print "2: i(%2d), p(%2d), q(%2d), j(%2d), k(%2d), n2(%2d) E2(%8.2f), " \
                                    % (i, i+k, i+k+1, j, k, n2, E2), branching
                            #
                        #
                    else:
                        if E2 < best_dGMk:
                            best_Mk = (i,j, [motif_J], E2, k)
                            best_dGMk = E2
                            if debug_localscan_for_M2:
                                n2x = [n2, ctp, motif_J.get_branches()]
                                print "2: i(%2d), p(%2d), q(%2d), j(%2d), k(%2d), n2(%2d) E2(%8.2f), " \
                                    % (i, i+k, i+k+1, j, k, n2, E2), branching
                            #
                        #
                    #
                else:
                    print "i(%2d), p(%2d), q(%2d), j(%2d), k(%2d), n12(%2d,%2d) E12(%8.2f, %8.2f) empty" \
                        % (i, i+k, i+k+1, j, k, n1, n2, E1, E2)
                    print "Something is wrong, either sector 1 or sector 2 should have a value"
                    sys.exit(1)
                #
            else:
                continue
            #
            
        #    
        lst = [(i + 1, j), (i+1,j-1), (i, j-1), (i,j)]
        for lstk in lst:
            p = lstk[0]; q = lstk[1]
            E3 = self.calc.dG[p][q]
            
            ctp3 = 'X'
            n3 = len(self.smap.glink[p][q].lg) # n3 > 0 if Lgroup has data
            if not n3 == 0:
                ctp3 = self.smap.glink[p][q].lg[0].motif[0].get_ctp()
            #
            if debug_localscan_for_M2:
                print "localscan_for_M: ij(%2d,%2d)[%8.2f][%s] V_HS(%8.2f)" \
                    % (p, q, E3, ctp3, self.V_HS)
            #
            if E3 < self.V_HS:
                if not E3 >= minFE:
                    continue
                #
                if debug_localscan_for_M2:
                    print "n3: ", n3, (p, q), ctp3
                #
                
                if not ctp3 == 'X': # was n3 > 0:
                    ctp  = self.smap.glink[p][q].lg[0].motif[0].get_ctp()
                    btp  = self.smap.glink[p][q].lg[0].motif[0].get_btp()
                    mjn  = self.smap.glink[p][q].lg[0].motif[0].get_branches()
                    base = self.smap.glink[p][q].lg[0].motif[0].get_base()
                    if debug_localscan_for_M2:
                        print "base: ", base
                    #
                    branching = []
                    branching = self.localscan_getBranches(p, q,   ctp3, branching)
                    # self.get_kref(i+k+1, j, E3, ctp)
                    ctpx = 'J'
                    btpx = '-'
                    if len(branching) > 1:
                        flag_M = True
                        ctpx = 'P'
                    elif i == p and q == j:
                        ctpx = ctp
                        btpx = btp
                        flag_I = True
                    else:
                        flag_I = True
                    #
                    motif_J = Motif(i, j, E3, ctpx, btpx, branching)
                    if debug_localscan_for_M2:
                        if ctp == 'S':
                            print "n3: "
                            print "boundaries: ", p, q
                            print motif_J.show_Motif()
                            print "stop at 2 in localscan_for_M"; sys.exit(0)
                        #
                    #
                    list_M += [(i,j, [motif_J], E3)]
                    if len(branching) == 1:
                        if E3 < best_dGIk:
                            k = fpp[(p,q)]
                            best_Ik = (i,j, [motif_J], E3, k)
                            best_dGIk = E3
                            if debug_localscan_for_M2:
                                print "3: p(%2d), q(%2d), n3(%2d) E3(%8.2f), " \
                                    % (p, q, n3, E3), branching
                            #
                        #
                    else:
                        if E3 < best_dGMk:
                            k = fpp[(p,q)]
                            best_Mk = (i,j, [motif_J], E3, k)
                            best_dGMk = E3
                            if debug_localscan_for_M2:
                                print "3: p(%2d), q(%2d), n3(%2d) E3(%8.2f), " \
                                    % (p, q, n3, E3), branching
                            #
                        #
                    #
                else:
                    continue
                #
            #
        #
        
        
        if debug_localscan_for_M1:
            location = "Trace: exiting localscan_for_M(%2d,%2d)" % (i, j)
            self.show_list_M(list_M, location)
            
            if flag_M:
                print best_dGMk, best_Mk[3]
                k = best_Mk[4]
                branching = best_Mk[2][0].get_branches()
                best_dGMk = best_Mk[3]
                tp1 = self.smap.glink[i    ][i+k].lg[0].motif[0].get_ctp()
                tp2 = self.smap.glink[i+k+1][j  ].lg[0].motif[0].get_ctp()
                print "best_dGMk->: ij(%2d,%2d)(%2d), (%2d,%2d){%s}|(%2d,%2d){%s}[%8.3f] " \
                    % (i, j, k, i, i+k, tp1, i+k+1, j, tp2, best_dGMk), branching
                #
                #if i == 0 and j == 14:
                #    sys.exit(0)
                #print "stop at 3 in localscan_for_M"; sys.exit(0)
            else:
                print "no hits for M/P in localscan_for_M"
            if len(best_Ik) > 0:
                print "best_Ik: ", len(best_Ik), best_Ik
                print "best_dGIk: ", best_dGIk
                print "best_dGIk, best_Ik[3]: ", best_dGIk, best_Ik[3]
                k = best_Ik[4]
                branching = best_Ik[2][0].get_branches()
                best_dGIk = best_Ik[3]
                tp1 = '-'; tp2 = '-'
                if k > 0:
                    if len(self.smap.glink[i    ][i+k].lg) > 0:
                        tp1 = self.smap.glink[i    ][i+k].lg[0].motif[0].get_ctp()
                    #
                    if len(self.smap.glink[i+k+1][j  ].lg) > 0:
                        tp2 = self.smap.glink[i+k+1][j  ].lg[0].motif[0].get_ctp()
                    #
                    print "best_dGIk->: ij(%2d,%2d)(%2d), (%2d,%2d){%s}|(%2d,%2d){%s}[%8.3f] " \
                        % (i, j, k, i, i+k, tp1, i+k+1, j, tp2, best_dGIk), branching
                    #
                else:
                    tp = best_Ik[2][0].ctp
                    p = ipp[k][0]; q = ipp[k][1]
                    print "best_dGIk->: ij(%2d,%2d), (%2d,%2d){%s}[%8.3f] " \
                        % (i, j, p, q, tp, best_dGIk), branching
                    #
                    
                #if i == 0 and j == 14: sys.exit(0)
                #print "stop at 4 in localscan_for_M"; sys.exit(0)
            else:
                print "no hits for I/J in localscan_for_M"
            #
            print "Threads(end localscan_for_M)-----------------"
            self.show_list_M(list_M, "add_localscan_for_M:")
            print "-----------------Threads(end localscan_for_M)"
            #print "stop at 5 in localscan_for_M"; sys.exit(0)
        #endif
        return list_M
    #                
    
    
    def show_list_M(self, list_M, location):
        # list_M -> [(i,j, [motif], V)]
        print location 
        print "list_M: --------"
        for lMv in list_M:
            print "(%2d,%2d) -> " % (lMv[0], lMv[1]), lMv[2][0].show_Motif(), lMv[3]
        #
        print "-------- :list_M"
    #
    
    
    
    
    # This operation removes redundant solutions from the list
    # obtained by localscan_for_M(i, j). The algorithm is derived from
    # the program reduce_list(t) that is listed at the end of this
    # package.
    def prune_list_M(self, list_M):
        debug_prune_list_M = False # True # 
        if debug_prune_list_M:
            print "prune_list_M():"
        #
        
        uniq_list_M = []
        n = len(list_M)
        for l in range(0, n - 1):
            lMl      = list_M[l]
            lMl_dG   = float( int(1000*lMl[3]))/1000.0
            lMl_base = lMl[2][0].get_base()
            flag_match = False
            
            for k in range(0, len(uniq_list_M)):
                lMk    = uniq_list_M[k]
                lMk_dG = float( int(1000*lMk[3]))/1000.0
                if lMl_dG == lMk_dG:
                    lMk_base = lMk[2][0].get_base()
                    
                    if len(lMl_base) == len(lMk_base):
                        if lMl_base == lMk_base:
                            if debug_prune_list_M:
                                print "match lM[2]: ", k, lMk_dG, lMk_base, lMl_base
                            #
                            flag_match = True
                            break
                        #
                    #
                #
            #
            if not flag_match:
                uniq_list_M += [lMl]
            #
        #

        
        if debug_prune_list_M:
            print "finished prune_list_M"
        #
        #print "stop at 0 in prune_list_M"; sys.exit(0)
        return uniq_list_M
    #
    
    def adjust_FEboundaries(self, bump_up):
        v = 10000*float(800)/float(self.N)
        wb_shift = int(v*100)/100
        self.wt_HS += bump_up
        self.V_HS = self.wt_HS * self.dGmin
        self.warning_boundary = 10000 + self.warning_boundary
        
        print "reset: upper limit to number of structures: %d" \
            % self.warning_boundary 
        print "reset: acceptance maximum free energy :     %12.3f / %12.3f" \
            % (self.V_HS, self.dGmin)
    #
    
    
    # find the minimum free energy
    def get_traces_top(self, hs, ndx, layer, flag_filter):
        #
        i_top = hs[0][0]; j_top = hs[0][1]
        ctp_top = hs[0][2]
        V_top   = hs[0][3]
        # print "V_top: %8.2f" % V_top
        
        if self.debug_traces_top:
            print "Enter get_traces_top(%d,%d)[%s] layer = %d" \
                % (i_top, j_top, self.smap.glink[i_top][j_top].lg[0].motif[0].get_ctp(), layer)
            #print "Stop at 0 in get_traces_top"; sys.exit(0)
        #endif
        
        if not (self.hv[i_top][j_top] > 0.0):
            print "get_traces_top(): enthalpy function hv is zero"
            print "                  at end point (%d,%d)!" % (i_top,j_top)
        #
        
        if self.debug_traces_top:
            print "making reference threads"
        #
        
            
        # 1) do a full scan using the mloop routine and record
        # whatever structures it pulls out from the scan.
        list_M = []
        list_M = self.localscan_for_M(i_top, j_top, list_M, V_top)
        if self.debug_traces_top:
            self.show_list_M(list_M, "get_traces_top: after localscan_for_M:")
            print "list_M length: ", len(list_M)
            #print "Stop at 1 in get_traces_top"; sys.exit(0)
        #endif

        # 2) now sort list_M so that the FE is ordered and organized 
        list_M = self.ins_SortHotSpot(list_M)
        if self.debug_traces_top:
            self.show_list_M(list_M, "get_traces_top: after ins_SortHotSpot:")
            print "list_M length: ", len(list_M)
            #print "Stop at 2 in get_traces_top"; sys.exit(0)
        #
        if len(list_M) > 1:
            list_M = self.prune_list_M(list_M)
            if self.debug_traces_top:
                self.show_list_M(list_M, "get_traces_top: after prune_list_M:")
                print "Stop at 3 in get_traces_top"; sys.exit(0)
            #
        #
        if self.debug_traces_top:
            print "add_localscan_for_M: enter"
        #
        self.add_localscan_for_M(i_top, j_top, list_M)
        if self.debug_traces_top:
            print "get_traces_top -> add_localscan_for_M"
            a = DispLThread(self.calc.N)
            a.disp_LThread(self.lt)
            print "finished add_localscan_for_M"
            print "Stop at 4 in get_traces_top"; sys.exit(0)
        #endif
        
        # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
        # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
        # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
        
        # 190128: I think this part is unnecessary because the task
        # should have already been done when I ran prune_list_M().
        
        use_prune_lt = False
        if use_prune_lt:
            if self.debug_traces_top:
                print "get_traces_top -> prune_lt"
            #
            self.lt = self.prune_lt(self.lt)
            if self.debug_traces_top:
                print "finished get_traces_top -> prune_lt: -----------"
                a = DispLThread(self.calc.N)
                a.disp_LThread(self.lt)
                print "----------- :finished get_traces_top -> prune_lt"
                #print "Stop at 5 in get_traces_top"; sys.exit(0)
            #
            if self.debug_traces_top:
                print "add_localscan_for_M: exit"
                #print "Stop at 6 in get_traces_top"; sys.exit(0)
            #
        #
        
        # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        # 
    #
    
    # removes redundant solutions from the list obtained by
    # localscan_for_M(i, j). This is a ratehr simple approach that
    # works on similarity based on lengths, free energy, and other
    # obvious kinds of tests.
    def prune_lt(self, lt):
        debug_prune_lt = False # True # 
        uniq_lt = []
        last_dG = 1000
        last_n  = -1
        for ltk in lt:
            #
            V = float( int(1000*ltk.dG))/1000.0
            if debug_prune_lt:
                print V, len(ltk.thread)
            if not (last_dG == V and last_n == len(ltk.thread)):
                last_dG = V; last_n = len(ltk.thread)
                uniq_lt += [ltk]
            #
        #
        if debug_prune_lt:
            print "finished prune_lt()"
            print "Stop at 0 in prune_lt"; sys.exit(0)
        #
        return uniq_lt
    #
    
    def add_localscan_for_M(self, i_top, j_top, list_M):
        # debug_add_localscan_for_M = self.debug_get_traces
        debug_add_localscan_for_M = False # True # 
        layer = 0
        k_ref = 0
        if debug_add_localscan_for_M:
            print "Enter: add_localscan_for_M(ij_top = (%d,%d))" % (i_top, j_top)
            self.show_list_M(list_M, "add_localscan_for_M:")
            #print "Stop at 0 in add_localscan_for_M"; sys.exit(0)
        #
        
        count = 0
        for lM in list_M:
            
            # Search through the constructed list of candidates
            
            # Here, we obviously have to add not just Link[0].motif[0]
            # anymore, but all items of the list. So this is the next
            # step in the process.
            
            i = lM[0]; j = lM[1]
            ctp      = lM[2][0].ctp
            btp      = lM[2][0].btp
            base     = lM[2][0].get_base()
            branches = lM[2][0].get_branches()
            pk       = lM[2][0].get_pks()
            wyspa    = lM[2][0].get_wyspa()
            V = lM[3]
            n = len(branches)
            
            # find the boundaries between the branches.
            ijo = branches[0]
            ijf = branches[n-1]
            ib = ijo[0]
            jb = ijf[1]
            
            
            if debug_add_localscan_for_M:
                print "lM: ", i, j, lM[2][0].show_Motif()
                print "new lM: n(%2d), ijt(%2d,%2d), ijr(%2d,%2d)[%s][%8.2f]" \
                    % (n, i_top, j_top, i,j, ctp, V), base, branches, pk, wyspa
                print "n(%2d), ijo = %s -> ib(%2d), ijf = %s -> jb(%2d)" % (n, ijo, ib, ijf, jb)
                
                for m in range(0,len(self.smap.glink[i][j].lg)):
                    print self.smap.glink[i][j].lg[m].motif[0].show_Motif()
                #
            #
            
            k_ref = 0
            for m in range(0,len(self.smap.glink[i][j].lg)):
                ctpx = self.smap.glink[i][j].lg[m].motif[0].get_ctp()
                btpx = self.smap.glink[i][j].lg[m].motif[0].get_btp()
                if  ctpx == ctp and btpx  == btp:
                    k_ref = m
                    break
                #
            #
            if debug_add_localscan_for_M:
                print "----------"
                print self.smap.glink[i][j].lg[k_ref].motif[0].show_Motif()
                
                if STOP_AT_FIND and ctp == 'S':
                    print "add_localscan_for_M(): planned exit when ctp = 'S'"
                    print "Stop at 1 in add_localscan_for_M"; sys.exit(0)
                #
                
                if ctp == 'P':
                    print "add_localscan_for_M(): planned exit when ctp = 'P'"
                    # print "Stop at 2 in add_localscan_for_M"; sys.exit(0)
                #
            #
            
            # now build a P-loop between [i_top, j_top] (inclusive)
            # based on search result obtained at bondaries [ib,jb]
            if debug_add_localscan_for_M:
                print "searching P: "
                if debug_add_localscan_for_M:
                    print "so-called searching P start: -----------------"
                    a = DispLThread(self.calc.N)
                    a.disp_LThread(self.lt)
                    print "----------------- :so-called searching P start"
                    #print "Stop at 3 in add_localscan_for_M"; sys.exit(0)
                    # xxx!!!
                #endif
            #
            self.lt += [LThread(self.N)] # create a new thread
            ndx = len(self.lt) - 1 # assign the particular LThread()
            if ctp == 'B':
                ctpB = self.smap.glink[i][j].lg[0].motif[0].get_ctp()
                btpB = self.smap.glink[i][j].lg[0].motif[0].get_btp()
                VpB = self.smap.glink[i][j].lg[0].Vij
                if debug_add_localscan_for_M:
                    print "B:              ij_B(%3d,%3d)[%s][%5s][%8.2f]" % (i, j, ctpB, btpB, VpB)  
                #
                self.get_traces(i, j, ctpB, VpB, ndx, layer+1)
                
            if ctp == 'S':
                # case for a stem (antiparallel or parallel)
                self.make_Stemtrace(i, j, ndx, k_ref, layer + 1, ctp, btp, "add_localscan_for_M")
            elif ctp == 'K':
                # case for a pseudoknot
                self.make_PKtrace(i, j, ndx, k_ref, layer + 1, ctp, btp, "add_localscan_for_M")
            elif ctp == 'W':
                # case for an island
                self.make_Islandtrace(i, j, ndx, k_ref, layer + 1, ctp, btp, "add_localscan_for_M")
            else:
                for lk in base:
                    # new ij
                    i_P = lk[0]; j_P = lk[1]
                    
                    ctpP = self.smap.glink[i_P][j_P].lg[0].motif[0].get_ctp()
                    btpP = self.smap.glink[i_P][j_P].lg[0].motif[0].get_btp()
                    VpP = self.smap.glink[i_P][j_P].lg[0].Vij
                    if debug_add_localscan_for_M:
                        print "lk: %10s, ij_P(%3d,%3d)[%s][%5s][%8.2f]" \
                            % (lk, i_P, j_P, ctpP, btpP, VpP)  
                    #
                    self.get_traces(i_P, j_P, ctpP, VpP, ndx, layer+1)
                #
            #
            if debug_add_localscan_for_M:
                print "add_localscan_for_M (after computing structure %d): -----------------" \
                    % count
                a = DispLThread(self.calc.N)
                a.disp_LThread(self.lt)
                print "----------------- :add_localscan_for_M (after computing structure %d)" \
                    % count
                if count == 2:
                    print "Stop at 4 in add_localscan_for_M"; sys.exit(0)
                #
                # xxx!!!190126
            #endif
            count += 1
        #
        self.lt = self.mergeSortTraces(self.lt)
        
    #
    
    
    
    # this is used to sort entries after obtaining the free energy for
    # various types of interactions ('B', 'I', 'M', 'J', or 'P'). It
    # is also used to sort the distribution of free energies when
    # obtaining the suboptimal structures
    def mergeSortTraces(self, alist):
        # print("Splitting ",alist)
        if len(alist)>1:
            mid = len(alist)//2
            lefthalf = alist[:mid]
            righthalf = alist[mid:]
            
            self.mergeSortTraces(lefthalf)
            self.mergeSortTraces(righthalf)
            
            i=0
            j=0
            k=0
            while i < len(lefthalf) and j < len(righthalf):
                if lefthalf[i].dG < righthalf[j].dG:
                    alist[k]=lefthalf[i]
                    i=i+1
                else:
                    alist[k]=righthalf[j]
                    j=j+1
                k=k+1
            #
            while i < len(lefthalf):
                alist[k]=lefthalf[i]
                i=i+1
                k=k+1
            #
            while j < len(righthalf):
                alist[k]=righthalf[j]
                j=j+1
                k=k+1
            # print("Merging ",alist)
        #
        return alist
    #    
    
    def make_Islandtrace(self, i, j, ndx, k_ref, layer, ctp, btp, iprog = "get_trace"):
        flag_debug = False # True # 
        if flag_debug:
            print "make_Islandtrace"
        #
        
        self.lt[ndx].add_lnode((i,j), 0.0, 'W', 'bgn')
        
        
        # the specific islands
        wyspa = self.smap.glink[i][j].lg[0].motif[0].get_wyspa()
        if flag_debug:
            print "wyspa: ", wyspa
        #
        for wk in wyspa:
            # some repetition, but who cares!
            i_W = wk[0]; j_W = wk[1]
            ###  hmm!!!, maybe this should be prefaced with
            
            #"if self.hv[i_W][j_W] > 0.0:" # <<<<<< ???
            
            # The function is not so accurate and doesn't zero!
            dGijW = self.calc.fe.calc_dG(i_W, j_W, self.hv[i_W][j_W], self.T, "make_Islandtrace")
            self.lt[ndx].add_lnode((i_W,j_W), dGijW, 'W', 'wyspa')
        #
        
        # the junctions between the islands
        jn_w = self.smap.glink[i][j].lg[0].motif[0].get_branches()
        
        if flag_debug:
            print "jn_w: ", jn_w
        #
        
        for lk in jn_w:
            # whether 1 or 100 elements, they can be processed this way.
            i_W = lk[0]; j_W = lk[1]     # new ij of this subdomain
            if len(self.smap.glink[i_W][j_W].lg) > 0:
                ctpW = self.smap.glink[i_W][j_W].lg[0].motif[0].get_ctp()
                VpW  = self.smap.glink[i_W][j_W].lg[0].Vij
                if flag_debug:
                    print "ij_W = (%d,%d)[%s], ndx(%d), layer(%d)" % (i_W, j_W, ctpW, ndx, layer)
                #
                if not (i == i_W and j == j_W):
                    # sometimes, there is nothing inside
                    self.get_traces(i_W, j_W, ctpW, VpW, ndx, layer+1)
                #
            #
            # if smap is empty, then this is all there is.
        #
        self.lt[ndx].add_lnode((i,j), 0.0, 'W', 'end')

        if flag_debug:
            n_wyspa = 5
            if len(wyspa) > n_wyspa:
                print "make_Islandtrace: found islands > %d" % n_wyspa
                print "Island results: from %s" % iprog
                for thr in self.lt[ndx].thread:
                    print thr.disp_lnode()
                if STOP_AT_FIND:
                    print "stop at end of make_Islandtrace"
                    sys.exit(0)
                #
            #
        #
        # #####12345
    #
    
    
    # need to add function for 'W'
    def make_PKtrace(self, i, j, ndx, k_ref, layer, ctp, btp, iprog = "get_trace"):
        flag_debug = False # True # 
        if flag_debug:
            print "make_PKtrace"
        #
        
        # root stem
        self.lt[ndx].add_lnode((i,j), 0.0, 'K', 'bgn')
        v = self.smap.glink[i][j].lg[k_ref].motif[0].get_base()
        i_R  = v[0][0]; j_R = v[0][1]
        ctpR = self.smap.glink[i_R][j_R].lg[k_ref].motif[0].get_ctp()
        btpR = self.smap.glink[i_R][j_R].lg[k_ref].motif[0].get_btp()
        VpR  = self.smap.glink[i_R][j_R].lg[k_ref].Vij
        
        if ctpR == 'I' or ctpR == 'M':
            dGij = self.calc.fe.calc_dG(i_R, j_R, self.hv[i_R][j_R], self.T, "make_PKtrace")
            self.lt[ndx].add_lnode((i_R,j_R), dGij, ctpR, btpR)
        #
        
        branching = self.smap.glink[i_R][j_R].lg[k_ref].motif[0].get_branches()
        if flag_debug:
            print "branching = ", branching
        for lk in branching:
            # whether 1 or 100 elements, they can be processed this way.
            i_M = lk[0]; j_M = lk[1]     # new ij of this subdomain
            if len(self.smap.glink[i_M][j_M].lg) > 0:
                ctpM = self.smap.glink[i_M][j_M].lg[0].motif[0].get_ctp()
                VpM  = self.smap.glink[i_M][j_M].lg[0].Vij
                if flag_debug:
                    print "ij_M = (%d,%d)[%s], ndx(%d), layer(%d)" % (i_M, j_M, ctpM, ndx, layer)
                #
                self.get_traces(i_M, j_M, ctpM, VpM, ndx, layer+1)
            else:
                if flag_debug:
                    print "ij_M = (%d,%d)[%s], ndx(%d), layer(%d): skipped!" % (i_M, j_M, ctpM, ndx, layer)
                #
            #
        #
            
        pks  = self.smap.glink[i][j].lg[k_ref].motif[0].get_pks()
        if flag_debug:
            print "pks: ", pks
        #
        
        btpx = ''
        for x in pks:
            btpx = x[2]
            self.lt[ndx].add_lnode((x[0][0],x[0][1]), x[1], 'l', x[2])
        #
        self.lt[ndx].add_lnode((i,j), 0.0, 'K', 'end')
        
        if flag_debug:
            print "pk results: from %s" % iprog
            for thr in self.lt[ndx].thread:
                print thr.disp_lnode()
            if btpx == 'sp':
                print "parallel"
                #print "stop at 0 in make_PKtrace"; sys.exit(0)
            #
        # #####12345pk
    #
    
    
    # need to add function for 'W'
    def make_Stemtrace(self, i, j, ndx, kref_t, layer, ctp, btp, iprog = "get_trace"):
        
        # ndx is the index for the particular thread of
        # self.lt[ndx].thread
        
        # kref_t the particular reference to the tail section of the
        # stem
        
        flag_debug = False # True # 
        if flag_debug:
            print "make_Stemtrace(called from %s)" % iprog
        #
        
        # root stem
        flag_pp = False
        self.lt[ndx].add_lnode((i,j), 0.0, 'S', 'bgn')
        ctpS = self.smap.glink[i][j].lg[kref_t].motif[0].get_ctp()
        btpS = self.smap.glink[i][j].lg[kref_t].motif[0].get_btp()
        VpS  = self.smap.glink[i][j].lg[kref_t].Vij
        branching = self.smap.glink[i][j].lg[kref_t].motif[0].get_branches()
        
        # establish the head and tail of the stem
        stmlen = len(branching) - 1
        i_t = i;          j_t = j
        i_h = i + stmlen; j_h = j - stmlen
        ihh = i_h + 1;    jhh = j_h - 1 # for parallel stems
        
        if flag_debug:
            print "antiparallel"
            print "ij_t: ", i_t, j_t
            print "ij_h: ", i_h, j_h
            print "motif(ij_t): ", self.smap.glink[i_t][j_t].lg[kref_t].motif[0].show_Motif()
        #
        
        # is it an antiparallel stem or a parallel stem?
        if not (btpS == 'sa' or btpS == 'c'):
            flag_pp = True
        
        kref_h = 0
        if not flag_pp:
            if flag_debug:
                print "antiparallel"
            #
            
            # scan up the antiparallel stem from the tail up to just
            # before the head
            
            #for k in range(0, len(branching)-1): # # before !!!xxx190116 !!!xxx190116
            for k in range(0, len(branching)-1):
                # whether 1 or 100 elements, they can be processed this way.
                i_S = branching[k][0]; j_S = branching[k][1]     # new ij of this subdomain
                # print "stem ij_S: ", i_S, j_S
                ctpS = 'S'
                bondtype = self.calc.fe.btype[i_S][j_S].btp
                # 190524 was self.calc.get_bondtype(self.hv[i_S][j_S])
                if bondtype == 's':
                    bondtype += 'a'
                else:
                    if not bondtype == 'sa':
                        # this decision must be enforced in the case of a stem!
                        bondtype = 'c'
                #
                dGij = self.calc.fe.calc_dG(i_S, j_S, self.hv[i_S][j_S], self.T, "make_Stemtrace 1")
                self.lt[ndx].add_lnode((i_S,j_S), dGij, 'S', bondtype)
            #
            
            # now add the head
            for m in range(0, len(self.smap.glink[i_h][j_h].lg)):
                # print "motif(ij_h): ", self.smap.glink[i_h][j_h].lg[m].motif[0].show_Motif()
                ctp = self.smap.glink[i_h][j_h].lg[m].motif[0].get_ctp()
                if ctp == 'W':
                    # chr11_111099793_111320552_res5kb.heat satisfied this condition
                    print "aa ij_t (tail): ", i_t, j_t
                    print "aa ij_h (head): ", i_h, j_h
                    print "found an island at the head of antiparallel stem."
                    if STOP_AT_FIND:
                        print "xx"
                        print "stop at 0 in make_Stemtrace"; sys.exit(0)
                    #
                #
                if ctp == 'B' or ctp == 'I' or ctp == 'M' or ctp == 'S' or ctp == 'W':
                    kref_h = m
                    break
                #
            #
            # now add the terminus (the head of the stem)
            ctpS = self.smap.glink[i_h][j_h].lg[kref_h].motif[0].get_ctp()
            btpS = self.smap.glink[i_h][j_h].lg[kref_h].motif[0].get_btp()
            VpS  = self.smap.glink[i_h][j_h].lg[kref_h].Vij
            
            # if you have antiparallel stem condition. All you need to
            # do is just jump to the next motif like in the old days
            # when this program started out.
            if flag_debug:
                print ctpS, btpS, VpS
                print "motif(ij_h): ", self.smap.glink[i_h][j_h].lg[kref_h].motif[0].show_Motif()
                print "aa: ij_h = (%d,%d)[%s], ndx(%d), layer(%d)" % (i_h, j_h, ctpS, ndx, layer)
            #
            
            self.get_traces(i_h, j_h, ctpS, VpS, ndx, layer+1)
            # print "xxx"; print "stop at 1 in make_Stemtrace"; sys.exit(0)
        else:
            if flag_debug:
                print "parallel"
            #
            
            # scan up the antiparallel stem from the tail up the head
            # (the whole stem without a cut off at the top)
            for k in range(0, len(branching)):
                # whether 1 or 100 elements, they can be processed this way.
                i_S = branching[k][0]; j_S = branching[k][1]     # new ij of this subdomain
                if flag_debug:
                    print "stem ij_S: ", i_S, j_S
                #
                ctpS = 'S'
                bondtype = self.calc.fe.btype[i_S][j_S].btp
                # 190524 was self.calc.get_bondtype(self.hv[i_S][j_S])
                if bondtype == 's':
                    bondtype += 'p'
                else:
                    if not bondtype == 'sp':
                        # this decision must be enforced in the case of a stem!
                        bondtype = 't'
                #
                dGij = self.calc.fe.calc_dG(i_S, j_S, self.hv[i_S][j_S], self.T, "make_Stemtrace 2")
                self.lt[ndx].add_lnode((i_S,j_S), dGij, 'S', bondtype)
            #
            
            # now add the head at (i_h + 1, j_h - 1)
            if len(self.smap.glink[ihh][jhh].lg) > 0:
                for m in range(0, len(self.smap.glink[ihh][jhh].lg)):
                    if flag_debug:
                        print "motif(ij_h): ", self.smap.glink[ihh][jhh].lg[m].motif[0].show_Motif()
                    #
                    ctp = self.smap.glink[ihh][jhh].lg[m].motif[0].get_ctp()
                    if ctp == 'J' or ctp == 'P' or ctp == 'K':
                        kref_h = m
                        break
                    #
            else:
                kref_h = -1
            #
            if flag_debug:
                if kref_h >= 0:
                    print "motif(ij_h): ", self.smap.glink[ihh][jhh].lg[kref_h].motif[0].show_Motif()
                    # technically, should verify that these exist, but
                    # presumably, we have already taken care of that
                    # in make_StemMotif().
                    ctph = self.smap.glink[ihh][jhh].lg[kref_h].motif[0].get_ctp()
                    btph = self.smap.glink[ihh][jhh].lg[kref_h].motif[0].get_btp()
                    Vph  = self.smap.glink[ihh][jhh].lg[kref_h].Vij
                    if flag_debug:
                        print "pp: ij_t = (%d,%d)" % (i_t, j_t)
                        print "pp: ij_h = (%d,%d)" % (i_h, j_h)
                        print "pp: ijhh = (%d,%d)[%s], ndx(%d), layer(%d)" % (ihh, jhh, ctph, ndx, layer)
                        if  i_t == 0 and j_t == 14:
                            print "planned exit when evaluating a parallel stem"
                            print "stop at 2 in make_Stemtrace"; sys.exit(0)
                        #
                    self.get_traces(ihh, jhh, ctph, Vph, ndx, layer+1)
                else:
                    if flag_debug:
                        print "motif(ij_h):  empty" 
        #    
        #print "stop at 3 in make_Stemtrace"; sys.exit(0)
        
        self.lt[ndx].add_lnode((i,j), 0.0, 'S', 'end')
        
        if flag_debug:
            print "stem results: from %s" % iprog
            for thr in self.lt[ndx].thread:
                print thr.disp_lnode()
            #
            
            #if flag_pp:
            #    print "stop when pp"
            #    print "stop at 3 in make_Stemtrace"; sys.exit(0)
        # #####12345stem
    #
    
    
    # find the minimum free energy
    def get_traces(self, i, j, ctp, V, ndx, layer):
        
        # ndx is the index for the particular thread of
        # self.lt[ndx].thread
        
        debug_get_traces = self.debug_get_traces
        if debug_get_traces and layer == 0:
            print "Enter: get_traces{ij(%d,%d), ctp(%s), V(%8.2f), ndx(%d), layer(%d)" \
                % (i, j, ctp, V, ndx, layer)
        #
        if layer > self.maxlayers:
            print "layer(%d) > maxlayer(%d)" % (layer, self.maxlayers)
            print "ERROR: something wrong in the recursion of get_traces!"
            print "stop at 0 in get_traces"; sys.exit(1)
        #
        
        
        #####
        flag_found = False
        kref_t = 0
        for m in range(0, len(self.smap.glink[i][j].lg)):
            ctpx   = self.smap.glink[i][j].lg[m].motif[0].get_ctp()
            Vx    = self.smap.glink[i][j].lg[m].Vij
            if debug_get_traces:
                basex = self.smap.glink[i][j].lg[m].motif[0].get_base()
                branchingx = self.smap.glink[i][j].lg[m].motif[0].get_branches()
                pkx   = self.smap.glink[i][j].lg[m].motif[0].get_pks()
                print "layer = %3d: %2d '%s' == ref(%s)?,  %8.3f <= ref(%8.3f)?" \
                    % (layer, m, ctpx, ctp, Vx, V),  basex, branchingx, pkx
                print "ctpx = ctp? ",(ctpx == ctp), ", Vx >= V? ", (Vx >= V)
            #
            if ctpx == ctp and Vx >= V:
                kref_t = m
                flag_found = True
                # before, everything was already unique, but now we
                # can two types of PKs and two types of stems due to
                # the parallel/antiparallel nature of chromatin (which
                # does not care so much about direction). So I have to
                # break after the first instance of this term being
                # satisfied.
                break
            #####
        if not flag_found:
            print "ERROR(Trace.get_traces()): didn't find a match for ", ctp
            print "stop at 1 in get_traces"; sys.exit(1)
        #
        
        
        # have to assemble the connection or termination. 
        if ctp == 'B' or ctp == 'I' or ctp == 'M':
            # these are the only structures that form links
            dGij = self.calc.fe.calc_dG(i, j, self.hv[i][j], self.T, "get_traces")
            # ctp is already defined, it should correspond to kref_t
            btp  = self.smap.glink[i][j].lg[kref_t].motif[0].get_btp()
            self.lt[ndx].add_lnode((i,j), dGij, ctp, btp)
            if debug_get_traces:
                print "get_traces(): (%2d,%2d)[%8.3f]" \
                    % (i, j, V), [(i,j)]
                
                print "bim-----------------"
                a = DispLThread(self.calc.N)
                a.disp_LThread(self.lt)
                # xxx!!!
                print "-----------------bim"
                #if not ctp == 'B': sys.exit(0)
                #print "stop at 2 in get_traces"; sys.exit(0)
                # xxx!!!
            #endif
            #
            # pseudoknot connections have to be handled in the
            # processing.
        #
        s = self.calc.fe.space(3*layer)
        if debug_get_traces:
            branching = self.smap.glink[i][j].lg[kref_t].motif[0].get_branches()
            print "%s[%s](%3d, %3d)[%8.3f][ndx = %d, layer = %d]: " \
                % (s, ctp, i,j, V, ndx, layer), branching
        #
        if ctp == 'M':
            if debug_get_traces:
                print "found a case of %s" % ctp
                print "%s------------------" % s
            #
            branching = self.smap.glink[i][j].lg[kref_t].motif[0].get_branches()
            for lk in branching:
                # new ij
                i_M = lk[0]; j_M = lk[1]
                ctpM = self.smap.glink[i_M][j_M].lg[0].motif[0].get_ctp()
                VpM = self.smap.glink[i_M][j_M].lg[0].Vij
                self.get_traces(i_M, j_M, ctpM, VpM, ndx, layer+1)
            #
            #
        elif ctp == 'P':
            if debug_get_traces:
                print "found a case of %s" % ctp
                print "%s------------------" % s
            #
            branching = self.smap.glink[i][j].lg[kref_t].motif[0].get_branches()
            for lk in branching:
                i_P = lk[0]; j_P = lk[1]
                ctpP = self.smap.glink[i_P][j_P].lg[0].motif[0].get_ctp()
                VpP = self.smap.glink[i_P][j_P].lg[0].Vij
                self.get_traces(i_P, j_P, ctpP, VpP, ndx, layer+1)
            #
            #
        elif ctp == 'I':
            if debug_get_traces:
                print "found a case of %s" % ctp
            #
            jn  = self.smap.glink[i][j].lg[kref_t].motif[0].get_branches()
            # new ij!!!
            i_I = jn[0][0]; j_I = jn[0][1]
            VpI, cpI, ctpI, pks, wyspa = self.calc.fe.filter_BIKMSW_from_glink(i_I,j_I)
            self.get_traces(i_I, j_I, ctpI, VpI, ndx, layer+1)
            
            if debug_get_traces:
                self.calc.show_smap_xy(i,j)
                print "Stop at ctp == I"
                print "q1-----------------"
                a = DispLThread(self.calc.N)
                a.disp_LThread(self.lt)
                # xxx!!!
                print "-----------------q1"
                # if i == 0 and j == 32: sys.exit(0)
                #print "stop at 3 in get_traces"; sys.exit(0) 
                # xxx!!!
            #endif
            
            #
        elif  ctp == 'J':
            if debug_get_traces:
                print "found a case of J"
            #
            jn  = self.smap.glink[i][j].lg[kref_t].motif[0].get_branches()
            # new ij!!!
            i_J = jn[0][0]; j_J = jn[0][1]
            # (i,j) -> (i_J, j_J)
            VpJ, cpJ, ctpJ, pks, wyspa = self.calc.filter_BIKMSW_from_glink(i_J,j_J)
            
            if VpJ == 1000.0:
                print "ERROR(get_traces): could not find a linkage at (%d,%d)" % (i_J, j_J)
                sys.exit(1)
            #
            self.get_traces(i_J, j_J, ctpJ, VpJ, ndx, layer+1)
            #
            
        elif ctp == 'S':
            
            if debug_get_traces:
                print "get_traces: S:"
            #
            
            btp  = self.smap.glink[i][j].lg[kref_t].motif[0].get_btp()
            self.make_Stemtrace(i, j, ndx, kref_t, layer+1, ctp, btp, "get_trace")
            
            #
        elif ctp == 'K':
            
            if debug_get_traces:
                print "get_traces: K:"
            #
            
            btp  = self.smap.glink[i][j].lg[kref_t].motif[0].get_btp()
            # print kref_t
            # print self.smap.glink[i][j].lg[kref_t].motif[0].show_Motif()
            
            #print "stop at 4 in get_traces"; sys.exit(0)
            self.make_PKtrace(i, j, ndx, kref_t, layer+1, ctp, btp, "get_trace")
            
            #
        elif ctp == 'W':
            
            if debug_get_traces:
                print "get_traces: W:"
            #
            btp  = self.smap.glink[i][j].lg[kref_t].motif[0].get_btp()
            self.make_Islandtrace(i, j, ndx, kref_t, layer+1, ctp, btp, "get_trace")
            if debug_get_traces:
            
                self.calc.show_smap_xy(i,j)
                print "Stop at ctp == W"
                print "v1-----------------"
                a = DispLThread(self.calc.N)
                a.disp_LThread(self.lt)
                # xxx!!!
                print "-----------------v1"
                # if i == 0 and j == 32: sys.exit(0)
                # xxx!!!
                #print "stop at 5 in get_traces"; sys.exit(0)
            #
            
            #
        else:
            if debug_get_traces:
                print "%s---" % s
            #
        #
    #
    
    
#




def main(cl):
    print "presently, this call doesn't really do anything..."
    print cl
    # presently doesn't have any tests
#




if __name__ == '__main__':
    # running the program
    main(sys.argv)

