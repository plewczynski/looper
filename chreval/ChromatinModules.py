#!/usr/bin/env python

"""@@@

Main Program:  ChromatinModules.py 

Classes:       SetUpBranchEntropy
               BranchEntropy

Author:        Wayne Dawson
creation date: mostly 2016 and up to March 2017 (separated from chreval 180709)
last update:   190119
version:       0


Purpose:

This module is intended to function specifically for chromatin. 

This unit is intended to be the central functioning module that can be
exchanged for other modules with minor alterations and allow
calculation of chromatin, RNA or proteins. The main objective is that
we can use the same main calculation machinery (Calculate) for
chromatin, RNA and proteins and we only have to interchange these
simple modules ChromatinModules for RNAmodules or ProteinModules to
accomplish this task.


Comments:

This module has been tested quite a bit now, so it is at least
basically functioning now.

"""



from math import exp
import sys
import os
import string

# main tool objects
from GetOpts import GetOpts

# Free energy parameters
from FreeEnergy   import FreeEnergy
from HeatMapTools import HeatMapTools

# Motif object representation
from Motif import Motif # a core object of these building programs
from Motif import Link
from Motif import LGroup
from Motif import Map # main map for all the FE and structures

from LoopRecords import Branch
from LoopRecords import MBLptr

# Other objects and tools
from FileTools   import getHeadExt
from ChPair      import ChPairData
from ChPair      import LThread2ChPair
from SimRNATools import SimRNAData

# for Vienna object representation
from Vienna      import Vstruct

# ################################################################
# ######################  Global constants  ######################
# ################################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# A hardwired infinity. Somewhat hate this, but it seems there is no
# way around this darn thing.
from Constants import INFINITY
from Constants import kB # [kcal/molK] (Boltzmann constant)

# coarse-grained resolution factor
from ChrConstants import seg_len
# for entropy calculations
from ChrConstants import xi # [bps]
from ChrConstants import lmbd # [bps]
from ChrConstants import gmm # dimensionless but related to D/2        
# constants: weights for the enthalpy terms 
from ChrConstants import febase # [kcal/mol]
from ChrConstants import feshift # (dimensionless, usually = 1)


# stem parameters
from ChrConstants import minStemLen      # [bp], minimum stem length (default is 1)
from ChrConstants import max_bp_gap
# loop parameters 
from ChrConstants import minLoopLen      # [nt], minimum loop length (default is 1)
from ChrConstants import dGMI_threshold  # [kcal/mol] threshold FE for M-/I-loops
# pseudoknot parameters
from ChrConstants import pk_scan_ahead   # [nt] the "hot lead" length for pk
from ChrConstants import dGpk_threshold  # [kcal/mol] threshold FE for pks

from ChrConstants import set_dangles
from ChrConstants import dG_range 

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ################################################################
# ################################################################


# ################################################################
# ###########  global functions/objects and constants  ###########
# ################################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# 1. control parameters for the program

PROGRAM      = "ChromModules.py"  # name of the program


# debugging options 
TEST    = False # Debugging in main()

# 2. special local global function

# debug cle_HloopE
DEBUG_cle_HloopE        = False # True # 
# debug trim_53pStems
DEBUG_T53S              = False # True # 
# debug cle_IloopE and tags 
DEBUG_CILE              = False # True # 
# debug cle_MloopEV
DEBUG_CMLEV             = False # True # 

# debugging traceback_mFE (used by pseudoknot building tools)
DEBUG_traceback_mFE     = False # True #


CHECK_ALL = False # True # 
# if you want to debug all parts simultaneously, then set CHECK_ALL to True.
def debug_all(fast_track):
    global DEBUG_cle_HloopE     
    global DEBUG_T53S           
    global DEBUG_CILE           
    global DEBUG_CMLEV          
    global DEBUG_traceback_mFE  
    
    DEBUG_cle_HloopE        = fast_track
    DEBUG_T53S              = fast_track
    DEBUG_CILE              = fast_track
    DEBUG_CMLEV             = fast_track
    DEBUG_traceback_mFE     = fast_track
#    

# debugging option (rarely used)
STOP_AT_FIND = False # True # 
# stops program when encouters condition


# This is not used so much now that I have introduced GetOpts.py,
# but it may still be useful at the very beginning of the program and
# possibly in other parts.

def usage():
    print "USAGE: %s {test0,test1}" % PROGRAM
    print "       default test0"
#

# 3. Other categories


# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ################################################################
# ################################################################


class SetUpBranchEntropy(object):
    def __init__(self, cseq = "cccccccccccccc"):
        self.source = "SetUpBranchEntropy"
        # This is a combination of SetUpFreeEnergy and
        # GenerateHeatMapTools
        
        self.f_heatmap = ["/home/dawson/python/chromatin/tests/test4.heat"]
        self.cseq      = cseq
        
        # input variables
        self.T        =  310.0
        self.N        = -1
        
        # constants: in entropy evaluation
        self.kB       = kB
        self.xi       = xi
        self.lmbd     = lmbd
        self.gmm      = gmm
        self.seg_len  = seg_len
        self.w        = self.seg_len*self.xi/(self.lmbd)**2
        
        
        
        # secondary structure stem parameters
        self.minStemLen     = minStemLen
        self.max_bp_gap     = max_bp_gap 
                
        
        # secondary structure loop parameters
        self.minLoopLen     = minLoopLen
        # minimum loop length (for chromatin it is 1, for RNA = 3)
        self.dGMI_threshold = dGMI_threshold # [kcal/mol]
        """threshold for MBL stability"""
        
        # pseudoknot parameters
        self.minPKloop      = 1              # minimum PK loop (2 nt)
        self.add_PK         = True           # include PK search
        self.scan_ahead     = pk_scan_ahead  # default 10
        self.dGpk_threshold = dGpk_threshold # default 2 (needs better definition!!!)
        
        self.Mg_binding     = False        
        self.dangles        = set_dangles
        self.dGrange        = dG_range
        
        # constants: weights for the enthalpy of binding
        self.dHbase  = febase
        self.dHshift = feshift
        
        # PET cluster weight
        self.add_PET_wt          = False
        self.add_PET_wt_to_edges = False
        self.PETwt               = 100.0
        
        # weights for selecting out PET clusters 
        self.CTCF_scale          = 100.0
        self.ctcf_tthresh        = 0.20*self.CTCF_scale
        self.ctcf_cthresh        = 0.40*self.CTCF_scale
        self.minLoopLen          = 1
        
        # This started when I was confronted with the somewhat strange
        # data I got from Nenski that had counts in almost every bin
        # and very large numbers.
        
        self.pssbl_ctcf          = {}
        self.edge_ctcf           = {}
        self.from_Nenski         = False
        
        # this re-scales the data by some fraction
        self.rescale_wt          = 1.0
        
        # other handling matters
        self.allowed_extns       = ["heat", "data"]
        self.set_GetOpts = True
    #
#
        

class BranchEntropy(FreeEnergy, HeatMapTools):
    def __init__(self, cl = SetUpBranchEntropy()):
        FreeEnergy.__init__(self, cl)   # inherit FreeEnergy
        HeatMapTools.__init__(self, cl) # inherit HeatMapTools
        
        if not cl.set_GetOpts:
            print "ERROR: options are not properly configured"
            sys.exit(1)
        #
        
        self.debug_BranchEntropy =  False
        """@
        
        181224: I think there is only one allowed file here, in the
        command line argument it expects only 1 item. However, this is
        one concern here. It seems like it would be possible to
        process more than one file in this way. If I remember
        correctly, multiple files were handled by a shell program.
        
        """
        self.flnm = ''
        
        
        if cl.source == "GetOpts":
            """@
            
            The "sources" are either SetUpBranchEntropy or GetOpts,
            where SetUp is used for testing. Presumably, SetUp has
            some different tests.
            
            This is where the program reads in the relevant
            information (the heatmap file, or presumably a sequence
            file, or whatever.
            
            """
            
            # assign the above variables: hv, ctcf_setv, N, btype,
            # all_ctcf. 
            self.flnm = cl.f_heatmap[0]
            self.N = self.assign_btypes(self.flnm)
            
        else:
            self.flnm = cl.f_heatmap[0]
            self.N = self.assign_btypes(self.flnm)
            
        #
        
        # build map layout
        self.smap       = Map(self.N)
        # build map layout
        # print "(3) N = ", self.N
        
        
        # allowed searches for ChromatinModules
        self.dSearchK    = { 'SB' : 0,
                             'IJ' : 1,
                             'MP' : 2 }
        
        self.dinvSearchK = {  0   : 'SB',
                              1   : 'IJ',
                              2   : 'MP' }
        
        self.dSearchTypes = { 'B' : 'SB',
                              'S' : 'SB',
                              'K' : 'SB', 
                              'R' : 'SB', 
                              'W' : 'SB', 
                              'I' : 'IJ',
                              'J' : 'IJ',
                              'M' : 'MP',
                              'P' : 'MP' }
        
        
        
        # Trace back variables
        self.opt_ss_seq = []
        self.maxlayers = self.N/2 + 10
        
        """@ 
        
        The maxlayers is used in traceback routines in both this
        package traceback_mFE() and in class Trace().
        
        In general, I found that maxlayers = 10 was
        sufficient. However, when using this program in very large
        sets of data, I discovered some cases where get_traces failed
        because maxlayer = N/6 was too small.
        
        Ultimately, it is not clear where the limit should be. In
        principle, there could be a solution where there are N/2
        bonds. Then, if you have 2 beads, N/2 = 1, maxlayers > 1
        really means you have a problem! Therefore, just to be safe, I
        set this to N/2 -- the absolute theoretical maximum.
        
        It seems kind of absurd that N/2 be required for a N > 100
        bead search, but because of the recursion, it is possible (in
        principle). The only issue I foresee is that I do not know the
        upper limit to the number of recursion steps before the
        program complains. In general, with multibranch looping (which
        is far more entropically favorable) this limit should never be
        necessary, but it is truly unclear how many levels are
        required, so this is what I am forced to conclude. Obviously,
        if the N/2 limit is exceeded, something is DEFINITELY wrong.
        
        """
        
        # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
        # -----------------------------------------------------
        self.iloop      = LGroup()
        self.V_iloop_mx = INFINITY
        self.n_iloop_mx = 21
        self.wt_iloop   = 0.20
        """@
        
        I am not really so sure that this method of search is
        particularly useful anymore. I tried using this when I was
        developing some speed up techniques with RNA structure
        prediction in vsfold6, and it appeared to be helpful in
        searching dsRNA structures. However, using it here, I found
        that it was not very effective. In general, the best speed up
        and most thorough search is generally derived from using
        searchForMBL, which retains both multiple branches and single
        branch results.

        """
        # -----------------------------------------------------
        # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        
    #
    
    
    def save_best_iloops(self, link):
        """@
        
        if the BIKMSW motif has the best FE of the list, then it is
        added to the list as a potential lookup item.
        
        The original purpose of this method was to provide a potential
        list of I-loops that would be good for testing new connection
        points at (i,j) because they were already used in the
        calculations. Since the new method works far different and
        doesn't really need this list, this probably should be
        decommissioned. I think the reason I consider leaving this
        concept around is because I will eventually upgrade things to
        double stranded RNA as well, where such a lookup list might
        prove useful. Therefore, I leave it around and still use it
        even though this would presently be applied only to single
        stranded RNA.
        
        """
        flag_debug = False
        n = len(self.iloop.lg)
        if flag_debug:
            loop = link.motif[0]
            print "save_best_iloops: dG(%2d,%2d)[%s] = %8.3f >= %8.3f, %d" \
                % (loop.i, loop.j, loop.get_ctp(), loop.Vij, self.V_iloop_mx, n)
        #
        if link.Vij < self.V_iloop_mx or n < self.n_iloop_mx:
            if flag_debug:
                print "saving"
            #
            self.iloop.add_link(link)
            self.smap.mergeSortLinks(self.iloop.lg) # sort the stack
            self.V_iloop_mx = self.wt_iloop * self.iloop.lg[0].Vij
            dv = LGroup()
            for k in range(min(n, self.n_iloop_mx-1), -1, -1):
                dv.add_link(self.iloop.lg[k])
            #
            self.iloop = dv
            if flag_debug:
                for dvk in self.iloop.lg:
                    print "(%2d,%2d) %8.3f" \
                        % (dvk.motif[0].base[0], dvk.motif[0].base[1], dvk.Vij)
                #
            #
        #        
    #
    
    
    
    # #######################################################
    # ##########  tools for checking the variables  #########
    # #######################################################
    # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    
    
    # this is mostly for viewing glink at a particular position (i,j)
    # when debugging the program.
    def show_smap_xy(self, i,j):
        print "ij(%d,%d):    " % (i, j)
        n = len(self.smap.glink[i][j].lg)
        print "lg:    ", n
        if n > 0:
            print "motif: ", len(self.smap.glink[i][j].lg[0].motif)
            for kx in range(0, len(self.smap.glink[i][j].lg)):
                print self.smap.glink[i][j].lg[kx].motif[0].show_Motif()
            #
        else:
            print "motif: not assigned"
        #
    #
    
    def show_smap_all(self, N, p, q):
        print "N from smap: %d, current pq = (%d,%d)" % (N, p, q)
        for q in range(0,N):
            for p in range(0, q-1):
                if len(self.smap.glink[p][q].lg) > 0:
                    n = len(self.smap.glink[p][q].lg)
                    nn = len(self.smap.glink[p][q].lg[0].motif)
                    print "smap.glink(%2d,%2d): " % (p, q), n, nn
                #
            #
        #
        sys.exit(0)
    #
    
    def stopWhenMatchFound(self, i, j, i_set, j_set, callpt = "not specified"):
        if i == i_set and j == j_set:
            print "planned stop when ij = (%d,%d):, called at [%s]" % (i, j, callpt)
            self.show_smap_xy(i,j)
            sys.exit(0)
        #
        
    #
    def stopWhenMatchFound_sfm(self, i, j, i_set, j_set, lmblh, callpt = "not specified"):
        if i == i_set and j == j_set:
            print "planned stop when ij = (%d,%d):, called at [%s]" % (i, j, callpt)
            print "smap:"
            self.show_smap_xy(i,j)
            print "lmblh:"
            self.show_lmblh(i, j, lmblh, callpt)
            sys.exit(0)
        #
        
    #
    
    def check_smap(self, i, j, i_set, j_set, flag_stop):
        if i == i_set and j == j_set:
            if not len(self.smap.glink[i][j].lg) > 0:
                print "ERROR: no values assigned for smap at (%2d,%2d)" % (i, j)
                sys.exit(1)
            print "smap: "
            s = ''
            for k in range(0, len(self.smap.glink[i][j].lg)):
                V  = self.smap.glink[i][j].lg[k].Vij
                for vvkn in self.smap.glink[i][j].lg[k].motif:
                    tp = vvkn.get_ctp() # conn type
                    bt = vvkn.get_btp() # bond type
                    bb = vvkn.get_base()
                    jn = vvkn.get_branches()
                    s += "[ij = (%2d, %2d)[%s_%s](%d): dG = %8.3f], %s -> %s\n" \
                         % (i, j, tp, bt, k, V, bb, jn)
            #
            print s
            
            show_iloop_data = False
            if len(self.iloop.lg) > 0 and show_iloop_data:
                print "iloop: "
                s = ''
                for ij_jn in self.iloop.lg:
                    V    = ij_jn.Vij
                    for ij_jnk in ij_jn.motif:
                        tp   = ij_jnk.get_ctp() # conn type
                        bt   = ij_jnk.get_btp() # bond type
                        jn   = ij_jnk.get_branches()
                        base = ij_jnk.get_base()
                        s += "[ij = (%2d, %2d)[%s_%s]: dG = %8.3f], %s -> %s\n" \
                             % (base[0][0], base[0][1], tp, bt, V, base, jn)
                print s
            #
            if flag_stop:
                sys.exit(0)
        #
        
    #


    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # #######################################################
    # ##########  tools for checking the variables  #########
    # #######################################################
    
    
    # ######################################################
    # ###############   loop display tools  ################
    # ######################################################
    # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    
    
    # (from RNAModules)
    def get_MBLptr(self, ph, qh, opt_iMBL = True, debug = False):
        """@
        
        Note(!!!!): Regardless of Whether opt_iMBL is True or False,
        get_MBLptr only calculates the pMBL part of the structure.
        
        The reason we must make the distinction about whether we have
        opt_iMBL True(or False) is because when it is an iMBL, the
        branching _must_ be taken from
        
        (ph+1, qh-1) <<< == the iMBL case
        
        not from 
        
        (ph, qh).    <<< == the pMBL case
        
        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        To add the base contribution at the base of the iMBL
        structure, to keep the form of the methods in line with
        cle_HloopE() and cle_IloopE(), this value is added when the
        cle_MloopEV() is called.
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        
        There is a little added overhead in doing things this way,
        but, to maintain the modularity of the program -- especially
        in using this both for sarabande.py and RThreads.py
        
        """
        # have to construct the unit to pass
        new_MBL    = MBLptr(ph, qh)
        phv = ph; qhv = qh
        if opt_iMBL:
            phv += 1; qhv -= 1
        #
        branchV = self.smap.glink[phv][qhv].lg[0].motif[0].get_branches()
        
        if len(branchV) == 0:
            ctp  = self.smap.glink[phv][qhv].lg[0].motif[0].get_ctp()
            print "ERROR(get_MBLptr()): called branching routines with "
            print "                     no branches! pqh = (%d,%d)[%s]"  % (ph, qh, ctp)
            sys.exit(1)
        #
        if opt_iMBL:
            if len(branchV) == 1:
                new_MBL.nm = 'I'
            else:
                new_MBL.nm = 'M'
            #
        else:
            if len(branchV) == 1:
                new_MBL.nm = 'J'
            else:
                new_MBL.nm = 'P'
            #
        #
        dG_pqhv = 0.0
        for l in range(0, len(branchV)): # sum the first part
            pv = branchV[l][0]; qv = branchV[l][1]
            dG_pqhv += self.smap.glink[pv][qv].lg[0].motif[0].get_Vij() 
            new_MBL.addBranch(pv, qv)
        #
        new_MBL.V = dG_pqhv
        
        if debug:
            print "get_MBLptr(%2d,%2d):" % (ph, qh)
            if opt_iMBL:
                print "stored dG(%2d,%2d): %8.2f   %s" \
                    % (ph+1, qh-1,
                       self.smap.glink[ph+1][qh-1].lg[0].motif[0].get_Vij(),
                       self.smap.glink[ph+1][qh-1].lg[0].motif[0].get_branches())
            else:
                print "stored dG(%2d,%2d): %8.2f   %s" \
                    % (ph,   qh,
                       self.smap.glink[ph][qh].lg[0].motif[0].get_Vij(),
                       self.smap.glink[ph][qh].lg[0].motif[0].get_branches())
            #
            print "constructed MBLptr:"
            print new_MBL
            print "Exit get_MBLptr"
        #
        return new_MBL
    #
    
    
    def showBranches(self, mbl, nm = "branch"):
        s = ''
        for i in range(0, len(nm)):
            s += ' '
        #
        s += "index    (  p,  q)       type          FE\n"
        for l in range(0, mbl.n):
            pv = mbl.Q[l].i; qv = mbl.Q[l].j
            ctpv = self.smap.glink[pv][qv].lg[0].motif[0].get_ctp()
            btpv = self.smap.glink[pv][qv].lg[0].motif[0].get_btp()
            Vijv = self.smap.glink[pv][qv].lg[0].motif[0].Vij
            s += "%s(%2d)      (%3d,%3d)     %s[%5s]     %8.2f\n"  \
                % (nm, l + 1, pv, qv,  ctpv, btpv, Vijv)
            #
        #
        return s
    #
    
    
    # (from RNAModules)
    def roughScan_MBL(self, mbl, # (class MBLptr)
                 opt_iMBL = False,
                 display  = False):
        
        """@
        
        This should not be used as the final result. It merely
        calculates the rough FE contribution from each of the
        branches. The final result should be combined with the FE
        obtained from trim_53pStems, which corrects for things like
        neighboring stems that (as a result) have no dangling base to
        contribute to the dangle FE. In short, this just dumbly adds
        the FE from the branches (_NOT_ the base of the stem) without
        any consideration about their arrangement of the branches
        relative to each other.
        
        Additionally, this can be used for displaying information
        about the arrangement of the branches and their type in the
        MBL (iMBL or pMBL).
        
        """
        VpMBL = 0.0
        for kv in range(0, len(mbl.Q)):
            pv = mbl.Q[kv].i; qv = mbl.Q[kv].j
            VpMBL += self.smap.glink[pv][qv].lg[0].motif[0].get_Vij()
        #
        
        if display: # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
            if opt_iMBL:
                print "roughScan_MBL(%2d,%2d) iMBL" % (mbl.i, mbl.j)
            else:
                print "roughScan_MBL(%2d,%2d) pMBL" % (mbl.i, mbl.j)
            #
            print "index      branch     tag       free energy"
            for kv in range(0, len(mbl.Q)):
                pv = mbl.Q[kv].i; qv = mbl.Q[kv].j
                ctpv = self.smap.glink[pv][qv].lg[0].motif[0].get_ctp()
                Vijv = self.smap.glink[pv][qv].lg[0].motif[0].get_Vij()
                print " %3d     (%3d,%3d)    %c       %8.2f" % (kv, pv, qv, ctpv, Vijv)
            #
            pv = mbl.i; qv = mbl.j
            if opt_iMBL:
                pv += 1; qv -= 1
            #
            print "-----      ------     ---       -----------"
            if len(self.smap.glink[pv][qv].lg) > 0:
                ctpv = self.smap.glink[pv][qv].lg[0].motif[0].get_ctp()
                Vijv = self.smap.glink[pv][qv].lg[0].motif[0].get_Vij()
                # may be useful for comparison, though they _should_ be
                # the same value. Maybe if they are not, something is
                # wrong.
                print "root     (%3d,%3d)    %c       %8.2f" % (mbl.i, mbl.j, ctpv, Vijv)
            else:
                print "root*    (%3d,%3d)    X         *empty*" % (mbl.i, mbl.j)
            #
            print "result                        %8.2f" % (VpMBL)
            print "opt_iMBL(%s), pqv(%2d,%2d), ij(%2d,%2d)" \
                % (opt_iMBL, pv, qv, mbl.i, mbl.j)
            print "Exiting roughScan_MBL()"
            #self.stopWhenMatchFound(mbl.i, mbl.j, 0, 6, "roughScan_MBL")
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        return VpMBL
    #
    
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # ######################################################
    # ###############   loop display tools  ################
    # ######################################################
    
    
    # ######################################################
    # ###############   loop building tools  ###############
    # ######################################################
    # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    
    
    
    def stem_dangle(self,
                    i,  j,   # boundaries (often [1 to N]) 
	            p,  q):   # 5'3' tail of the stem
        """@
        
        This is presently a blank, but it is intended to work with the
        concepts of dangling bonds on I-loops and MBLs.
        
        Currently just a blank (for chromatin, this is sufficient).
        Normally, it would look at the local sequence and then decide
        what the residual FE is of the dangling regions. However, with
        chromatin, this is not really something we know presently.
        
        """
        return 0.0
    #
    
    
    def trim_53pStems(self,
                      i, j,       ## int:  5' to 3' of [ip]MBL
                      mbl,        ## class MBPptr:
                      smap,       ## the general mapping 
                      opt_iMBL):  ## True if iMBL, False if pMBL
        """
        This is presently a blank, but it is intended to work with the
        concepts of dangling bonds on I-loops and MBLs.
        """
        # currently just a blank (for chromatin, this is sufficient).
        return mbl.V
    #
    
    
    
    #########################################
    ###  multi-branch loop building tool  ###
    #########################################
    
    
    def cle_MloopEV(self,
		    i,  j,      ## iMBL terminal coordinates 
		    mbl,        ## iMBL branches
                    smap,       ## main map
                    opt_iMBL):  ## True if iMBL, False if pMBL
        """
        This has some working components. It is intended to work with the
        concepts of dangling bonds on I-loops and MBLs found in
        vsfold5.
        """
        # currently, we just make this a blank for chromatin. We
        # presently employ rules based upon the heatmap, so we don't
        # really have to consider whether interactions are allowed or
        # or disallowed or that they require some specific
        # parameterization.
        
        # ptype -> btype[i][j].pair (For chromatin, this would
        # distinguish between CTCFs and singletons
        
        # ptype -> btype[i][j].btp (Perhaps it would consider whether
        # it is a parallel-strand interaction or a antiparallel-strand
        # interaction.)
        
        # get_RNA_seq() -> seq
        
        dG = 0.0
        if opt_iMBL:
            dG += self.calc_dG(i, j, self.hv[i][j], self.T, "cle_MloopEV")
            # dG = TdS(i, j, self.T)
        else:
            dG += 0.0
        #
        
        Vij = 0.0
        for brk in mbl.Q:
            pv = brk.i; qv = brk.j
            # print "pqv: ", pv, qv
            ctp = smap.glink[pv][qv].lg[0].motif[0].ctp
            if ctp == 'B' or ctp == 'S' or ctp == 'R' or ctp == 'K' or ctp == 'W':
                Vij += smap.glink[pv][qv].lg[0].motif[0].Vij
                if self.debug_BranchEntropy:
                    ctpv = smap.glink[pv][qv].lg[0].motif[0].ctp
                    btpv = smap.glink[pv][qv].lg[0].motif[0].btp
                    Vijv = smap.glink[pv][qv].lg[0].motif[0].Vij
                    print "ij{(%2d,%2d), pqv(%2d,%2d)}[%s][%5s][%8.2f]  Vtotal = %8.2f" \
                        % (i, j, pv, qv, ctpv, btpv, Vijv, Vij)
                    #
                #endif
            #
        #
        # print Vij
        dG += Vij
        dG  = self.trim_53pStems(i, j, mbl, smap, opt_iMBL)
        dG += self.cle_loopInt(i, j, j-i+1)
        if self.debug_BranchEntropy:
            print "cle_MloopEV(result): dG = %8.2f" % dG
        #
        return dG
    #
    
    
    def cle_loopInt(self, i, j, ilen):
        """
        This is presently a blank, but it is intended to work with the
        concepts of dangling bonds on I-loops and MBLs.
        """
        # currently just a blank (for chromatin, this is sufficient).
        return 0.0
    #
    
    
    
    #####################################
    ###  internal loop building tool  ###
    #####################################
    
    
    def cle_IloopE(self,
                   i,  j,      ## tail of I-loop 
		   p,  q,      ## head of I-loop 
		   mbl,        ## MBLptr 
                   smap,       ## smap (motif)
                   opt_iMBL):  ## iMBL=True, pMBL = False
        """
        I think this may be redundant. It is working. It is intended to
        work with the concepts of dangling bonds on I-loops found in
        vsfold5.
        """
        # currently, we just make this a blank for chromatin. We
        # presently employ rules based upon the heatmap, so we don't
        # really have to consider whether interactions are allowed or
        # or disallowed or that they require some specific
        # parameterization.
        
        # ptype -> btype[i][j].pair (For chromatin, this would
        # distinguish between CTCFs and singletons
        
        # ptype -> btype[i][j].btp (Perhaps it would consider whether
        # it is a parallel-strand interaction or a antiparallel-strand
        # interaction.)
        
        # get_RNA_seq() -> seq
        
        dG = 0.0
        if opt_iMBL:
            dG += self.calc_dG(i, j, self.hv[i][j], self.T, "cle_IloopE")
            # dG = TdS(i, j, self.T)
        else:
            dG += 0.0
        #
        
        Vij = 0.0
        ctp = smap.glink[p][q].lg[0].motif[0].ctp
        if ctp == 'B' or ctp == 'S' or ctp == 'R' or ctp == 'K' or ctp == 'W':
            Vij += smap.glink[p][q].lg[0].motif[0].Vij
            if self.debug_BranchEntropy:
                ctpv = smap.glink[p][q].lg[0].motif[0].ctp
                btpv = smap.glink[p][q].lg[0].motif[0].btp
                Vijv = smap.glink[p][q].lg[0].motif[0].Vij
                print "ij{(%2d,%2d), pqv(%2d,%2d)}[%s][%5s][%8.2f]  Vtotal = %8.2f" \
                    % (i, j, pv, qv, ctpv, btpv, Vijv, Vij)
                #
            #endif
        else:
            print "ERROR: I-loop should link to a object of type B,S,R,K or W"
            print "       position(i(%d) <= p(%d) < q(%d) <= j(%d)" % (i, p, q, j)
            sys.exit(1)
        #
        dG += Vij
        if self.debug_BranchEntropy:
            print "cle_IloopE(result): dG = %8.2f" % dG
        #
        return dG
    #
    
    
    ####################################
    ###  hairpin loop building tool  ###
    ####################################
    
    def cle_HloopE(self,
                   i,  j,  ## tail of I-loop 
		   mbl):   ## MBLptr 
        """@
        
        This is a basic function for H loops. It is working. It is
        intended to work with the concepts of dangling bonds on
        H-loops found in vsfold5.
        
        currently, we just make this a blank for chromatin. We
        presently employ rules based upon the heatmap, so we don't
        really have to consider whether interactions are allowed or or
        disallowed or that they require some specific
        parameterization.
        
        ptype -> btype[i][j].pair (For chromatin, this would
        distinguish between CTCFs and singletons
        
        ptype -> btype[i][j].btp (Perhaps it would consider whether it
        is a parallel-strand interaction or a antiparallel-strand
        interaction.)

        """
        
        # get_RNA_seq() -> seq
        dG = self.calc_dG(i, j, self.hv[i][j], self.T, "cle_HloopE")
        # dG = self.TdS(i, j, self.T)
        return dG
    #
    
    
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # ######################################################
    # ###############   loop building tools  ###############
    # ######################################################
    
    # #######################################################
    # #############  stem & PK filtering tools  #############
    # #######################################################
    # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    
    
    def filter_BIKMSW_from_glink(self, k, l):
        """@ 
        
        The top of the list need not be an 'M', 'K', 'I' or 'B'!!!  So
        this makes sure that the structure that is selected is one of
        these. It is still used find_ifStem(), though it is for the
        parallel stem search where I have shut that part off
        presently.  The purpose is to find a structure with a single
        closing point at (k,l), not a pMBL or a J-loop.
        
        I start to realize that the purpose of this thing was not
        exactly what I thought. Originally, it was used with I-loop
        searches, but recently, I got rid of those in favor of
        searchForMBL (which does a far more complete job of
        searching). I'm not so sure this should be used where it is
        with the parallel stem. Mainly there, I want something
        attached on top, whatever it is.
        
        """
        DEBUG_filter_glink = False # True # 
        #
        V   = INFINITY
        cp  = []
        ctp = ''
        pks = []
        wps = []
        if DEBUG_filter_glink:
            print "filter_BIKMSW_from_glink -------------"
        for m in range(0, len(self.smap.glink[k][l].lg)):
            ctpx = self.smap.glink[k][l].lg[m].motif[0].get_ctp()
            if DEBUG_filter_glink:
                print "(%2d,%2d)[%d]: %8.3f [%s]" \
                    % (k,l, m, self.smap.glink[k][l].lg[m].Vij, ctpx)
            if ctpx == 'B' or ctpx == 'I' or ctpx == 'M':
                # 160601wkd: M and K are also closed at the ends!!!!
                if self.smap.glink[k][l].lg[m].Vij < V:
                    for vv in self.smap.glink[k][l].lg[m].motif:
                        cp +=  vv.get_base() # original: cp = (k,l)
                    V  = self.smap.glink[k][l].lg[m].Vij
                    ctp = ctpx
                #
            elif ctpx == 'S':
                # this could be mixed together with the first group
                if self.smap.glink[k][l].lg[m].Vij < V:
                    for vv in self.smap.glink[k][l].lg[m].motif:
                        cp +=  vv.get_base() # original: cp = (k,l)
                    V  = self.smap.glink[k][l].lg[m].Vij
                    ctp = ctpx
                    
            elif  ctpx == 'K':
                if self.smap.glink[k][l].lg[m].Vij < V:
                    for vv in self.smap.glink[k][l].lg[m].motif:
                        cp +=  vv.get_base() # original: cp = (k,l)
                    V  = self.smap.glink[k][l].lg[m].Vij
                    ctp = ctpx
                    pks = self.smap.glink[k][l].lg[m].motif[0].get_pks()
                #
            elif ctpx == 'W':
                if self.smap.glink[k][l].lg[m].Vij < V:
                    for vv in self.smap.glink[k][l].lg[m].motif:
                        cp +=  vv.get_base() # original: cp = (k,l)
                    V  = self.smap.glink[k][l].lg[m].Vij
                    ctp = ctpx
                    wps = self.smap.glink[k][l].lg[m].motif[0].get_wyspa()
                    if DEBUG_filter_glink:
                        self.smap.glink[k][l].lg[m].motif[0].show_Motif()
                        #sys.exit(0)
                        # xxx!!!
                    #
                #
                #
        if DEBUG_filter_glink:
            print ".", cp
        #
        return V, cp, ctp, pks, wps 
    #

    # (from RNAModules)
    def filter_KSW_from_glink(self, k, l):
        """@ 
        
        The top of the list need not be an 'S', 'K' or 'W'!!!  So this
        makes sure that the structure that is selected is the best one
        of the list or None. When this procedure is called, a K, S, or
        W should be the top item. The main reason for scanning through
        the list is just in case it hasn't already be sorted. The
        purpose is to find a structure with a single closing point at
        (k,l) that is not a pMBL or a J-loop or a single base pair.
        
        This tool is used by find_best_PK().

        """
        DEBUG_filter_glink = False # True # 
        #
        
        V   = INFINITY
        Ob  = None
        if DEBUG_filter_glink:
            print "filter_KSW_from_glink -------------"
        #endif
        for m in range(0, len(self.smap.glink[k][l].lg)):
            ctpx = self.smap.glink[k][l].lg[m].motif[0].get_ctp()
            if DEBUG_filter_glink:
                print "(%2d,%2d)[%d]: %8.3f [%s]" \
                    % (k,l, m, self.smap.glink[k][l].lg[m].Vij, ctpx)
            elif ctpx == 'S':
                # this could be mixed together with the first group
                if self.smap.glink[k][l].lg[m].Vij < V:
                    Ob = self.smap.glink[k][l].lg[m].motif[0].Ob
                    V  = self.smap.glink[k][l].lg[m].Vij
            elif  ctpx == 'K':
                if self.smap.glink[k][l].lg[m].Vij < V:
                    Ob = self.smap.glink[k][l].lg[m].motif[0].Ob
                    V  = self.smap.glink[k][l].lg[m].Vij
            elif ctpx == 'W':
                if self.smap.glink[k][l].lg[m].Vij < V:
                    Ob = self.smap.glink[k][l].lg[m].motif[0].Ob
                    V  = self.smap.glink[k][l].lg[m].Vij
                    if DEBUG_filter_glink:
                        self.smap.glink[k][l].lg[m].motif[0].show_Motif()
                        #sys.exit(0)
                        # xxx!!!
                    #endif
                #
            #
        #
        if DEBUG_filter_glink:
            print ".", Ob
        #
        return Ob
    #
    
    
    
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # #######################################################
    # #############  stem & PK filtering tools  #############
    # #######################################################
    
    



    # ######################################################
    # ###############   stem building tools  ###############
    # ######################################################
    # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    def calc_ddGlocal_Stem(self,
                           i_t, j_t,   # stem tail (needed for checking)
                           i_h, j_h,   # stem head (needed for checking)
                           num_bp,     # the number of bp in the stem
                           flag_debug = False):  # option to debug
        """@
        
        Calculate the difference in the effective stem FE vs the free
        strand contribution.
        
        For chromatin, the Kuhn length is not changed and the free
        strand free energy is zero. This is here to maintain
        compatability with other modules.
        
        """
        L_eff = float((i_h - i_t + 1) + (j_t - j_h +1))/2.0
        xi_stm = self.xi
        dGloc_stm = 0.0
        dGloc_fs  = 0.0
        
        return xi_stm, L_eff, dGloc_stm, dGloc_fs
    #
    
    
    
    # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
    def find_branch_weight(self, ih, jh, branches, opt_iMBL = True):
        """@
        
        Used in RTrace to calculate the connecting junction FE between
        the stem and its neighboring object
        
        """
        dGh = 0.0
        if len(branches) > 1:
            new_MBL = MBLptr(ih, jh)
            new_MBL.nm = ctp
            new_MBL.dr = btp
            new_MBL.V  = 0.0
            new_MBL.pushBranchlist(branches, True)
            # dangling bond corrections
            dGpqh_d53p = self.trim_53pStems(ih, jh,   # 5' to 3' iMBL closing pt
                                            new_MBL,  # mbl branch info
                                            opt_iMBL) # iMBL or pMBL? B.Cs
            # loop closing FE
            MLloopEV   = self.cle_MloopEV(ih,  jh,    # iMBL terminal coordinates 
                                          new_MBL,  # iMBL branches
                                          opt_iMBL) # True if iMBL, False if pMBL
            dGh = dGpqh_d53p + MLloopEV
        elif len(branches) == 1:
            p = branches[0][0]; q = branches[0][1]
            xi_stm = self.xi
            dGh    = self.cle_IloopE(ih, jh,  # tail of I-loop 
       	                                     p, q,    # head of I-loop 
                                             xi_stm)  # [nt] Kuhn length
        else:
            dGh = self.cle_HloopE(ih, jh)  ## tail of H-loop
            
        #
        
        return dGh
    #
    # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    
    
    
    def new_make_aaStemMotif(self, stem, fr_zones, flag_debug = False):
        """copied from RNAModules for update in ChromatinModules"""
        flag_monitor = False # True # 
        """@
        
        In building stems with a variable Kuhn length, we have to
        basically compute them individually. This method accepts a
        list of base pairs in a stem (they must have some consecuitive
        pattern -- e.g., (0,12),(1,11),(3,10), (4,9) -- but otherwise,
        it is treated as a list).
        
        stem -> list [(0,20), (1,19), (2,18) ... ]
        
        """
        
        # ###########################################################
        # Now we have all the information about what stem we need to
        # build. Now we have to build the stem in part.
        # ###########################################################
        
        # build the current antiparallel stem at (i_t,j_t) = (i,j) -> (i_h,j_h)
        stmlen_full = len(stem)  # integer
        if stmlen_full == 0:
            print "ERROR: stems absolutely must contain at least one base pair"
            print "       this stem is empty!"
            sys.exit(1)
        #
        
        n = stmlen_full - 1
        i  = stem[0][0]; j  = stem[0][1]
        ph = stem[n][0]; qh = stem[n][1]
        
        xi_stm, L_eff, dGloc_stm, dGloc_fs \
            = self.calc_ddGlocal_Stem(i,  j,       # stem tail (needed for checking)
                                      ph, qh,      # stem head (needed for checking)
                                      stmlen_full) # the number of bp in the stem
        
        ddG_l = dGloc_stm - dGloc_fs
        
        # build the free energy for the stem
        pairs_aa = []
        dGaa_stem = 0.0
        dGijpmk_aa = 0.0
        
        for k in range(0, n):
            dGijpmk_aa = 0.0
            iaa    = stem[k  ][0]; jaa    = stem[k  ][1]
            iaa_p1 = stem[k+1][0]; jaa_m1 = stem[k+1][1]
            if not (iaa_p1 - iaa == 1 and jaa - jaa_m1 == 1):
                # add connected I-loop between the effective stems
                
                dGijpmk_aa = self.cle_IloopE(iaa,    jaa,    # tail of I-loop 
       	                                     iaa_p1, jaa_m1, # head of I-loop 
                                             xi_stm)         # [nt] Kuhn length
                
                if flag_debug:
                    print "iloop(k=%2d): ijaa(%2d,%2d)::ijaa_pm1(%2d,%2d)[%s]" \
                        % (k, iaa, jaa, iaa_p1, jaa_m1, self.btype[iaa][jaa].f2nt)
                    print "dGijpmk_aa: %8.2f" % (dGijpmk_aa)
                    #sys.exit(0)
                #
                
            else:
                dGijpmk_aa  = self.calc_pair_dG(iaa, jaa, xi_stm)
                if flag_debug:
                    print "stem (k=%2d): ijaa(%2d,%2d)::ijaa_pm1(%2d,%2d)[%s]" \
                        % (k, iaa, jaa, iaa_p1, jaa_m1, self.btype[iaa][jaa].f2nt)
                    print "dGijpmk_aa: %8.2f" % (dGijpmk_aa)
                    #sys.exit(0)
                #
                
            #
            
            dGaa_stem += dGijpmk_aa
            
            
            if flag_monitor and dGijpmk_aa >= 0.0:
                print "Note!!! (%d,%d): fe.dGp = %8.2f, ==>> dGijpmk_aa = %8.2f <<==" \
                    % (iaa, jaa, self.btype[iaa][jaa].dGp, dGijpmk_aa)
            #
            
            pairs_aa += [((iaa, jaa), dGijpmk_aa, "sa")]
            # there is no other option than anti-parallel than "sa"
            
            # pairs_aa: in effect, stem is expanded to contain more
            # information in this step
            
        #
        
        # testRNAseq3
        #self.stopWhenMatchFound(i, j, 0,  5, "new_make_aaStemMotif 1")
        #self.stopWhenMatchFound(i, j,  2, 45, "new_make_aaStemMotif 1")
        
        # testRNAseq3a
        #self.stopWhenMatchFound(i, j, 0, 13, "new_make_aaStemMotif 1")
        
        # testRNAseq3b
        #self.stopWhenMatchFound(i, j, 0, 19, "new_make_aaStemMotif 1")
        
        # testRNAseq4
        #self.stopWhenMatchFound(i, j, 5, 18, "new_make_aaStemMotif 1")
        #self.stopWhenMatchFound(i, j, 5, 48, "new_make_aaStemMotif 1")
        
        
        MLclose   = 0.0
        dGpqh_V   = 0.0
        dGpqh_cap = 0.0
        ctp       = ''
        btp       = ''
        branchlist = []
        
        tMLclose   = 0.0
        tdGpqh_V   = 0.0
        tdGpqh_cap = 0.0
        tctp       = ''
        tbtp       = ''
        tpairs_aa = []
        for v in pairs_aa:
            tpairs_aa += [v]
        #
        
        
            
        cap_aa = self.find_branching(ph, qh, [], 'sa', flag_debug)
        MLclose    = cap_aa[0]
        dGpqh_V    = cap_aa[1]
        dGpqh_cap  = cap_aa[2]
        ctp        = cap_aa[3]
        btp        = cap_aa[4]
        branchlist = cap_aa[5]
        pairs_aa += [((ph, qh), tMLclose, "sa")]
            
        """
         cap_aa = [MLclose,         #  0 closing entropy correction at (ph,qh)
                   dGpqh_V,         #  1 FE from adjoining connections to (ph, qh)
                   dGpqh_cap,       #  2 total FE of the cap region
                   ctp,             #  3 connection type
                   btp,             #  4 bond type (thinking ahead for parallel stems)
                   branchlist]      #  5 branches at (ph,qh)
        
        """
        
        
        if flag_debug:
            print "pairs_aa: ", tpairs_aa
            print self.rnaseq
        #
        
        dGaa_full = dGaa_stem + dGpqh_cap + ddG_l
        
        stem_aa = [(i, j),          #  0
                   dGaa_full,       #  1 full FE from (i,j)
                   pairs_aa,        #  2 individual bp pair weights
                   xi,              #  3 stem Kuhn length
                   dGaa_stem,       #  4 pair contribution to Stem 
                   dGloc_stm,       #  5 local stem entropy corrections
                   dGloc_fs, #  6 local free strand entropy corrections
                   MLclose,         #  7 closing entropy correction at (ph,qh)
                   dGpqh_V,         #  8 FE from adjoining connection to (ph, qh)
                   branchlist]      #  9 branches at (ph,qh)
        
        if flag_debug:
            
            V = 0.0
            for k in range(0, len(pairs_aa)):
                V += pairs_aa[k][1]
                ik = pairs_aa[k][0][0]; jk = pairs_aa[k][0][1]
                print "(%2d,%2d)[%8.2f]" % (ik, jk, pairs_aa[k][1]) 
            #
            
            i_t = i;   j_t = j
            i_h = ph;  j_h = qh
            
            print "dGstem = %8.2f" % V
            dGpqh_TdS = self.TdS(ph, qh, xi_stm) # used for reference
            print "aa pq_h(%2d,%2d)[dGpqh_V]:         %8.2f" % (i_h, j_h, dGpqh_V)
            print "aa pq_h(%2d,%2d)[MLclose]:         %8.2f" % (i_h, j_h, MLclose)
            print "--------------------------------    -------"
            print "aa pq_h(%2d,%2d)[dGpqh_cap]:       %8.2f" % (i_h, j_h, dGpqh_cap)
            print "********************************    *******"
            print "aa ij_t(%2d,%2d)[dGaa_stem]:       %8.2f" % (i_t, j_t, dGaa_stem)
            print "aa ij_t(%2d,%2d)[dGloc_stm]:       %8.2f" % (i_t, j_t, dGloc_stm)
            print "aa pq_h(%2d,%2d)[dGloc_fs]:        %8.2f" % (i_t, j_t, -dGloc_fs)
            print "aa pq_h(%2d,%2d)[dGpqh_cap]:       %8.2f" % (i_h, j_h, dGpqh_cap)
            
            print "--------------------------------    -------"
            print "aa ij_t(%2d,%2d)[dGaa_full]:       %8.2f" % (i_t, j_t, dGaa_full)
            
            print "reference (global entropy at pqh):"
            print "aa pq_h(%2d,%2d)[dGpqh_TdS]:       %8.2f" % (i_h, j_h, dGpqh_TdS)
            print stem_aa
            print "Exiting new_make_aaStemMotif()"
            
            """
            if ctpaa == 'I' or ctpaa == 'M':
                print "planned exit"
                sys.exit(0)
            #
            sys.exit(0)
            """
            
        #
        
        dGij = dGaa_full
        
        stm = Stem(stem)
        
        stm.Vij = dGij
        for ss in pairs_aa:
            stm.dGp += [ss[1]]
            # each stem should be of one singlar type, so there can be
            # no mixing of parallel stems here. Hence, btp is sufficient.
        #
        
        stm.xi     = xi_stm
        stm.dG_l   = dGloc_stm
        stm.dG_fs  = dGloc_fs
        stm.dGloop = MLclose
        stm.Vpqh   = dGpqh_V
        stm.branching = []
        for k in range(0, len(branchlist)):
            stm.branching += [branchlist[k]]
        #
        
        
        link_S = Link(stm)
        if flag_debug:
            print link_S.motif[0].show_Motif()
        #
        
        self.smap.glink[i][j].add_link(link_S)
        
        if flag_debug:
            print "XXXXXXXXXX"
            self.show_StemMotif(i, j, stm)
            print "Exit new_make_aaStemMotif{(%d,%d)}" % (i,j)
            print self.smap.glink[i][j].lg[0].motif[0].show_Motif()
            #if i == 1 and 30 == j:  sys.exit(0)
            #if i == 0 and 11 == j:  sys.exit(0)
            #sys.exit(0)
        #
        
        # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        
        # Assuming an antiparallel connection, we are done
        # here. However, parallel stems require a bit more attention.
        
        # testRNAseq2
        #self.stopWhenMatchFound(i, j, 0, 30, "new_make_aaStemMotif, end of")
        
        # testRNAseq2b
        #self.stopWhenMatchFound(i, j,15, 20, "new_make_aaStemMotif, end of")
        #self.stopWhenMatchFound(i, j, 0, 36, "new_make_aaStemMotif, end of")
        
        # testRNAseq3
        #self.stopWhenMatchFound(i, j, 2, 45, "new_make_aaStemMotif, end of")
        
        # testRNAseq4
        #self.stopWhenMatchFound(i, j, 5, 18, "new_make_aaStemMotif, end of")
        #self.stopWhenMatchFound(i, j,17, 36, "new_make_aaStemMotif, end of")
        #self.stopWhenMatchFound(i, j,14, 39, "new_make_aaStemMotif, end of")
        #self.stopWhenMatchFound(i, j,11, 42, "new_make_aaStemMotif, end of")
        #self.stopWhenMatchFound(i, j, 8, 45, "new_make_aaStemMotif, end of")
        #self.stopWhenMatchFound(i, j, 5, 48, "new_make_aaStemMotif, end of")
        
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        return stem_aa
    
    #
    
    
    def show_StemMotif(self, i, j, stm):
        print "Stem:"
        print stm.disp_Stem()
        
        V = 0.0
        for k in range(0, len(stm.stem)-1):
            V += stm.dGp[k]
            print " (%2d,%2d)   %8.2f" % (stm.stem[k].i, stm.stem[k].j, stm.dGp[k])
        #
        
        dGij = stm.get_Vij()
        V += stm.dG_l - stm.dG_fs
        V += stm.Vpqh + stm.dGloop
        print "V_tot = %8.2f, dGloop = %8.2f" % (V, stm.dGloop)
        
        btp = "sa"
        print "general: ", (i,j), dGij, stm.stem, btp, xi
                
        stmlen = len(stm.stem)
        i_t = i;        j_t = j
        i_h = stm.ih;   j_h = stm.jh
        ihh = i_h + 1;  jhh = j_h - 1
        
        
        """@
        
        The computation is (i_t,j_t) -- the tail -- and (i_h,j_h) --
        the head. From the diagram below, perhaps one can understood
        that the entry points and exit points are in the same position
        on the chain.
        
        
           antiparallel        
             _____                
            /     \            
            |     |            
            |     |            
            |     |            
             \   /             
         i_h  \_/ j_h          
              |_|              
         i_t  |_| j_t
        
        """
        
        print "new_make_aaStemMotif(): anti-parallel stem"
        print "aa ij_t (tail): ", i_t, j_t
        print "aa ij_h (head): ", i_h, j_h
        print "aa ijhh (jnct): ", ihh, jhh, ", dijhh: ", (jhh - ihh)
        print "region does not have a defined terminus"
        # !!!xxx
        
        #if i_t == 5 and j_t == 48:
        #if i_t == 5 and j_t == 13:            
        if i_t == 2 and j_t == 45:
            print "btp(%d,%d) = %s" % (i, j, btp)
            for m in range(0, len(self.smap.glink[i_t][j_t].lg)):
                print self.smap.glink[i_t][j_t].lg[m].motif[0].show_Motif()
            #
            
            print "get_bondtype: ", self.btype[i][j].btp
            #sys.exit(0)
        # 
        #sys.exit(0)
    #
    
    def find_branching(self, ph, qh, fr_zones, btp, flag_debug = False):
        """@
        
        (ph, qh): the closing point at the "head" of the stem
        
        fr_zones: a list of free "zones" where we are to search for
                  branches. This is needed particularly when
                  evaluating pseudoknots where the linkage and root
                  stems do not have their branching starting from (ph,
                  qh) but are divided into zones between different
                  parts of the loop. If fr_zones = [], then the free
                  zone isassumed to be the same as (ph,qh).
        
        """
        branchlist = []
        if len(fr_zones) == 0:
            """@ 

            If len(fr_zones) == 0, then the requested operation can
            only be on a standard stem. In that case, we want to be
            sure that we don't suck up some nearest neighboring bp at
            (ph+1, qh-1) because that would presume an I-loop with
            
            (ph,qh) -> (ph+1, qh-1), 
            
            which is just another bp in a _contiguous_ stem. That maes
            no sense to evaluate it as a different stem. Indeed, the
            results would be erroneous. Therefore, in this step, we
            make sure that the stem connects to something other than
            part of the same contiguous stem.
            
            """
            br = self.filter_XJP_from_glink(ph+1, qh-1)
            if len(br) > 0:
                branchlist = br
            #
            
            
        else: #  len(fr_zones) > 0:
            for zn in fr_zones:
                pz = zn[0]; qz = zn[1]
                br = self.smap.glink[pz+1][qz-1].lg[0].motif[0].get_branches()
                if len(br) > 0:
                    branchlist += br
                #
            #
        #
        
        ctp = 'B'
        if len(branchlist) == 1:
            ctp = 'I'
        elif len(branchlist) > 1:
            ctp = 'M'
        #
        
        if flag_debug:
            print "stem terminates at pqh(%2d,%2d)" % (ph, qh)
            print "scan zone:   %s" % (fr_zones)
            print "branchlist:  %s" % (branchlist)
            #print "find_branching 1: exit"; sys.exit(0)
        #endif
        
        """@
        
        Compute the free energy of the terminal motif.  In the case of
        antiparallel stems, they are joined by an antiparallel
        connector like 'B', 'I', 'K', or 'M'. Therefore, the _last
        item_ on the list should NOT be computed from the free energy
        of a closing loop, as this will result in double counting the
        closing region.
        
        """
        
        dGpqh_cap = 0.0
        MLclose   = 0.0
        dGpqh_V   = 0.0
        
        if len(branchlist) == 0:
            # terminate at (ph,qh) with B loop
            ctp = 'B'
            MLclose   = self.cle_HloopE(ph, qh)  ## tail of H-loop
            dGpqh_cap = MLclose
            
            if flag_debug:
                print "stem head forms an H-loop:"
                print "dGpqh_V    = %8.2f" % (dGpqh_V)
                print "MLclose    = %8.2f" % (MLclose)
                print "--------------    -------"
                print "dGpqh_cap  = %8.2f" % (dGpqh_cap)
            #
        elif len(branchlist) == 1:
            # terminate at (ph,qh) with I loop
            ctp = 'I'
            ptk = branchlist[0][0]; qtk = branchlist[0][1]
            MLclose  = self.cle_IloopE(ph, qh,     # tail of I-loop 
       	                               ptk, qtk,   # head of I-loop 
                                       self.xi_fs) # [nt] Kuhn length
            dGpqh_V   = self.smap.glink[ptk][qtk].lg[0].motif[0].get_Vij()
            dGpqh_cap = dGpqh_V + MLclose
            branchlist = [(ptk, qtk)]
            if flag_debug:
                print "stem head forms an I-loop:"
                print "dGpqh_V    = %8.2f" % (dGpqh_V)
                print "MLclose    = %8.2f" % (MLclose)
                print "--------------    -------"
                print "dGpqh_cap  = %8.2f" % (dGpqh_cap)
            #
            
        else: # len(branchlist) > 1
            # terminate at (ph,qh) with M loop
            opt_iMBL = True
            
            ctp = 'M'
            # have to construct the unit to pass
            new_MBL = self.build_MBLptr(ph, qh, branchlist, flag_debug)
            
            # branching FE
            dGpqh_V    = new_MBL.V
            dGpqh_d53p = self.trim_53pStems(ph, qh,   # 5' to 3' iMBL closing pt
                                            new_MBL,  # mbl branch info
                                            opt_iMBL) # iMBL or pMBL? B.Cs
            # loop closing FE
            MLloopEV   = self.cle_MloopEV(ph, qh,    # iMBL terminal coordinates 
	                                  new_MBL,   # iMBL branches
                                          opt_iMBL)  # True if iMBL, False if pMBL
            MLclose   = dGpqh_d53p + MLloopEV
            dGpqh_cap = dGpqh_V + MLclose 
            if flag_debug:
                print "stem head forms an MBL:"
                print "dGpqh_d53p = %8.2f" % (dGpqh_d53p)
                print "MLloopEV   = %8.2f" % (MLloopEV)
                print "--------------    -------"
                print "MLclose    = %8.2f" % (MLclose)
                print "dGpqh_V    = %8.2f" % (dGpqh_V)
                print "--------------    -------"
                print "dGpqh_cap  = %8.2f" % (dGpqh_cap)
            #
            #sys.exit(0)
            
        #
        
        cap_aa = [MLclose,         #  0 closing entropy correction at (ph,qh)
                  dGpqh_V,         #  1 FE from adjoining connections to (ph, qh)
                  dGpqh_cap,       #  2 total FE of the cap region
                  ctp,             #  3 connection type
                  btp,             #  4 bond type (thinking ahead for parallel stems)
                  branchlist]      #  5 branches at (ph,qh)
        
        return cap_aa
    #
    
    
    
    
    def make_StemMotif(self, stem, debug = False):
        flag_debug = debug
        i = stem[0][0]; j = stem[0][1]; dGij = stem[1]
        if flag_debug:
            print "Enter make_StemMotif{(%d,%d), dGij=%8.2f}" % (i,j, dGij)
            print "stem: ", stem
        #
        
        join = []
        btp_special = ''
        btp_general = ''
        for ss in stem[2]:
            join += [ss[0]]
            if ss[2] == 'c' or ss[2] == 'l' or ss[2] == 'r' or ss[2] == 't':
                btp_special = ss[2] # btype
            else:
                btp_general = ss[2] # btype
        #
        btp = btp_general
        if len(btp_special) > 0:
            if flag_debug:
                print "make_StemMotif()"
                print "special: ", (i,j), dGij, join, btp_special
            #
            btp = btp_special
        else:
            if flag_debug:
                print "general: ", (i,j), dGij, join, btp_special
                print "btp_general: ", btp_general
                print "btp_special: ", btp_special
            #
        #
        i_t = i; j_t = j
        stmlen = len(join)
        i_h = i + stmlen - 1; j_h = j - stmlen + 1
        ihh = i_h + 1;        jhh = j_h - 1
        if btp_general == "sp":
            i_h = i + stmlen; j_h = j + stmlen
            ihh = i_h - 1;    jhh = j_h - 1
        #

        """@
        
        I know, it is kind of strange, but whether the stem is
        parallel or antiparallel, the computation of (i_t,j_t) -- the
        tail -- and (i_h,j_h) -- the head -- is the same.
        
        From the diagram below, perhaps one can understood that the
        entry points and exit points are in the same position on the
        chain, it is just that the order of the pairing arrangement
        is changed.
        
        
           antiparallel                parallel
             _____                
            /     \                            ______
            |     |                    j_h    /
            |     |                   __._ _./ j_t
            |     |                  / _|_|_|___
             \   /               ____|/i_t i_h  \
         i_h  \_/ j_h                |           |
              |_|                    |___________|
         i_t  |_| j_t   
                        
        
        the main difference is that the head of the antiparallel
        structure can encorporate a motif of type 'B', 'I' or 'M', but
        the parallel stem motif cannot. This is because the closing
        points in the parallel stem are at disparate positions
        (separated by a distance (i_h - i_t) or (j_t - j_h). On the
        other hand, the closing point (i_h, j_h) on the antiparallel
        stem is exactly a point where an 'M' motif can attach.
        
        """
        
        if flag_debug:
            print "make_StemMotif(): parallel stem"
            print "pp ij_t (tail): ", i_t, j_t
            print "pp ij_h (head): ", i_h, j_h
            print "pp ijhh (jnct): ", ihh, jhh, ", dijhh: ", (jhh - ihh)
            print "region does not have a defined terminus"
            # !!!xxx
            if i_t == 1 and j_t == 30:
            #if i_t == 5 and j_t == 13:            
                print "btp(%d,%d) = %s" % (i, j, btp)
                for m in range(0, len(self.smap.glink[i_t][j_t].lg)):
                    print self.smap.glink[i_t][j_t].lg[m].motif[0].show_Motif()
                #
                print "get_bondtype: ", self.btype[i][j].btp
                # was: self.get_bondtype(self.btype[i][j].wt)
                # self.hv[i][j] -> self.btype[i][j].wt
                print "all_ctcf:     ", self.all_ctcf
                #sys.exit(0)
                #
            # sys.exit(0)
        #
        
        link_S = Link(i_t, j_t, dGij, 'S', btp, join)
        if flag_debug:
            print link_S.motif[0].show_Motif()
        #
        self.smap.glink[i_t][j_t].add_link(link_S)
        
        if flag_debug:
            print "Exit make_StemMotif{(%d,%d)}" % (i,j)
            print self.smap.glink[i_t][j_t].lg[0].motif[0].show_Motif()
            #if i_t == 1 and 30 == j_t:  sys.exit(0)
            #if i_t == 5 and 13 == j_t:  sys.exit(0)
            
        #
        # Assuming an antiparallel connection, we are done
        # here. However, parallel stems require a bit more attention.
        
    #
    
    
    # sort the branching results from find_ifStem() and or the pks
    # results from find_best_PK()
    def ins_sort_stempp(self, pp):
        for i in range(1,len(pp)):    
            j = i                    
            while j > 0 and pp[j][0] < pp[j-1][0]: 
                pp[j], pp[j-1] = pp[j-1], pp[j] # syntactic sugar: swap the items
                j=j-1 
            #
        #
        return pp
    #
    
    def lookupCap(self, i_h, j_h, dG_h, flag_debug = False):

        flag_debug = False # True #
        has_cap = False
        ctp_h = 'X'
        btp_h = '-'
        jn_h = []
        i_x = i_h + 1; j_x = j_h - 1
        
        #      *  *
        #    *      *
        # i_x *    * j_x
        #  i_h |--| j_h
        #      |--|
        #      |--|
        
        n_x = len(self.smap.glink[i_x][j_x].lg)
        n_h = len(self.smap.glink[i_h][j_h].lg)
        if flag_debug:
            print "Enter lookupCap{ij(%d,%d), dG_h(%8.2f)}" % (i_h, j_h, dG_h)
            print "lg x: ", n_x
            print "lg h: ", n_h
            print "btype(%d,%d).pair = %d" % (i_h, j_h, self.btype[i_h][j_h].pair)
        #
        
        if n_x > 0 and n_h > 0:
            # 190129: In part 1, we used searchForMBL to determine the
            # best I- and M-loop. However, I am concerned that we can
            # end up doing adding on adding if we don't redo this
            # search from under the requirement that the structure
            # have boundaries p and q that satisfy the requirement
            
            # ({pq}|i_h < p < q < j_h-1) && !(p==i_h+1 && q==j_h-1).
            
            # Additionally, the binding point (i_h+1,j_h-1) should be
            # such that V(i_h+1,j_h-1) < 0.0. All this information is
            # obtained in searchForMBL.
            
            ctp_h = self.smap.glink[i_h][j_h].lg[0].motif[0].get_ctp()
            ctp_x = self.smap.glink[i_x][j_x].lg[0].motif[0].get_ctp()
            if flag_debug:
                print "ctp_h: ", ctp_h, ", ctp_x: ", ctp_x 
            #
            
            if ctp_x == 'J' or ctp_x == 'I':
                # unfortunately, we have to extract from ij_x because
                # otherwise we will keep adding more and more I to the
                # outcome
                Vij_h = self.smap.glink[i_x][j_x].lg[0].Vij
                jn_h  = self.smap.glink[i_x][j_x].lg[0].motif[0].get_branches()
                if not (jn_h[0][0] == i_x and jn_h[0][1] == j_x):
                    if Vij_h < 0.0:
                        dG_h += Vij_h
                        ctp_h = 'I'
                        btp_h = 's'
                        has_cap = True
                        if flag_debug:
                            print "lookupCap, J: ", ctp_x
                            print "Vij_h: %8.2f, (%s)" % (Vij_h, ctp_h)
                            print "len lg: ", len(self.smap.glink[i_x][j_x].lg)
                            print self.smap.glink[i_h][j_h].lg[0].motif[0].show_Motif()
                            print self.smap.glink[i_x][j_x].lg[0].motif[0].show_Motif()
                            print "dG_h: ", dG_h, ", jn_h: ", jn_h
                            # sys.exit(0)
                        #
                        
                    #
                    
                #
                
            elif ctp_x == 'P' or ctp_x == 'M':
                # unfortunately, we have to extract from ij_x because
                # otherwise we will keep adding more and more M to the
                # outcome
                Vij_h = self.smap.glink[i_x][j_x].lg[0].Vij
                jn_h  = self.smap.glink[i_x][j_x].lg[0].motif[0].get_branches()
                if Vij_h < 0.0:
                    dG_h += Vij_h
                    ctp_h = 'M'
                    btp_h = 's'
                    has_cap = True
                    if flag_debug:
                        print "lookupCap, P: ", ctp_x
                        print "Vij_h: %8.2f, (%s)" % (Vij_h, ctp_h)
                        print "len lg: ", len(self.smap.glink[i_x][j_x].lg)
                        print self.smap.glink[i_h][j_h].lg[0].motif[0].show_Motif()
                        print self.smap.glink[i_x][j_x].lg[0].motif[0].show_Motif()
                        print "dG_h: ", dG_h, ", jn_h: ", jn_h
                        # sys.exit(0)
                    #
                    
                #
                
            elif ctp_h == 'B': # i.e., if (i_h,j_h) is already type 'I'
                Vij_h = self.smap.glink[i_h][j_h].lg[0].Vij
                btp_h = self.smap.glink[i_h][j_h].lg[0].motif[0].get_btp()
                jn_h  = self.smap.glink[i_h][j_h].lg[0].motif[0].get_branches()
                if flag_debug:
                    print "B: btp_h: ", btp_h, ", jn_h: ", jn_h
                #
                
                if Vij_h < 0.0:
                    dG_h += Vij_h
                    ctp_h = 'B'
                    has_cap = True
                    if flag_debug:
                        print "lookupCap, B: ", ctp_h
                        print "Vij_h: %8.2f, (%s)" % (Vij_h, ctp_h)
                        print "len lg: ", len(self.smap.glink[i_h][j_h].lg)
                        print self.smap.glink[i_h][j_h].lg[0].motif[0].show_Motif()
                        print "dG_h: ", dG_h, ", jn_h: ", jn_h
                        # sys.exit(0)
                    #
                        
                #
                
                """@
                
                I don't understand what I was trying to do below, even
                in principle. It should __always__ be the case that
                
                jn_h ==>  (i_h,j_h) 
                
                Likewise, if it were a parallel stem that we were
                looking at, then the order should be such that
                
                (i_h,j_h), jn_x | i_h < i_x and j_h < j_x
                
                It is conceivable that the later was the original
                intent, but what I have below would not achieve it! I
                am completely at a loss.
                
                As far as I can tell, the changes above have not
                changed the result, and that should have been the case
                anyway.
                
                """
                
                # VVVVVVVVVVVVVVV  original code  VVVVVVVVVVVVVVVV
                """
                if ((jn_h[0][0] < i_h) and (jn_h[0][1] < j_h)) and \
                   not ((jn_h[0][0] == i_x) and (jn_h[0][1] == j_x)):
                    if Vij_h < 0.0:
                        dG_h += Vij_h
                        ctp_h = 'I'
                        has_cap = True
                        if flag_debug:
                            print "lookupCap, B: ", ctp_h
                            print "Vij_h: %8.2f, (%s)" % (Vij_h, ctp_h)
                            print "len lg: ", len(self.smap.glink[i_h][j_h].lg)
                            print self.smap.glink[i_h][j_h].lg[0].motif[0].show_Motif()
                            print "dG_h: ", dG_h, ", jn_h: ", jn_h
                            # sys.exit(0)
                        #
                        
                    #
                    
                #
                """
                # AAAAAAAAAAAAAAA  original code  AAAAAAAAAAAAAAAA
                
            #
            
        #
        
        return has_cap, ctp_h, btp_h, dG_h, jn_h
    #
    
    def is_connected_aaStem(self, slen1, ph1, qh1, slen2, pt2, qt2):
        """@
        
        For chromatin, this is basically just an empty function to
        passify the torturous scanning routines in
        Vienna2TreeNode. Other systems like RNA, require all this
        parsing to establish what is a connected stem, but at least
        for 5 kbp resolution chromatin, this is completely irrelevant.

        """
        return False
    #
    
    
    # see if there is a Motif of type Stem at (i,j)
    def find_ifStem(self, i, j, dGij, btp, DEBUG_find_ifStem):
        pairs_aa = []
        stem_aa = []
        pairs_pp = []
        stem_pp = []
        
        if DEBUG_find_ifStem:
            print "Enter find_ifStem{ij=(%d,%d), dGij=%8.2f, btp(%s)}" % (i, j, dGij, btp)
        #
        
        
        # search for an antiparallel connection
        flag_pairs_aa = True
        k = 1
        while flag_pairs_aa:
            # i+k -?- j-k 
            #  ..........
            # i+2 -?- j-2
            # i+1 -?- j-1
            #   i --- j
            
            iaa = i + k; jaa = j - k
            dGp_aa = self.btype[iaa][jaa].dGp # self.hv[iaa][jaa]
            if DEBUG_find_ifStem:
                print "(%d,%d): fe.hv = " % (iaa, jaa), self.btype[iaa][jaa].dGp, \
                    "btype.wt = ", dGp_aa
            #
            if dGp_aa < 0.0 and (jaa - iaa) > (1 + self.minLoopLen):
                btpx = self.btype[iaa][jaa].btp # was self.get_bondtype(dGp_aa)
                if btpx == 's':
                    btpx = 'sa'
                else:
                    # have to enforce this, even if the prediction is
                    # 't' from get_bondtype()!
                    btpx = 'c'
                #
                
                dGijaa  = dGp_aa
                # 190524 was self.calc_dG(iaa,  jaa, dGp_aa, self.T, "find_ifStem 1")
                
                
                pairs_aa += [((iaa, jaa), dGijaa, btpx)]
                if DEBUG_find_ifStem:
                    print "ijaa: ", iaa, jaa
                k += 1
            else:
                flag_pairs_aa = False
            #
        #
        if len(pairs_aa) > 0:
            btpx = btp
            if not (btp == 's' or btp == 'sa'):
                btpx = 'c'
            else:
                btpx = 'sa'
            #
            pairs_aa = [((i, j), dGij, btpx)] + pairs_aa
            
            if DEBUG_find_ifStem:
                print "pairs_aa: ", pairs_aa
            #
            # Compute the total free energy of the motif.  In the case
            # of antiparallel stems, they are joined by an
            # antiparallel connector like 'B', 'I', 'K', 'M' or
            # 'W'. Therefore, the _last item_ on the list should NOT
            # be computed from the free energy of a closing loop, as
            # this will result in double counting the closing region.
            dGaa = 0.0
            for k in range(0, len(pairs_aa)):
                dGaa += pairs_aa[k][1]
                if DEBUG_find_ifStem:
                    print pairs_aa[k][0], pairs_aa[k][1]
                #
            #
            kstm = len(pairs_aa)
            iaa = i + kstm -1 ; jaa = j - kstm + 1
            if DEBUG_find_ifStem:
                print "dGaa: %8.2f" % dGaa
                print "ijaa: (%2d,%2d)" % (iaa, jaa)
            #
            has_cap, ctpaa, btpaa, dGaa, jnaa \
                = self.lookupCap(iaa, jaa, dGaa, DEBUG_find_ifStem)
            
            stem_aa = [(i, j), dGaa, pairs_aa]
            if DEBUG_find_ifStem:
                i_t = i;   j_t = j
                i_h = iaa; j_h = jaa
                print "aa ij_t (tail): ", i_t, j_t, dGaa
                print "aa ij_h (tail): ", i_h, j_h, ctpaa, btpaa 
                print "aa ijaa:        ", iaa, jaa
                print stem_aa
                # if ctpaa == 'I' or ctpaa == 'M':
                #     print "planned exit"
                #     sys.exit(0)
                # #
            #
        #
        
        use_parallel_option = False # True # 
        if not use_parallel_option:
            if DEBUG_find_ifStem:
                print "Exit find_ifStem{ij=(%d,%d)}" % (i, j)
                print pairs_aa, pairs_pp
                # sys.exit(0)
            #

            """@
            
            190122: So, although in essence, it only returns stem_aa
            
            I opted for this because the parallel stem issue
            introduces a great deal of complication to the
            lookupBranches() and searchForMBL(). In particular, I
            noticed clashes when I was trying to fith the structure
            chr1_10751707_10970928_res5kb.heat in the test
            examples. This was particularly exemplary around
            [(25,31),(26,30)]ap vs [(26,30),(27,31),(28,32)]pp, and
            [(35,42),(36,41)]ap vs [(36,41),(37,42),(38,43)]pp.
            
            To remove those clashes, I would have to create a mask
            for the region involving parallel stems so that the
            antiparallel stems do not clash. There is even more
            complication because fitting the parallel stem starts at
            the base and works rightward. So the base is
            ".(.....).." and the next point is ".([....)].". I would
            have to supply a description where the (i_t,j_t) and
            (i_h,j_h) are described as a block [i_t - j_h] to
            specify the area encompassed by the parallel stem
            properly. If multiple parallel stems are involved, the
            grid would have to be expanded to [i_t_1 - j_h_n]. It may be possible 
            
            So it seems that the best workaround is to just presume a
            complex parallel k-type pseudoknot. For example,
            testI-2Sp.ss/testI-2Sp.heat:
            
            when "use_parallel_option = True"
            
              total free energy:  -13.837
              number of pairs:  9
              (  2,  20)[ 0][K-bgn  ] dGij_B =    0.000
              (  2,  13)[ 0][S-bgn  ] dGij_B =    0.000
              (  2,  13)[ 0][S-sp   ] dGij_B =   -2.767
              (  3,  14)[ 0][S-sp   ] dGij_B =   -2.767
              (  4,  15)[ 0][S-sp   ] dGij_B =   -2.767
              (  2,  13)[ 0][S-end  ] dGij_B =    0.000
              (  9,  20)[ 0][l-sp   ] dGij_B =   -2.767
              (  8,  19)[ 0][l-sp   ] dGij_B =   -2.767
              (  7,  18)[ 0][l-sp   ] dGij_B =   -2.767
              (  2,  20)[ 0][K-end  ] dGij_B =    0.000
  
              > 00001    dG =  -16.604   p =   0.98910511
              ccxxxccxxxcccyyyccyyycc
              ..<AB..EDC...>ab..edc..
            
            when "use_parallel_option = False"
            
              total free energy:  -16.604
              number of pairs:  8
               (  2,  20)[ 0][K-bgn  ] dGij_B =    0.000
               (  2,  13)[ 0][B-s    ] dGij_B =   -2.767
               (  9,  20)[ 0][l-sp   ] dGij_B =   -2.767
               (  8,  19)[ 0][l-sp   ] dGij_B =   -2.767
               (  7,  18)[ 0][l-sp   ] dGij_B =   -2.767
               (  4,  15)[ 0][l-sp   ] dGij_B =   -2.767
               (  3,  14)[ 0][l-sp   ] dGij_B =   -2.767
               (  2,  20)[ 0][K-end  ] dGij_B =    0.000
             
              > 00001    dG =  -16.604   p =   0.98910511
              ccxxxccxxxcccyyyccyyycc
              ..(DC..BA<...)dc..ba>..

            so the output is essentially the same, but the way the
            structure is represented of course quite different.
            
            For the sake of the stability of the program and avoiding
            as much additional complication as possible, this latter
            approach appears to make better sense presently.
            
            so much of having a more elegant style of depiction and so
            forth.
            
            """
            return stem_aa, stem_pp
        #
            
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # search for an parallel connection
        flag_pairs_pp = True
        k = 1
        i_t = i;     j_t = j
        if DEBUG_find_ifStem:
            print "ijpp: (%2d,%2d)" % (i, j)
        #
        while flag_pairs_pp:
            # we only search backwards with this because we go over
            # old solutions.
            ipp = i   + k; jpp = j   + k
            if jpp == self.N:
                flag_pairs_pp = False
                break
            #
            dGp_pp = self.btype[ipp][jpp].dGp # self.hv[ipp][jpp]
            if dGp_pp > 0.0 and (jpp < self.N):
                btpx = self.btype[ipp][jpp].btp # was self.get_bondtype(dGp_pp)
                if btpx == 's':
                    btpx = 'sp'
                else:
                    # have to enforce this, even if the prediction is
                    # 'c' from get_bondtype()!
                    btpx = 't'
                #
                dGijpp  = dGp_pp
                # 190524 was self.calc_dG(ipp, jpp, dGp_pp, self.T, "find_ifStem 2")
                pairs_pp += [((ipp, jpp), dGijpp, btpx)]
                if DEBUG_find_ifStem:
                    print "ijpp: (%2d,%2d)" % (ipp, jpp)
                #
                k += 1
            else:
                flag_pairs_pp = False
            #
        #
        stmlenpp = k 
        
        if len(pairs_pp) > 0:
            btpx = btp
            if not (btp == 's' or btp == 'sp'):
                btpx = 't'
            else:
                btpx = 'sp'
            #
            pairs_pp = [((i, j), dGij, btpx)] + pairs_pp
            
            dGpp = 0.0
            for k in range(0, len(pairs_pp)):
                dGpp += pairs_pp[k][1]
                if DEBUG_find_ifStem:
                    print pairs_pp[k][0], pairs_pp[k][1]
                #
            #
            kstmpp = len(pairs_pp)
            ipp = i + kstmpp; jpp = j - kstmpp
            if DEBUG_find_ifStem:
                print "dGpp: %8.2f" % dGpp
                print "ijpp: (%2d,%2d)" % (ipp, jpp)
            #
            i_t = i;          j_t = j
            i_h = i + kstmpp; j_h = j + kstmpp
            ihh = i_h - 1;    jhh = j_h - 1
            
            dGh1pp = 0.0; 
            
            if j_t - ipp - 1 > 0:
                if not self.smap.glink[ipp][j_t - 1].lg[0].motif[0].get_ctp() == 'X':
                    # 190121: If anything is attached at (i_h,j_h),
                    # then it must be a structure that would be found
                    # at (ipp,jpp); i.e., at (i_h + 1, j_t - 1). It is
                    # a little trickier than the case for anti-
                    # parallel stems because the boundaries are ipp
                    # and j_t - 1. To understand why this is so,
                    # consider the following parallel stem with
                    # internal structure...
                    
                    #         ipp   j_t-1
                    #         v       v
                    # ...ABCDE..((..)).abcde...
                    #    ^    ^        ^    ^
                    #    |    |        |    |
                    #   i_t   ipp     j_t   jpp
                    
                    # So it is a bit messy.
                    
                    dGh1pp, cp, ctp, pks, wps = self.filter_BIKMSW_from_glink(ipp, j_t - 1)
                    if DEBUG_find_ifStem:
                        print "dGh1pp: %8.2f, (%s)" % (dGh1pp, ctp), cp, pks, wps
                    #
                #
            #
            
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # 190121: We run the risk of problems if, for example, we
            # have some further parallelism: for example,,
            
            # ...ABCDEF...GHIJKL....abcdef...ghijkl...
            #    ^     ^            ^     ^
            #    |     |            |     |
            #   i_t   ipp          j_t   jpp
            
            # I don't really have anything to manage that directly.
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if DEBUG_find_ifStem:
                # I go to greater effort because I wanted to check
                # python's sort vs a clearly known sort of exactly the
                # particular article. It seems that python sort finds
                # the same result, but I don't know what it is doing
                # or why, so I worry.
                print "pairs_pp organization"
                print "before:         ", pairs_pp
                pairs_pp = self.ins_sort_stempp(pairs_pp)
                print "after ins_sort: ", pairs_pp
                print "----------------"
                print "ij:          ", i, j
                print "ij_t (tail): ", i_t, j_t
                print "ij_h (head): ", i_h, j_h
                print "ijhh:        ", ihh, jhh
                
                # this was for a particular problem, but it may be
                # needed again.
                if i_t == 2 and j_t == 13:
                    print "btp(%d,%d) = %s" % (i, j, btp)
                    for m in range(0, len(self.smap.glink[i_t][j_t].lg)):
                        print self.smap.glink[i_t][j_t].lg[m].motif[0].show_Motif()
                    print "get_bondtype: ", self.btype[i][j].btp
                    # was self.get_bondtype(self.btype[i][j].wt)
                    # self.hv[i][j] -> self.btype[i][j].wt
                    print "all_ctcf:         ", self.all_ctcf
                    #sys.exit(0)
                #
            #
            
            pairs_pp = self.ins_sort_stempp(pairs_pp)
            # This sort operation appears to also be done identically
            # using the intrinsic python function on lists
            # 'pairs_pp.sort()'. However, it is not clear what exactly
            # the sort() function is doing in python, whereas I know
            # exactly what my sort function is doing.
            
            
            # appears to sort in ascending order when confronted with this.
            dGpp = 0.0
            for pp in pairs_pp:
                if DEBUG_find_ifStem:
                    print pp[1]
                #
                dGpp += pp[1]
            #
            
            # unlike the antiparallel stem, this one has to be ground
            # all the way through.
            
            #
            stem_pp = [(i_t, j_t), dGpp, pairs_pp]
            if DEBUG_find_ifStem:
                print "pp ij_t (tail): ", i_t, j_t, dGpp
                print "pp ij_h (head): ", i_h, j_h              
                print stem_pp
                if len(pairs_pp) > 1:
                    print "pp planned exit"
                    # sys.exit(0)
                #
            #
        #
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        if DEBUG_find_ifStem:
            print "Exit find_ifStem{ij=(%d,%d)}" % (i, j)
            print pairs_aa, pairs_pp
            # sys.exit(0)
        #
        return stem_aa, stem_pp
    #
    
    def searchForPKs(self, iz, jz, debug_searchForPKs = False):
        if debug_searchForPKs:
            print "Enter searchForPKs{ijz = (%d,%d)}" % (iz, jz)
        #
        ijzlist  = [] # scan window where we search for PKs
        ijrtlist = [] # the domains in the largest window region
        edgebase = jz
        if iz == 0:
            """####
            
            We are scanning for PKs over the entire sequence
            length. Happens once at the end of scanning over i and
            before incrementing to the next j.
            
            """
            
            ctp      = self.smap.glink[iz][jz].lg[0].motif[0].get_ctp()
            btp      = self.smap.glink[iz][jz].lg[0].motif[0].get_btp()
            V        = self.smap.glink[iz][jz].lg[0].Vij
            ijrtlist = self.smap.glink[iz][jz].lg[0].motif[0].get_branches()
            # root Pair list: this is the window in which we search
            # for potential PK structures and it differs from the
            # region (iz,jz) when len(ijrtlist) > 1.
            
            if debug_searchForPKs:
                print "ijz = (%2d,%2d)" % (iz, jz)
                print "ctp(%s), btp(%s), V(%8.2f), ijrtlist = %s" \
                    % (ctp, btp, V, ijrtlist)
            #
            
            if len(ijrtlist) == 1 and not ijrtlist[0][0] == 0:
                """###
                
                Means this has already been check and doesn't need to
                be checked again. For example, the following case.
                
                ...((((..[[..))))...]]
                
                On the other hand, because we scan from i = j- loop to
                0, this is the last position that is being
                checked. Therefore, the following case must be
                checked.
                
                ((((....))))***********
                            ^         ^
                             hot lead
                """
                return 0
            #
            
            """@@@
            
            Now we have to set up a list of (iz,jz).  If there is more
            than one domain over the span iz(=0) to jz, then (iz,jz)
            itself is a variable that contains the boundaries of the
            domains between 0 and jz. For example, suppose we have the
            following ijrtlist
            
            ijrtlist = [(2, 19), (23, 31), (35, 46)]
            
            Then our ijzlist will be the following
            
            ijzlist = [(2, 46), (23, 46), (35, 46)]
            
            """
            
            n = len(ijrtlist) - 1
            ijzlist = []
            for k in range(0, len(ijrtlist)):
                ijzlist += [(ijrtlist[k][0], ijrtlist[n][1])]
            #
            
            
            if debug_searchForPKs:
                print "ijzlist: ", ijzlist
                
                for zz in ijzlist:
                    
                    dGmin = self.traceback_mFE(zz[0], zz[1], 0, debug_searchForPKs)
                    sstr = string.join(self.opt_ss_seq, '')
                    print sstr
                    print dGmin
                #
                print "searchForPKs"
                # if jz == 46: print "xx pk exit"; sys.exit(0)
            
            #
            # return 0
        else:
            ijzlist  = [(iz, jz)] # root domain span
            ijrtlist = [(iz, jz)] # search window of root domain
        #
        if debug_searchForPKs:
            print "ijzlist:  ", ijzlist
            print "ijrtlist: ", ijrtlist
        #
        
        for k in range(0, len(ijzlist)):
            irt = ijrtlist[k][0]; jrt = ijrtlist[k][1]
            iz  = ijzlist[k][0];  jz  = ijzlist[k][1]
            
            pkaa, pkpp = self.find_best_PK(k, ijzlist, ijrtlist, \
                                           edgebase, debug_searchForPKs)
            dGbest = self.smap.glink[iz][jz].lg[0].Vij # !!!!!
            
            # anti-parallel stem results
            if not len(pkaa) < 1:
                i_pk = pkaa[0][0]; j_pk = pkaa[0][1]; dGpk = pkaa[1]
                domains = []
                for kd in range(k, len(ijrtlist)):
                    domains += [(ijrtlist[kd][0], ijrtlist[kd][1])]
                #
                link_K = Link(i_pk, j_pk, dGpk, 'K', 'sa', domains, pkaa[2])
                self.smap.glink[i_pk][j_pk].add_link(link_K)
                
                if debug_searchForPKs:
                    L = self.scan_ahead
                    if edgebase + L >= self.N:
                        L = (self.N - 1) - edgebase
                    #
                    print "pkaa: ", pkaa
                    print "ijrt  = (%3d,%3d): " % (irt,  jrt)
                    print "ijdz  = (%3d,%3d): " % (iz,   jz)
                    print "ij_pk = (%3d,%3d): " % (i_pk, j_pk)
                    print "edgebase(%d) + L(%d): %d" % (edgebase, L, edgebase + L)
                    
                    print "dGpk(%8.2f)   dGbest(%8.2f)   ddGpk(%8.2f)" \
                        % (dGpk, dGbest, (dGpk - dGbest))
                #
            #
            
            # parallel stem results
            if not len(pkpp) < 1:
                i_pk = pkpp[0][0]; j_pk = pkpp[0][1]; dGpk = pkpp[1]
                domains = []
                for kd in range(k, len(ijrtlist)):
                    domains += [(ijrtlist[kd][0], ijrtlist[kd][1])]
                #
                link_K = Link(i_pk, j_pk, dGpk, 'K', 'sp', domains, pkpp[2])
                self.smap.glink[i_pk][j_pk].add_link(link_K)
                
                if debug_searchForPKs:
                    L = self.scan_ahead
                    if edgebase + L >= self.N:
                        L = (self.N - 1) - edgebase
                    #
                    print "pkpp: ", pkpp
                    print "ijrt  = (%3d,%3d): " % (irt,  jrt)
                    print "ijdz  = (%3d,%3d): " % (iz,   jz)
                    print "ij_pk = (%3d,%3d): " % (i_pk, j_pk)
                    print "edgebase(%d) + L(%d): %d" % (edgebase, L, edgebase + L)
                    print "dGpk(%8.2f)   dGbest(%8.2f)   ddGpk(%8.2f)" \
                        % (dGpk, dGbest, (dGpk - dGbest))
                #
            #
            
            # self.stopWhenMatchFound(idz, jdz, 125, 134, "after pk search")
            
            if debug_searchForPKs and (len(pkaa) > 0 or len(pkpp) > 0):
                # still checking this
                print ">>> after find_best_PK()"
                if len(pkaa) > 0: 
                    print "pkaa: pk stem length %2d; " % len(pkaa[2]), pkaa 
                #
                if len(pkpp) > 0:
                    print "pkpp: pk stem length %2d; " % len(pkpp[2]), pkpp 
                #
            #
        #
        if STOP_AT_FIND:
            print "planned stop"
            sys.exit(0)
        #
        
        return 0
    #
    
    def is_inside_domain(self, branches, i_pk, j_pk):
        show = False
        #if len(branches) > 1:
        #    show = True
        ##
            
        is_inside = False
        for bb in branches:
            ir = bb[0]; jr = bb[1]
            if ir < i_pk and i_pk < jr:
                if show:
                    print "ir(%d) < i_pk(%d) < jr(%d) < j_pk(%d)" % (ir, i_pk, jr, j_pk)
                    print "branches: ", branches
                #
                is_inside = True
                break
            #
        #
        return is_inside
    #
    
    def find_best_PK(self, k, ijzlist, ijrtlist, edgebase, debug = False):
        irt = ijrtlist[k][0]; jrt = ijrtlist[k][1] # scan window for PK
        iz  = ijzlist[k][0];  jz  = ijzlist[k][1]  # region between irt and edgebase
        
        flag_debug_PK = debug 
        if flag_debug_PK:
            print "Enter find_best_PK{ijz(%3d,%3d), ijrt(%3d,%3d), edgebase(%d)}" \
                % (iz, jz, irt, jrt, edgebase)
        #
        L = self.scan_ahead
        if edgebase + L >= self.N:
            L = self.N - edgebase - 1
        #

        branches = None
        dGmin = self.traceback_mFE(iz, jz, 0, flag_debug_PK)
        if self.smap.glink[iz][jz].lg[0].motif[0].get_ctp() == 'W':
            # print "found W at ijz(%d,%d)" % (iz,jz)
            branches =   self.smap.glink[iz+1][jz-1].lg[0].motif[0].get_base()
        else:
            branches =   self.smap.glink[iz][jz].lg[0].motif[0].get_base()
        #
        """
        if len(branches) > 1:
            print "found a region with %d branches" % len(branches), branches
            sys.exit(0)
        #
        """
        
        # it would probably be faster to use a modified version of
        # traceback to look for interaction points, but it is also a
        # bit more complicated.
        ssv = self.opt_ss_seq
        # provide a list of open locations for pk binding
        
        if flag_debug_PK:
            print "find_best_PK{ijrt(%d,%d), leading edge(L=%d):  %d}" \
                % (irt, jrt, L, edgebase+L)
            print string.join(self.opt_ss_seq, '')
        #
        
        best_dG = INFINITY
        best_ijp = ()
        
        best_Spklistaa = []
        last_iaa = irt
        best_Spklistpp = []
        last_ipp = jrt
        flag_find = False
        best_dGaa = 1000
        best_dGpp = 1000
        
        
        # first, scan to find the best attachment point
        for jp in range(edgebase+L, edgebase, -1):
            if flag_debug_PK:
                print "jp = %2d" % jp
            #
            
            """@
            
            this finds the best hit on each scan
            
            we scan across (irt,jrt) over the range edgebase to
            edgebase + L (in reverse order) and look for the best
            result on that scan """
            
            best_ddGaa = 1000
            best_ddGpp = 1000
            best_daa = ()
            best_dpp = ()
            
            # find best PK along iaa for given jp
            for ip in range(irt+1, jrt):
                iaa = ip
                ipp = jrt - iaa + irt
                # print "ijpp = (%3d,%3d), ijaa = (%3d,%3d)" % (ipp, jp, iaa, jp)
                
                """
                
                #############################
                first the anti-parallel case:
                #############################
                """
                
                # print "anti-parallel case:"
                
                dGp_f = self.btype[iaa][jp].dGp # dGp_f = dGp Forward
                
                # should choose something significant!!!
                if (dGp_f < self.dGpk_threshold and
                    self.is_inside_domain(branches, iaa, jp)):
                    # must satisfy the threshold and must be inside a domain
                    
                    dGijp = dGp_f
                    """190524 was
                    
                    dGijp = self.calc_dG(iaa, jp,
                                            dGp_f,
                                            self.T,
                                            "find_best_PK ap search")
                    """
                    # print "aa: dGijp = %8.2f" % dGijp
                    
                    if ssv[iaa] == '.' and ssv[jp] == '.':
                        
                        """@ 
                        
                        Currently, we just look for any slot where we
                        can cram in a link. I'm not sure this is the
                        best strategy, but anyway, that is what this
                        is doing. I think we should probably test and
                        prefer contiguous linkage in general and only
                        permit single cases when investigating systems
                        like chromatin -- where the flexibility is
                        such that this is possible.  """
                        
                        dGaa_cur = 0.0
                        for aa in best_Spklistaa:
                            dGaa_cur += aa[1]
                        #
                        
                        # the best grab is the one with the best
                        # FE. Again, I don't know this is the best
                        # kind of search, but this is the strategy.
                        
                        # find _a_ best on for a given jrt
                        if last_iaa < iaa and dGijp < best_ddGaa:
                            # print "1(aa). iaa, jp: ", iaa, jp
                            best_ddGaa = dGijp
                            best_daa   = (iaa, jp)
                            #
                        elif dGijp < best_dGaa and dGijp < dGaa_cur:
                            # print "2(aa). iaa, jp: ", iaa, jp
                            best_ddGaa = dGijp
                            best_daa   = (iaa, jp)
                            best_Spklistaa = []
                            
                        #
                    #
                #
                
                # #############################
                # now the parallel case:
                # #############################
                
                # print "parallel case:"
                
                dGp_b = self.btype[ipp][jp].dGp # dGp_b = dGp Backward
                
                # should choose something significant!!!

                if (dGp_b < self.dGpk_threshold and
                    self.is_inside_domain(branches, ipp, jp)):
                    # must satisfy the threshold and must be inside the domain
                    
                    dGijp = dGp_b
                    """190524 was

                    dGijp = self.calc_dG(ipp, jp,
                                            dGp_b,
                                            self.T,
                                            "find_best_PK pp search")
                    """
                    # print "pp: dGijp = %8.2f" % dGijp
                    if ssv[ipp] == '.' and ssv[jp] == '.':
                        # currently, we just look for any slot where
                        # we can cram a link. (See above comments on
                        # the ap case)
                        dGpp_cur = 0.0
                        for pp in best_Spklistpp:
                            dGpp_cur += pp[1]
                        #
                        
                        # the best grab is the one with the best
                        # FE. Again, I don't know this is the best
                        # kind of search, but this is the strategy.
                        
                        # find _a_ best on for a given jrt
                        if last_ipp > ipp and dGijp < best_ddGpp:
                            # print "1(pp). ijp, jp: ", ipp, jp
                            best_ddGpp = dGijp
                            best_dpp   = (ipp, jp)
                            #
                        elif dGijp < best_dGpp and dGijp < dGpp_cur:
                            # print "2(pp). ijp, jp: ", ipp, jp
                            best_ddGpp = dGijp
                            best_dpp   = (ipp, jp)
                            best_Spklistpp = []
                            
                        #
                    #
                #
            #

            # anti-parallel result
            if best_ddGaa < 0.0:
                if last_iaa < iaa:
                    i_aa = best_daa[0]; j_aa = best_daa[1]
                    btp = self.btype[i_aa][j_aa].btp
                    # was self.get_bondtype(self.btype[i_aa][j_aa].wt) 
                    if btp == 's' or btp == 'sp':
                        #print "aa pk: s or sp: ", i_aa, j_aa, btp 
                        #sys.exit(0)
                        btp = 'sa'
                        
                    else:
                        """@
                        
                        We found this in the anti-parallel
                        case. Therefore, even when
                        
                        hv[i][j] > self.htools.ctcf_tthresh 
                        
                        regardless of whether the bondtype is 't' or
                        'c', we have to force the situation (btp =
                        'c') because this is a case wherein the stem
                        is an antiparallel stem
                        
                        """
                        print "aa pk: not s or sp: ", i_aa, j_aa, btp
                        #sys.exit(0)
                        btp = 'c'
                    #

                    best_Spklistaa += [(best_daa, best_ddGaa, btp)] 
                    last_iaa      = i_aa
                    flag_find = True
                    if flag_debug_PK:
                        print "last_iaa; ", last_iaa
                        print "best_daa(%2d,%2d), best_ddGaa(%8.2f), dGijp(%8.2f)" \
                            % (best_daa[0], best_daa[1], best_ddGaa, dGijp)
                        print "best_Spklistaa: ", best_Spklistaa
                        print "ij:  ", irt, jrt
                    #
                #
            #
            
            # parallel result
            if best_ddGpp < 0.0:
                if last_ipp > ipp:
                    i_pp = best_dpp[0]; j_pp = best_dpp[1]
                    btp = self.btype[i_pp][j_pp].btp
                    # was self.get_bondtype(self.btype[i_pp][j_pp].wt)
                    if btp == 's' or btp == 'sa':
                        #print "pp pk: s or sa: ", i_pp, j_pp, btp 
                        #sys.exit(0)
                        btp = 'sp'
                    else:
                        """@
                        
                        We found this in the parallel case. Therefore,
                        even when
                        
                        hv[i][j] > self.htools.ctcf_tthresh 
                        
                        regardless of whether the bondtype is 't' or
                        'c', we have to force the situation (btp =
                        't') because this is a case wherein the stem
                        is a parallel stem.
                        
                        """
                        print "pp pk: not s or sa: ", i_pp, j_pp, btp 
                        # sys.exit(0)
                        btp = 't'
                    #
                    
                    
                    best_Spklistpp += [(best_dpp, best_ddGpp, btp)] 
                    last_ipp = i_pp
                    flag_find = True
                    if flag_debug_PK:
                        print "last_ipp; ", last_ipp
                        print "best_dpp(%2d,%2d), best_ddGpp(%8.2f), dGijp(%8.2f)" \
                            % (best_dpp[0], best_dpp[1], best_ddGpp, dGijp)
                        print "best_Spklistpp: ", best_Spklistpp
                        print "ijrt:  ", irt, jrt
                    #
                #
            #
        #
        if flag_debug_PK:
            print "find_best_PK(): finished search"
            print "best_Spklistaa: ", best_Spklistaa
            print "best_Spklistpp: ", best_Spklistpp
        #
        
        pklist_aa = []
        pklist_pp = []
        if flag_find:
            
            # if we are here, it means that we found at least one pk
            # connect.
            dGaa = INFINITY
            dGpp = INFINITY
            if len(best_Spklistaa) > 0:
                ip = irt; jp = 0
                max_j = jp
                
                dGaa = 0.0
                for aa in best_Spklistaa:
                    if aa[0][1] > max_j:
                        max_j = aa[0][1]
                    dGaa += aa[1]
                #
                jp = max_j
                dGaa += dGmin
                if flag_debug_PK:
                    print "aa, ijp: ", ip, jp
                #
                pklist_aa = [(ip, jp), dGaa, best_Spklistaa]
                
            #
            if len(best_Spklistpp) > 0:
                ip = irt; jp = 0
                max_j = jp
                
                dGpp = 0.0
                for pp in best_Spklistpp:
                    if pp[0][1] > max_j:
                        max_j = pp[0][1]
                    dGpp += pp[1]
                #
                jp = max_j
                dGpp += dGmin
                if flag_debug_PK:
                    print "pp, ijp: ", ip, jp
                #
                # to avoid duplication of the same structure found both ways
                if not (dGaa == dGpp and len(best_Spklistaa) == len(best_Spklistpp)):
                    pklist_pp = [(ip, jp), dGpp, best_Spklistpp]
            #
            if flag_debug_PK:
                print "aa:  ", pklist_aa
                print "pp:  ", pklist_pp
                #
                nlen = 4
                if len(pklist_aa) > nlen or len(pklist_pp) > nlen:
                    print "find_best_PK(): found more than %d items in pklist_aa/pp" % nlen
                    
                    if STOP_AT_FIND:
                        print "planned stop"
                        sys.exit(0)
                    #
                #
            #
        #
        if flag_debug_PK:
            print "Exiting find_best_PK(%d,%d)" % (irt, jrt)
            print "pklist_aa: ", pklist_aa
            print "pklist_pp: ", pklist_pp
            print "==============================="
            # sys.exsxit(0)
        #endif
        return pklist_aa, pklist_pp
    #
    
    
    
    
    
    # ######################################################
    # ############   pseudoknot building tools  ############
    # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    
    
    # this is currently used in display in the various trace back
    # functions: traceback_mFE, get_traces, get_traces_top, etc.
    def space(self, n):
        s = ''
        for i in range(0,n):
            s += ' '
        #
        return s
    #
    
    # traces out the structure with the minimum free energy
    def traceback_mFE(self, i, j, layer, show_structure = False):
        flag_debug = DEBUG_traceback_mFE
        if flag_debug:
            print "traceback_mFE: ij = (%d,%d), layer=%d" % (i, j, layer)
        #
        if layer > self.maxlayers:
            print "ERROR: something wrong in the recursion of traceback_mFE!"
            sys.exit(1)
        #
        
        
        # construct a secondary structure sequence
        if layer == 0:
            self.opt_ss_seq = []
            for k in range(0, self.N):
                self.opt_ss_seq += ['.']
            #
        #
        
        #if not len(self.smap.glink[i][j].lg) > 0:
        #    return 0.0 # absolutely empty
        #
        
        ctp  = self.smap.glink[i][j].lg[0].motif[0].get_ctp()
        btp  = self.smap.glink[i][j].lg[0].motif[0].get_btp()
        V    = self.smap.glink[i][j].lg[0].Vij
        join = self.smap.glink[i][j].lg[0].motif[0].get_branches()
        if flag_debug:
            print ctp, btp, V, join
        #
        s = self.space(3*layer)
        if show_structure:
            print "%s[%s](%3d, %3d)[%8.3f]: " % (s, ctp, i,j, V), join
        #
        if ctp == 'M' or ctp == 'P':
            if show_structure:
                print "%s------------------" % s
            #
            if ctp == 'M':
                self.opt_ss_seq[i] = '('
                self.opt_ss_seq[j] = ')'
            #
            for lk in join:
                i_M = lk[0]; j_M = lk[1]
                self.traceback_mFE(i_M, j_M, layer+1, show_structure)
            #
        elif ctp == 'I' or ctp == 'J':
            if ctp == 'I':
                self.opt_ss_seq[i] = '('
                self.opt_ss_seq[j] = ')'
            #
            jn = self.smap.glink[i][j].lg[0].motif[0].get_branches()
            i_I = jn[0][0]; j_I = jn[0][1]
            self.traceback_mFE(i_I, j_I, layer+1, show_structure)
            #
        elif ctp == 'S':
            flag_debug_S = False # True # 
            if flag_debug_S:
                print "inside traceback_mFE():"
                print "stem Motif: ", self.smap.glink[i][j].lg[0].motif[0].show_Motif()
            #
            
            btp = self.smap.glink[i][j].lg[0].motif[0].get_btp()
            jn = self.smap.glink[i][j].lg[0].motif[0].get_branches()
            
            for k in range(0, len(jn)):
                i_h = jn[k][0]; j_h = jn[k][1]
                if flag_debug_S:
                    print "ij_h: ", i_h, j_h # head of Stem
                #
                self.opt_ss_seq[i_h] = '('
                self.opt_ss_seq[j_h] = ')'
            if flag_debug_S:
                print "final ij_h: ", i_h, j_h
            #
            if btp == 'c' or btp == 'sa':
                self.traceback_mFE(i_h, j_h, layer+1, show_structure)
            else:
                i_h = jn[k][0]; j_h = jn[k][1]
                if flag_debug_S:
                    print "ij_h: ", i_h, j_h # head of Stem
                    print "traceback_mFE(): evaluating parallel stem"
                    print self.smap.glink[i][j].lg[0].motif[0].show_Motif()
                #
                
                stmlen = len(jn) - 1
                i_t = i;          j_t = j
                i_h = i + stmlen; j_h = j - stmlen
                if flag_debug_S:
                    print "pp ij_t (setup): ", i_t, j_t
                    print "pp ij_h (setup): ", i_h, j_h
                #
                
                if len(self.smap.glink[i_h][j_h].lg) > 0:
                    if not self.smap.glink[i_h][j_h].lg[0].motif[0].get_ctp() == 'B':
                        if len(self.smap.glink[i_h+1][j_h-1].lg) > 0:
                            self.traceback_mFE(i_h+1, j_h-1, layer+1, show_structure)
                        #
                    #
                #
                if flag_debug_S:
                    print "traceback_mFE(): stem results " 
                    print ''.join(self.opt_ss_seq)
                    #sys.exit(0)
                #
                
            #
        elif ctp == 'K':
            flag_debug_PK = False # True # 
            if flag_debug_PK:
                print "tracelink_mFE(), K:"
            #
            # the root stem
            vR    = self.smap.glink[i][j].lg[0].motif[0].get_base()
            if flag_debug_PK:
                print vR
            #
            i_R  = vR[0][0]; j_R = vR[0][1] 
            join = self.smap.glink[i_R][j_R].lg[0].motif[0].get_branches()
            self.opt_ss_seq[i_R] = '('
            self.opt_ss_seq[j_R] = ')'
            
            # trace up the root
            for lk in join:
                # whether 1 or 100 elements, they can be processed this way.
                i_M = lk[0]; j_M = lk[1]     # new ij of this subdomain
                self.traceback_mFE(i_M, j_M, layer+1, show_structure)
            #
            
            # now compute the PK
            pks  = self.smap.glink[i][j].lg[0].motif[0].get_pks()
            if flag_debug_PK:
                print "pks: ", pks
            #
            for x in pks:
                i_K = x[0][0]; j_K = x[0][1]
                if flag_debug_PK:
                    print "ijK: ", i_K, j_K
                self.opt_ss_seq[i_K] = '['
                self.opt_ss_seq[j_K] = ']'
                if show_structure:
                    print "%s[%s](%3d, %3d)[%8.3f]: " % (s, "l", i_K,j_K, x[1])
                #
            #
            if show_structure:
                print "%s---" % s
            if flag_debug_PK:
                print "pk results: from traceback_mFE()" 
                print ''.join(self.opt_ss_seq)
                # sys.exit(0)
            #
        #
        elif ctp == 'W':
            flag_debug_W = False # True # 
            # island
            wyspa = self.smap.glink[i][j].lg[0].motif[0].get_wyspa()
            if flag_debug_W:
                print "wyspa: ", wyspa
            #
            kv = 0; kmx = len(wyspa)-1
            for wk in wyspa:
                # some repetition, but who cares!
                if kv == 0:
                    self.opt_ss_seq[wk[0]] = '{'
                    self.opt_ss_seq[wk[1]] = '}'
                else:
                    self.opt_ss_seq[wk[0]] = '|'
                    self.opt_ss_seq[wk[1]] = '|'
            #
            
            # island connections
            jn_w = self.smap.glink[i][j].lg[0].motif[0].get_branches()
            if flag_debug_W:
                print "jn_w: ", jn_w
            #
            if not (jn_w[0][0] == i and jn_w[0][1] == j):
                for lk in jn_w:
                    # whether 1 or 100 elements, they can be processed this way.
                    i_W = lk[0]; j_W = lk[1]     # new ij of this subdomain
                    if len(self.smap.glink[i_W][j_W].lg) > 0:
                        self.traceback_mFE(i_W, j_W, layer+1, show_structure)
                    #
                #
            #
            
            if flag_debug_W:
                print "island results: from traceback_mFE()" 
                print ''.join(self.opt_ss_seq)
                
                nw = 5
                if len(wyspa) > nw:
                    print "traceback_mFE(): found more than %d CTCF units, they exist" % nw
                    print "ij: ", i, j
                    if STOP_AT_FIND:
                        sys.exit(0)
                    #
                #
            #
            
        else:
            self.opt_ss_seq[i] = '('
            self.opt_ss_seq[j] = ')'
            if show_structure:
                print "%s---" % s
            #
        #
        return V 
    #
    
    def find_ctcf_islands(self, i, j, best_dG, DEBUG_find_ctcf_islands = False):
        
        if DEBUG_find_ctcf_islands:
            print "find_ctcf_islands(%d,%d), best_dG = %8.2f" % (i,j, best_dG)
        #
        
        """160915wkd: CURRENT STRATEGY
        
        If I compute every case, I have i2,i3,... iN possible binding
        points. This translates into N individual, and N!/(r!(N-r)!)
        subsets of tuples: the "nCr" function on the Casio. This is a
        lot of cases as N grows!
        
        We have claimed already that we form these CTCF islands, we
        show them as being a domain and their frequency is much
        greater than the singletons.
        
        Therefore, for the moment, I will work with the heuristic that
        any region that predicts a group of CTCF interactions at point
        i or point j and within (i,j) has suffient free energy to bind
        to the respective i or j all together.
        
        Since these singletons can extend from either the i-side or
        the j-side of the _largest_ domain, and both can contribute
        significantly to the stability of the domain, I think the only
        thing we can do is just find all the CTCF connections (k,l)
        that satisfy i == k OR j == l and split up the structure into
        these sub-regions. In general, the CTCF connections are
        significantly stronger than any of the singleton connections.
        
        If we really want to know all the possible free energy results
        from all possible patterns (in all the various tuples and
        alone), then we can add some additional functionality later.
        
        At least for the relatively simple cases, this will work with
        no problems. It is also conceivable that we can test the
        solutions directly at some point. At any rate, since the free
        energy of binding the CTCFs is typically so much greater than
        the singletons, I think it is sufficient to simply ASSume that
        they all bind, presently.
        
        """
        
        keylist = self.all_ctcf.keys()
        if DEBUG_find_ctcf_islands:
            print "keylist: ", keylist
        #
        vrlps = [] # oVeRLaPS
        zone  = []
        joinW = []
        wyspa = []
        islands = []
        flag_found = False
        
        # first go through the keylist and find common contacts
        for h in keylist:
            if DEBUG_find_ctcf_islands:
                print "keylist h: i(=%2d) <= h[0](=%2d) and h[1](=%2d) <= j(=%2d)" \
                    % (i, h[0], h[1], j)
            #
            if h[0] >= i and h[1] <= j and not h == (i,j):
                if h[0] == i:
                    vrlps += [h[1]]
                    wyspa += [h]
                    flag_found = True
                elif h[1] == j:
                    vrlps += [h[0]]
                    wyspa += [h]
                    flag_found = True
                #
            #
        #
        
        # if there is are no other points inside, there is no meaning
        # in continuing this search
        if not flag_found:
            dGijh  = self.btype[i][j].dGp # self.hv[i][j]

            """190524 was 

            dGijh  = self.calc_dG(i,  j,
                                     self.btype[i][j].wt, # self.hv[i][j]
                                     self.T,
                                     "find_ctcf_islands 1")
            """
            
            dGh = 0.0
            if len(self.smap.glink[i +1][j-1].lg) > 0:
                dGh += self.smap.glink[i+1][j-1].lg[0].motif[0].Vij
            #
            dGijh += dGh
            joinW += [(i,j)]
            if j - i - 3 > 0:
                try:
                    ctp = self.smap.glink[i+1][j-1].lg[0].motif[0].get_ctp()
                    mjn = self.smap.glink[i+1][j-1].lg[0].motif[0].get_branches()
                except  (NameError, IndexError, AttributeError) as error:
                    print "ij(%2d,%2d) -> ip1jm1(%2d,%2d) failed" % (i,j, i+1, j-1)
                    #print self.smap.glink[i+1][j-1].lg[0].motif[0].show_Motif()
                    #print ctp
                    #print mjn
                    sys.exit(1)
                #
                if ctp == 'X':
                    joinW = [(i,j)]
                elif ctp == 'J' or ctp == 'P':
                    joinW = mjn
                else:
                    joinW = [(i+1, j-1)]
                #
            #    
            islands += [([(i,j)], joinW, dGijh)]
            return islands
        #
        
        wyspa += [(i,j)]
        
        
        # "list(set(vrlps)).sort()" sorts and eliminates redundant
        # items and it doeesn't matter if it is empty
        vrlps = list(set(vrlps))
        vrlps.sort()
        #
        
        if DEBUG_find_ctcf_islands:
            print "vrlps: ", vrlps
        #
        
        
        # set up the specific zones where calculation will occur.
        if len(vrlps) > 0:
            last_i = i
            for iw in vrlps:
                zone += [(last_i,iw)]
                last_i = iw
            #
            zone += [(last_i, j)]
        #
        
        if DEBUG_find_ctcf_islands:
            print "zone: ", zone
        #
        
        
        dGijh  = self.btype[i][j].dGp # self.hv[i][j]
        """190524 was
        
        dGijh  = self.calc_dG(i,  j,
                                 self.btype[i][j].wt, # self.hv[i][j]
                                 self.T,
                                 "find_ctcf_islands 1")
        """
        
        if DEBUG_find_ctcf_islands:
            print "dGijh(W) = %8.2f, dGbest = %8.2f" % (dGijh, best_dG)
        #
        dGh = INFINITY
        
        """@
        
        now go through the list of zones and compute the singleton
        contributions
        
                          common point (can happen)
                               V
                                 ..................        
                               |           ........|
                               v          v        v
             |..................          |        |       free energy
             |........         |          |        |
             v       v         v          |        |
           ..|.......|.........|..........|........|...
             i     ijh1      ijh2       ijh3       j
              ------    ------    -------    -----
        
                 <-------  loop regions ------>
        
        naturally, if there are interactions between (jh1,jh2),
        (jh2,jh3), etc., these would also show up in the growth of the
        structure
        
        """
        
        if len(zone) > 0:
            # you would have to calculate L and R separately anyway
            dGh = 0.0 
            for ijh in zone:
                ih = ijh[0]; jh = ijh[1]
                if DEBUG_find_ctcf_islands:
                    print "ijh              i(%3d) <- ih(%3d) <- jh(%3d) <- j(%3d): " \
                        % (i, ih, jh, j)
                #
                # case where we have an internal CTCF between the 
                if self.all_ctcf.has_key((ih,jh)) and i < ih and jh < j:
                    # need to save this result for the next step
                    wyspa += [(ih,jh)]
                    if DEBUG_find_ctcf_islands:
                        print "internal island: i(%3d) <  ih(%3d) <  jh(%3d) <  j(%3d)" \
                            % (i, ih, jh, j)
                    #
                #
                if (jh - ih) > 2: # 2 -> ih + 1 = jh - 1
                    #
                    if len(self.smap.glink[ih +1][jh-1].lg) > 0:
                        dGh += self.smap.glink[ih+1][jh-1].lg[0].motif[0].Vij
                    #
                    joinW += [(ih+1,jh-1)]

                    """@
                    
                    160923wkd: Even if there is nothing inside, this
                    has to be saved so that other parts of the program
                    know what to do with the domains. The problem was
                    that I had two islands
                    
                      given: wyspa  = [(38, 43), (40, 43)]
                         =>   zone  = [(38, 40), (40, 43)]
                         =>  joinW -> [(39, 39), (41, 42)]
                    
                    Although the insides are obviously empty, if you
                    don't indicate this relationship for joinW, then
                    the program just keeps cycling through (38,43) as
                    the reference until the recursion limits kill it.

                    However, recently, I found that the cases like
                    joinW = [(39,39)] are also a problem. so I
                    introduced the condition that jh - ih > 2.

                    """
                else:
                    sx  = "find_ctcf_islands; warning ctcf ij(%2d,%2d)" % (ih, jh)
                    sx += " --> proximal ligation (j(%d) - i(%d) = %d)" \
                        % (jh, ih, jh - ih)
                    print sx
                    # sys.exit(0)
                #
            #
        #
        wyspa.sort()
        if len(wyspa) > 1:
            for ijW in wyspa:
                i_W = ijW[0]; j_W = ijW[1]
                
                dGh += self.btype[i_W][j_W].dGp # self.hv[i_W][j_W]
                
                """190524 was
                
                dGh += self.calc_dG(i_W, j_W,
                                    self.btype[i_W][j_W].wt, # self.hv[i_W][j_W]
                                    self.T,
                                    "find_ctcf_islands 2")
                """
            #
        #
        
        if DEBUG_find_ctcf_islands:
            print "wyspa: ", wyspa
            print "    dGh(%8.2f) vs best_dG(%8.2f)" % (dGh, best_dG)
            if dGh < best_dG:
                print "    -- island generates more negative free energy"
            else:
                print "    -- regular structure is more stable"
        #
        
        if len(vrlps) > 0:
            joinW.sort()
            islands += [(wyspa, joinW, dGh)]
            flag_found = True
        #
        
        if flag_found:
            if DEBUG_find_ctcf_islands:
                print "full island:  zone(", zone, ")"
                for ff in islands:
                    print "             wyspa(", ff[0], ")"
                    print "             joinW(", ff[1], ")"
                    print "             dG        %8.3f" % ff[2]
                    print "             dG_best   %8.3f" % best_dG
                    print "             ddG       %8.3f" % (ff[2] - best_dG)
                print "found %d CTCF island" % len(wyspa)
            #
            
            
            # still looking for examples
            if STOP_AT_FIND:
                nzones = 4
                if len(zone) > nzones:
                    print "number of CTCFs > %d, stopping the program for reference" \
                        % nzones
                    sys.exit(0)
                #
            #
            #
        #
        #
        
        return islands
    #
    
    
#



def test0():
    # 
    fe = BranchEntropy() ## default settings
    Tds = fe.TdS(0, 1, fe.T)
    
    print "weight    enthalpy       entropy     free energy"
    print "         [kcal/mol]     [kcal/mol]    [kcal/mol]"
    
    for dw in range(0,10,1):
        dh = fe.dH(float(dw))
        print "%4d     %8.3f      %8.3f       %8.3f" % (dw, dh, Tds, (dh + Tds))
    #
#


def test1(cl):
    from LThreadBuilder     import LThreadBuilder
    from Vienna2TreeNode    import Vienna2TreeNode
    from TreeNode           import TreeNode2Motif
    
    # Evidently, Vienna() still has a few problems presently because
    # it cannot convert ".ABCD.......abcd...." properly.
    
    
    # regular secondary structure
    #          0         10        20        30        40        50        60        70        80        90
    #          |         |         |         |         |         |         |         |         |         |
    #ss_seq = "..(((((((..(((((...(((..((((..(...).))))..(((.(.....).)))..)))..))))).(((....)))..))))))).."
    
    # parallel stems and mixtures
    #          0         10        20        30        40        50        60        70        80        90
    #          |         |         |         |         |         |         |         |         |         |
    #ss_seq  = ".ABCD....EFG...efg.....HIJ...hij....abcd.."
    
    # with PK contacts
    #          0         10        20        30        40        50        60        70        80        90
    #          |         |         |         |         |         |         |         |         |         |
    #ss_seq = ".(((.(((((((.[[...)))))..((((.]]...))))..))))).....([)]([)].......(..)..((..)).(..)..........................."
    
    # with CTCF contacts  
    #          0         10        20        30        40        50        60        70        80        90
    #          |         |         |         |         |         |         |         |         |         |
    ss_seq  = "(.(((((......))))).((((.....)))).)"
    #            (2, 17)          (19,31)
    
    be = BranchEntropy()
    
    vs = Vstruct()
    print cl
    if len(cl) > 1:
        ss_seq = cl[1]
    #
    try:
        print "main: input structure sequence:"
        print ss_seq
    except(UnboundLocalError):
        print "ERROR, ss_seq is not assigned -- you idiot!"
        usage()
        sys.exit(1)
    #
    vs.parse_fullDotBracketStructure(ss_seq, True)
    be.add_hv(vs)
    be.add_btype(vs)
    
    # print "planned exit after running convert_CTCFstruct"; sys.exit(0);

    v2t = Vienna2TreeNode(vs)
    v2t.vienna2tree()
    print "main:"


    t2m = TreeNode2Motif(v2t)
    t2m.visit(v2t.genTree)
    print v2t.MPlist
    t2m.post_graftMP()

    tf = LThreadBuilder(v2t)
    tf.visit(v2t.genTree)
    print "LThread notation: "
    tf.disp_lt()
    
    # print "planned exit after running vienna2tree"; sys.exit(0);
    t2m = TreeNode2Motif(v2t)
    t2m.visit(v2t.genTree)
    print v2t.MPlist
    t2m.post_graftMP()
    print len(t2m.smap.glink)
    
    
    mbl_M = MBLptr(0,33)
    mbl_M.addBranch(2,17)
    mbl_M.addBranch(19,31)
    mbl_M.n = len(mbl_M.Q)
    mbl_M.V = be.cle_MloopEV(0,33, mbl_M, l2m.smap, True)
    print mbl_M.V
    
    #i = 2; j = 17
    i = 0; j = 33
    mbl_H = MBLptr(i,j)
    mbl_H.addBranch(i,j)
    mbl_H.n = len(mbl_H.Q)
    mbl_H.V = be.cle_HloopE(i, j, mbl_H)
    print mbl_H.V
    
    i = 0; j = 33
    p = 2; q = 17
    mbl_I = MBLptr(i,j)
    mbl_I.addBranch(p,q)
    mbl_I.n = len(mbl_H.Q)
    mbl_I.V = be.cle_IloopE(i, j, p, q, mbl_H, l2m.smap, True)
    print mbl_H.V
#   



def main(cl):
    print cl
    m = ''
    if len(cl) == 1:
        # default
        test0()
    else:
        m = cl[1]
    #
    if m == "test0":
        test0()
    elif m == "test1":
        v = [cl[0]]
        if len(cl) == 3:
            v += [cl[2]]
        #
        print v
        test1(v)
    else:
        print "%s is ignored" % cl[1]
        usage()
    #
#

if __name__ == '__main__':
    # running the program
    main(sys.argv)
    
