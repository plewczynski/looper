#!/usr/bin/env python3

"""@

Main Program:  Calculate.py 

Classes:       Calculate
               TestCalculate

Author:        Wayne Dawson
creation date: mostly 2016 and up to March 2017 
last update:   200210 (merged all current changes Kopernik and YuriV, upgraded to python3)
version:       0

Purpose:

This is perhaps the most essential part of the program, the
generalized engine for doing the dynamic programming algorithm.


Comments:

190403: Even though is is probably the most important single program
in the directory other than Thread, and the main program that I
developed in Poland in The Laboratory for Structural and Functional
Genomics, the introduction here had very little writting of this part.

"""

from math import exp
import sys
import os
import string

# main tool objects
from GetOpts import GetOpts

# Free energy parameters
# from FreeEnergy import FreeEnergy
from ChromatinModules import BranchEntropy

# Motif object representation
from Motif import Motif # a core object of these building programs
from Motif import Link
from Motif import LGroup
from Motif import Map # main map for all the FE and structures
#from Motif import Branch

from LoopRecords import Branch
from LoopRecords import MBLptr
from LoopRecords import MBLHandle
from LoopRecords import show_MBLHandle

# Other objects and tools
from FileTools   import getHeadExt
from ChPair      import ChPairData
from ChPair      import LThread2ChPair
from Chromatin2SimRNA import SimRNARestraints #from SimRNATools import SimRNAData
from BasicTools  import initialize_matrix


# for Vienna object representation
from Vienna          import Vstruct
from LThreadBuilder  import LThreadBuilder
from Vienna2TreeNode import Vienna2TreeNode
from TreeNode        import TreeNode2Motif

# ################################################################
# ######################  Global constants  ######################
# ################################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

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

PROGRAM      = "Calculate.py"  # name of the program

# Debugging in main()
SHOWMAIN     = False # True #

# debugging settings for routines in Calculate()


# debugging minFE 
DEBUG_minFE             = False # True #
# debugging traceback_mFE
DEBUG_traceback_mFE     = False # True #
# debugging searchForMBL 
DEBUG_searchForMBL      = False # True # 
# debugging lookupBranches 
DEBUG_lookupBranches    = False # True # 
# debugging find_ifStem 
DEBUG_find_ifStem       = False # True # 
# debugging make_StemMotif 
DEBUG_make_StemMotif    = False # True # 
# debugging find_best_PK 
DEBUG_find_best_PK      = False # True # 
# debugging find_ctcf_islands 
DEBUG_find_ctcf_islands = False # True # 

CHECK_ALL = False # True # 
# if you want to debug all parts simultaneously, then set CHECK_ALL to True.
def debug_all(fast_track):
    
    global DEBUG_minFE             
    global DEBUG_traceback_mFE     
    global DEBUG_searchForMBL      
    global DEBUG_lookupBranches    
    global DEBUG_find_ifStem       
    global DEBUG_make_StemMotif    
    global DEBUG_find_best_PK      
    global DEBUG_find_ctcf_islands 
    
    DEBUG_minFE             = fast_track
    DEBUG_traceback_mFE     = fast_track
    DEBUG_searchForMBL      = fast_track
    DEBUG_lookupBranches    = fast_track
    DEBUG_find_ifStem       = fast_track
    DEBUG_make_StemMotif    = fast_track
    DEBUG_find_best_PK      = fast_track
    DEBUG_find_ctcf_islands = fast_track
#    

# debugging option (rarely used)
STOP_AT_FIND = False # True # 
# stops program when encouters condition


# 2. special local global function
dconnection = { 'B' : True,
                'S' : True,
                'I' : True,                
                'M' : True,                
                'K' : True,                
                'R' : True,                
                'W' : True,
                'w' : True }                

# This is not used so much now that I have introduced GetOpts.py,
# but it may still be useful at the very beginning of the program and
# possibly in other parts.


def usage():
    print ("USAGE: %s -f file.heat " % PROGRAM)
#

def kickOn(i, i_set, j, j_set):
    flag_on = False
    if i_set <= i and j <= j_set:
        flag_on = True
    #
    return flag_on
#


# 3. tests
TEST0 = False # True #
TEST1 = True # False # 

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ################################################################
# ################################################################


####################################################################
######################  MINIMUM FREE ENERGY  #######################
####################################################################

class TestCalculate(object):
    def __init__(self, seqs = ["(.(((((......))))).((((.....)))).)"]):
        self.source = "TestCalculate"
        self.set_GetOpts = True
        
        self.f_heatmap = ["faketest1.heat"] # elements are a list now
        self.T         = 310.0 # temperature is undefined!
        
        # output option
        self.p_all_1D  = False
        # In general, I don't see that I particularly even want to see
        # more than the first 10 structures in these created
        # directories. For a large calculation, it is easy to produce
        # tens of thousands of structures. For quite large structure,
        # the program can easily generate one hundred thousand
        # structures. Therefore, to avoid an explosion of files, the
        # default is to print a maximum of 50 structures. This is
        # primarily so that the program (or the computer) doesn't
        # crash due to the creation of so many files.
        
        # It might be better to have a size limit.
        
        
        self.dGrange = 10.0 # kcal/mol 
        # The search for suboptimal structures will range between
        # dGmin and dGmin + dGrange.
        
        # settings related to the PET clusters
        
        self.PETwt = 100.0  # PET wt scale
        self.N = len(seqs[0])
        # print ("N = ", self.N)
    #
#


# Calculate() is the main machine of this program
class Calculate(object):
    def __init__(self, cl = TestCalculate()):
        flag_debug =  False
        if not cl.set_GetOpts:
            print ("ERROR(Calculate): options are not properly configured")
            sys.exit(1)
        #
        
        self.flnm = cl.f_heatmap[0]
        
        if cl.source == "TestCalculate":
            print ("using test module, not the real thing!")
            self.fe      = BranchEntropy()
        else:
            self.fe      = BranchEntropy(cl)
        #
        
        self.T       = cl.T # temperature is undefined!
        print ("sequence length N = ", self.fe.N)
        
        # output option
        self.p_all_1D      = cl.p_all_1D
        """@ 
        
        In general, I don't see that I particularly even want to see
        more than the first 10 structures in these created
        directories. For a large calculation, it is easy to produce
        tens of thousands of structures. For quite large structure,
        the program can easily generate one hundred thousand
        structures. Therefore, to avoid an explosion of files, the
        default is to print a maximum of 50 structures. This is
        primarily so that the program (or the computer) doesn't crash
        due to the creation of so many files.
        
        It might be better to have a size limit.
        
        """
        
        self.dGrange = cl.dGrange
        # The search for suboptimal structures will range between
        # dGmin and dGmin + dGrange.
        
        # settings related to the PET clusters
        
        self.PET_wt0 = cl.PETwt                  # PET wt scale
        
        # set up the heat map "enthalpy"
        self.ctcf_setv = self.fe.ctcf_setv
        if cl.source == "TestCalculate":
            print ("(2) N = ", cl.N)
            self.N  = cl.N
        else:
            self.N  = self.fe.N
        #
        
        
        if flag_debug:
            print ("finished Calculate constructor")
        #
        
    #
    
    
    # ###################################################################
    # ###############  New I-loop and MBL search method   ###############
    # ###################################################################
    # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
    # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
    
    """@
    
    
              .|__|.
            .  p  q .
           .        .
            .      .
             . __ . 
             i|  |j
    
       i ..... (p, q) .... j
    
    
    
    from searchForMBL, the argument is always of the following form:
    (i, j, i, i + k) or (i, j, k + 1, j). Therefore, one of the
    boundaries is always either i or j in this call. The thing is, we
    are looking at a "window" (p,q). We don't care about anything else
    except that small window between (p,q) inclusive. Nevertheless, we
    sometimes need to know the window in which this unit
    fits. Therefore, i,j is listed separate from (p, q). Just remember
    that (p,q) has one i or on j and i<=p<q<=j.
    
    """
    def lookupBranches(self,
                       i,  j,    # boundaries of search
                       p,  q,    # sector of search
                       mblplnk): # class MBLptr
        
        global DEBUG_lookupBranches
        
        # Q is an indefinite description MBL-like structure buffer 
        found_a_Q = False # found a value for Q
        dangleFE = 0.0
        
        """@
        
        _effectively_ reset mblplnk; doesn't remove prior values, but
        it erases the tag indicating how many there are and what the
        free energy of the content was in the previous fill. Since we
        rely on the entry 0 for the total FE and the number of
        branches, this is effectively similar to "erasing" a file on a
        HD. By erasing all markers about the file, the space is made
        open for writing again, or ignoring; whichever is needed at
        the time.
        
        """
        
        if DEBUG_lookupBranches:
            print ("enter lookupBranches{ij(%d,%d), pq(%d,%d)}" % (i, j, p, q))
            print ("      btype{pq(%2d,%2d)}.pair = %d" % (p, q, self.fe.btype[p][q].pair))
        #endif
        
        if (q - p <= 0) or (j - i <= 0):
            # hopefully, shouldn't happen, but if it does, there is
            # something seriously wrong with the call.
            print ("ERROR lookupBranches{ij(%d,%d), pq(%d,%d)}" % (i, j, p, q))
            if q - p <= 0:
                print ("p(%d) > q(%d)!!!!" % (p, q))
            if j - i <= 0:
                print ("p(%d) > q(%d)!!!!" % (i, j))
            sys.exit(0)
        #
        
        mctp = 'X'; mbtp = '-'; mV   = INFINITY; mjn  = [];  mbs  = []
        sctp = 'X'; sbtp = '-'; sV   = INFINITY; sjn  = [];  sbs  = []
        
        """@
        
        181212(190725): Unlike the vsfold version, here we have the
        possibility of degenerate solutions. These are contained in
        objects from class Link. Therefore, eventually, this will have
        to be constructed so that it processes more than one
        link. Because of the multiple states of similar free energy
        and a rather simple evaluation function, chromatin, in
        particular, may require more consideration than RNA or
        proteins would. On the other hand, it is generally rather
        difficult to get degenerate solutions in general (except for
        contrived examples) and even this likelihood is very small for
        typical RNA and protein sequences.
        
        At the current juncture, since I am just coming back to these
        modules after a long haitus, I will focus on making the
        primary link work; i.e., self.fe.smap.glink[p][q].lg[0]
        
        With this chromatin program, the nature of the processing unit
        has changed drastically. From within this much higher object
        of class Calculate(), we now have an object smap from class
        Map() that contains all objects (lg) of class LGroup that
        contain the main structurl features at any position (i,j). The
        object lg from class LGroup serves as a container for the
        object glink from class Link(). Each lg object is a container
        that stores structures possessing the same _minimum_ free
        energy but with different structural configurations. Then
        lg[k].glink is a container for the actual objects motif from
        class Motif(). The motifs contain the actual connections
        between different structural motifs in smap. These can be
        assembled as various fragments and only require the expected
        pattern in the neighborhood of (i,j).
        
        After having worked with this for some time, I am convinced
        that the better approach is to store the degenerate structures
        within Link and essentially make link specify the Motifs. So,
        in fact, lg[0] is the only case, even though I am still
        parsing lg[0].motif.  Just as in the case of vsfold, the main
        search is directed to the minimum free energy; hence, only the
        topmost structure is tracked down in this initial
        search. Therefore, in this stage, we basically look only at
        lg[0].motif[0] in this section.

        There is a small wrinckle. if the solution is 'S', this
        information is stored without any dangling interactions, but
        the others are not. This is because a solution for 'S' must be
        tacked onto a stem, so we don't want to double count false
        free energies. On the other hand, this dangle contribution can
        be significant enough that at the immediate position (i,j),
        the better solution is I, M etc, yet when this is joined with
        the stem, the 'S' solution is better.
        
        """
        
        
        m_mtf = None
        N = self.fe.smap.N
        # self.show_smap_all(N, p, q)
        
        if len(self.fe.smap.glink[p][q].lg) == 0:
            # print ("lg(%d,%d) = 0" % (p, q))
            mblplnk.V = self.fe.cle_HloopE(p, q, mblplnk)
            mblplnk.addBranch(p, q)
            mblplnk.nm = 'X'
            mblplnk.dr = '-'
            # print (mblplnk)
            return mblplnk
        #
        
        if self.fe.smap.glink[p][q].lg[0].motif[0].get_ctp() == 'X':       
            # print ("lg(%d,%d).motif.ctp = X" % (p, q))
            mblplnk.V = self.fe.cle_HloopE(p, q, mblplnk)
            mblplnk.addBranch(p, q)
            mblplnk.nm = 'X'
            mblplnk.dr = '-'
            # print (mblplnk)
            return mblplnk
        #
        
        # This is a rather goofy approach, but ultimately it expresses
        # the contents of self.fe.smap.glink[p][q].lg[0].motif[0]
        for mtfk in self.fe.smap.glink[p][q].lg[0].motif:
            m_mtf = mtfk
            if mtfk.get_ctp() == 'M' or mtfk.get_ctp() == 'P':
                mctp = mtfk.get_ctp()
                mbtp = mtfk.get_btp()
                mV   = mtfk.get_Vij()
                mjn  = mtfk.get_branches()
                mbs  = mtfk.get_base()
            elif mtfk.get_ctp() == 'B':
                # these are in many ways the same as S but still I am
                # concerned that they may need to be treated different
                # from 'S'
                mctp = mtfk.get_ctp()
                mbtp = mtfk.get_btp()
                mV   = mtfk.get_Vij()
                mjn  = [(p, q)] # mtfk.get_branches()
                mbs  = mtfk.get_base()
            elif mtfk.get_ctp() == 'I' or mtfk.get_ctp() == 'J':
                # these are in many ways the same as S but still I am
                # concerned that they may need to be treated different
                # from 'S'
                mctp = mtfk.get_ctp()
                mbtp = mtfk.get_btp()
                mV   = mtfk.get_Vij()
                mjn  = mtfk.get_branches()
                mbs  = mtfk.get_base()
            elif mtfk.get_ctp() == 'R' or mtfk.get_ctp() == 'K' or mtfk.get_ctp() == 'W':
                # right now, I have defined these the same as 'M', but
                # I am still mulling over whether I should treat them
                # different from 'S'
                mctp = mtfk.get_ctp()
                mbtp = mtfk.get_btp()
                mV   = mtfk.get_Vij()
                mjn  = mtfk.get_branches()
                mbs  = mtfk.get_base()
            elif mtfk.get_ctp() == 'S':
                # presently, the distinction is whether it is a stem
                # or something else.
                sctp = mtfk.get_ctp()
                sbtp = mtfk.get_btp()
                sV   = mtfk.get_Vij()
                sjn  = mtfk.get_branches()
                sbs  = mtfk.get_base()
            else:
                print ("ERROR: unrecognized structure type (%s)" % mtfk.get_ctp())
                sys.exit(1)
            # in the future, this should be continued to assign throughout
        #
        
        if DEBUG_lookupBranches:
            print ("s/mctp[%s,%s], s/mbtp[%s,%s]: " % (sctp, mctp, sbtp, mbtp))
        #
        
        ss = ''
        if sctp == 'S':
            if N <= 0:
                print ("lookupBranches N = %d\n" % N)
                sys.exit(1)
            #
            
            # ss contains stem info, mm contains mbl
            # optimizations. However, ss is without the dangles
            
            dangleFE = self.fe.stem_dangle(0,  N,      # boundaries (often [1 to N]) 
                                           p,  q)      # 5'3' tail of the stem
            
            if sV <= (mV - dangleFE): 
                ss = 'S'
            #
            if DEBUG_lookupBranches:
                print ("pq(%d,%d), compare cc[%8.2f] <= {fML[%8.2f] - dangleFE(%8.2f) = %8.2f}" \
                    % (p, q, sV, mV, dangleFE, mV - dangleFE))
                if sV <= (mV - dangleFE): 
                    print ("[%s](%8.2f) vs %8.2f: stem wins the day\n" % (ss, sV, mV))
                else:
                    print ("[%s](%8.2f) vs %8.2f: mbl is best choice\n" % (ss, mV, sV))
                #
            #
            
        #
        
        if ss == 'S': #  found stem tail at (p,q) 
            # best structure at pq was found to be a stem
            if DEBUG_lookupBranches:
                print ("pq(%d,%d): ss[%s][dG(%8.2f)], mm[%s][dG(%8.2f)]" \
                    % (p, q, sctp, sV, mctp, mV))
            #
            best_cS = INFINITY;
            flag_pass_stem = True;
            
            if sV == INFINITY:
                # this is just to make sure that cc has not been used somewhere!!!
                print ("ERROR: lookupBranches{ij(%d,%d), pq(%d,%d)}" % (i, j, p, q))
                print ("       (%2d,%2d): ss[%s][dG(%8.2f)], mm[%s][dG(%8.2f)]" \
                    % (p, q, sctp, sV, mctp, mV))
                print ("has not been assigned!")
                sys.exit(1)
            #
            
            
            # 180609: maybe should be testing sV. So far, the program
            # has not stopped here, though I haven't pushed it so much.
            if mV > 0:
                if DEBUG_lookupBranches:
                    print ("fML[(%d,%d)](%8.2f) > 0" % (p, q, mV))
                #
                
                if mV >= INFINITY:  # this is less likely now.
                    # Personally, I think we should reject stems that
                    # have a positive FE. They really shouldn't form
                    # unless you are looking for "possible" states.
                    slen = len(sjn)
                    L_stm = slen
                    if DEBUG_lookupBranches: # lookupBranch
                        i_xi = 4
                        # this is the more usual situation, where the Kuhn length
                        # is long and the stem is just a nubbin. We want to
                        # discourage these sorts of stems.
                        s  = "lookupBranches: i(%d) <= p (%d) < q(%d) <= j(%d)\n" \
                            % (i, p, q, j)
                        s += "pq(%3d,%3d)[%s][%5s][V=%8.2f], MBLbase = %s\n" \
                            % (p, q, mctp, mbtp, mV, mbs)
                        s += "pq(%3d,%3d)[%s][%5s][V=%8.2f],   sbase = %s\n" \
                            % (p, q, sctp, sbtp, sV, sbs)
                        s += "slen(%d), i_xi(%d), L_stm = %d,   sjn -> %s" \
                            % (slen, i_xi, L_stm,  sjn)
                        print (s)
                    #endif
                    #print ("planned exit")
                    #sys.exit(0);
                #
            #
            
            # 180424: NOTE!!! Here, we should be looking at fML[pq] because
            # this element contains the dangle information and the stem bend
            # cle weight. The other elements c[pq].dG is merely the stem
            # energy (without cle stem bend or dangle) and likewise FPK[pq],
            # the dangle information is added with the function stem_dangle()
            # and terminal issues corrected when the result passes through
            # trim_53pStems(). Anyway, for the purpose of using searchForMBL,
            # we should use fML[pq] because it contains the dangle
            # corrections and the stem corrections for cle.
            
            best_cS = sV + dangleFE; # since dangle typically improves the FE
            
            
            if flag_pass_stem:
                mblplnk.nm   = sctp
                mblplnk.dr   = sbtp
                mblplnk.V    = best_cS
                mblplnk.addBranch(p,q) 
                mblplnk.n    = len(mblplnk.Q)
                found_a_Q = True;
                if DEBUG_lookupBranches:
                    s  = "lookupBranches: ij(%d,%d), pq(%d,%d)[%s]: " \
                        % (i, j, p, q, mblplnk.nm)
                    s += "best_cS = %8.2f (index = %d)\n" % (mblplnk.V, mblplnk.n) 
                    s += self.showBranches(mblplnk, "mblplnk")
                    print (s)
                #endif
            
            else:
                if DEBUG_lookupBranches:  
                    # ss_name
                    print ("stem(%d,%d)[%s] skipped!" \
                        % (p,   q,   self.fe.smap.glink[p  ][q  ].lg[0].motif[0].get_btp()))
                    print ("stem(%d,%d)[%s]" \
                        % (p+1, q-1, self.fe.smap.glink[p+1][q-1].lg[0].motif[0].get_btp()))
                #endif    
            #
            # sys.exit(0)
            
        elif mctp == 'B' or mctp == 'I' or mctp == 'J':
            
            # 180425: this is essentially the same as 'M' and 'P' but
            # it involves only one branch. Also, 'J' was essentially
            # the same as an pMBL, but here it means a pMBL with only
            # one branch. This nonsense is mostly preserved because of
            # the fact that the program still has issues with I-loops
            # and M-loops, but in principle, they should be treated
            # almost exactly the same.
            best_cQ = INFINITY;
            
            if mV < best_cQ:
                best_cQ = mV
            
            else:
                print ("lookupBranches: warning (%d,%d)[Q] not assigned properly" % (p, q))
                print ("                fML = %8.2f >= %8.2f" % (mV, INFINITY));
                sys.exit(1)
            #
            
            #brnv  = self.fe.smap.glink[p][q].lg[0].motif[0].get_branches()
            mblplnk.nm = self.fe.smap.glink[p][q].lg[0].motif[0].get_ctp()
            mblplnk.dr = self.fe.smap.glink[p][q].lg[0].motif[0].get_btp()
            mblplnk.V  = best_cQ
            # note: for PKs, mbs is the base of the PK, when we detect
            # an I-loop, we must include the base, if we detect a
            # J-loop, we don't include the base, ... or should I say,
            # the base of the J-loop is the branches of the J-loop,
            # the base of the I-loop is (p,q).
            if mctp == 'I':
                mblplnk.pushBranchlist(mbs, True)  # the set of branches
                if DEBUG_lookupBranches:
                    print (self.fe.smap.glink[p][q].lg[0].motif[0].show_Motif())
                #
                
            else:
                mblplnk.pushBranchlist(mjn, True)  # the set of branches
            #
            
            mblplnk.n  = len(mblplnk.Q) # should be 1
            found_a_Q  = True;
            
            
            if DEBUG_lookupBranches:
                s  = "%s: lookupBranches: ij(%d,%d), pq(%d,%d)[%s]: " \
                     % (mctp, i, j, p, q, mblplnk.nm)
                s += "best_cI = %8.2f (index = %d)\n" % (mblplnk.V, mblplnk.n) 
                s += self.showBranches(mblplnk, "mblplnk")
                print (s)
                print ("mbs:  ", mbs)
                print ("mjn:  ", mjn)
                print ("mctp: ", mctp)
                # if p == 52 and q == 63:
                #     sys.exit(0)
                # #
                # if p == 52 and q == 59:
                #     sys.exit(0)
                # #
            #endif
            
        elif mctp == 'K' or mctp == 'R' or mctp == 'W': 
            # also search for pseudoknots
            # best structure is a pseudoknot ('K' or 'R' type)
            best_cX = INFINITY
            
            if mV < best_cX: 
                best_cX = mV
            
            else:
                print ("lookupBranches: warning (%d,%d)[K] not assigned properly" % (p, q))
                sys.exit(1);
            #
            
            mblplnk.nm = mctp
            mblplnk.dr = mbtp
            mblplnk.V  = best_cX
            mblplnk.pushBranchlist([(p,q)], True)  # the set of branches
            mblplnk.n  = len(mblplnk.Q) # should be 1
            found_a_Q = True;
            
            if DEBUG_lookupBranches:
                s  =  "lookupBranches: ij(%d,%d), pq(%d,%d)[%s]: " \
                     % (i, j, p, q, mblplnk.nm)
                s += "best_c%s = %8.2f (index = %d)\n" % (mctp, mblplnk.V, mblplnk.n) 
                s += self.showBranches(mblplnk, "mblplnk")
                print (s)
                print (self.fe.smap.glink[p][q].lg[0].motif[0].show_Motif())
                print ("mbs: ", mbs)
                print ("mjn: ", mjn)
                #sys.exit(0)
            # endif
            
        elif mctp == 'M' or mctp == 'P': 
            
            # best structure at pq was found to be a pMBL
            best_cM = INFINITY;
            if mV < best_cM:
                best_cM = mV;
            #
            else: 
                print ("lookupBranches: warning iMBL at (%d,%d) not assigned properly" \
                    % (p, q))
                print ("                fML = %8.2f >= %8.2f" \
                    % (mV, INFINITY))
                sys.exit(1);
            #
            
            mblplnk.nm = mctp
            mblplnk.dr = mbtp
            mblplnk.V  = best_cM
            # note: When we detect
            # an M-loop, we must include the base, if we detect a
            # P-loop, we don't include the base, ... or should I say,
            # the base of the P-loop are the branches of the P-loop,
            # the base of the M-loop is (p,q).
            if mctp == 'M':
                mblplnk.pushBranchlist(mbs, True)  # the set of branches
                if DEBUG_lookupBranches:
                    print (self.fe.smap.glink[p][q].lg[0].motif[0].show_Motif())
                #
                
            else:
                mblplnk.pushBranchlist(mjn, True)  # the set of branches
            #
            
            mblplnk.n  = len(mblplnk.Q) # should be 1
            found_a_Q  = True;
            
            if DEBUG_lookupBranches:
                s  =  "lookupBranches: ij(%d,%d), pq(%d,%d)[%s]: " \
                     % (i, j, p, q, mblplnk.nm)
                s += "best_cM = %d (index = %8.2f)\n" % (mblplnk.V, mblplnk.n) 
                s += self.showBranches(mblplnk, "mblplnk")
                print (s)
                print ("mbs:  ", mbs)
                print ("mjn:  ", mjn)
                print ("mctp: ", mctp)
                if False and mctp == 'M':
                    print ("exit if answer is M")
                    sys.exit(0)
                #
            #endif
            
        #
        
        
        # check if all values are > 0, means no potentially stable
        # structure found anywhere
        if DEBUG_lookupBranches:
            s = ''
            if not found_a_Q:
                s += "lookupBranches(ij(%2d,%2d),pq(%2d,%2d)) -> no I-loops or stems" \
                     % (i, j, p, q)
            #
            else:
                if mblplnk.n == 1: 
                    s += "lookupBranches(ij(%2d,%2d),pq(%2d,%2d))[%c][E=%8.2f] " \
                         % (i, j, p, q, mblplnk.nm, mblplnk.V)
                    s += "-> found %d branch\n" % mblplnk.n
                else:
                    s += "lookupBranches(ij(%2d,%2d),pq(%2d,%2d))[%c][E=%8.2f] " \
                         % (i, j, p, q, mblplnk.nm, mblplnk.V)
                    s += "-> found %d branchs\n" % mblplnk.n
                #
                s += self.showBranches(mblplnk, "mblplnk")
            #
            print (s)
        #endif  
        return mblplnk;
    #
    
    
    
    
    
    """180413: 
    
    This is the new approach to doing mbls (this module was introduced
    to chromatin package 1807x)
    
    One thing that occurs to me now is that, since ss_iMBL(i,j) will
    now incorporate both the Iloop solution and the iMBL solution,
    saving the best result, we would lose this set of information. for
    optimal strutures, this is irrelevant, but for suboptiml
    structures, we may lose information. It might be good to have both
    solutions, particularly if the difference between the iMBL and the
    Iloop are very small. So it might be better to run this mbl search
    with two variables; one for Iloops and one for iMBLs,
    
    """
    
    def searchForMBL(self,
                     i, j,          # search region
                     lmblh,         # list of type MBLHandle()
                     opt_iMBL):     # iMBL=1/pMBL=0
        global DEBUG_searchForMBL
        
        # mbl pointer link
        dG_lb = lmblh[0].dG_lb
        
        if DEBUG_searchForMBL: 
            print ("\n!!!Enter searchForMBL(): from %3d to %3d" % (i, j))
        #endif
        
        # //int mnHairPin = 2*i_xi_min + 1; 
        mnHairPin = self.fe.minLoopLen;
        
        """181211: 
        
        Now that this module serves the function of searching for both
        a single branch and multiple branches, I really have to search
        all the way across the region with no cutoff. Originally,
        because there were restrictions on multibranch loop sizes and
        I required two branches, this matter was not an issue, but it
        is an issue here. We have to scan all the way across as a
        result. So it probably should not be named minimum hair pin
        (mnHairPin) anymore, but historically, its role was to help
        minimize the loop searching region and it made sense to call
        it that.
        """
        
        """180413(Friday!): 
        
        mnHairPin: Originally, I built these functions with a pedantic
        insistence that an MBL _must_ contain at least two stems or it
        could not be called an _multi-branch_ loop. I think this
        restriction is not at all necessary. An MBL means there is _at
        least one_ branch.
        
        Therefore, the idea here is that an MBL must contain at least
        one hairpins, and in vsfold, the sequence length of one
        hairpin must be at least 3*i_xi_min [bp]. So _surely_,
        2*i_xi_min + 1 is simply too small to allow even one hairpin
        to form. I think it would even be safe to say 4*i_xi_min +
        2. The "2" is there because the scan is always from i+1 to
        j-1.
        
        """
        
        if (j - i) < mnHairPin:
            # I know it looks absurd when mnHairPin = 1, but this is
            # because we have to fill smap.glink with something.
            
            """181215: 
            
            comments added
            
            ######################################################
            In this design, this should be the first time that
            position (i,j) has been looked at, so I don't think we
            have to check this at all. In general, the thing that will
            have to be added to this method is a way for handling
            degeneracy, but that most certainly does not matter here.
            ########################################################

            """
            
            new_MBL    = MBLptr(i,j)  # default tags 
            new_MBL.V  = self.fe.cle_HloopE(i, j, new_MBL)
            new_MBL.addBranch(i, j)
            new_MBL.nm = 'X'
            new_MBL.dr = '-'
            #
                
            """@
            
            btp is always 'sa' in this case, since 'sp' and 'sa' are
            indistinguishable for a single bond in a single loop in
            the _current_ problems. If this comes to proteins, this
            would have to be decided.
            
            Admittedly, the job is already done, but both vsfold and
            this new approach expect something to be returned, so we
            return something. Actually, this routine, in general, is
            looking up the best I-loop and M-loop, so it is still only
            one of several possibilities, at this stage. For example,
            it could still be that a PK wins. This routine will figure
            out that a stem is the best solution [(i,j), (i+1,j-1)] or
            [(i,j), (i+1,j+1)]. Still, a PK or an island could compete
            with this structure. So we pass the best secondary
            structure solution, then compare it with PK and island,
            and make a final decision. As a result, I think maybe in
            addtion to Vij and mjn, we need to return ctp and btp.
            
            """
            """181217:
            
            vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
            I decided that since the best way to pass the variable in
            the latter part of this routine is as the object new_MBL
            (from class MBLptr), we should do the same here. The final
            processing after exiting this routine is as follows:
            
               link_H = Link(i, j, Vij, 'B', 'sa', mjn)
               self.fe.smap.glink[i][j].add_link(link_H)
            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            
            """
            
            k = self.fe.dSearchK['SB']
            lmblh[k].update_mbl(new_MBL, i, j, DEBUG_searchForMBL, "first juncture")
            return lmblh  # ... no viable solution found
        
        #
        if DEBUG_searchForMBL:
            print ("Evaluate [%d to %d]" % (i, j))
        #
        
        
        # These are all variables that need to be available outside of
        # the mbl search loop that will start subsequently.
        
        p  = i; q  = j;  # pq = ssIndex(p,q);
        
        new_cP = INFINITY;
        new_iMBL_dG = INFINITY
        lmblh[0].best_iMBL_dG = INFINITY
        
        # temporary storage objects 
        mblplnk1 = MBLptr(i, j)  
        mblplnk2 = MBLptr(i, j)  
        new_MBL  = MBLptr(i, j)
        
        
        Mloop_wt = 0.0;
        if opt_iMBL:

            """180421: 

            For the purposes I am using this tool, I cannot think of
            an instance where this should be applied for this
            function. It is basically a holdover from the early days
            of vsfold5 where it was used for both iMBLs and pMBLs.
            
            At any rate, because the "closing point" of this function
            is (i,j), we only have to do this task once, here.
            
            %%%% 051021wd (updated 180102wkd)
            
            The operations below attempt to distinguish between
            structures in a pMBL and structures in an iMBL.  There
            should be differences because the branches in a pMBL are
            simply optimizing their FE, but the iMBL must also close
            around in a loop.  I still think that the cost is not only
            in the formation and size of the loop, but there are
            probably additional costs that set a lower bound for the
            smallest allowed loop size. Otherwise, the best iMBL is
            the short loop; a minimum loop size.  Nevertheless, it
            seems that exactly how to weight these things is a bit
            more complicated than I expected.

            """
            
            Mloop_wt = self.fe.cle_MloopEV(i, j, new_MBL, self.fe.smap, opt_iMBL)
        #
        
        lmblh[0].Mloop_wt = Mloop_wt
        
        
        no_V1 = False
        no_V2 = False
        no_Vf = False
        found_one = False
        
        """@
        
        Scans everything from 
        
        from {V(i,i+1),V(i+2,j)} to {V(i,j-2),V(j-1,j)}. 
        
        Of course V(i,i+1) and V(j-1,j) are obviously not taken
        seriously and are effectively thrown out. The point is that we
        achieve most of the coverage. The remaining missing indices
        
        {V(i,j-1),V(i+1,j-1),V(i+1,j),V(i,j)} 
        
        are calculated in the next step. Of course V(i,j) would be a
        stem or a PK if it is the optimal object.
        
        """
        
        # !!!VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
        on_ij = False
        # on_ij True/False is the search region (i,j), (i+1,j), (i+1,j-1), (i,j-1)
        
        # 190130: this should be removed, eventually
        # !!!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        
        for k in range(1, j - i - 1):
            # (int k = 1; k <= j - i - 2; k++) { } # (j - i - 2 + 1)
            no_V1 = False
            no_V2 = False
            no_Vf = False

            """@
            
            Here, rather than consider whether there are two or more
            branches, one or more branches are saarched for. I don't
            know why it took me so long to realize it, but whether one
            or more branches, the behavior should not really diffr all
            that much. We look at both sector (i,i+k) and (i+k+1,j)
            essentially simultaneously.  

            """
            
            # ############################################
            # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
            if i+k > j:
                # This problem really _shouldn't_ occur as far as I can tell;
                # however, for the moment, I think it is better to have this
                # trap here than not.
                print ("problems: i = %d, k = %d, j = %d" % (i, k, j))
                sys.exit(1)
            #
            # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            # ############################################
            
            p = i+k;
            q = i+k+1;
            mblplnk1.resetBranch()
            if DEBUG_searchForMBL:       
                print ("MBL scanning(i=%d,p=%d,q=%d,j=%d)" % (i, p, q, j))
            #
            
            no_V1 = False
            no_V2 = False
            if p - i > mnHairPin:
                mblplnk1 = self.lookupBranches(i, j, i, p, mblplnk1)
                if mblplnk1.nm == 'X':
                    no_V1 = True
                    mblplnk1.n = 0;
                #
                
            else:
                no_V1 = True;
            #
            
            mblplnk2.resetBranch()
            if j - q > mnHairPin:
                mblplnk2 = self.lookupBranches(i, j, q, j, mblplnk2)
                if mblplnk2.nm == 'X':
                    no_V2 = True
                    mblplnk2.n = 0;
                #
                
            else: 
                no_V2 = True
            #
            
            
            """@ 
            
            There are arguments both ways, but my opinion is that if
            we are "minimizing", a positive branch is basically a
            disallowed branch. Certainly if the branch is more
            expensive than a hairpin (and especially if it is more
            expensive than free strand), then it should not be
            selected in the first place, even if such a thing is
            found.
            
            """
            
            E1 = mblplnk1.V
            E2 = mblplnk2.V
            
            # for all __branched__ structures
            if (not no_V1) and (not no_V2):
                if (E1 > 0) or (E2 > 0):
                    
                    """180621: 
                    
                    In summing the free energies of the branches, it
                    makes no sense to add a very unstable branch to a
                    more stable one and call it an "answer". In fact,
                    because the FE is added, it might produce spurious
                    results where the only reason the "suboptimal"
                    solution was selected (the next step) was because
                    it counterbalanced a more stable branch. It may
                    even amount to rejecting an optimal branch in
                    favor of a suboptimal branch. If both branches
                    have a positive FE, then the less expensive one is
                    chosen.
                    
                    """
                    
                    if (E1 <= E2):
                        no_V2 = True
                        mblplnk2.Q = []
                        mblplnk2.n = 0;
                    else:
                        no_V1 = True
                        mblplnk1.Q = []
                        mblplnk1.n = 0;
                    #
                #
            #
            
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
            # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
            # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
            
            """180621; 
            
            for selecting suboptimal structures.
            
            In the vsfold5 version, selection of suboptimal structure
            was done by going back through this routine and gradually
            bouncing dG_lb (dG lower bound) up from the minimum free
            energy.
            
            It looks to me like this aspect of the method is not
            needed anymore. The way I build suboptimal structure in
            this newer program is I work backwards through the
            pushdown stack with the fragments and use fragment
            assembly to build the suboptimals based on the results of
            localscan_for_M().
            
            Nevertheless, I have left this code here anyway. There is
            a good chance that this code will eventually be removed.
            Additionally, further down in the method, perhaps the test
            "new_iMBL_dG < dG_lb" can also be removed.
            
            """
            
            # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
            # possible to turn this off. It is not used for anything
            if (not no_V1) and (not no_V2): # both are found
                if E1 + E2 < dG_lb + 0.5:
                    if (E1 <= E2) and (E1 > dG_lb):
                        no_V2 = True # turn off V2
                        mblplnk2.n = 0;
                        mblplnk2.Q = []
                        
                    elif (E2 < E1) and (E2 > dG_lb):
                        no_V1 = True  # turn off V1
                        mblplnk1.n = 0;
                        mblplnk1.Q = []
                        
                    else:
                        no_V1 = True  # turn off V1
                        no_V2 = True  # turn off V2
                        mblplnk1.n = 0;
                        mblplnk1.Q = []
                        mblplnk2.n = 0;
                        mblplnk2.Q = []
                        # both might have to be turned off
                    #
                    
                #
                
                if DEBUG_searchForMBL:
                    print ("[E1(%8.2f) + E2(%8.2f) > dG_lb(%8.2f)]?  => no_V1(%d), no_V2(%d)" \
                        % (E1, E2, dG_lb, no_V1, no_V2))
                #endif
                
            elif (not no_V1) or (not no_V2):
                if E1 <= E2:
                    if E1 < dG_lb + 0.5:
                        no_V1 = True
                        mblplnk1.n = 0;
                        mblplnk1.Q = []
                    #
                    
                    if DEBUG_searchForMBL:       
                        print ("[E1(%8.2f) > dG_lb(%8.2f)]? => no_V1(%d)" \
                            % (E1, dG_lb, no_V1))
                    #endif
                    
                else:
                    if E2 < dG_lb + 0.5:
                        no_V2 = True
                        mblplnk2.n = 0;
                        mblplnk2.Q = []
                    #
                    
                    if DEBUG_searchForMBL:
                        print ("[E2(%8.2f) > dG_lb(%8.2f)]? => no_V2(%d)" \
                            % (E2, dG_lb, no_V2))
                    #endif
                    
                #
                
            #
            
            # possible to turn this off. It is not used for anything
            # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            
            # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            if no_V1 and no_V2:
                if DEBUG_searchForMBL:
                    print ("no_V12: %d, %d; skip" % (no_V1, no_V2))
                #endif
                continue # nothing found over (i,j)
            # vvvvvvvvvvvvvvvvvv
            else: 
                if DEBUG_searchForMBL:
                    print ("no_V12: %d, %d; -> next stage" % (no_V1, no_V2))
                #endif
            # ^^^^^^^^^^^^^^^^^^
            
            """190103:
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            I still wonder if it would be far better to run
            trim_53pStems() here instead of further down. There are
            some costs for certain types of MBL arrangements which
            might be missed in the current scheme. On the other hand,
            some of this should be compensated in the current design
            because we will look for the best I-loop, the best M-loop
            and maybe also add PKs and Islands (where it applies).
            
            The original vsfold tried to move the stems around when
            there were conflicts ... perhaps this is another way to
            approach this problem.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            """
            
            new_MBL.resetBranch()
            if  (not no_V1) and (not no_V2):
                if DEBUG_searchForMBL:
                    print ("both: no_V12 = %d, %d" % (no_V1, no_V2))
                #endif
                
                # branches exist for both sectors
                new_MBL.resetBranch()
                # construct the branching structure from the obtained data
                new_MBL.n = mblplnk1.n + mblplnk2.n # formerly branchV =
                # new_MBL =>  branch(0, 1, 2, ..., n1-1)
                #           + branch(n1+0, n1+1, ... n2-1)
                
                # => bramnch(0 ... n2-1)
                new_MBL.V = mblplnk1.V + mblplnk2.V
                
                shiftV1 = mblplnk1.n
                for l in range(0, shiftV1): # sum the first part
                    new_MBL.pushBranch(mblplnk1.getBranch(l))
                #
                
                shiftV2 = mblplnk2.n
                for l in range(0, shiftV2): # add
                    new_MBL.pushBranch(mblplnk2.getBranch(l))
                #
                
                new_MBL.n = len(new_MBL.Q)
                        
                if DEBUG_searchForMBL:
                    print ("sectors 1 and 2: branches: %d, E = %8.2f" \
                        % (new_MBL.n, new_MBL.V))
                #endif
                
            elif (not no_V1) and no_V2:
                # sector 1 has a viable branch
                if DEBUG_searchForMBL:
                    print ("sector 1: no_V12 = %d, %d" % (no_V1, no_V2))
                #endif
                
                # construct the branching structure from the obtained data
                new_MBL.n = mblplnk1.n  # formerly branchV = new_mbl
                new_MBL.V = mblplnk1.V
                
                for l in range(0, mblplnk1.n):
                    new_MBL.pushBranch(mblplnk1.getBranch(l))
                #
                
                new_MBL.n = len(new_MBL.Q)
                if DEBUG_searchForMBL:       
                    print ("sector 1: branches: %d, E = %8.2f" % (new_MBL.n, new_MBL.V))
                #endif
                
            else: # no_V1 and (not no_V2):
                # sector 2 has a viable branch
                if DEBUG_searchForMBL:
                    print ("sector 2: no_V12 = %d, %d" % (no_V1, no_V2))
                #endif
                
                # construct the branching structure from the obtained data
                new_MBL.n = mblplnk2.n  # formerly branchV = new_mbl
                new_MBL.V = mblplnk2.V
                for l in range(0, mblplnk2.n):
                    # for (int l = 1; l <= mblplnk2.Q[0].i; l++) {
                    new_MBL.pushBranch(mblplnk2.getBranch(l))
                #
                
                new_MBL.n = len(new_MBL.Q)
                if DEBUG_searchForMBL:       
                    print ("sector 2: branches: %d, E = %8.2f" % (new_MBL.n, new_MBL.V))
                #endif
                
            #
            
            
            if new_MBL.n > 0:
                # Admittedly, it would be nice not to have to do this
                # rigamoral all the time, but unfortunately, we have
                # to. Every time we go through here, we have something
                # unique about the input structure.
                
                # currently doesn't exist 
                new_cP = self.fe.trim_53pStems(i, j,      # 5' to 3' [i/p]MBL bounds 
                                               new_MBL,   # mbl branch info
                                               self.fe.smap, # the general mapping
                                               opt_iMBL)  # iMBL or pMBL? B.Cs 
                
            else:
                print ("ERROR searchForMBL(%d,%d): really shouldn't be here." % (i, j))
                print ("      This operation is after the search claimed at least")
                print ("      one branch to be present. So why would the number of")
                print ("      branches now be %d? Something is seriously wrong." % (new_MBL.n))
                sys.exit(1);
            #
            
            new_MBL.V = new_cP
            # for chromatin, trim_53pStem() doesn't actually change
            # anything. However, for RNA, it could change slightly.
            new_iMBL_dG = new_cP + Mloop_wt;
            # if Mloop_wt is non-zero it is still a single value because
            # it reflects closing at (i,j).
            
            if new_MBL.n > 1:
                new_MBL.nm = 'P'
                new_MBL.dr = '-'
            else:
                new_MBL.nm = 'J'
                new_MBL.dr = '-'
            #
            
            
            
            # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
            lmblh = self.check_MBLHandle(new_MBL, i, j, lmblh, on_ij, DEBUG_searchForMBL, "v1+v2")
            # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            
            if DEBUG_searchForMBL:
                E1 = mblplnk1.V
                E2 = mblplnk2.V
                pv = -1; qv = -1
                
                print ("MLclose(%8.2f) + new_cP(%8.2f) = %8.2f"  \
                    % (Mloop_wt, new_cP, Mloop_wt+new_cP))
                print ("opt_iMBL = %d" % (opt_iMBL))
                if opt_iMBL:
                    print ("NOTE: iMBL weight ADDED: Mloop_wt(%d,%d) = %8.2f" \
                        % (i, j, Mloop_wt))
                    print ("      new_cP(%8.2f) + Mloop_wt(%8.2f) = %8.2f, best_iMBL_dG = %8.2f" \
                        % (new_cP, Mloop_wt, new_cP + Mloop_wt, best_iMBL_dG))
                #
                print ("E1(%3d,%3d) = %8.2f; E2(%3d,%3d) = %8.2f; new_iMBL_dG = %8.2f; best_iMBL_dG(%d) = %8.2f" \
                    % (i, i+k, E1, i+k+1, j, E2, new_iMBL_dG, 0, lmblh[0].best_iMBL_dG))
                print ("                                                                        best_iMBL_dG(%d) = %8.2f" \
                    % (1, lmblh[1].best_iMBL_dG))
                print ("                                                                        best_iMBL_dG(%d) = %8.2f" \
                    % (2, lmblh[2].best_iMBL_dG))
                
                print (self.showBranches(new_MBL, "new_MBL"))
            #endif
            
            
        # end for (int k = 1; k <= j - i - 2; k++)
        
        if DEBUG_searchForMBL:  
            print ("finished part1 .... now part2")
        #
        
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # for k in range(1, j - i - 2): # (int k = 1; k <= j - i - 2; k++) {
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        
        # the final set of cases to test to find who is boss
        
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        branch  = []
        branch += [Branch(i,   j-1)]
        branch += [Branch(i+1, j-1)]
        branch += [Branch(i+1, j)]
        #branch += [Branch(i,   j)]  # just a hairpin
        
        # we don't need to look at (i,j) because that question is
        # addressed in the next part of the operation that follows
        # searchForMBL(i,j).
        
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for k in range(0, len(branch)):
            mblplnk1.resetBranch()
            new_MBL.resetBranch()
            if DEBUG_searchForMBL:
                print ("MBL scanning(pq=(%d, %d))" \
                    % (branch[k].i, branch[k].j))
            # endif
            
            no_Vf = False
            p = branch[k].i; q = branch[k].j
            # !!!VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
            if p == i and q == j:
                on_ij = True
            #
            
            # 190130: this should be removed, eventually
            # !!!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            if q - p > mnHairPin: 
                mblplnk1 = self.lookupBranches(i, j,
                                               p, q,
                                               mblplnk1)
                if mblplnk1.nm == 'X':
                    no_Vf = True
                    mblplnk1.n = 0;
                #
                
                
            else:
                no_Vf = True
            #
            
            if DEBUG_searchForMBL:       
                print ("no_Vf(%d,%d) = %d, mblplnk1.nm = %s" \
                    % (p, q, no_Vf, mblplnk1.nm))
            #endif
            
            if not no_Vf:
                
                new_MBL.n  = mblplnk1.n  # number of branches
                new_MBL.V  = mblplnk1.V  # total free energy
                btp = '-'
                if len(self.fe.smap.glink[i][j].lg) > 0:
                    btp  = self.fe.smap.glink[i][j].lg[0].motif[0].get_btp()
                #
                
                
                if new_MBL.n > 1:
                    new_MBL.nm = 'P'
                    new_MBL.dr = '-'
                else:
                    new_MBL.nm = 'J'
                    new_MBL.dr = '-'
                #
                
                # !!!VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
                # !!!VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
                
                # I think that this part (between the VVVVV and AAAAA)
                # can be removed now. 190130
                if on_ij and (new_MBL.nm == 'K' or new_MBL.nm == 'R'):
                    # 190128: The order of this is important. We have
                    # to assign new_MBL.nm before we can test it
                    # here. Because we are skimming the topmost
                    # optimal configurtion, This step has to come
                    # after the above operations.
                    print (self.fe.smap.glink[i][j].lg[0].motif[0].show_Motif())
                    print ("pk on_ij")
                    # sys.exit(0)
                    break
                #
                
                # I think with the current arrangement of searchForMBL
                # in the analysis sequence, this test is not
                # necessary. The main feature to check for that is not
                # immediately known and accounted for should be only K
                # (and R).
                if on_ij and (new_MBL.nm == 'B' or
                              new_MBL.nm == 'S' or
                              new_MBL.nm == 'W'):
                    # shouldn't happen!!!???!!!
                    print ("found a case of BSW, nm = ", new_MBL.nm, (i, j))
                    sys.exit(1)
                    break
                #
                # !!!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                # !!!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                
                for l in range(0, mblplnk1.n):
                    new_MBL.pushBranch(mblplnk1.getBranch(l))
                #
                if DEBUG_searchForMBL:
                    ct = mblplnk1.nm
                    bt = mblplnk1.dr
                    print ("p2: new_MBL(%d,%d)[%s][%5s] = %8.2f" % (p, q, ct, bt, mblplnk1.V))
                    print (self.showBranches(new_MBL, "new_MBL[pt2]"))
                #endif
                
                
                # 180421wkd: It seems that this routine was designed more
                # flexible than I expepcted. The method trim_53pStems is
                # able to handle new_MBL.n == 1 without any modification
                # at all!
                new_cP = self.fe.trim_53pStems(i, j,      # 5' to 3' [i/p]MBL bounds 
                                               new_MBL,   # mbl branch info
                                               self.fe.smap, # the general mapping
                                               opt_iMBL)  # iMBL or pMBL? B.Cs
                new_MBL.V = new_cP
                new_iMBL_dG   = new_cP + Mloop_wt
                
                
                # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
                lmblh = self.check_MBLHandle(new_MBL, p, q, lmblh, on_ij, DEBUG_searchForMBL, "v1")
                # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            # ##  finish if not no_Vf: { ... }
            
        # ## finish for (int k = 0; k < n_branch; k++) 
        
        del branch
        
        
        has_branches = False
        # print ("1")
        for gid in range(0, len(lmblh)):
            
            if len(lmblh[gid].mbls) > 0:
                # print ("assign gid(%d){ij(%d,%d)}" % (gid, i, j))
                for kv in range(0, len(lmblh[gid].mbls)):
                    has_branches = True
                    if DEBUG_searchForMBL:
                        p, q = lmblh[gid].mbls[kv].getBranchEnds()
                        print ("found mblh[%d]{ij(%d,%d),pq(%d,%d)}:" % (kv, i, j, p, q))
                    #endif
                #
                if DEBUG_searchForMBL:
                    show_MBLHandle(lmblh[gid], self.fe.smap)
                #endif
            else:
                if DEBUG_searchForMBL:
                    print ("empty gid(%d){ij(%d,%d)}" % (gid, i,j))
                #endif
            #
                
        #
        
        
        # summarize the results of the search 
        
        # 180413(indeed, Friday the 13th): I am changing this function to
        # find both I-loops and M-loops. I think this strategy will make
        # this subprogram an essential tool for doing the speedup. The
        # other part is simplifying the I-loop technology
        if not has_branches:
            # print ("2")
            if not self.fe.btype[i][j].pair > 0:
                k = self.fe.dSearchK['SB']
                lmblh[k].mbls += [MBLptr(i,j)]
                lmblh[k].mbls[0].V = self.fe.cle_HloopE(i, j, lmblh[0].mbls[0])
                lmblh[k].mbls[0].addBranch(i, j)
                lmblh[k].mbls[0].nm = 'X'   # fake motif
                lmblh[k].mbls[0].dr = '-'
                
                if DEBUG_searchForMBL:
                    print ("final assignment ij(%d,%d) is %s" % (i, j, lmblh[k].mbls[0].nm))
                    #sys.exit(0)
                #endif
                
            else:
                if DEBUG_searchForMBL:
                    print ("location has a pair ij(%d,%d)" % (i, j))
                #endif
            #
        #
        
        if DEBUG_searchForMBL:
            for gid in range(0, len(lmblh)):
                for kv in range(0, len(lmblh[gid].mbls)):
                    print ("lmblh[%d].mbls[%d]" % (gid, kv))
                    print (lmblh[gid].mbls[kv])
                #
            #
        #endif
        
        # chr10_3894786_4781825_res5kb.heat
        #self.fe.stopWhenMatchFound(i, j,   0,   5, "end of searchForMBL()")
        # result: S, sa, -1.32 (0,5) .. (2,3)
        #self.fe.stopWhenMatchFound(i, j,   0,   6, "end of searchForMBL()")
        #self.fe.stopWhenMatchFound(i, j,   0,   7, "end of searchForMBL()")
        #self.fe.stopWhenMatchFound(i, j,   0,   8, "end of searchForMBL()")
        #self.fe.stopWhenMatchFound(i, j,   4,   9, "end of searchForMBL()")
        #self.fe.stopWhenMatchFound(i, j,   0,   9, "end of searchForMBL()")
        # self.fe.stopWhenMatchFound_sfm(i, j, 10, 99, lmblh, "end of searchForMBL()")
        
        # print (self.fe.btype[i][j])
        # testB-1S.heat
        #self.fe.stopWhenMatchFound(i, j, 4, 13, "end of searchForMBL()")
        #self.fe.stopWhenMatchFound(i, j, 3, 14, "end of searchForMBL()")
        #self.fe.stopWhenMatchFound(i, j, 2, 15, "end of searchForMBL()")
        
        # testI-2S.heat
        #self.fe.stopWhenMatchFound(i, j, 5, 24, "end of searchForMBL()")
        #self.fe.stopWhenMatchFound(i, j, 4, 25, "end of searchForMBL()")
        #self.fe.stopWhenMatchFound(i, j, 3, 26, "end of searchForMBL()")
        
        # test1.heat
        #self.fe.stopWhenMatchFound(i, j, 3, 12, "end of searchForMBL()")
        #self.fe.stopWhenMatchFound(i, j, 1, 14, "end of searchForMBL()")
        #self.fe.stopWhenMatchFound(i, j, 1, 19, "end of searchForMBL()")
        #self.fe.stopWhenMatchFound(i, j, 0, 24, "end of searchForMBL()")
        
        # test4.heat
        #self.fe.stopWhenMatchFound(i, j, 5, 25, "end of searchForMBL()")
        #self.fe.stopWhenMatchFound(i, j, 3, 28, "end of searchForMBL()")
        #self.fe.stopWhenMatchFound(i, j, 2, 29, "end of searchForMBL()")
        #self.fe.stopWhenMatchFound(i, j, 1, 30, "end of searchForMBL()")
        #self.fe.stopWhenMatchFound(i, j, 0, 32, "end of searchForMBL()")
        del new_MBL
        del mblplnk1
        del mblplnk2
        
        return lmblh  
    #
    
    
    
    
    
    def show_lmblh(self, i, j, lmblh, callpt = "not specified"):
        for gid in range(0, len(lmblh)):
            
            if len(lmblh[gid].mbls) > 0:
                # print ("assign gid(%d){ij(%d,%d)}" % (gid, i, j))
                for kv in range(0, len(lmblh[gid].mbls)):
                    p, q = lmblh[gid].mbls[kv].getBranchEnds()
                    nm = lmblh[gid].mbls[kv].nm
                    dr = lmblh[gid].mbls[kv].dr
                    V  = lmblh[gid].mbls[kv].V
                    print ("found mblh[%d]{ij(%d,%d),pq(%d,%d)[%s][%5s][%8.2f]}:" \
                        % (kv, i, j, p, q, nm, dr, V))
                    #endif
                #
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                # 190108: for the moment, we will leave this here
                # because mbls[0] is the only structure
                # produced. However, this variable is not needed and
                # it is becoming increasingly problematical.
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                show_MBLHandle(lmblh[gid], self.fe.smap)
            else:
                if DEBUG_searchForMBL:
                    print ("empty gid(%d){ij(%d,%d)}" % (gid, i,j))
                #endif
            #
                
        #
    #
    
    
    
    def check_MBLHandle(self, new_MBL, p, q, lmblh, on_ij, debug = True, lbl = "x"):
        i = new_MBL.i; j = new_MBL.j
        if i == p and j == q and on_ij:
            # on_ij True/False is the search region (i,j), (i+1,j), (i+1,j-1), (i,j-1)
            if self.fe.btype[i][j].pair > 0:
                lmblh[self.fe.dSearchK['SB']].update_mbl(new_MBL, p, q, debug, lbl)
            #
            
        #
        
        if new_MBL.n == 1 and not on_ij:
            lmblh[self.fe.dSearchK['IJ']].update_mbl(new_MBL, p, q, debug, lbl)
        #
        
        if new_MBL.n > 1:
            lmblh[self.fe.dSearchK['MP']].update_mbl(new_MBL, p, q, debug, lbl)
        #
        
        return lmblh
    #
    
    
    
    
    # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    # ###################################################################
    # ###############  New I-loop and MBL search method   ###############
    # ###################################################################
    
    
    
    
    
    
    
    # Main calculation script for finding the minimum free energy and
    # (in the future) suboptimal structures.
    def minFE(self, T):
        global DEBUG_find_best_PK
        global DEBUG_searchForMBL
        global DEBUG_lookupBranches
        
        debug_all(CHECK_ALL)        
        opt_iMBL = False
        
        self.T = T
        
        self.fe.T = T
        for j in range(1, self.N):
            for i in range(j-1, -1, -1):
                #DEBUG_find_best_PK   = kickOn(i, 9, j, 14)
                #DEBUG_find_ifStem    = kickOn(i, 9, j, 14)
                
                #DEBUG_searchForMBL   = kickOn(i, 45, j, 65)
                #DEBUG_lookupBranches = kickOn(i, 45, j, 65)
                
                # build search weighting handle
                lmblh = []
                ndx = 0
                for pp in self.fe.dSearchK.keys():
                    if DEBUG_minFE:
                        print ("building assigning %d, %s" % (self.fe.dSearchK[pp], pp))
                    #
                    mblh = MBLHandle(ndx, i, j, -1.0e6) # tester
                    # 90103: Unfortunately, we have to pass this
                    # shit (fe and smap) for various reasons that
                    # are unavoidable. The object smap probably
                    # can be defind as None if no debugging is
                    # required.
                    
                    lmblh += [mblh]
                    ndx += 1
                # 
                
                lmblh = self.searchForMBL(i, j,     # search region
                                          lmblh,    # MBLHandle()
                                          opt_iMBL) # opt_iMBL (True -> iMBL)
                # to get branchlist as (i1,j1), (i2,j2) etc. ==>  mbl.getBranchlist()
                
                # 190101: I see now that I probably can store the
                # degenerate solutions even right here! Therefore,
                # I don't really need glink. What is important is
                # having a consistent way to construct the data.
                link_X = None
                for gid in range(0, len(lmblh)):
                    if len(lmblh[gid].mbls) == 0:
                        continue
                    #
                    
                    flag_link_initialized = False
                    for kv in range(0, len(lmblh[gid].mbls)):
                        if not flag_link_initialized:
                            # first on is always the case
                            ctp_X  = lmblh[gid].mbls[kv].nm
                            btp_X  = lmblh[gid].mbls[kv].dr
                            cp_X   = lmblh[gid].mbls[kv].getBranchlist()
                            V_X    = lmblh[gid].mbls[kv].V
                            if (ctp_X == 'I' or ctp_X == 'M') and btp_X == '-':
                                print ("1, found btp_X = %s for I-loop" % btp_X)
                                sys.exit(0)
                            #
                            
                            link_X = Link(i, j,
                                          V_X,
                                          ctp_X,
                                          btp_X,
                                          cp_X)
                            flag_link_initialized = True
                            
                        else:
                            # if there are structures of a degenerate FE
                            V_X   = lmblh[gid].mbls[kv].V
                            cp_X  = lmblh[gid].mbls[kv].getBranchlist()
                            ctp_X = lmblh[gid].mbls[kv].nm
                            btp_X = lmblh[gid].mbls[kv].dr
                            if (ctp_X == 'I' or ctp_X == 'M') and btp_X == '-':
                                print ("2, found btp_X = %s for I-loop" % btp_X)
                                sys.exit(0)
                            #
                            
                            link_X.add_Motif(i, j, V_X, ctp_X, btp_X, cp_X)
                            self.fe.save_best_iloops(link_X)
                            
                        #
                        
                        if DEBUG_minFE:
                            print ("assigned %s(%s): (%2d,%2d)[%8.2f]" \
                                % (ctp_X, btp_X, i, j, V_X), cp_X)
                        #endif
                        
                    #|endfor kv
                    
                    self.fe.smap.glink[i][j].add_link(link_X)
                    
                #|endfor gid
                
                self.fe.smap.mergeSortLinks(self.fe.smap.glink[i][j].lg)
                
                
                if self.fe.btype[i][j].pair > 0: # if hvij > 0.0:
                    link_Q = None
                    flag_link_initialized = False
                    bondtype = self.fe.btype[i][j].btp
                    # was: self.get_bondtype(self.fe.btype[i][j].wt)
                    # self.fe.hv[i][j] -> self.fe.btype[i][j].wt
                    # 190524 was dGij = self.fe.calc_dG(i, j, self.fe.btype[i][j].wt, T)
                    dGij = self.fe.btype[i][j].dGp
                    # self.fe.hv[i][j] -> self.fe.btype[i][j].wt
                    
                    if DEBUG_minFE:
                        # print (self.fe.btype[i][j].wt)
                        # print (self.fe.hv[i][j])
                        print ("assigned B: (%2d,%2d)[%8.3f]" % (i, j, dGij), [(i,j)])
                        # sys.exit(0)
                    #
                    
                    # ######  since 'B' _exists_, must be recorded  ######
                    
                    if bondtype == 's' or bondtype == 'sa' or bondtype == 'sp': 
                        if dGij < self.fe.dGMI_threshold:
                            has_cap, ctpB, btpB, dGB, jnB \
                                = self.fe.lookupCap(i, j, dGij, DEBUG_minFE)
                            if has_cap:
                                if DEBUG_minFE:
                                    print ("Q: btpB = ", btpB, ", B:", has_cap)
                                #
                                
                                link_X = Link(i, j, dGB, ctpB, btpB, jnB)
                                self.fe.smap.glink[i][j].add_link(link_X)
                                if len(jnB) == 1:
                                    self.fe.save_best_iloops(link_X)
                                #
                                
                            #
                            
                        #
                        
                        link_Q = Link(i, j, dGij, 'B', bondtype, [(i,j)])
                        self.fe.smap.glink[i][j].add_link(link_Q)
                        self.fe.save_best_iloops(link_Q)
                        if DEBUG_minFE:
                            self.show_smap_xy(i,j)
                        #
                        
                        flag_link_initialized = True
                        
                    elif (bondtype == 'c' or
                          bondtype == 't' or
                          bondtype == 'r' or
                          bondtype == 'l'):
                        link_Q = Link(i, j, dGij, 'W', bondtype, [(i,j)])
                        self.fe.smap.glink[i][j].add_link(link_Q)
                        self.fe.save_best_iloops(link_Q)
                        flag_link_initialized = True
                        
                    else:
                        link_Q = Link(i, j, dGij, 'X', '-', [(i,j)])
                        self.fe.smap.glink[i][j].add_link(link_Q)
                        self.fe.save_best_iloops(link_Q)
                        print ("called X in H loop ij=(%d,%d)" % (i,j))
                        sys.exit(0)
                        
                    #
                    
                    flag_link_initialized = True
                    
                    stem_aa, stem_pp = self.fe.find_ifStem(i, j, dGij, bondtype,
                                                           DEBUG_find_ifStem)
                    # print (stem_aa, len(stem_aa), stem_pp, len(stem_pp))
                    if len(stem_aa) > 0:
                        self.fe.make_StemMotif(stem_aa, DEBUG_make_StemMotif)
                        # !!!XXX
                    #
                    
                    if len(stem_pp) > 0:
                        self.fe.make_StemMotif(stem_pp, DEBUG_make_StemMotif)
                    #
                    
                    # /home/yuriv/cent/loops_CCDs/paper_examples_try161012/190514/test1.heat
                    #self.fe.stopWhenMatchFound(i, j, 3, 6, "after find_ifStem")
                    
                    # testB-1Sp.heat
                    # self.fe.stopWhenMatchFound(i, j, 2, 13, "after find_ifStem")
                    
                    # chr1_10751707_10970928_res5kb.heat
                    #self.fe.stopWhenMatchFound(i, j, 25, 41, "after find_ifStem")
                    #self.fe.stopWhenMatchFound(i, j, 36, 41, "after find_ifStem")
                    #self.fe.stopWhenMatchFound(i, j, 35, 42, "after find_ifStem")
                    #self.fe.stopWhenMatchFound(i, j, 0, 35, "after find_ifStem")
                    
                    
                    if (i,j) in self.fe.all_ctcf:
                        dGbest = self.fe.smap.glink[i][j].lg[0].Vij
                        # if it doesn't have the key, then forget it!
                        island = self.fe.find_ctcf_islands(i, j, dGbest,
                                                           DEBUG_find_ctcf_islands)
                        
                        # print (dGbest, island)
                        if len(island) > 0:
                            
                            # if there's nothing there, then don't bother with it!
                            for ff in island:
                                wyspa = ff[0]; join = ff[1]; dGW = ff[2] 
                                link_W = Link(i, j, dGW, 'W', 'wyspa', join, [], wyspa)
                                self.fe.smap.glink[i][j].add_link(link_W)
                            #|endfor
                            
                            # self.traceback_mFE(i,j, 0, True)
                            
                        #
                        
                        #self.fe.stopWhenMatchFound(i, j, 0, 32, "in find_ctcf_islands")
                        #print ("island: planned exit")
                        #sys.exit(0)
                    #
                    
                    #self.fe.stopWhenMatchFound(i, j, 1, 14, "after find_ctcf_islands")
                #
                
                # To do the pseudoknot part, we have to sort the current data
                self.fe.smap.mergeSortLinks(self.fe.smap.glink[i][j].lg)
                
                # pseudoknot and CTCF-island solutions
                ctp = self.fe.smap.glink[i][j].lg[0].motif[0].get_ctp()
                
                check_PK = False
                if ctp in dconnection:
                    """@
                    
                    Have to decide if we will check for PKs or not.
                    Originally, the test was
                    
                    if self.fe.btype[i][j].pair > 0:
                    
                    and before that, it was 
                    
                    if hvij > 0.0:
                    
                    so it always has had some problems.
                    """
                    check_PK = True 
                    if ctp == 'M' or ctp == 'I':
                        dGijH = self.fe.btype[i][j].dGp
                        
                        """190524 was
                        
                        dGijH = self.fe.calc_dG(i, j,
                                                self.fe.btype[i][j].wt,
                                                T) 
                        """
                        
                        if dGijH >= 0.0:
                            # the basic structure at ij should at
                            # least have a favorable FE at the
                            # closing position to make a PK a
                            # feasible possibility
                            if DEBUG_find_best_PK:
                                print (dGijH)
                                print ("check for PK: skipping (%d,%d)[%s]" % (i, j, ctp))
                                #self.fe.stopWhenMatchFound(i, j, i, j, "test_PK")
                                # sys.exit(0)
                            #
                            
                            check_PK = False
                            
                        #
                        
                    #
                    
                    if check_PK:
                        #print ("call searchForPKs(%d,%d):" % (i, j))
                        self.fe.searchForPKs(i, j, DEBUG_find_best_PK)
                        check_PK = False
                    #
                    
                #
                
                # for longer sequences, we also simply check each time
                # we go through this cyle at i = 0.
                if i == 0 and j-i > 15:
                    if DEBUG_find_best_PK:
                        print ("searching the full span for PK")
                    #
                    
                    self.fe.searchForPKs(i, j, DEBUG_find_best_PK)
                #
                
                # /home/yuriv/cent/loops_CCDs/paper_examples_try161012/190514/test1.heat
                # self.fe.stopWhenMatchFound(i, j, 3, 6, "after searchForPKs")
                # self.fe.stopWhenMatchFound(i, j, 3, 7, "after searchForPKs")
                # /home/yuriv/cent/loops_CCDs/paper_examples_try161012/190514/check3.heat
                # self.fe.stopWhenMatchFound(i, j, 3, 8, "after searchForPKs")
                
                
                
                self.fe.smap.mergeSortLinks(self.fe.smap.glink[i][j].lg)
                self.fe.dG[i][j] = self.fe.smap.glink[i][j].lg[0].motif[0].Vij # !!!!!
                
                # corePKwithbranch.heat
                #self.fe.stopWhenMatchFound(i, j,   0,  46, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j,  11,  52, "finishing minFE loop at pos ij")
                
                # chr10_3894786_4781825_res5kb.heat
                #self.fe.stopWhenMatchFound(i, j,   0,   5, "finishing minFE loop at pos ij")
                # result: S, sa, -1.32 (0,5) .. (2,3)
                #self.fe.stopWhenMatchFound(i, j,   0,   6, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j,   0,   7, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j,   0,   8, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j,   4,   9, "finishing minFE loop at pos ij")
                
                #self.fe.stopWhenMatchFound(i, j,   2,   7, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j,   2,   8, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j,   2,   9, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j,   0,   9, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j,   6,  16, "finishing minFE loop at pos ij")
                
                #self.fe.stopWhenMatchFound(i, j,  30,  36, "finishing minFE loop at pos ij")
                #self.show_smap_xy(38,41)
                #self.fe.stopWhenMatchFound(i, j,  37,  42, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j,  38,  41, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j,  37,  43, "finishing minFE loop at pos ij")
                # result: S, sa, -2.11 (6,16) ... (10,12)
                
                #self.fe.stopWhenMatchFound(i, j,  52,  59, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j,  51,  59, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j,  50,  59, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j,  46,  60, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j,  51,  61, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j,  48,  62, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j,  47,  62, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j,  52,  63, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j,  45,  64, "finishing minFE loop at pos ij")
                
                #self.fe.stopWhenMatchFound(i, j,  45,  65, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j,  10,  99, "finishing minFE loop at pos ij")
                
                #self.fe.stopWhenMatchFound(i, j, 125, 132, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j, 125, 134, "finishing minFE loop at pos ij")
                
                
                # chr1_10751707_10970928_res5kb.heat
                #self.fe.stopWhenMatchFound(i, j, 0, 8, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j,  0, 42, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j, 27, 42, "finishing minFE loop at pos ij")
                
                # chr10_106061350_106114340_res5kb.heat
                #self.fe.stopWhenMatchFound(i, j, 2,  8, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j, 1,  9, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j, 0, 10, "finishing minFE loop at pos ij")
                
                
                # testP-3B.heat
                #self.fe.stopWhenMatchFound(i, j, 2, 35, "finishing minFE loop at pos ij")
                
                
                # testB-1S.heat
                #self.fe.stopWhenMatchFound(i, j, 4, 13, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j, 3, 14, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j, 2, 15, "finishing minFE loop at pos ij")
                
                # testI-2S.heat
                #self.fe.stopWhenMatchFound(i, j, 5, 24, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j, 4, 25, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j, 3, 26, "finishing minFE loop at pos ij")
                
                # test1.heat
                #self.fe.stopWhenMatchFound(i, j, 3, 11, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j, 1, 14, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j, 1, 19, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j, 0, 24, "finishing minFE loop at pos ij")
                # test4.heat
                #self.fe.stopWhenMatchFound(i, j, 5, 25, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j, 3, 28, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j, 2, 29, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j, 1, 30, "finishing minFE loop at pos ij")
                #self.fe.stopWhenMatchFound(i, j, 0, 32, "finishing minFE loop at pos ij")

                # chr10_1462921_1758490_res5kb.heat
                #self.fe.stopWhenMatchFound(i, j, 9, 15, "finishing minFE loop at pos ij")
                
                
            #|endfor i in range(j-1, -1, -1):
            
        #|endfor j in range(1, self.N):
        
        if DEBUG_minFE:
            print ("finished minFE")
            # sys.exit(0)
        #endif
        
        return self.fe.dG, self.fe.smap
    #
    
#


def test0():
    
    # Evidently, Vienna() still has a few problems presently because
    # it cannot convert ".ABCD.......abcd...." properly.
    
    
    # a set structure
    #          0         10        20        30        40        50        60        70        80        90
    #          |         |         |         |         |         |         |         |         |         |
    ss_seq  = "(.(((((......))))).((((.....)))).)"
    #            (2, 17)          (19,31)
    
    be = BranchEntropy()
    
    vs = Vstruct()
    vs.convert_CTCFstruct(ss_seq, True)
    be.add_hv(vs)
    be.add_btype(vs)
    
    # print ("planned exit after running convert_CTCFstruct"); sys.exit(0);
    v2t = Vienna2TreeNode(vs)
    v2t.vienna2tree()
    print ("main:")
    tf = LThreadBuilder(v2t)
    tf.visit(v2t.genTree)
    print ("LThread notation: ")
    tf.disp_lt()
    
    
    # print ("planned exit after running vienna2thread"); sys.exit(0);
    l2m = LThread2Motif(v2t.lt.sqlen)
    l2m.thread2motif(v2t.lt)
    print (len(l2m.smap.glink))
    
    
    mbl_M = MBLptr(0,33)
    mbl_M.addBranch(2,17)
    mbl_M.addBranch(19,31)
    mbl_M.n = len(mbl_M.Q)
    mbl_M.V = be.cle_MloopEV(0,33, mbl_M, l2m.smap, True)
    print (mbl_M.V)
    
    #i = 2; j = 17
    i = 0; j = 33
    mbl_H = MBLptr(i,j)
    mbl_H.addBranch(i,j)
    mbl_H.n = len(mbl_H.Q)
    mbl_H.V = be.cle_HloopE(i, j, mbl_H)
    print (mbl_H.V)
    
    i = 0; j = 33
    p = 2; q = 17
    mbl_I = MBLptr(i,j)
    mbl_I.addBranch(p,q)
    mbl_I.n = len(mbl_H.Q)
    mbl_I.V = be.cle_IloopE(i, j, p, q, mbl_I, l2m.smap, True)
    print (mbl_I.V)
#   





def test1():
    
    # Evidently, Vienna() still has a few problems presently because
    # it cannot convert ".ABCD.......abcd...." properly.
    
    
    # a set structure
    #          0         10        20        30        40        50        60        70        80        90
    #          |         |         |         |         |         |         |         |         |         |
    ss_seq  = ["(.(((((......))))).((((.....)))).)"]
    #            (2, 17)          (19,31)
    
    clclt = Calculate()
    
    vs = Vstruct()
    vs.parse_fullDotBracketStructure(ss_seq[0], True)
    clclt.fe.add_hv(vs)
    clclt.fe.add_btype(vs)
    
    # print ("planned exit after running convert_CTCFstruct"); sys.exit(0);
    v2t = Vienna2TreeNode(vs)
    v2t.vienna2tree()
    print ("main:")
    tf = LThreadBuilder(v2t)
    tf.visit(v2t.genTree)
    print ("LThread notation: ")
    tf.disp_lt()

    t2m = TreeNode2Motif(v2t)
    t2m.visit(v2t.genTree)
    print (v2t.MPlist)
    t2m.post_graftMP()
    
    print (len(t2m.smap.glink))
    
    
    mbl_M = MBLptr(0,33)
    mbl_M.addBranch(2,17)
    mbl_M.addBranch(19,31)
    mbl_M.n = len(mbl_M.Q)
    mbl_M.V = clclt.fe.cle_MloopEV(0,33, mbl_M, t2m.smap, True)
    print (mbl_M.V)
    
    #i = 2; j = 17
    i = 0; j = 33
    mbl_H = MBLptr(i,j)
    mbl_H.addBranch(i,j)
    mbl_H.n = len(mbl_H.Q)
    mbl_H.V = clclt.fe.cle_HloopE(i, j, mbl_H)
    print (mbl_H.V)
    
    i = 0; j = 33
    p = 2; q = 17
    mbl_I = MBLptr(i,j)
    mbl_I.addBranch(p,q)
    mbl_I.n = len(mbl_H.Q)
    mbl_I.V = clclt.fe.cle_IloopE(i, j, p, q, mbl_I, t2m.smap, True)
    print (mbl_I.V)
    
    i = 2; j = 31
    
    clclt.smap = t2m.smap

    ndx = 0
    lmblh = [] # build objects of class MBLHandle() for processing
    for pp in clclt.fe.dSearchK.keys():
        print ("building assigning %d, %s" % (clclt.fe.dSearchK[pp], pp))
        
        # build search weighting handle
        mblh = MBLHandle(ndx, i, j, -1.0e6) # main tester
        lmblh += [mblh]
        ndx += 1
        # 90103: Unfortunately, we have to pass all this shit (fe and
        # smap) for various reasons that are unavoidable. The object
        # smap probably can be defind as None if no debugging is
        # required.
    #|endfor
    
    lmblh = clclt.searchForMBL(i, j,     # search region
                              lmblh,    # list of type MBLHandle()
                              False)    # iMBL=1/pMBL=0
    for pp in clclt.fe.dSearchK.keys():
        print ("results %d, %s" % (clclt.fe.dSearchK[pp], pp))
        if len(lmblh[clclt.fe.dSearchK[pp]].mbls) > 0:
            show_MBLHandle(lmblh[clclt.fe.dSearchK[pp]], clclt.smap)
        else:
            print ("nothing found")
        #
        
    #|endfor
    
#   



def main(cl):
    # presently doesn't do anything
    print (cl)
#


if __name__ == '__main__':
    if TEST0:
        test0()
    elif TEST1:
        test1()
    else:
        main(sys.argv)
    #
#
