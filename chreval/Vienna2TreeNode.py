#!/usr/bin/env python

"""@@@

Main Module:   Vienna2TreeNode.py 

Objects:       LThread2DotBracket 
               Vienna2TreeNode

Functions:     getStem 
               makeStemList 
               dispStemList 
               ins_sort_StemList 
               matchTupleList2VsList

Author:        Wayne Dawson
creation date: 170126 (originally part of Threads.py & RThreads.py) 
last update:   190718
version:       0.1

Purpose:

This module is currently used only by "sarabande.py" to analyze RNA
structure. Because there are some differences in the analysis routines
of the chromatin, RNA and protein based program, I have had to make
these separate.

Comments:

190705:

   This is actually divided in two parts: one is the LThread data set
   and the other is the TreeNode analysis of the meta-structures. I
   think these can and should be separated.

   1. The LThread data is rather independent of the data sets. It is
      very similar to Pair data and I am thinking that eventually, I
      may merge these two data sets, as I don't much like the
      proliferation of data sets that I have created so far.

   2. The TreeNode data structures are more domain specific, so they
      probably cannot be so easily merged under a single processing
      strategy, but maybe.

   The idea situation is that both the LThread and TreeNode modules
   can be operated rather universally. I'm not sure this can be
   achieved, but I am still thinking the matter over.

earlier comments (updated 190705)

   The _intention_ was that this code would be shared among a variety of
   programs and modules.

   The current objective is to always produce human readable output,
   regardless of whether it comes from Calculate (an actual
   computation of a structure) or from an input structural sequence
   (dot bracket format or otherwise).

   The other goal is that the TreeNode programs can created
   intermediate test structures so that I can verify whether modules
   are working. So the design should work both from bottom up with the
   prediction programs and top down with the structure input programs.

"""

import sys
import string
from copy import deepcopy


# #######################
# Thermodynamic constants
# #######################

# A hardwired infinity. Somewhat hate this, but it seems there is no
# way around this darn thing.
from Constants import INFINITY
from Constants import kB   # [kcal/molK] (Boltzmann constant)

# coarse-grained resolution factor (for RNA, _always_ should be 1)
#from RConstants import seg_len
# constants for entropy calculations
#from RConstants import xi    # [bp] default stem Kuhn length
#from RConstants import xi_fs # [nt] free strand Kuhn length
#from RConstants import lmbd  # [nt] bead-to-bead bond distance
#from RConstants import gmm   # dimensionless
#from RConstants import T37C  # [K] absolute temp at 37 C.
#from RConstants import T_0C  # [K] absolute temp at  0 C.
"""@

gamma is the dimensionless self avoiding walk parameter. Its value
represents roughly 2*gmm ~ D, where D is the dimensionality of the
system.

"""
# #########################
# general constants for RNA
# #########################
#from RConstants import minStemLen      # minimum stem length
#from RConstants import max_bp_gap      # max gap between bps
#from RConstants import minLoopLen      # minimum loop length
#from RConstants import pk_scan_ahead   # hot lead length for pk
#from RConstants import dGpk_threshold  # threshold FE for pks
#from RConstants import dGMI_threshold  # threshold FE for M-/I-loops
#from RConstants import set_dangles     # dangle parameter (always = 2)
#from RConstants import dG_range        # FE range in suboptimal structures


"""for Vienna Pair object representation"""
from Pair import find_this_ij
from Pair import is_ap_stem
from Pair import is_pp_stem
from Pair import Pair
from Pair import SortPair
from Pair import vsPair2list

from Vienna    import Vstruct      # 1D structure analysis tool
from Constants import sysDefLabels # system RNA, Chromatin

# LThread object representation
from LThread import LNode
from LThread import DispLThread
from LThread import LThread
from LThread import LThread2Vienna

"""Motif (helping make all these things unified)"""
from Motif import MBL
from Motif import XLoop
from Motif import Stem
from Motif import PseudoKnot
from Motif import ArchMP
from Motif import MultiPair
from Motif import disp_Motif
from Motif import copyStem
from Motif import ins_sort_StemList

"""for Motif data representation"""
from Motif      import Motif # a core object of these building programs
from Motif      import Link
from Motif      import LGroup
from Motif      import Map # main map for all the FE and structures



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
from BasicTools import sortPairListWRT_n 
from BasicTools import copyList
from BasicTools import tuple2List
from BasicTools import roundoff
from BasicTools import initialize_matrix

"""Tree structure building/analysis tools"""
from NaryTree   import Node
from NaryTree   import NodeAnalysis
from NaryTree   import TreeBuilder

from LThreadBuilder  import LThreadBuilder
from TreeNode        import TreeNode2DotBracket
from TreeNode        import TreeNode2Motif



# ################################################################
# ############  Local-global functions and constants  ############
# ################################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# 1. control parameters tests etc in the program

# Test   function
#  0     insert into a list
#  1     parser for simple structures
#  2     convert a dot bracket structure representation to Motif

TEST = 2


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




class LThread2DotBracket(object):
    
    def __init__(self, N, fe):
        self.N      = N
        self.vSeq   = []
        self.fe     = fe
        self.sseq   = self.N*'c'
        
    #
    
    
    def getLThread2DotBracket(self, lt, structure_layout = 1, is_chromatin = True):
        """@
        
        This method generates the 1D meta-structure representation
        from an input from the module and class LThread (specifically,
        a list of objects of class LNodes).
        
        lt: input should generally be from the list "thread" in class
            LThread. It is a list of objects of class LNode.
        
        structure_layout: a formatting option for the output. 
           0: only structure (no sequence)
           1: both sequence and structure
        
        is_chromatin: chromatin does not have a real sequence like RNA
           and proteins. Therefore, a pseudosequence has to be
           generated
        
        """
        
        debug_getLThread2DotBracket = False # True # 
        if debug_getLThread2DotBracket:
            print "Enter getLThread2DotBracket"
        #
        
        # obtain viable vs Pair list from lt data
        vf = LThread2Vienna()
        vf.lt2vs(lt)
        
        # process data in Vstruct
        vsp = Vstruct()
        vsp.init_vs_from_lt(vf)
        
        # print self.vSeq 
        # build Tree diagram image
        v2tp = Vienna2TreeNode(vsp, self.fe, self.sseq)
        v2tp.vienna2tree()
        if debug_getLThread2DotBracket:
            print v2tp.dispTreeStruct(v2tp.ssTree,  "Secondary Structure Tree")
            print v2tp.dispTreeStruct(v2tp.genTree, "Full Structure Tree")
        #
        # print "0: ", v2tp.vseq
        
        # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
        tn2db = TreeNode2DotBracket(v2tp, is_chromatin)
        #print tn2db.seqv
        if not is_chromatin:
            tn2db.seqv = []
            for vv in self.vSeq:
                tn2db.seqv += [vv]
            #
        #
        
        tn2db.vseq = tn2db.makeFinal(self.vSeq)
        # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        
        #print "1:" , self.vSeq
        tn2db.make_dotbracket(is_chromatin)
        self.vSeq = tn2db.vSeq
        self.vstr = tn2db.vstr
        #print "2:" , self.vSeq
        
        s = ''
        if structure_layout == 0:
            s = tn2db.vstr + '\n' # this will only return the ss string
            
        elif structure_layout == 1:
            s = self.vSeq + '\n' + self.vstr + '\n'
        
        else:
            s  = self.vSeq + '\n'
            s += self.vstr + '\n'
        #
        
        return s
    #
#




# ###############################################################
# ###############   VARIOUS ASSORTMENT OF TOOLS   ###############
# ###############     USED BY Vienna2TreeNode     ###############
# ###############################################################

# these tools are heavily used by the various TreeNode type programs.

def getStem(ib, jb, stemlist, show = False):
    stem = None
    found = False
    for ssk in stemlist:
        if ib == ssk[0].i and ssk[0].j == jb:
            stem = Stem(ssk)
            found = True
            break
        #
    #
    if not found:
        print "ERROR: getStem(%d,%d)" % (ib, jb)
        sys.exit(1)
    if show:
        print stem.disp_Stem()
    #
    return stem
#



def makeStemList(ib, jb, stemlist):
    stems = []
    for ssk in stemlist:
        i = ssk[0].i; j = ssk[0].j 
        if ib <= i and j <= jb:
            stems += [Stem(ssk)]
        #
    #
    return stems
#

def dispStemList(snode, name):
    s = "** list of %s:" % name
    for ss in snode:
        s += '\n'
        s += disp_Motif(ss) # ss.disp_Stem(), ss.disp_PseudoKnot(), etc.
    #
    return s
#
    
    


def matchTupleList2VsList(tuplelist, Xlist):
    """
    Build a the vienna listing (Xlist) from a list of tuples
    (tuplelist): match (i,j) -> (i,j)[dir]
    """
    debug_matchTupleList2VsList = False
    if debug_matchTupleList2VsList:
        print "matchTupleList2VsList: "
    #
    newlist = []
    
    for tlk in tuplelist:
        i_r = tlk[0]; j_r = tlk[1]
        if debug_matchTupleList2VsList:
            print "ij_r = (%2d,%2d)" % (i_r, j_r)
        #
        flag_cont = True
        kr = 0
        while kr < len(Xlist) and flag_cont:
            i_p = Xlist[kr].i; j_p = Xlist[kr].j
            if debug_matchTupleList2VsList:
                print "ij_p = (%2d,%2d)" % (i_p, j_p)
            #
            if i_r == i_p and j_r == j_p:
                newlist += [Xlist[kr]]
                flag_cont = False
            else:
                kr += 1
            #
        #
    #
    return newlist
#


class Vienna2TreeNode(SortPair):
    """@@@
    
    Originally, this was called Vienna2LThread. I thought that it
    would work to use Threads as an intermediate to build objects of
    the Motif class. However, as I tried to make it work, it became
    clear that Threads is too poor a way to express this
    information. Perhaps it could be done that way, but the labor to
    convert it to something useful for Motif and making secondary
    structure diagrams would still be required like what is done
    here. So finally I abandoned that approach in favor of what is
    done now.
    
    Furthermore, I saw that I would have to do this conversion of the
    data to produce the structure diagrams anyway. So it make no sense
    to seriously pursue actually generating some output of class
    LThreads and then try to use it for other things. I will have to
    process it like I have here anyway.
    
    I think eventually, I add a few small modules that will generate a
    real LThread output; however, I conclude that it will be symbolic
    for the most part.
    
    """
    def __init__(self, vs, fe = None, seq = ""):
        # #######################  Vienna Data   #######################
        self.vs        = vs          # for employing functions
        self.system    = vs.system   # what system? RNA? Chromatin? 
        self.N         = vs.N        # sequence length
        self.vstr      = vs.vstr     # the 1D structure
        self.vseq      = vs.vseq     # the 1D sequence
        self.vsBPlist  = vs.BPlist   # the decomposed 1D ss structure
        self.vsBProots = vs.BProots  # the roots for the 1D ss structure
        self.vsPKlist  = vs.PKlist   # the decomposed 1D apPK structure
        self.vsPKroots = vs.PKroots  # the roots for the 1D apPK structure
        self.vsMPlist  = vs.MPlist   # the decomposed 1D CTCF structure
        # #####################   engine driver   #####################
        self.fe = fe # RNAModules or ChromatinModules
        
        if self.fe == None:
            
            # 1. Is there a sequence and is it the same length as the
            # structure
            
            # print "Module not set:"
            if self.vseq == "":
                print "beads in sequence not set"
                pair_nm = sysDefLabels[self.system]
                if seq == "":
                    # no information, hence, chose default
                    self.vseq = len(self.vstr)*pair_nm
                    
                else:
                    if len(seq) == len(self.vstr):
                        self.vseq = seq
                        
                    else:
                        # information inaccurate, hence, chose default
                        self.vseq = len(self.vstr)*pair_nm
                    #
                    
                #
                
                #print "1"
                
            else:
                if len(seq) == len(self.vstr):
                    self.vseq = seq
                else:
                    # information inaccurate, hence, chose default
                    pair_nm   = sysDefLabels[self.system]
                    self.vseq = len(self.vstr)*pair_nm
                #
                #print "2"
                
            #
            # print "vseq: ", self.vseq
            
            # 2. Choose the system under study and set the proper
            # modules to operate it. 
            
            if   self.system == "RNA":
                """contains RNA free energy parameters"""
                from RNAModules import BranchEntropy 
                from RNAModules import SetUpBranchEntropy 
                self.fe = BranchEntropy(SetUpBranchEntropy(self.vseq))
                #print self.fe.g_kBT() # just testing
                
            elif self.system == "Chromatin":
                """contains Chromatin free energy parameters"""
                from ChromatinModules import BranchEntropy 
                from ChromatinModules import SetUpBranchEntropy 
                self.fe = BranchEntropy(SetUpBranchEntropy(self.vseq))
                #print self.fe.g_kBT() # just testing
                #print self.fe.sseq
                
            else:
                print "ERROR: unrecognized system type (%s)" % self.system
                print "       allowed names "
                for snm in sysDefLabels.keys():
                    print snm
                #
                sys.exit(1)
            #
            #print "Vienna2TreeNode constructor: stop here"; sys.exit(0)
        #
        
        # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        # ############  holding variables at various stages ############
        # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        
        # these variables hold essential tail information that can
        # be used to figure out the stem structures from the secondary
        # structure inputs.
        self.apBPtails = [] # holds all the anti-parallel ss stem tails  
        self.ppBPtails = [] # holds all the parallel ss stem tails  
        self.apPKtails = [] # holds all the anti-parallel pk stem tails  
        self.ppPKtails = [] # holds all the parallel pk stem tails  
        
        # these are the next step in the process ... an intermediate
        # list that contains a group of pairs that form a
        # stem. Presently, the stem itself is a list, but this will be
        # used later to build the actual object of class Stem.
        self.ssaptails = [] # ** list of objects of class Stem
        self.pkaptails = [] # ** list of objects of class Stem
        self.sspptails = [] # *
        self.pkpptails = [] # *
        self.pkstructs = [] # holds objects of type PseudoKnot
        
        # ** used extensively
        # *  used a little bit
        
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # ############  holding variables at various stages ############
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        
        
        # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        # ###############  main results from this class  ###############
        # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        self.ssStems   = [] # final list of ss Stems (both a and p)
        self.genStems  = [] # final list PKs and Stems
        self.MPlist    = [] # the decomposed 1D CTCF structure
        self.ssTree    = None # final ss Tree
        self.genTree   = None # final general Tree (minus MultiPair)
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # ###############  main results from this class  ###############
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        
        self.max_bp_gap = fe.max_bp_gap # 
        
        # ######  program flags  ########
        # Presently, these tags initialized_X have very limited
        # utility. The tag "self.initialize_v2t" is already set to
        # True when this constructor is run.
        
        # =================================
        self.initialized_bplist = False # not really used
        # =================================
        self.initialize_v2t   = True
    #endInit
    
    def vienna2tree(self):
        if not self.initialize_v2t:
            print "ERROR: Vienna2TreeNode() is not initialized!"
            sys.exit(1)
        #
        
        debug_vienna2tree = False # True # 
        if debug_vienna2tree:
            self.disp_ViennaLists("Vienna data on entering vienna2tree")
            # print "stop 0 in vienna2tree", sys.exit(0)
        #
        
        self.build_BPstruct()
        
        self.build_MPstruct()
    #endMethod
    
    def disp_ViennaLists(self, sttmnt = "in Thread from Vienna"):
        print sttmnt
        print "structural sequence"
        self.vs.print_vseq()
        self.vs.print_vstr()
        print "secondary structure"
        self.vs.print_Xlist_n(self.vsBPlist)
        print "secondary structure roots"
        self.vs.print_Xlist_n(self.vsBProots)
        
        if len(self.vsPKlist) > 0:
            print "pseudoknot linkages: "
            self.vs.print_Xlist_n(self.vsPKlist)
            print "pseudoknot linkage roots"
            self.vs.print_Xlist_n(self.vsPKroots)
        #
        
        if len(self.vsMPlist) > 0:
            print "MP connects"
            self.vs.print_Xlist_n(self.vsMPlist)
        #
        
        if len(self.vsMPlist) > 0:
            islands = []
            for cl in self.vsMPlist:
                islands += [self.vs.expand_island(cl)]
            #endfor
            
            print "MP breakdown:"
            kk = 1
            for island_k in islands:
                print "island(%d):" % kk
                for ii in island_k:
                    print ii.disp_Pair()
                #endfor
                
                kk += 1
            #endfor
            
        #
    #endMethod
    
    
    def build_BPstruct(self):
        debug_build_BPstruct = False # True # 
        if debug_build_BPstruct:
            print "Enter build_BPstruct()"
        #
        
        # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
        """
        I suspect that there is a way to economize this step, but
        presently, I am not completely certain how.
        
        It seems that find_apStem_tails() and find_ppStem_tails()
        serve two purposes and are therefore called twice. The purpose
        _here_ is to actually look up stem tails as markers. This is
        needed because we are looking for pseudoknots (PKs), so we
        need to find their locations relative to the structures they
        are attached to. So here, we actually need to simply "find"
        the stem tails.
        
        When build_BPstruct() is completed, then we need to go through
        and assign the "bgn" and "end" terms to the corresponding
        stems in the structure. After this whole thing is completed,
        we need to mark the new stems we found in the structure to
        help the next level of abstraction (where we build these
        things into motifs) to build the corresponding motif
        correctly.
        
        Anyway, if it seems like find_apStem_tails and
        find_ppStem_tails are used a lot, well ... they are.
        
        """
        # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        
        # here, we just look up the tails, nothing else!
        
        # build the first level of bp list
        self.apBPtails = self.find_apStem_tails(self.vsBPlist, False)
        if debug_build_BPstruct:
            print "apBPtails:  ", self.apBPtails
            #print "stop at 0 in build_BPstruct"; sys.exit(0)
        #
        
        self.ppBPtails = self.find_ppStem_tails(self.vsBPlist, False)
        if debug_build_BPstruct:
            print "ppBPtails:  ", self.ppBPtails
            #print "stop at 1 in build_BPstruct";  sys.exit(0)
        #
        
        
        apPKtails = self.find_apStem_tails(self.vsPKlist, False)
        ppPKtails = self.find_ppStem_tails(self.vsPKlist, False)
        
        pkroots = []
        for vv in apPKtails:
            pkroots += [vv]
        #endfor
        
        for vv in ppPKtails:
            pkroots += [vv]
        #endfor
        
        pkroots = sortPairListWRT_n(pkroots, 0)
        
        # Update self.vsPKroots
        del self.vsPKroots
        self.vsPKroots = matchTupleList2VsList(pkroots, self.vsPKlist)
        
        if debug_build_BPstruct:
            print "apPKtails:  ", apPKtails
            print "ppPKtails:  ", ppPKtails
            print "pkroots:    ", pkroots
            #print "stop at 2 in build_BPstruct"; sys.exit(0)
        #
        
        # anti-parallel stems
        ss_apstemlist, pk_apstemlist = self.findStemsIn_apList(apPKtails)
        #sys.exit(0)
        # parallel stems
        ss_ppstemlist, pk_ppstemlist = self.findStemsIn_ppList(ppPKtails)
        #sys.exit(0)
        if debug_build_BPstruct:
            print "ap ssstems:"
            for k in range(0, len(ss_apstemlist)):
                print "%2d: %s" % (k, self.disp_stem(ss_apstemlist[k]))
            #
            print "ap pkstems: "
            for k in range(0, len(pk_apstemlist)):
                print "%2d: %s" % (k, self.disp_stem(pk_apstemlist[k]))
            #
            print "pp ssstems:"
            for k in range(0, len(ss_ppstemlist)):
                print "%2d: %s" % (k, self.disp_stem(ss_ppstemlist[k]))
            #
            print "pp pkstems: "
            for k in range(0, len(pk_ppstemlist)):
                print "%2d: %s" % (k, self.disp_stem(pk_ppstemlist[k]))
            #
            #print "stop at 3 in build_BPstruct"; sys.exit(0)
        #
        
        # This next step is long and probably could be broken down
        # into smaller methods, but I had to figure out some way to
        # get this to work, and so this is what I presently have come
        # up with.
        rt_Kpk, rt_ppKpk, rt_Rpk \
            = self.buildTreeDmns(ss_apstemlist, # anti-parallel
                                 ss_ppstemlist, # parallel
                                 pk_apstemlist, # anti-parallel
                                 pk_ppstemlist) # parallel   
        if debug_build_BPstruct:
            print "rt_Kpk:   ", rt_Kpk
            print "rt_ppKpk: ", rt_ppKpk
            print "rt_Rpk:   ", rt_Rpk
            print "BProots"
            for vv in self.vsBProots:
                print vv.disp_Pair()
            #endfor
            
            print "PKroots"
            for vv in self.vsPKroots:
                print vv.disp_Pair()
            #endfor
            
            print self.dispTreeStruct(self.ssTree,  "Secondary Structure Tree")
            print self.dispTreeStruct(self.genTree, "Full Structure Tree")
            #print "stop at 4 (exit from) build_BPstruct"; sys.exit(0)
        #
    #endMethod
    
    
    def find_apStem_tails(self, Xlist, assign_lt = True):
        """
        This goes through the entire cycle looking for all anti-parallel
        stems and parallel stems selecting out the anti-parallel
        stems. So, unlike find_apStem_head() that looks for the head
        of the stem, this routine looks though the list and each time
        it encouters something that looks like a stem, it search for
        the next stem head. So this scans the whole list, it is not
        just a search of a given stem.
        
        test example
         0         10        20        30        40        50        60        70        80        90
         |         |         |         |         |         |         |         |         |         | 
        "{((((((.A.AAAAA......))))).((((........a..aaaaa....))))..).|..([)]([)].ABC.abc.}....{.....}"
        
        """
        
        debug_find_apStem_tail = False # True # 
        if debug_find_apStem_tail:
            print "Enter find_apStem_tail"
        #
        
        tails = []
        kvs = 0      # vienna index
        flag_cont = True
        n = len(Xlist)
        if kvs >= n:
            flag_cont = False
        #
        
        while flag_cont:
            if kvs >= n:
                flag_cont = False
            else:
                if Xlist[kvs].v == 'p' or Xlist[kvs].v == '-':
                    kvs += 1
                    
                else: 
                    if debug_find_apStem_tail:
                        print "ap: Xlist tail index: ", kvs, (Xlist[kvs].get_ssPair())
                    #
                    kvso = kvs
                    # Here, kvs and kvso are the tail of some new stem
                    i, j, v, c = Xlist[kvs].get_ssPair()
                    # now find the head of that stem at kvso
                    kvs = self.find_apStem_head(kvs, Xlist)
                    tails += [(i,j)]
                    kvs +=1
                #
                
            #
            
        #endwhile
        
        if debug_find_apStem_tail:
            print "exit find_apStem_tail: index = ", kvs
            print tails
            #print "planned exit at end of find_apStem_tail(): "; sys.exit(0)
        #
        
        return tails
    #endMethod
    
    def find_apStem_head(self, kvs, Xlist):
        """
        This works up the channel at a given stem and finds the head
        of the stem. It is very primative so the criteria for
        getting to the head of the stem is that the contiguity is
        lost. So it still needs broader definitions for stem in
        general.
        """
        debug_find_apStem_head = False # True # 
        if debug_find_apStem_head:
            print "enter find_apStem_head"
        #
        
        flag_cont = True
        if kvs >= len(Xlist) - 1:
            flag_cont = False
        #
        
        while flag_cont:
            if debug_find_apStem_head:
                print "ap head check: ", \
                    kvs,     (Xlist[kvs  ].get_ssPair()), \
                    (kvs+1), (Xlist[kvs+1].get_ssPair())
            #
            i_0, j_0, v_0, c_0 = Xlist[kvs  ].get_ssPair()
            i_1, j_1, v_1, c_1 = Xlist[kvs+1].get_ssPair()
            if kvs + 1 == len(Xlist)-1:
                flag_cont = False
            if i_0 + 1 == i_1 and j_0 - 1 == j_1:
                kvs+=1
            else:
                flag_cont = False
            #
            
        #endWhile
        
        if debug_find_apStem_head:
            print "exit find_apStem_head: index = ", kvs
        #
        
        return kvs
    #endMethod
    
    
    def find_ppStem_tails(self, Xlist, lt_assign = True):
        """
        This goes through the entire cycle looking for all
        anti-parallel and parallel stem selecting out the parallel
        stems. So, unlike find_ppStem_head() that looks for the head
        of the stem, this routine looks though the whole list and
        each time it encouters something that looks like a stem, it
        search for the next stem head. So this scans the whole list,
        it is not just a search of a given stem.
        
        test example
         0         10        20        30        40        50        60        70        80        90
         |         |         |         |         |         |         |         |         |         | 
        "{((((((.A.AAAAA......))))).((((........a..aaaaa....))))..).|..([)]([)].ABC.abc.}....{.....}"
        """
        
        debug_find_ppStem_tails = False # True # 
        if debug_find_ppStem_tails:
            if lt_assign:
                print "enter find_ppStem_tail, lt_assign = True"
            else:
                print "enter find_ppStem_tail, lt_assign = False"
            #
            
        #
        
        tails = []
        kvs = 0          # vienna index for Xlist
        flag_cont = True
        n = len(Xlist)
        while flag_cont:
            if kvs >= n:
                flag_cont = False
            else:
                if Xlist[kvs].v == 'a'  or Xlist[kvs].v == '-':
                    kvs += 1
                else:
                    if debug_find_ppStem_tails:
                        print "pp: Xlist tail index: ", kvs, (Xlist[kvs].get_ssPair()) 
                    #
                    kvso = kvs
                    i, j, v, c = Xlist[kvs].get_ssPair()
                    kvs = self.find_ppStem_head(kvs, Xlist)
                    tails += [(i,j)]
                    kvs +=1
                #
                
            #
            
        #endWhile
        
        if debug_find_ppStem_tails:
            print "exit find_ppStem_tail: index = ", kvs
        #
        
        return tails
    #endMethod
    
    def find_ppStem_head(self, kvs, Xlist):
        """
        This works up the channel at a given stem and finds the head
        of the parallel stem. It is very primative so the criteria
        for getting to the head of the stem is that the contiguity
        is lost. So it still needs broader definitions for stem in
        general.
        """
        debug_find_ppStem_head = False # True # 
        if debug_find_ppStem_head:
            print "enter find_ppStem_head"
        #
        
        flag_cont = True
        if kvs >= len(Xlist) - 1:
            flag_cont = False
        #
        
        while flag_cont:
            if kvs >= len(Xlist) - 1:
                flag_cont = False
            else:
                if debug_find_ppStem_head:
                    print "pp head check: ", \
                        kvs,     (Xlist[kvs  ].get_ssPair()),\
                        (kvs+1), (Xlist[kvs+1].get_ssPair())
                #
                i_0, j_0, v_0, c_0 = Xlist[kvs  ].get_ssPair()
                i_1, j_1, v_1, c_1 = Xlist[kvs+1].get_ssPair()
                if i_0 + 1 == i_1 and j_0 + 1 == j_1:
                    kvs+=1
                else:
                    flag_cont = False
                #
                
            #
            
        #endWhile
        
        if debug_find_ppStem_head:
            print "exit find_ppStem_head: index = ", kvs
        #
        
        return kvs
    #endMethod
    
    
    
    
    def findStemsIn_apList(self, apPKtails):
        """
        Find the anti-parallel root stems in self.vsBPlist that forms the
        base of the PK found in self.vsPKlist
        
        This routine is for __anti-parallel__ stems. It takes the
        given list of pk stem tails (apPKtails) and reduces the list
        down to the unique set of root domains. In multilayered MBLs,
        this helps to find the base of the stem forming the hierarchy.
        
        This proved to be useful in PK calculations because it finds
        the base of the PK-stem(s) and removes the other components in
        the track of the PK-stem(s). In principle, the PK can overlap
        at multiple places on the chain in its current design. I
        suppose this might be an issue where I really do need to
        include an "I-s" at every juncture of the PK-stem, but
        presently, I think these sort of bizzare examples other than a
        single semi-contiuous stem (which we can deal with) is
        sufficient.
        
        """
        debug_findStemsIn_apList  = False # True # 
        debug_findStemsIn_apListx = False # True # 
        stop_at_end = False # True # 
        
        # x = show details of calculations done on the different pairs
        
        flag_cont = True
        if debug_findStemsIn_apList:
            print "Enter findStemsIn_apList: "
            print "initial apPKtails: ", apPKtails
            #print "stop at 0 in findStemsIn_apList"; sys.exit(0)
        #
        
        """
        First reduce the contents of self.apBPtails list to express only
        the root domains in the secondary structure.
                
        Because the search producing apPKtails is a rather crude scan
        for any anti-parallel (ap) stems, it seems that pks can be
        contained in the list of tails. Here, we try to remove the
        pseudoknots, where we define the root stem as the one on the
        left (going 5' to 3' from left to right). If there is evidence
        of a tail from any linkage stem on the right, then we write
        that into pkTailList and remove the item from ssTailList.
        
        """
        pkapTailList = copyList(apPKtails)
        ssapTailList = tuple2List(copyList(self.apBPtails))
        
        
        if debug_findStemsIn_apList:
            print "before doing this additional sorting into ss and pk stems:"
            print "0: ssapTailList: ", ssapTailList
            print "0: pkapTailList: ", pkapTailList
            #print "stop at 1 in findStemsIn_apList"; sys.exit(0)
        #
        
        kr = 0
        while kr < (len(ssapTailList)-1):
            ikr0 = ssapTailList[kr  ][0]; jkr0 = ssapTailList[kr  ][1]
            ikr1 = ssapTailList[kr+1][0]; jkr1 = ssapTailList[kr+1][1]
            if debug_findStemsIn_apListx:
                print "ikr0(%2d) < ikr1(%2d) < jkr0(%2d) < jkr1(%2d)" \
                    % (ikr0, ikr1, jkr0, jkr1)
            #
            
            if ikr0 < ikr1 and jkr0 < jkr1 and ikr1 < jkr0:
                if debug_findStemsIn_apListx:
                    print "rearranging to pklist", ssapTailList[kr+1]
                #
                
                pkapTailList += [(ikr1, jkr1)]
                del ssapTailList[kr+1]
                
            else:
                kr += 1
            #
            
        #endWhile
        
        pkapTailList = sortPairListWRT_n(pkapTailList, 0)
        
        # build a vslist (sstails) from self.vsBPlist
        self.ssaptails =  matchTupleList2VsList(ssapTailList, self.vsBPlist)
        # build a vslist (pktails) from self.vsPKlist and self.vsBPlist
        self.pkaptails =  matchTupleList2VsList(pkapTailList, self.vsPKlist)
        self.pkaptails += matchTupleList2VsList(pkapTailList, self.vsBPlist)
        # second component (self.vsBPlist and self.vsPKlist) are used as
        # reference to pkapTailList and written to pktails.
        
        self.pkaptails = self.sortvsList(self.pkaptails, 'i')
        
        if debug_findStemsIn_apList:
            print "BPlist:"
            for vv in self.vsBPlist:
                print vv.disp_Pair()
            #endfor
            
            print "BProots:"
            for vv in self.vsBProots:
                print vv.disp_Pair()
            #endfor
            
            print "ssaptails: "
            for vv in self.ssaptails:
                print vv.disp_Pair()
            #endfor
            
            print "pkaptails"
            for vv in self.pkaptails:
                print vv.disp_Pair()
            #endfor
            
            #print "stop at 2 in findStemsIn_apList";  sys.exit(0)
        #
        
        
        """
        In this next exercise, we play pktails against ssapTailList to
        identify the root stem-tail of the PK regions. For core PKs,
        the root will correspond to a stem. For an extended PK, the PK
        root will consist of the first stem on the left that opens the
        critical PK and the last stem on the right side that closes
        it.
        
        For example, the following structure (without the PK overlap)
        has two domains: (1,36) and (39, 74). Two linkage stems
        connect these two domains at (10,49) and (25, 64). Therefore,
        the extended PK domain extends between 1 and 74. In this part
        of the routine, ssapTailList will be modified so that (1,36) and
        (39,74) are conjoined to produce a domain (1,74).
        
                  0         10        20        30        40        50        60        70    
                  |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .
        ss_seq = ".(((.(((..AA..)))...(((..BB...))).)))..(((.(((..aa..)))...(((..bb...))).)))"
        initial apBPtails:  [(1, 36), (5, 16), (20, 32), (39, 74), (43, 54), (58,70)]
                             -------                     --------
        
        initial pkapTailList:  [(10, 49), (25, 64)]
        updated ssapTailList: [(1,74)]
        
        In the next example, the result will be that the overall
        domain (id1,jd4); i.e., encompassing four domains where the
        first and fourth domain indices are saved in this step; i.e.,
        the extended pseudoknot runs over (2,45) with the domains
        involving the pseudoknot located at (2,13) and
        (33,45). However, the domains between [2, 45] are (2,13),
        (17,21), (25,29) and (33,45).
        
        0         10        20        30        40        50
        |    .    |    .    |    .    |    .    |    .    |
        ..(((.[[[..)))...(...)...(...)...(((.]]]...)))..,
          ---- 1 -----   - 2 -   - 3 -   ---- 4 ------
          ^          ^   ^   ^   ^   ^   ^           ^
         id1        jd1 id2 jd2 id3 jd3 id4         jd4
        
        """
        
        for pktail in self.pkaptails: 
            ipk = pktail.i; jpk = pktail.j
            if debug_findStemsIn_apListx:
                print "fsiapl: ijpk = (%2d,%2d)" % (ipk, jpk)
            #
            
            kr = 0
            while kr < (len(ssapTailList)-1):
                ikr = ssapTailList[kr][0]; jkr = ssapTailList[kr][1]
                if debug_findStemsIn_apListx:
                    print "ijkr = (%2d,%2d)[%2d]" % (ikr, jkr, kr)
                #
                
                if ipk == ikr and jpk == jkr:
                    if debug_findStemsIn_apListx:
                        print "delete ssapTailList(ks=%d) : " % ks, ssapTailList[ks]
                        print "I think this job should already have been done, "
                        print "so I should stop to investigate. Presently, it has"
                        print "not stopped, so maybe I don't need to consider this"
                        sys.exit(0)
                    #
                    
                    del ssapTailList[ks]
                #
                
                kr += 1
            #endWhile
            
        #endfor
        
        
        
        if debug_findStemsIn_apList:
            print "build stems from the current information:"
        #
        
        ss_apstemlist  = self.build_apStem(self.ssaptails, self.vsBPlist)
        ss_apstemlist  = self.ins_sort_stemlist(ss_apstemlist)
        
        pk_apstemlist  = self.build_apStem(self.pkaptails, self.vsPKlist)
        pk_apstemlist += self.build_apStem(self.pkaptails, self.vsBPlist)
        
        
        #print pk_apstemlist
        #sys.exit(0)
        pk_apstemlist  = self.ins_sort_stemlist(pk_apstemlist)
        if debug_findStemsIn_apList:
            print "BPlist:"
            for vv in self.vsBPlist:
                print vv.disp_Pair()
            #endfor
            
            print "BProots:"
            for vv in self.vsBProots:
                print vv.disp_Pair()
            #endfor
            
            print "ap ssstem:"
            for k in range(0, len(ss_apstemlist)):
                print "%2d: %s" % (k, self.disp_stem(ss_apstemlist[k]))
            #endfor
            
            print "ap pkstem: "
            for k in range(0, len(pk_apstemlist)):
                print "%2d: %s" % (k, self.disp_stem(pk_apstemlist[k]))
            #
            
            #print "stop at 0 in findStemsIn_apList"; sys.exit(0)
            if stop_at_end:
                print "stop at the end of findStemsIn_aplist"; sys.exit(0)
            #endfor
            
        #
        
        return ss_apstemlist, pk_apstemlist
        ##################################
    #endMethod


    def buildTreeDmns(self,
                      ss_apstemlist,
                      ss_ppstemlist,
                      pk_apstemlist,
                      pk_ppstemlist):
        
        """#### 
        
        Builds a ss tree and a full general tree with the exception of
        the MP interactions. 
        
        For the MP, I am basically trying to work around the matter. For free
        energy calculations, This is easily done because we do not
        have to work depth first but can simply add the individual
        contributions together from each LNode in LThread. However, it
        is a bit more tricky with TreeNode2Motif, because the
        structure _must_ be built depth first. I think for any
        arbitrary horrendous mess, the best we can do is evaluate the
        tree nodes first and back door introduce the corrected free
        energy for the MP interactions. The other possibility might be
        to build a tree with the MPs, and the fit smaller trees inside
        of the MPs. It also depends partly on the nature of the MP
        binding. If it is weak, it is easier to do the former. If it
        is strong, like the ctcfs it was originally designed to
        conquer, then probably the latter is better.
        
        """
        
        debug_buildTreeDmns  = False # True # 
        debug_buildTreeDmnsx = False # True # 
        stop_at_end          = False # True # 
        
        
        if debug_buildTreeDmns:
            print "Enter buildTreeDmns()"
        #
        
        ss_stemlist = []
        for ssk in ss_ppstemlist:
            ss_stemlist += [ssk]
        #endfor
        
        for ssk in ss_apstemlist:
            ss_stemlist += [ssk]
        #endfor
        
        # sort stems
        ss_stemlist = self.ins_sort_stemlist(ss_stemlist) # must sort!!!
        
        pk_stemlist = []
        for ssk in pk_ppstemlist:
            pk_stemlist += [ssk]
        #endfor
        
        for ssk in pk_apstemlist:
            pk_stemlist += [ssk]
        #endfor
        
        # sort stems
        pk_stemlist = self.ins_sort_stemlist(pk_stemlist) # must sort!!!
        
        if debug_buildTreeDmns:
            print "----------------------------"
            print "ss_stemlist(after sorting): "
            for vv in ss_stemlist:
                print self.disp_stem(vv)
            #endfor
            
            print "pk_stemlist(after sorting): "
            for vv in pk_stemlist:
                print self.disp_stem(vv)
            #endfor
            
            print "----------------------------"
            #print "stop at 0 in buildTreeDmns"; sys.exit(0)
        #
        
        
        
        if debug_buildTreeDmns:
            print "make ssStems and pkStems and build ss Tree:"
        #
        
        self.ssStems = makeStemList(0, self.N, ss_stemlist)
        self.pkStems = makeStemList(0, self.N, pk_stemlist)
        # pkStems is not really used at all!
        
        # basic tree to be used to make structural assessments.  Note
        # that this is not the final tree, it is an intermediate tree.
        self.ssTree  = self.buildTreeStruct(self.ssStems, 0, self.N-1)
        """ 
        NOTE!!! N-1 
        
        .... Whereas the length of the list is N, the list runs from 0
        to N-1, __NOT__ 1 to N. Because future versions of this tool
        may permit inclusion of MultiPair interactions, the boundaries
        (ib,jb) <=> (0, N-1) are the actual positions on the sequence,
        not merely the end of the sequence. On subsequences, either
        all of them would have to be entered as jb+1, or we do that
        here for the root stem. I think it is better to have
        buildTreeStruct assume that the boundaries (ib,jb) are actual
        indices of actual positions _on_ the sequence, not convenient
        notation conventions.
        
        """
        # print self.ssTree
        
        # #######################################
        # Deal with the parallel stem case:
        # #######################################
        
        # because the parallel cases are rather special and don't
        # involve so many very complex, intricately entwined
        # structures, we address the parallel stems first.
        ssppStems = makeStemList(0, self.N, ss_ppstemlist)
        pkppStems = makeStemList(0, self.N, pk_ppstemlist)
        # unlike self.pkStems, pkppStems is actually used.
        
        if debug_buildTreeDmns:
            print dispStemList(self.ssStems, "regular stems")
            print dispStemList(self.pkStems, "pseudoknot stems")
            print self.dispTreeStruct(self.ssTree, "secondary structure tree")
            print "subset of parallel stems:"
            print "len: ssppStems = %2d, pkppStems = %2d" \
                % (len(ssppStems), len(pkppStems))
            print dispStemList(ssppStems, "root    part of pp type stems")
            print dispStemList(pkppStems, "linkage part of pp type stems")
            #print "stop at 1 in buildTreeDmns"; sys.exit(0)
        #
        
        if debug_buildTreeDmns:
            print "** parallel stem case"
        #
        
        dpkpp = {} # keep a record of whether a connection was found
        for kpk in range(0, len(pkppStems)):
            ipk = pkppStems[kpk].it; jpk = pkppStems[kpk].jt
            if debug_buildTreeDmnsx:
                print "btd: ijpk(%2d,%2d)" % (ipk, jpk)
            #
            
            dpkpp.update({(ipk, jpk) : False })            
        #endfor
        
        rt_ppKpk = [] # potential extended PKs
        for krt in range(0, len(ssppStems)):
            if debug_buildTreeDmnsx:
                print "ssppStems(%2d): %s" % (krt, ssppStems[krt])
            #
            
            irt = ssppStems[krt].it; jrt = ssppStems[krt].jt
            newpk = [[(irt, jrt)], []] 
            
            kpk = 0
            while kpk < len(pkppStems):
                ipk = pkppStems[kpk].it; jpk = pkppStems[kpk].jt
                if debug_buildTreeDmnsx:
                    print "ijrt(%2d,%2d), ijpk(%2d,%2d)" % (irt, jrt, ipk, jpk)
                #
                
                if irt < ipk and ipk < jrt and jrt < jpk:
                    if debug_buildTreeDmnsx:
                        print "use rt(%2d,%2d): -> pk(%2d,%2d)" % (irt, jrt, ipk, jpk)
                    #
                    
                    newpk[1] += [(ipk, jpk)]
                    if dpkpp.has_key((ipk, jpk)):
                        dpkpp[(ipk, jpk)] = True
                    else:
                        dpkpp.update({(ipk, jpk) : True})
                    #
                    
                else:
                    
                    if debug_buildTreeDmnsx:
                        if (jrt < ipk) or (jpk < irt):
                            if jrt < ipk:
                                print "irt(%d) < jrt(%d) < ipk(%d) < jpk(%d)" \
                                    % (irt, jrt, ipk, jpk)
                            elif jpk < irt:
                                print "ipk(%d) < jpk(%d) < irt(%d) < jrt(%d)" \
                                    % (ipk, jpk, irt, jrt)
                            else:
                                print "here xx"
                                sys.exit(1)
                            #
                            
                        else:
                            if irt < ipk:
                                print "irt(%d) < ipk(%d) < jpk(%d) < jrt(%d)" \
                                    % (irt, ipk, jpk, jrt)
                            elif ipk < irt:
                                print "ipk(%d) < irt(%d) < jrt(%d) < jpk(%d)" \
                                    % (ipk, irt, jrt, jpk)
                            else:
                                print "here yy"
                                sys.exit(1)
                            #
                            
                        #
                        
                    #endif
                    
                #
                kpk += 1
                # kpk always increments because this is just
                # searching, not searching and deleting.
            #endWhile
            
            rt_ppKpk += [newpk]
            
        #endfor
        
        k = 0
        while k < len(rt_ppKpk):
            if len(rt_ppKpk[k][1]) == 0:
                if debug_buildTreeDmns:
                    print "delete rt_ppKpk[%d]: %s" % (k, rt_ppKpk[k])
                #
                
                del rt_ppKpk[k]
            else:
                k += 1
            #
            
        #endWhile
        
        
        if debug_buildTreeDmns:
            print "(pp) rt_ppKpk:  ", rt_ppKpk
            print "dpkpp:          ", dpkpp
            print "pkppStems:      ", pkppStems
            print "self.pkpptails: ", self.pkpptails
            #print "stop at 2 in buildTreeDmns"; sys.exit(0)
        #
        
        # still processing the parallel stem cases
        
        self.pkstructs = []
        for k in range(0, len(rt_ppKpk)):
            if debug_buildTreeDmnsx:
                print "root: ", rt_ppKpk[k][0], ", linkage: ", rt_ppKpk[k][1]
            #
            
            """@@@
            
            Parallel stems are a rather special case because the
            progression for the PK is always towards the right. In the
            most general case, a PK can be either on the left or the
            right-hand side of the root stem. Here, we don't have to
            worry about that, so we just order the stems so they are
            always right leaning. In principle, this should already be
            ordered because the parallel stems were initially ordered.
            
            """
            
            pktails = sortPairListWRT_n(rt_ppKpk[k][1], 1) # order by j (index 1). 
            
            if debug_buildTreeDmnsx:
                print "rt_ppKpk[k]:     ", rt_ppKpk[k]
                print "rt_ppKpk[k][0]:  ", rt_ppKpk[k][0]
                print "rt_ppKpk[k][1]:  ", rt_ppKpk[k][1]
                print "pktails(sorted): ", pktails
            #
            
            pklinks = []
            for pkx in pktails:
                pL = pkx[0]; qL = pkx[1]
                if debug_buildTreeDmnsx:
                    print "pqL = (%2d,%2d)" % (pL, qL)
                #
                
                pklinks += [getStem(pL, qL, pk_ppstemlist)]
            #endFor
            
            
            if rt_ppKpk[k][0][0][0] < pktails[0][0]:
                n = len(pklinks) - 1
                m = len(pklinks[n].stem) - 1
                iz = rt_ppKpk[k][0][0][0]; jz = pklinks[n].stem[m].j
                
                
                # This looks a bit strange because position "it" is
                # derived from the root stem, position "jt" is derived
                # from the linkage stem.
                if debug_buildTreeDmnsx:
                    print "rt_ppKpk[%d][0][0][0] = %d, pktails[0][0] = %d" \
                        % (k, rt_ppKpk[k][0][0][0], pktails[0][0])
                    print "len(pklinks): ", len(pklinks)
                    print "pklinks:                    %s" % (pklinks)
                    print "pklinks[%2d]:                %s" % (n, pklinks[n].stem)
                    print "pklinks[n=%2d].stem[m=%2d].j: %s" % (n, m, pklinks[n].stem[m].j)
                    #print rt_ppKpk[k]
                    print "rt_ppKpk[%2d][0] = %s" % (k, rt_ppKpk[k][0])
                    #nv = len(rt_ppKpk[k][0]) - 1
                    #print rt_ppKpk[k][0][0]
                    #print rt_ppKpk[k][0][nv]
                    print "ijz: (%2d,%2d)" % (iz, jz)
                    #print "stop at 3 in buildTreeDmns"; sys.exit(0)
                #
                
                newpk = PseudoKnot(iz, jz, 'K')
                newpk.roottmp  = rt_ppKpk[k][0]
                newpk.linkages = pklinks
                self.pkstructs += [newpk]
                
            else:
                n = len(rt_ppKpk[k][0]) - 1
                iz = pklinks[0].stem[0].i; jz = rt_ppKpk[k][0][n][1] 
                
                
                # This looks a bit strange because position ".it" is
                # derived from the root stem, position ".jt" is derived
                # from the linkage stem.
                if debug_buildTreeDmnsx: 
                    print "rt_ppKpk[%d][0][0][0] = %d, pktails[0][0] = %d" \
                        % (k, rt_ppKpk[k][0][0][0], pktails[0][0])
                    print "len(pklinks): ", len(pklinks)
                    print "pklinks:                    %s" % (pklinks)
                    print "pklinks[0].stem[0]:              %s" % (pklinks[0].stem[0])
                    #print rt_ppKpk[k]
                    print "rt_ppKpk[%2d][0] = %s" % (k, rt_ppKpk[k][0])
                    print "rt_ppKpk[%2d][0][%2d] = %s" % (k, n, rt_ppKpk[k][0][n])
                    #nv = len(rt_ppKpk[k][0]) - 1
                    #print rt_ppKpk[k][0][0]
                    #print rt_ppKpk[k][0][nv]
                    print "ijz: (%2d,%2d)" % (iz, jz)
                    #print "stop at 3 in buildTreeDmns"; sys.exit(0)
                #
                
                newpk = PseudoKnot(iz, jz, 'K')
                newpk.roottmp  = rt_ppKpk[k][0]
                newpk.linkages = pklinks
                self.pkstructs += [newpk]
                #sys.exit(0)
        #endFor
        
        if debug_buildTreeDmnsx:
            print "-----"
            print "pkstructs: ", self.pkstructs
            print "-----"
            #print "stop at 4 in buildTreeDmns"; sys.exit(0)
        #
        
        # #######################################
        # Deal with the antiparallel stem case:
        # #######################################
        
                
        """
        Here, we need to update the search list to include parallel stem
        items that were not picked up in the above operation, as
        indicated by the assignment False in dpkpp. These need to be
        combined (in the next few steps) and processed in subsequent
        steps.
        
        This section could probably use some consolidation, but I have
        no time to do that right now. Currently, the tree searching
        algorithms use the class Pair -> (i,j), so they index slightly
        different from class Stem -> (it,jt). As a result, all this
        rigmarole here is because I need to construct an updated list
        PK stems using class Stem. It is probably not so hard to
        change self.scanForRdmns and self.is_extendedPK to receive
        type Stem, but presently, I just felt better to do it this
        way.
        
        """
        
        if debug_buildTreeDmns:
            print "update the pk stem tails with the actual references for the stems"
        #
        
        udpktails = []
        for dd in dpkpp.keys():
            # print dd
            if not dpkpp[dd]:
                i = dd[0]; j = dd[1]
                for kv in range(0, len(self.pkpptails)):
                    it = self.pkpptails[kv].i; jt = self.pkpptails[kv].j 
                    if i == it and j == jt:
                        udpktails += [ self.pkpptails[kv] ]
                        if debug_buildTreeDmns:
                            print "%s --> %s" % (dd, self.pkpptails[kv])
                        #
                        
                    #
                    
                #endFor
                
            #
            
        #endFor
        
        
        for vv in self.pkaptails:
            udpktails += [vv]
        #
        
        udpktails = self.sortvsList(udpktails, 'i')
        
        kt = 0
        while kt < len(udpktails):
            ipkt = udpktails[kt].i; jpkt = udpktails[kt].j
            flag_found = False
            for kr in range(0, len(self.pkStems)):
                ipkr = self.pkStems[kr].it; jpkr = self.pkStems[kr].jt
                
                if debug_buildTreeDmnsx:
                    print "ijpkr(%2d,%2d) < > ijpkt(%2d,%2d)" % (ipkr, jpkr, ipkt, jpkt)
                #
                
                if ipkr == ipkt and jpkt == jpkr:
                    flag_found = True
                    break
                #
                
            #endFor
            
            if not flag_found:
                if debug_buildTreeDmnsx:
                    print "delete: ", udpktails[kt]
                #
                
                del udpktails[kt]
            else:
                kt += 1
            #
            
        #endWhile
        
        
        # Now, back to the matter at hand!!!!
        
        
        if debug_buildTreeDmns:
            print "begin search for extended PKs"
            
            print "udpktails: (this should be the full list of viable candidates)"
            for vv in udpktails:
                print vv.disp_Pair()
            #endFor
            
            print "self.pkaptails: (a subset list)"
            for vv in self.pkaptails:
                print vv.disp_Pair()
            #endFor
            
            print "self.pkpptails: (another subset list)"
            for vv in self.pkpptails:
                print vv.disp_Pair()
            #endFor
            
            #print "stop at 5 in buildTreeDmns"; sys.exit(0)
        #
        
        # #####################################
        # #####  search for extended PKs  #####
        # #####################################
        
        rt_Rpk = [] # potential extended PKs
        #for pktail in self.pkaptails: 
        for pktail in udpktails:
            ipk = pktail.i; jpk = pktail.j
            iR, jR = self.scanForRdmns(self.ssTree, pktail)
            if debug_buildTreeDmnsx:
                print "ijR = (%2d,%2d)" % (iR, jR)
            #
            
            if jR > 0:
                if self.is_extendedPK(iR, jR, pktail):
                    dmns = self.findDmnsInRegion(self.ssTree, iR, jR)
                    if len(dmns) == 1:
                        print "only one domain!"
                        sys.exit(0) # should not end up here?
                    #
                    
                    if debug_buildTreeDmnsx:
                        print "is extended"
                        print dmns
                    #
                    
                    if dmns == None:
                        print "ERROR: after findDmnsInRegion{tree, ijR(%d,%d)" % (iR, jR)
                        print "       obtained domains was \"None\"!" 
                        sys.exit(1)
                    #
                    
                    if len(dmns) == 0:
                        print "ERROR: after findDmnsInRegion{tree, ijR(%d,%d)" % (iR, jR)
                        print "       obtained domains was \empty\"!" 
                        sys.exit(1)
                    #
                        
                    rt_Rpk += [[dmns, [(ipk, jpk)]]]
                    if debug_buildTreeDmnsx:
                        print "ijpk = (%2d,%2d), ijR = (%2d,%2d)" % (ipk, jpk, iR, jR)
                        print "rt_Rpk = ", rt_Rpk[len(rt_Rpk)-1]
                    #
                    
                #
                
            #
            
        #endFor
        
        if debug_buildTreeDmnsx:
            print "initial rt_Rpk: ", rt_Rpk
            #print "stop at 6 in buildTreeDmns"; sys.exit(0)
        #
        
        # compact the results if necessary.
        
        # Presently, the rt_Rpk is a long list that can be rudundant.
        kr = 0
        while kr < (len(rt_Rpk)-1):
            iRr = rt_Rpk[kr][0][0][0]; jRr = rt_Rpk[kr][0][0][1]
            if debug_buildTreeDmnsx:
                print "ijRr(%2d,%2d)", iRr, jRr
            #
            
            ks = kr + 1
            flag_del = False
            while ks < len(rt_Rpk):
                iRs = rt_Rpk[ks][0][0][0]; jRs = rt_Rpk[ks][0][0][1]
                if debug_buildTreeDmnsx:
                    print "ijRs = ", iRs, jRs
                #
                
                if iRr == iRs and jRr == jRs:
                    if debug_buildTreeDmnsx:
                        print "delete ks = ", ks
                        print "rt_Rpk(kr = %2d, %s), rt_Rpk(kr = %2d, %s)" \
                            % (kr, rt_Rpk[kr][1], ks, rt_Rpk[ks][1])
                    #
                    
                    rt_Rpk[kr][1] += rt_Rpk[ks][1]
                    del rt_Rpk[ks]
                    flag_del = True
                    break
                
                else:
                    ks += 1
                #
                
            #endWhile
            
            if not flag_del:
                kr += 1
            #
            
        #endWhile
        
        if debug_buildTreeDmns:
            print "revised rt_Rpk: ", rt_Rpk
            #sys,exit(0)
        #
        
        
        # #################################
        # #####  search for core PKs  #####
        # #################################
        
        if debug_buildTreeDmns:
            print "begin search for core PKs"
        #
        
        rt_Kpk = [] # potential core PKs
        #for pktail in self.pkaptails: 
        for pktail in udpktails:
            ilk = pktail.i; jlk = pktail.j
            idr, jdr = self.scanForKdmns(self.ssTree, pktail)
            if debug_buildTreeDmnsx:
                print "ijlk = (%2d,%2d)" % (ilk, jlk)
                print "ijdr = (%2d,%2d)" % (idr, jdr)
            #
            
            if jdr > 0:
                testcore = self.is_corePK(pktail)
                #print "testcore = ", testcore
                if testcore:
                    iK = jK = 0
                    if ilk < idr:
                        iK = ilk; jK = jdr
                    else:
                        iK = idr; jK = jlk
                    #
                    
                    dmns = self.findDmnsInRegion(self.ssTree, iK, jK)
                    
                    rt_Kpk += [[dmns, [(ilk, jlk)]]]
                    if debug_buildTreeDmnsx:
                        print "ijlk = (%2d,%2d), ijdr = (%2d,%2d)" % (ilk, jlk, idr, jdr)
                        print "domain region: %s" % dmns
                    #
                    
                #
                
            #
            
        #endFor
        
        if debug_buildTreeDmnsx:
            print "initial rt_Kpk: ", rt_Kpk
            #print "stop at 7 in buildTreeDmns"; sys.exit(0)
        #
        
        
        # compact the results if necessary
        kr = 0
        while kr < (len(rt_Kpk)-1):
            n_dr = len(rt_Kpk[kr][0]) - 1
            idr = rt_Kpk[kr][0][0][0]; jdr = rt_Kpk[kr][0][n_dr][1]
            if debug_buildTreeDmnsx:
                print "rt_Kpk: ", rt_Kpk[kr][0], ", n_dr = ", n_dr
                print "ijdr = ", idr, jdr
                # print "stop at 7.1: compact the results"; sys.exit(0)
            #
            
            ks = kr + 1
            flag_del = False
            while ks < len(rt_Kpk):
                n_ds = len(rt_Kpk[ks][0]) - 1
                ids = rt_Kpk[ks][0][0][0]; jds = rt_Kpk[ks][0][n_ds][1]
                if debug_buildTreeDmnsx:
                    print "ijds = ", ids, jds
                #
                
                if idr == ids and jdr == jds:
                    if debug_buildTreeDmnsx:
                        print "delete ks = ", ks
                        print "rt_Kpk(kr = %2d, %s), rt_Kpk(kr = %2d, %s)" \
                            % (kr, rt_Kpk[kr][1], ks, rt_Kpk[ks][1])
                        #print "stop at 7.2: compact the results"; sys.exit(0)
                    #
                    
                    rt_Kpk[kr][1] += rt_Kpk[ks][1]
                    del rt_Kpk[ks]
                    flag_del = True
                    break
                else:
                    ks += 1
                #
                
            #
            if not flag_del:
                kr += 1
            #
            
        #endWhile
        
        if debug_buildTreeDmns:
            print "revised rt_Kpk: ", rt_Kpk
            # print "stop at 7.3: finished compact the results"; sys,exit(0)
        #
        
        # ######################################################################
        # #####  Deal with the resulting mixture of core and extended PKs  #####
        # ######################################################################
        
        # Now we look for overlaps between the different the R and K
        # pks and adjusting these points
        
        kR = 0 #  R => Extended PK
        while kR < len(rt_Rpk):
            nR = len(rt_Rpk[kR][0]) - 1
            iR = rt_Rpk[kR][0][0][0]; jR = rt_Rpk[kR][0][nR][1]
            if debug_buildTreeDmnsx:
                print "ijR = (%2d,%2d)" % (iR, jR),  ", rt_Rpk: ", rt_Rpk[kR]
            #
            
            kK = 0 # K => Core PK
            remove_kR = False
            while kK < len(rt_Kpk):
                nK = len(rt_Kpk[kK][0]) - 1
                iz = rt_Kpk[kK][0][0][0]; jz = rt_Kpk[kK][0][nK][1]
                if debug_buildTreeDmnsx:
                    print "ijz = (%2d,%2d)" % (iz, jz), ", rt_Kpk: ", rt_Kpk[kK]
                #
                
                """
                              |----------  linkage  ---------|
                              v                              v
                ....((((......[[[....)))).....(((.....)))..]]]...
                    |         ^         ^               |    ^
                    id        pL        jd              |    qL
                    |------ root -------|               |
                    |--------------  full region  ------|
                    iz                                  jz
                """
                
                if iz <=  iR and jR <= jz:
                    if iz == iR and jR == jz:
                        """Qxxxxxx
                        
                        I think this case is not possible because we
                        cannot have both a core PK and a true extended
                        PK both occupying exactly the same sites. Of
                        course, there could be two thermodynamically
                        stable structures of equal probability that
                        could satisfy this condition, but they would
                        simultaneously in the result here.
                        
                        """
                        print "ftd 01 reset iz(%2d) == iR(%2d) == jR(%2d) < jz(%2d)" \
                            % (iz, iR, jR, jz)
                        print "rt_Kpk[%2d][0] = %s, rt_Kpk[%2d][0] = %s" \
                            % (kK, rt_Kpk[kK][0], kR, rt_Rpk[kR][0])
                        print "rt_Kpk[%2d][1] = %s, rt_Kpk[%2d][1] = %s" \
                            % (kK, rt_Kpk[kK][1], kR, rt_Rpk[kR][1])
                        print "part 01 ftd. Probably should not happen, please check ";
                        sys.exit(0)
                    else:
                        """#@@@@@
                        Captures structures such as the following: 
                                  0         10        20        30        40        50        60        70        80        90
                                  |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |
                        ss_seq = ".....(((....AAA....)))...(((....[[[...)))...(((...]]]....)))....aaa...."
                                       |      |        |   |--------- extended PK -----------|      |
                                       |--------root---|---|-------------core PK ------------|------|
                                       id ----| ------ jd  iR                               jR      |
                                       iz ----|--------|----------------------------------- jz      |
                                              pL --------------------- linkage -------------------- qL
                        
                        So this is an example of an embedded extended
                        PK. The extended PK from 25 to 599, the core
                        PK that embeds this extends from 5 to 66
                        
                        we can ignore this case because the PK is inside the other PK
                        
                        """
                        print "02 ijR = (%2d, %2d) is embedded in ijz = (%2d,%2d)" \
                            % (iR, jR, iz, jz)
                    #
                    
                elif iR < iz and jz < jR:
                    if debug_buildTreeDmnsx:
                        print "4 ijz = (%2d, %2d) internal to ijR = (%2d,%2d)" \
                            % (iz, jz, iR, jR)
                    #
                    
                elif iR <  iz and jz  <= jR:
                    """#@@@@@
                    Captures structures such as the following: 
                              0         10        20        30        40        50        60        70        80        90
                              |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |
                    ss_seq = ".((....AA..BB..CC..))..DD....aa..bb..cc...EE...FF..GG..dd.....ee..ff..gg." # structure 2 K-type only
                    ss_seq = ".((....AA..BB..CC..))..((....aa..bb..cc...EE...FF..GG..)).....ee..ff..gg." # structure 2 K-type only
                               |                     |-------------------- core PK ---|--------------|
                               |--------------------- -extended PK -------------------|              |
                               iR                    |                               jR              |
                                                     iz=id                           jz=jd           |
                                                     pL ------------  linkage  -------|------------- qL
                    
                    ijz ==>    |++++++++++++++++++++++++++++++++++++++++++++++++++++++|    
                    
                    The core PK extends from 23 to 71
                    The extended PK extends from 1 to 56
                    The full mixed PK extends from 1 to 72
                    #@@@@@
                    """
                    if debug_buildTreeDmnsx:
                        print "1 reset iR(%2d) < iz(%2d) < jz(%2d) <= jR(%2d)" \
                            % (iR, iz, jz, jR)
                        print "rt_Kpk[%2d][0] = %s, rt_Kpk[%2d][0] = %s" \
                            % (kK, rt_Kpk[kK][0], kR, rt_Rpk[kR][0])
                        print "rt_Kpk[%2d][1] = %s, rt_Kpk[%2d][1] = %s" \
                            % (kK, rt_Kpk[kK][1], kR, rt_Rpk[kR][1])
                    #
                    
                    rt_Kpk[kK][0]  = rt_Rpk[kR][0] 
                    # include all domains of ijR in the root box ([0]) of core PK
                    rt_Kpk[kK][1] += rt_Rpk[kR][1]
                    # add the PKs from R into K
                    del rt_Rpk[kR]
                    remove_kR = True
                    # print "buildTreeDmns: stopped in mixture section 1"; sys.exit(0)
                    
                elif iR <= iz and jz <  jR:
                    
                    """#@@@@@
                    Captures structures such as the following: 
                              0         10        20        30        40        50        60        70        80        90
                              |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |
                    ss_seq = ".......AA..BB..CC......((....aa..bb..cc...DD...EE..FF..))..GG...dd..ee..ff...gg.." # structure 2 K-type only
                    ss_seq = ".......AA..BB..CC......((....aa..bb..cc...DD...EE..FF..))..((...dd..ee..ff...)).." # structure 2 K-type only
                                     |---------------|-- core PK ---------------------|
                                     |               |-----------------extended PK ---|---------------------|
                                     pL              iR                               qL                    jR
                                                     iz=id ---------------------------jz=jd
                    
                    ijz ==>                          |++++++++++++++++++++++++++++++++++++++++++++++++++++++|    
                    
                    The core PK extends from 7 to 56
                    The extended PK extends from 23 to 78
                    The full mixed PK extends from 7 to 78
                    #@@@@@
                    """
                    
                    if debug_buildTreeDmnsx:
                        print "2 reset iR(%2d) <= iz(%2d) < jz(%2d) < jR(%2d)" \
                            % (iR, iz, jz, jR)
                        print "rt_Kpk[%2d][0] = %s, rt_Kpk[%2d][0] = %s" \
                            % (kK, rt_Kpk[kK][0], kR, rt_Rpk[kR][0])
                        print "rt_Kpk[%2d][1] = %s, rt_Kpk[%2d][1] = %s" \
                            % (kK, rt_Kpk[kK][1], kR, rt_Rpk[kR][1])
                    #
                    
                    rt_Kpk[kK][0]  = rt_Rpk[kR][0]
                    # include all domains of ijR in the root box ([0]) of core PK
                    rt_Kpk[kK][1] += rt_Rpk[kR][1]
                    # add the PKs from R into K
                    del rt_Rpk[kR]
                    remove_kR = True
                    # print "buildTreeDmns: stopped in mixture section 2"; sys.exit(0)
                #
                
                kK += 1
                # here also, we just scroll through rt_Kpk without
                # deleting anything. The only thing that can be
                # deleted is rt_Rpk.
                
            #endWhile
            
            if not remove_kR:
                kR += 1
            #
            
        #endWhile
        
        if debug_buildTreeDmns:
            print "final rt_Kpk: ", rt_Kpk
            print "final rt_Rpk: ", rt_Rpk
            #print "pklinks:      ", pklinks
            print "pk_stemlist:  ", pk_stemlist
            print "next, build R-type PK"
            #print "stop at 8 in buildTreeDmns"; sys.exit(0)
        #
        
        # ##################################################################
        # #####   Build the linkage stems of the various R and K PKs   #####
        # ##################################################################
        
        for pkR in rt_Rpk:
            
            n = len(pkR[0])
            iR = pkR[0][0][0]; jR = pkR[0][n-1][1]
            pklinks = sortPairListWRT_n(pkR[1], 1) # order by j (index 1)
            if debug_buildTreeDmnsx:
                #print pkR
                #print pkR[0]
                #print pkR[0][0]
                #print pkR[0][n-1]
                print "ijR = (%d, %d)" % (iR, jR)
            #
            
            newpk = PseudoKnot(iR, jR, 'R')
            newpk.roottmp = pkR[0]
            for pkx in pklinks:
                pL = pkx[0]; qL = pkx[1]
                newpk.linkages += [getStem(pL, qL, pk_stemlist)]
            #endFor
            
            self.pkstructs += [newpk]
        #endFor
        
        if debug_buildTreeDmnsx:
            print "current pkstructs after running R"
            for pkk in self.pkstructs:
                print pkk.disp_PseudoKnot()
            #
            print "next, build K-type PK"
            #print "stop at 9 in buildTreeDmns"; sys.exit(0)
        #
        
        
        for pkK in rt_Kpk:
            
            n_rt = len(pkK[0]) - 1
            iz = pkK[0][0][0]; jz = pkK[0][n_rt][1]
            n_pk = len(pkK[1]) - 1
            #print "pkK0: ", pkK[0]
            #print "pkK1: ", pkK[1]
            pklinks = sortPairListWRT_n(pkK[1], 1) # order by j (index 1)
            
            # need to find the minimum pL and maximum qL for the next
            # step where we build the PK
            pLmin = self.N
            qLmax = 0
            for pkK1k in pklinks:
                # print pkK1k
                pL   = pkK1k[0]; qL = pkK1k[1]
                if pL < pLmin:
                    pLmin = pL
                #
                
                if qL > qLmax:
                    qLmax = qL
                #
                
            #endFor
            
            if debug_buildTreeDmnsx:
                # print pkK
                # print pkK[0]
                # print pkK[0][0]
                # print pkK[0][n_rt-1]
                print "ijz = (%2d, %2d), {pLmin(%2d), qLmax(%2d)}" \
                    % (iz, jz, pLmin, qLmax)
                print "pklinks: ", pklinks
                
            #
            
            newpk = None
            
            iK = self.N; jK = 0
            if pLmin < iz and qLmax < jz and iz < qLmax:
                #print "1"
                iK = pLmin; jK = jz
                newpk = PseudoKnot(iK, jK, 'K')
                
            elif iz < pLmin and jz < qLmax and  pLmin < jz:
                #print "2"
                iK = iz; jK = qLmax
                newpk = PseudoKnot(iK, jK, 'K')
                
            else:
                print "ERROR(buildTreeDmns near 10): somehow, core PK is not properly defined"
                print "       iz(%2d)   < pLmin(%2d) < jz(%2d)   < qLmax(%2d)" \
                    % (iz,   pLmin, jz,   qLmax)
                print "<or>   pLmin(%2d) < iz(%2d)   < qLmax(%2d) < jz(%2d)"   \
                    % (pLmin, iz,   qLmax, jz)
                self.disp_ViennaLists("buildTreeDmns near 10")
                sys.exit(1)
            #
            
            #print "pkK[0]: ", pkK[0]
            newpk.roottmp = pkK[0]
            for pkx in pklinks:
                pL = pkx[0]; qL = pkx[1]
                newpk.linkages += [getStem(pL, qL, pk_stemlist)]
            #
            
            self.pkstructs += [newpk]
        #endFor
        
        if debug_buildTreeDmns:
            print "current pkstructs after running K"
            for pkk in self.pkstructs:
                print pkk.disp_PseudoKnot()
            #
            print "next, adding stems and setting up Stems for pruning"
            #print "stop at 10 in buildTreeDmns"; sys.exit(0)
        #
        
        # debug_buildTreeDmnsx = True
        
        kss = 0
        while kss < len(self.ssStems):
            stmk = self.ssStems[kss]
            iss = stmk.it; jss = stmk.jt
            flag_del = False
            for kpk in range(0, len(self.pkstructs)):
                ipk = self.pkstructs[kpk].it; jpk = self.pkstructs[kpk].jt
                if (ipk == iss and jss <= jpk) or (ipk <= iss and jss == jpk):
                    
                    
                    if len(self.pkstructs[kpk].rootstems) > 0:
                        
                        """### 
                        
                        I had to add this because for both types of
                        PK, I may pass through here more than once
                        when there are multiple stems in the extended
                        PK or the core PK subregions. if I don't skip
                        when rootstems is already assigned, it will do
                        this search twice and may inadvertently add
                        multiple roots to this part.
                        
                        """
                        if debug_buildTreeDmnsx:
                            print self.pkstructs[kpk].pktype 
                            print "skipped because rootstem assigned. len ", \
                                len(self.pkstructs[kpk].rootstems)
                            #print "stop at 11 in buildTreeDmns"; sys.exit(0)
                        #
                        
                        continue
                    #
                    
                    
                    if debug_buildTreeDmnsx:
                        print "got to here"
                        print "self.pkstructs[kpk].pktype: ", self.pkstructs[kpk].pktype
                        print "ipk(%d) <= iss(%d) and jss(%d) == jpk(%d)" \
                        % (ipk, iss, jss, jpk)
                        if (ipk <= iss and jss == jpk):
                            print "pk linkage on the left hand side"
                        elif (ipk == iss and jss <= jpk):
                            print "pk linkage on the right hand side"
                        #
                        
                        print "ijpk(%2d,%2d), ijss(%2d,%2d)" % (ipk,jpk, iss, jss)
                        # print "del ijss region: ", self.pkstructs[kpk].roottmp
                    #
                    
                    for kv in range(0, len(self.pkstructs[kpk].roottmp)):
                        idv = self.pkstructs[kpk].roottmp[kv][0];
                        jdv = self.pkstructs[kpk].roottmp[kv][1];
                        if debug_buildTreeDmnsx:
                            print "ijdv: (%2d,%2d)" \
                                % (idv, jdv), ", roottmp: ", self.pkstructs[kpk].roottmp[kv]
                        #
                        
                        for lv in range(0, len(self.ssStems)):
                            issv = self.ssStems[lv].it; jssv = self.ssStems[lv].jt
                            if issv == idv and jssv == jdv:
                                if debug_buildTreeDmnsx:
                                    print "marking ijssv = (%2d,%2d) <=> ijdv = (%2d,%2d)" \
                                        % (issv, jssv, idv, jdv)
                                #
                                
                                stem = copyStem(self.ssStems[lv].stem)
                                self.pkstructs[kpk].rootstems += [stem]
                                self.ssStems[lv].mark = True
                            #
                            
                            if debug_buildTreeDmnsx:
                                if len(self.pkstructs[kpk].rootstems) > 0:
                                    print "ijssv = (%2d,%2d) vs ijdv = (%2d,%2d)" \
                                        % (issv, jssv, idv, jdv)
                                    if issv == idv and jssv == jdv:
                                        print "pkstructs[%2d].rootstem: %s" \
                                            % (kpk, self.pkstructs[kpk].rootstems)
                                    #
                                    
                                #
                                
                            #
                            
                        #endfor
                        
                    #endfor
                            
                    
                    flag_del = True
                    kss += 1
                    break
                #
                
            #endfor
            
            if not flag_del:
                kss += 1
            #
            
        #endwhile
        
        self.ssStems = ins_sort_StemList(self.ssStems)
        if debug_buildTreeDmns:
            print "---------------------------------------------------"
            print dispStemList(self.ssStems, "current secondary structure")
            print "---------------------------------------------------"
            print "current pkstructs pruning"
            for pkk in self.pkstructs:
                print pkk.disp_PseudoKnot()
            #
            
            #print "stop at 12 in buildTreeDmns"; sys.exit(0)
        #
        
        self.genStems = []
        for k in range(0, len(self.ssStems)):
            # print self.ssStems[k]
            stem = copyStem(self.ssStems[k].stem)
            if self.ssStems[k].mark:
                stem.mark = True
            #
            
            self.genStems += [stem]
        #endfor
        
        for pks in self.pkstructs:
            self.genStems += [pks]
        #endfor
        
        self.genStems = ins_sort_StemList(self.genStems)
        self.ssStems = ins_sort_StemList(self.ssStems)
        if debug_buildTreeDmnsx:
            print "---------------------------------------------------"
            print dispStemList(self.genStems, "general stem structure (before pruning)")
            print "---------------------------------------------------"
            #print "stop at 13 in buildTreeDmns"; sys.exit(0)
            print "Next, pruning the Stems and PKs"
        #
        
        k = 0
        while k < len(self.genStems):
            typenm = type(self.genStems[k]).__name__
            i = self.genStems[k].it; j = self.genStems[k].jt
            if debug_buildTreeDmnsx:
                print "%s(%2d,%2d) " % (typenm, i, j)
            #
            
            if typenm == "Stem":
                if self.genStems[k].mark:
                    if debug_buildTreeDmnsx:
                        print "delete ", self.genStems[k]
                    #
                    
                    del self.genStems[k]
                else:
                    k += 1
                #
                
            else:
                k += 1
            #
            
        #endwhile
        
        if debug_buildTreeDmns:
            print "---------------------------------------------------"
            print dispStemList(self.genStems, "general stem structure (after pruning)")
            print "---------------------------------------------------"
            
            for astm in self.genStems:
                print disp_Motif(astm)
            #endFor
            
            #print "stop at 14 in buildTreeDmns"; sys.exit(0)
        #
        
        """
        ***********************************************************
        ###########################################################
        !!!!!!!!!!!!!!!!!!!!!!!!  FINALLY  !!!!!!!!!!!!!!!!!!!!!!!!
        ###########################################################
        ***********************************************************
       
        ...not to be underwhelming after all of this ordeal
        """
        
        # ###########################################################
        self.genTree = self.buildTreeStruct(self.genStems, 0, self.N-1)
        """NOTE!!! N-1 
        
        .... Whereas the length of the list is N, the list runs from 0
        to N-1, __NOT__ 1 to N. Because future versions of this tool
        may permit inclusion of MultiPair interactions, the boundaries
        (ib,jb) <=> (0, N-1) are the actual positions on the sequence,
        not merely the end of the sequence. On subsequences, either
        all of them would have to be entered as jb+1, or we do that
        here for the root stem. I think it is better to have
        buildTreeStruct assume that the boundaries (ib,jb) are actual
        indices of actual positions _on_ the sequence, not convenient
        notation conventions.
        
        """
        # ###########################################################
        # print genTree
        self.initialized_bplist = True # almost no meaning!!!
        
        
        if debug_buildTreeDmns:
            print "final results from buildTreeDmns:"
            print "BProots"
            for vv in self.vsBProots:
                print vv.disp_Pair()
            #endFor
            
            print "rt_Kpk:    "
            for rt_Kpk_k in rt_Kpk:
                print rt_Kpk_k
            #endFor
            
            print "rt_Rpk:    "
            for rt_Rpk_k in rt_Rpk:
                print rt_Rpk_k
            #endFor
            
            print "rt_ppKpk:    "
            for rt_ppKpk_k in rt_ppKpk:
                print rt_ppKpk_k
            #endFor
            
            print "sstails:   "
            for sst in self.ssaptails:
                print sst.disp_Pair()
            #endFor
            
            print "pktails:   "
            for pkt in self.pkaptails:
                print pkt.disp_Pair()
            #endFor
            
            print self.dispTreeStruct(self.ssTree,  "Secondary Structure Tree")
            print self.dispTreeStruct(self.genTree, "Full Structure Tree")
            print "-----"
            
            if 0: #debug_buildTreeDmnsx:
                print "planned exit at end of buildTreeDmns: "
                #print "stop at 15 in buildTreeDmns"; sys.exit(0)
            #
            
        #
        
        """
        Examples of above output for the cases below:
        
                  0         10        20        30        40        50        60        70        80        90
                  |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |
        ss_seq = ".((....[[......))..]].."                       # example 1
        ss_seq = ".((....AA......))..((.....BB....))...bb..aa.." # example 2
        ss_seq = ".((....AA...))..BB....aa....bb."               # example 3
        ss_seq = ".((....AA..BB..CC..))..((....aa..bb..cc...DD...EE..FF..)).....dd..ee..ff." # example 4
        ss_seq = ".((....AA..BB..CC..))..DD....aa..bb..cc...EE...FF..GG..dd.....ee..ff..gg." # example 4 (alternative)
        
        example 1 produces the following output display below:
        rt_Kpk:
        ((1, 16), (7, 20))
        rt_Rpk: (empty)
        sstails:
        (   1,  16)[a]
        pktails:
        (   7,  20)[a]
        
        example 2 produces the following output display below:
        rt_Kpk:
        ((1, 16), (7, 42))
        ((19, 33), (26, 38))
        rt_Rpk: (empty)
        sstails:
        (   1,  16)[a]
        (  19,  33)[a]
        pktails:
        (   7,  42)[a]
        (  26,  38)[a]
        
        example 3 produces the following output display below:
        rt_Kpk: (empty)
        rt_Rpk:
        ((1, 29), (7, 23))
        ((19, 33), (26, 38))
        sstails:
        (   1,  13)[a]
        (  16,  29)[a]
        pktails:
        (   7,  23)[a]
        
        example 4 produces the following output display below:
        rt_Kpk:
        ((1, 56), (42, 63))
        ((1, 56), (47, 67))
        ((1, 56), (51, 71))
        rt_Rpk:
        ((1, 56), (7, 30))
        ((1, 56), (11, 34))
        ((1, 56), (15, 38))
        sstails:
        (   1,  20)[a]
        (  23,  56)[a]
        pktails:
        (   7,  30)[a]
        (  11,  34)[a]
        (  15,  38)[a]
        (  42,  63)[a]
        (  47,  67)[a]
        (  51,  71)[a]
        
        
        This pairs down repetition of the same PK. In the above
        action, it selects _all_ cases where the pk linkage
        satisfies the conditions. For example
        
                  0         10        20        30        40        50        60        70        80        90        100       110       120
                  |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |
        ss_seq = "...(((..((((..(...[[[..).))))..(((.(..]]].....).)))..)))...(((..((((..(...AAA..).))))..(((.(..BBB..aaa...).)))..)))..bbb"
        
        here, rt_Rpk = [((3, 55), (18, 40)), ((8, 46), (18, 40)), ((59, 105), (74, 101))]
        
        clearly, the pk (18,40) is expressed twice but the correct
        point from which to extract the extended PK is from (8, 50)
        rather than from (3, 55) or (14, 46). The entry ((3, 55),
        (18, 40)) is removed in this next step.
        
        rt_Kpk: 
        ((59, 114), (94, 119))
        rt_Rpk:
        ((8, 50), (18, 40))
        ((64, 109), (74, 101))
        sstails:
        (   3,  55)[a]
        (   8,  28)[a]
        (  14,  23)[a]
        (  31,  50)[a]
        (  35,  46)[a]
        (  59, 114)[a]
        (  64,  84)[a]
        (  70,  79)[a]
        (  87, 109)[a]
        (  91, 105)[a]
        pktails:
        (  18,  40)[a]
        (  74, 101)[a]
        (  94, 119)[a]
        
        """
        
        
        if stop_at_end:
            print "planned stop at the end of buildTreeDmns()";
            # sys.exit(0)
        #
        
        return rt_Kpk, rt_ppKpk, rt_Rpk
    #endMethod
    
    
    
    # ########################################################
    # #####   PseudoKnot tests using TreeNode approach   #####
    # ########################################################
    
    
    def is_corePK(self, lnkg):
        debug_is_corePK = False # True # 
        
        
        flag_adjoining_region = False
        ilk = lnkg.i; jlk = lnkg.j
        contact_l = 0
        contact_r = 0
        contact_point = {}
        for kr in range(0, len(self.ssaptails)):
            ikr = self.ssaptails[kr].i; jkr = self.ssaptails[kr].j
            if ilk == ikr and jlk == jkr:
                if debug_is_corePK:
                    print "skip: (%2d,%2d)" % (ikr, jkr)
                #
                
                continue
            #
            
            if ilk < ikr and jlk < jkr and ikr < jlk:
                
                # ...[[...((((......]]...))))...
                #    ^    ^          ^      ^
                #   ilk  ikr        jlk    jkr
                #              (right side)
                
                if debug_is_corePK:
                    print "adjoining ilk(%2d) < ikr(%2d) < jlk(%2d) < jkr(%2d)" \
                        % (ilk, ikr, jlk, jkr)
                #
                
                contact_r += 1
                contact_point.update({(ikr,jkr) : (ilk,jlk)})
                
            elif ikr < ilk and jkr < jlk and ilk < jkr:
                
                # ...((((...[[.....)))).....]]...
                #    ^      ^         ^      ^
                #   ikr    ilk       jkr    jlk
                #       (left side)
                
                if debug_is_corePK:
                    print "adjoining ikr(%2d) < ilk(%2d) < jkr(%2d) < jlk(%2d)" \
                        % (ikr, ilk, jkr, jlk)
                #
                
                contact_l += 1
                contact_point.update({(ikr,jkr) : (ilk,jlk)})
                
            #
            
            if flag_adjoining_region:
                break
            #
            
        #endFor
        if debug_is_corePK:
            print "contact_lr: ", contact_l, contact_r
            print "contact_point: ", contact_point
        #
        
        if contact_l > 0 and contact_r > 0:
            flag_adjoining_region = True
        #
        
        return not flag_adjoining_region
    #endMethod
    
    
    def is_extendedPK(self, iR, jR, pk):
        
        debug_is_extendedPK = False # True # 
        
        
        ipk = pk.i; jpk = pk.j
        if debug_is_extendedPK:
            print "ijpk = (%2d,%2d), ijR(%2d,%2d)" % (ipk, jpk, iR, jR)
            for ssk in self.ssaptails:
                print ssk
            #
            
        #
        
        flag_is_extendedPK  = False # True # 
        kr = 0
        while kr < (len(self.ssaptails)-1):
            ikr = self.ssaptails[kr].i; jkr = self.ssaptails[kr].j
            if debug_is_extendedPK:
                print "ie: ijkr = (%2d,%2d)[%2d]" % (ikr, jkr, kr)
            #
            
            if ipk == ikr and jpk == jkr:
                if debug_is_extendedPK:
                    print "ie: kr=%2d: pk(%2d,%2d) <=> sstail(%2d,%2d), skip" \
                        % (kr, ipk, jpk, ikr, jkr)
                #
                
                kr += 1
                continue
            elif not (iR <= ikr and jkr <= jR):
                if debug_is_extendedPK:
                    print "ie: kr=%2d: not (iR(%2d) <= ikr(%2d) < jkr(%2d) <= jR(%2d), skip" \
                        % (kr, iR, ikr, jkr, jR)
                #
                
                kr += 1
                continue
            elif (ikr < ipk and jkr < jpk and ipk < jkr):
                # We assume an ordered list, so we will encounter
                # (ikr,jkr) first and subsequently, we will encounter
                # (iks, jks). 
                ks = kr + 1
                ks_cont = True
                while ks < len(self.ssaptails):
                    if ks > 40:
                        print "recursion problems: ikr(%d) < ipk(%d) < jkr(%d) < jpk(%d)" \
                            % (ikr, ipk, jkr, jpk)
                        #
                        
                        sys.exit(1)
                    #
                    
                    iks = self.ssaptails[ks].i; jks = self.ssaptails[ks].j
                    if debug_is_extendedPK:
                        print "ie: ijks = (%2d,%2d)[%2d]" % (iks, jks, ks)
                    #
                    
                    if ipk == iks and jpk == jks:
                        if debug_is_extendedPK:
                            print "ie: ks=%2d: pk(%2d,%2d) <=> sstail(%2d,%2d), skip" \
                                % (ks, ipk, jpk, ikr, jkr)
                        #
                        
                        ks += 1
                        continue
                    
                    elif not (iR <= ikr and jkr <= jR):
                        if debug_is_extendedPK:
                            print "ie: ks=%2d: not (iR(%2d) <= ikr(%2d) < jkr(%2d) <= jR(%2d), skip" \
                                % (ks, iR, ikr, jkr, jR)
                        #
                        
                        ks += 1
                        continue
                    
                    elif (ipk < iks and jpk < jks and iks < jpk):
                        if debug_is_extendedPK:
                            print "found ss(%2d,%2d)[kr=%d], ss(%2d,%2d)[ks=%d]: -> ijR(%d, %d)" \
                                % (self.ssaptails[kr].i, self.ssaptails[kr].j, kr, 
                                   self.ssaptails[ks].i, self.ssaptails[ks].j, ks, 
                                   self.ssaptails[kr].i, self.ssaptails[ks].j) 
                        #
                        """
                        Extended PK: expand the region to include the
                        two root stems the bound the extended PK. For
                        example,
                        
                        ..((((....[[...))))...((((...]]..))))..
                          ^               ^   ^             ^
                         ikr            jkr iks           jks
                        
                          ^                                 ^
                          iR                                jR
                        """
                        flag_is_extendedPK = True
                        ks_cont = False
                        ks += 1
                        break
                    
                    elif jpk < iks:
                        # ordered list!!!!
                        if debug_is_extendedPK:
                            print "no additional branch found"
                        #
                        flag_is_extendedPK = False
                        ks_cont = False
                        break
                        
                    else:
                        ks += 1
                    #
                    
                #endwhile
                if flag_is_extendedPK:
                    break
                #
                
                kr += 1
            #
            
            kr += 1
        #endwhile
        
        if debug_is_extendedPK:
            print "Exiting is_extendedPK"
            #print "stop at 0 (return from) in is_extended"; sys.exit(0)
        #
        
        return flag_is_extendedPK
    #endMethod
    
    
    
    def scanForRdmns(self, tree, pk):
        if tree.name == None:
            # really nothing there! ... if we got to here, the
            # structure calculation very likely has some serious
            # problems!
            return 0, 0
        else:
            return self.scanForRdmnsUtil(tree, pk)
        #
        
    #endMethod
    
    def scanForRdmnsUtil(self, curr, pk):
        debug_scanForRdmns = False # True # 
        if curr == None:
            # no information! ... node could be empty
            return 0, 0
        
        else:
            iR = curr.name.it; jR = curr.name.jt
            n = len(curr.children)
            if n > 0:
                iR = curr.children[0].name.it; jR = curr.children[n-1].name.jt
                if iR < pk.i and pk.j < jR:
                    iRndx = -1
                    jRndx = -1
                    for kr in range(0, len(curr.children)):
                        cc = curr.children[kr].name
                        ir = cc.it; jr = cc.jt
                        if ir < pk.i and pk.i < jr:
                            iRndx = kr
                        elif ir < pk.j and pk.j < jr:
                            jRndx = kr
                        #
                        
                    #
                    
                    if debug_scanForRdmns:
                        print "ijRndx: ", iRndx, jRndx
                    #
                    
                    if iRndx >= 0 and iRndx < jRndx:
                        iR = curr.children[iRndx].name.it
                        jR = curr.children[jRndx].name.jt
                    else:
                        for k in range(0, len(curr.children)):
                            iR_t, jR_t = self.scanForRdmnsUtil(curr.children[k], pk)
                            if jR_t > 0:
                                iR = iR_t; jR = jR_t
                            #
                            
                        #
                        
                    #
                    
                else:
                    iR = 0; jR = 0       
                #
                
            else:
                iR = 0; jR = 0       
            #
            
            return iR, jR
        #
        
    #endMethod
    
    def scanForKdmns(self, tree, pk):
        if tree.name == None:
            print "scanForKdmns: root node is empty"
            # really nothing there! ... if we got to here, the
            # structure calculation very likely has some serious
            # problems!
            return 0, 0
        else:
            return self.scanForKdmnsUtil(tree, pk)
        #
        
    #endMethod
    
    def scanForKdmnsUtil(self, curr, pk):
        debug_scanForKdmns = False # True # 
        
        if curr == None:
            #print "n"
            # no information! ... node could be empty
            return 0, 0
        
        else:
            id_r = curr.name.it; jd_r = curr.name.jt
            ib  = curr.name.ih; jb  = curr.name.jh
            if ib == 0 and jb == self.N-1:
                ib = -1; jb = self.N
            #
            
            if debug_scanForKdmns:
                print "ijd_r(%2d,%2d), ijb(%2d,%2d) ijpk(%2d,%2d)" \
                    % (id_r, jd_r, ib, jb, pk.i, pk.j)
                print "pk.i(%d) < id_r(%d) < pk.j(%d) < jd_r(%d)" % (pk.i, id_r, pk.j, jd_r)
            #endif
            
            if ((id_r < pk.i and jd_r < pk.j and pk.i < jd_r) or
                (pk.i < id_r and pk.j < jd_r and id_r < pk.j)):
                if debug_scanForKdmns:
                    if (id_r < pk.i and jd_r < pk.j and pk.i < jd_r):
                        print "id_r(%2d) < pk.i(%2d) < jd_r(%2d) < pk.j(%2d)" \
                            % (id_r, pk.i, jd_r, pk.j)
                    else:
                        print "pk.i(%2d) < id_r(%2d) < pk.j(%2d) < jd_r(%2d)" \
                            % (pk.i, id_r, pk.j, jd_r)
                    #
                    
                #
                
                return id_r, jd_r
            
            elif ib < pk.i and pk.j < jb:
                if debug_scanForKdmns:
                    print "ib(%2d) < pk.i(%2d) < pk.j(%2d) < jb(%2d)" \
                        % (ib, pk.i, pk.j, jb)
                #
                
                n = len(curr.children)
                if n > 0:
                    for k in range(0, len(curr.children)):
                        id_t, jd_t = self.scanForKdmnsUtil(curr.children[k], pk)
                        if jd_t > 0:
                            id_r = id_t; jd_r = jd_t
                        #
                        
                    #
                    
                else:
                    id_r = 0; jd_r = 0       
                #
                
            else:
                id_r = 0; jd_r = 0       
            #
            
            return id_r, jd_r
        #
        
    #endMethod
    
    
    def findDmnsInRegion(self, tree, iK, jK):
        if tree.name == None:
            # really nothing there! ... if we got to here, the
            # structure calculation very likely has some serious
            # problems!
            return []
        else:
            return self.findDmnsInRegionUtil(tree, iK, jK)
        #
        
    #endMethod
    
    def findDmnsInRegionUtil(self, curr, iK, jK):
        debug_findDmnsInRegion = False # True # 
        if debug_findDmnsInRegion:
            print "Enter findDmnsInRegionUtil()"
        #
        
        if curr == None:
            # no information! ... node could be empty
            return []
        
        else:
            ib = curr.name.it; jb = curr.name.jt
            if debug_findDmnsInRegion:
                print "ib(%2d) <= iK(%2d) < jK(%2d) <= jb(%2d)" % (ib, iK, jK, jb)
            #
            
            if ib <= iK and jK <= jb:
                
                n = len(curr.children)
                dmns = []
                if debug_findDmnsInRegion:
                    print "n = ", n
                #
                
                for k in range(0, n):
                    iV = curr.children[k].name.it; jV =  curr.children[k].name.jt
                    if debug_findDmnsInRegion:
                        print "ib(%2d) <= iV(%2d) < jV(%2d) <= jb(%2d)" \
                            % (ib, iV, jV, jb)
                    #
                    
                    # we encounter the boundary of the core PK
                    if iK <= iV and jV <= jK:
                        if debug_findDmnsInRegion:
                            print "ijV = (%d,%d)" % (iV, jV)
                        #
                        
                        dmns += [(iV, jV)]
                    #
                    
                #endfor
                
                if debug_findDmnsInRegion:
                    print "dmns: ", dmns
                #
                
                # print dmns, len(dmns)
                if len(dmns) == 0:
                    # print "entered len(dmns) = 0"
                    if n > 0:
                        for k in range(0, len(curr.children)):
                            ibx = curr.children[k].name.it; jbx = curr.children[k].name.jt
                            if debug_findDmnsInRegion:
                                # print "vv: ", curr.children[k].name
                                print "ib(%2d) <= iK(%2d) < jK(%2d) <= jb(%2d)" \
                                    % (ibx, iK, jK, jbx)
                            #
                            
                            if ibx <= iK and jK <= jbx:
                                dmns = self.findDmnsInRegionUtil(curr.children[k], iK, jK)
                            #
                            
                        #endfor
                        
                    #
                    
                #
                
                # print "xx: ", dmns
                return dmns
            
            else:
                return []
            #
            
        #
        
    #endMethod
    
    # #################################################
    # #####    General Tree Building Utilities    #####
    # #################################################
    
    def dispTreeStruct(self, tree, name = "tree diagram:"):
        s = "%s:\n" % name
        na = NodeAnalysis()
        na.structLayout(tree)
        s += na.disp_structLayout()
        return s
    #endMethod
    
    def buildTreeStruct(self, stems, ib, jb):
        debug_bts = False # True # 
        tb = TreeBuilder()
        root = Pair()
        root.put_ssPair(ib, jb, 'X', 'base')
        # the base of the structure is the root of the tree
        base = Stem([root]) # <<<====
        base.vtype = "root"
        base.name  = "Base"
        tb.tree.name = base
        # print tb.tree
        k, tb.tree = self.buildTreeStructUtil(tb.tree, stems, 0, 1, [1, 0], ib, jb)
        if debug_bts:
            print "k = %d"
            print self.dispTreeStruct(tb.tree, "check results")
            #print "stop at 0 (end of) in buildTreeStruct"; sys.exit(0)
        #
        
        return tb.tree
    #endMethod
    
    def increment_child(self, set_arr):
        ch_arr = []
        for cc in set_arr:
            ch_arr += [cc]
        #
        ch_arr[0] += 1
        ch_arr += [0]
        return ch_arr
    #endMethod
    
    def buildTreeStructUtil(self, node, stems, kstm, level, set_arr, ib, jb):
        debug_btsU = False # True # 
        if debug_btsU:
            print "buildTreeStructUtil{kstm(%2d), level(%2d), set_arr(%s), ijb(%2d,%2d)" \
                % (kstm, level, set_arr, ib, jb)
        #
        
        if level > 30:
            print "ERROR: !! something must be wrong in the recursion"
            print "kstm = %d, level = %d, set_arr %s" % (kstm, level, set_arr)
            sys.exit(1)
        #
        
        if kstm >= len(stems):
            if debug_btsU:
                print "kstm(%d) + 1 = %d is out of range" % (kstm, kstm+1)
            #
            
            return kstm, node
        #
        
        it1 = stems[kstm].it; jt1 = stems[kstm].jt
        if debug_btsU:
            print "ib(%2d) < it1(%2d) < jt1(%2d) < jb(%2d), level = %d" \
                % (ib, it1, jt1, jb, level)
        #
        
        if ib <= it1 and jt1 <= jb:
            ch_arr = self.increment_child(set_arr)
            new_level = level + 1
            node_type = type(stems[kstm]).__name__
            if debug_btsU:
                print "set_arr(%s) -> ch_arr(%s)" % (set_arr, ch_arr)
                print "%s(%2d,%2d) -> " % (node_type, stems[kstm].it, stems[kstm].jt)
            #
            
            
            node.children = [Node(stems[kstm])]
            kstm, node.children[set_arr[level]] \
                = self.buildTreeStructUtil(node.children[set_arr[level]],
                                           stems,
                                           kstm + 1,
                                           new_level,
                                           ch_arr,
                                           it1, jt1)
            if debug_btsU:
                print "return from buildTreeStructUtil: level = %d, kstm = %d" \
                    % (level, kstm) 
                print "node = %s" % (node)
            #
            
            flag_cont = True
            while kstm < len(stems) and flag_cont:
                set_arr[level] += 1
                it1 = stems[kstm].it; jt1 = stems[kstm].jt
                if debug_btsU:
                    print "update set_arr(%s)" % (set_arr)
                    print "while kstm(%2d), level(%2d), set_arr(%s)" \
                        % (kstm, level, set_arr)
                    print "                 ijt1 = (%2d,%2d), ijb(%2d,%2d)" \
                        % (it1, jt1, ib, jb)
                #
                
                if not (ib <= it1 and jt1 <= jb):
                    if debug_btsU: print "break"
                    flag_cont = False
                    
                else:
                    ch_arr = self.increment_child(set_arr)
                    #ch_arr[new_level] += 1
                    if debug_btsU:
                        print "                              ch_arr(%s)" % (ch_arr)
                    #
                    
                    node.children += [Node(stems[kstm])]
                    kstm, node.children[set_arr[level]] \
                        = self.buildTreeStructUtil(node.children[set_arr[level]],
                                                   stems,
                                                   kstm + 1,
                                                   new_level,
                                                   ch_arr,
                                                   it1, jt1)
                #
                
            #
            
            if debug_btsU:
                print "return from buildTreeStructUtil: kstm = %2d, level = %2d"  \
                    % (kstm, level)
                # if level == 1: print "planned exit";  print node; sys.exit(0)
            #
            
            return kstm, node
        
        else:
            if debug_btsU:
                print "return from level %d, kstm = %d" % (level, kstm)
            #
            
            return kstm, node
        
        #
        
        if debug_btsU:
            print "got to here!!! xxxx"
            print "ijb = ", (ib, jb), ", tree = ", tree
            print "back kstm = %d, level = %d, ijt1 = (%2d,%2d)" \
                % (kstm, level, it1, jt1)
            print "tree:   ", node
            #print "stop at 0 (return from) in buildTreeStructUtil"; sys.exit(0)
        #
        
        return kstm, node
    #endMethod
    
    
    # ###################################vvvvvvvvvvvvvvvvvvv
    # ####  Stem construction tools  ####vvvvvvvvvvvvvvvvvvv
    # ###################################vvvvvvvvvvvvvvvvvvv
    
    """
    !!!!!!!!!!!!!!
    #### NOTE ####
    !!!!!!!!!!!!!!
    
    The stem in this part of the code is STILL A LIST!!!!
    
    later it will be converted to the object of Stem class.
    """    
    
    # sort stemlist according to i.
    def ins_sort_stemlist(self, stemlist):
        for i in range(1,len(stemlist)):    
            j = i                    
            while j > 0 and stemlist[j][0].i < stemlist[j-1][0].i: 
                stemlist[j], stemlist[j-1] = stemlist[j-1], stemlist[j]
                # syntactic sugar: swap the items
                j=j-1 
            #endwhile
            
        #endfor
        
        return stemlist
    #endMethod
    
    
    def disp_stem(self, stem):
        s = "len =  %2d, " % len(stem)
        
        for k in range(0, len(stem)-1):
            s += "%s, " % stem[k].disp_Pair()
        #endfor
        
        s += "%s" % stem[len(stem)-1].disp_Pair()
        return s
    #endMethod
    
    def build_apStem(self, tailtags, Xlist):
        """
        This is a support method for build_apStem(). Anti-parallel stems
        have a property that follows the rule
        
        stem = [(i, j), (i+1, j-1), (i+2, j-2) .... ]
        
        in short, (i,j) -> (i + ki, j - kj), where ki and kj > 0 and
        generally, ki = kj = 1.
        
        
        The important aspect of this part of the code is that we
        define the criteria (typically rather simple) for what
        constitutes an "anti-parallel stem". This comes from my days
        with RNA where I came to understand that a stem need not be a
        contiguous segment. For example, the following parallel stem
        is surely contiguous.
        
        ....((((......))))....
        
        However, what can we say about the following stem?
        
        ...(..(..(..(......)...)...)...)...
        
        
        It would surely depend on the type of interactions, which
        should be module sepecific. Nevertheless, perhaps my intuition
        would be that it is weakly bound and perhaps the pairs Aa, Bb,
        Cc and Dd are largely independent.
        
        How about this example. 
        
        ...((((((..........((((....))))...........))))))....
        
        Here, I think the two stems should be considered independent.
        
        How about this example?
        
        ...(.(.(.(......).).).)...
        
        or this one
        
        ...((((((..((((....))))..))))))....
        
        
        Its regularity, closeness and symmetry would tend to pursuade
        me that it is effectively a single stem, not four separated
        and isolated stems.
        
        In the very first external call to this general program, the
        scan of the sequence only identifies stems with a contiguous
        set of base pairs. In this method, some simple rules for
        parsing for effective stems is introduced to help resolve
        these issues.
        
        For a concrete example on some dsRNA, consider this case. 
        
        GGGGGGGGGuuuCCCCCCCCCC
        |||||||||   ||||||||||
        CCCCCCCCCuuuGGGGGGGGGG
        
        I would probably call this effectively a single stem, not two
        different stems. 
        
        This method is intended to help merge these types of cases for
        ssRNA. A later method would handle dsRNA.
        
        """
        
        
        debug_build_apStem = False # True # 
        if debug_build_apStem:
            print "Enter build_apStem:"
            print "tailtags: ", tailtags
            print "Xlist:    ", Xlist
            #sys.exit(0)
        #
        
        stemlist = []
        for sst in tailtags:
            if debug_build_apStem:
                print "sst = ", sst
            #
            
            for k_apr in range(0, len(Xlist)):
                ib = Xlist[k_apr].i; jb = Xlist[k_apr].j
                if Xlist[k_apr].v == 'p':
                    # parallel stem!
                    continue
                #
                
                if ib == sst.i and jb == sst.j:
                    if debug_build_apStem:
                        print "ijb = ", (ib, jb)
                    #
                    
                    stem = self.scanFor_apStem(Xlist, k_apr)
                    if debug_build_apStem:
                        print "build_apStem result:"
                        print stem, len(stem)
                        # sys.exit(0)
                    #endif
                    
                    stemlist += [stem]
                    break
                
                #
                
            #endfor
            
        #
        
        if debug_build_apStem:
            for k in range(0, len(stemlist)):
                print k, stemlist[k]
            #
            
        #
        
        kr = 0
        while kr < (len(stemlist)):
            if debug_build_apStem:
                print kr, stemlist
                print "kr(%2d), stemlist: %s" % (kr, stemlist[kr])
            #
            
            tail = 0
            i_tr = stemlist[kr][tail].i; j_tr = stemlist[kr][tail].j
            head = len(stemlist[kr])-1
            i_hr = stemlist[kr][head].i; j_hr = stemlist[kr][head].j
            if debug_build_apStem:
                print "stem{(%2d,%2d):(%2d,%2d)}" % (i_tr, j_tr, i_hr, j_hr)
            #
            
            kt = kr + 1
            while kt < len(stemlist):
                i_tt = stemlist[kt][tail].i; j_tt = stemlist[kt][tail].j
                if debug_build_apStem:
                    print "build_apStem, test (k=%2d)(%2d,%2d)" % (kt, i_tt, j_tt)
                #
                
                if (i_tr <= i_tt and
                    i_tt <= i_hr and 
                    j_hr <= j_tt and
                    i_tt <= j_tr):
                    if debug_build_apStem:
                        print "del %d, %s" % (kt, stemlist[kt])
                    #
                    del stemlist[kt]
                else:
                    kt += 1
                #
                
            #endwhile
            
            if debug_build_apStem:
                print "kr(%d) -> %d" % (kr, kr + 1)
            #
            
            kr += 1
        #endwhile
        
        if debug_build_apStem:
            print "results from build_apStem:"
            for slk in stemlist:
                print slk
            #endfor
            
            #print "stop at 0 (end of) build_apStem"; sys.exit(0)
        #
        
        return stemlist
    #endMethod
    
    
    
    
    
    def scanFor_apStem(self, Xlist, k_init):
        """
        This is a support method for build_apStem(). Anti-parallel stems
        have a property that follows the rule
        
        stem = [(i, j), (i+1, j-1), (i+2, j-2) .... ]
        
        in short, (i,j) -> (i + ki, j - kj), where ki and kj > 0 and
        generally, ki = kj = 1.
        """
        
        stem_len = 1
        i_t = Xlist[k_init].i; j_t = Xlist[k_init].j
        
        debug_scanFor_apStem = False # True # 
        if debug_scanFor_apStem:
            print "Enter scanFor_apStem(%2d,%2d):" % (i_t, j_t)
            print Xlist
        #
        
        i_prv = i_t; j_prv = j_t # prv = previous
        cnStemSegs = [] # proposed connected stem segments (list of lists)  
        stem_n = [Xlist[k_init]]
        kr = k_init
        flag_cont = True
        gap = 2 # allowed gap in [bp]
        nx = len(Xlist)-1
        while flag_cont and kr < nx:
            kr += 1
            if Xlist[kr].v == 'p':
                # parallel stem!
                continue
            #
            
            i_nxt = Xlist[kr].i; j_nxt = Xlist[kr].j # nxt = next
            
            
            if debug_scanFor_apStem:
                print "kr = %2d, ij_t(%2d,%2d), ij_nxt(%2d,%2d)" \
                    % (kr, i_t, j_t, i_nxt, j_nxt)
                
                print "i_nxt(%2d) - i_prv(%2d) <= gap(%2d)" % (i_nxt, i_prv, gap)
                print "j_prv(%2d) - j_nxt(%2d) <= gap(%2d)" % (j_prv, j_nxt, gap)
            #
            
            
            if i_nxt - i_prv == 1 and j_prv - j_nxt == 1:
                stem_n += [Xlist[kr]]
                stem_len += 1
                gap = stem_len / 2 + 1 # <== integer!!!!
                if gap > self.max_bp_gap:
                    gap = self.max_bp_gap
                #
                
                if kr == nx:
                    
                    cnStemSegs += [deepcopy(stem_n)]
                    # close out all stem_n and save in cnStemSegs
                    stem_n = []      # reset stem_n 
                    flag_cont = False
                #
                
            else:
                
                cnStemSegs += [deepcopy(stem_n)] # close out stem_n and save in cnStemSegs
                stem_n = []      # reset stem_n 
                
                # now test if the next stem_n (if exists) is also
                # satisfied by the current criteria
                if i_nxt - i_prv <= gap and j_prv - j_nxt <= gap:
                    if i_t < i_nxt and j_nxt < j_t: 
                        stem_n += [Xlist[kr]]
                        stem_len += 1
                        gap = stem_len / 2 + 1 # <== integer!!!!
                        if gap > self.max_bp_gap:
                            gap = self.max_bp_gap
                        #
                        
                        if debug_scanFor_apStem:
                            print "gap  = ", gap
                            print "slen = ", stem_len
                        #
                        
                    #
                    
                else:
                    flag_cont = False
                #
                
            #
            
            i_prv = i_nxt; j_prv = j_nxt
        #
        
        
        # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
        # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
        # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
        """@
        
        Here we finalize the true effective length of the stem
        according to the rules used in RNAModules(). This ensures that
        the estimations are consistent between the two independent
        ways of tackling the problem (the predictions obtained from
        Calculate vs the reverse calculation using this thread
        engine).
        
        """
        stem = []
        stem_lenf = 0
        if len(cnStemSegs) > 1:
            """@
            
            Whereas normally we might want to use the free energy to
            decide if we will end matters with the previous stem or
            combine the previous stem with this next stem (and any
            additional cnStemSegs that satisfied the rules), here we
            claim that this _is_ the linkage structure
            (fact). Therefore, we simply say, "you asked for it, you
            got it" in the later parts where the free energy is
            assigned.
            
            """
            
            n = len(cnStemSegs)
            stem_p = deepcopy(cnStemSegs[0])
            if debug_scanFor_apStem:
                print "stem_p", stem_p
            #
            
            n_p = len(stem_p) - 1
            ph_p = stem_p[n_p].i; qh_p = stem_p[n_p].j;
            slen_p = len(stem_p)
            
            stem_lenf += slen_p
            stem += stem_p 
            for k in range(1, n):
                stem_n = deepcopy(cnStemSegs[k])
                if debug_scanFor_apStem:
                    print "stem_n", stem_n
                #
                
                slen_n = len(stem_n)
                pt_n = stem_n[0].i; qt_n = stem_n[0].j;
                
                is_connect = self.fe.is_connected_aaStem(slen_p, ph_p, qh_p,
                                                         slen_n, pt_n, qt_n)
                
                if is_connect:
                    stem += stem_n
                    stem_lenf += slen_n
                    
                    if k < n - 1:
                        stem_p = deepcopy(cnStemSegs[k])
                        if debug_scanFor_apStem:
                            print "stem_p", stem_p
                        #
                        
                        slen_p = len(stem_p)
                        n_p  = slen_p - 1
                        ph_p = stem_p[n_p].i; qh_p = stem_p[n_p].j;
                    else:
                        break
                    #
                    
                else:
                    break
                #
                
            #
            
            # print "stem_lenf:", stem_lenf
            stem_len = stem_lenf
        elif len(cnStemSegs) > 0:
            stem = deepcopy(cnStemSegs[0])
        else:
            stem = [Xlist[k_init]]
        #
        
        # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        
        
        if debug_scanFor_apStem:
            print "stem_len:     ", stem_len
            print "gap:          ", gap
            print "cnStemSegs:   ", cnStemSegs
            print "stem (final): ", stem
            #self.fe.stopWhenMatchFound(i_t, j_t, 17, 24, "scanFor_apStem")
            
            #print "planned exit, end of scanFor_apStem"; sys.exit(0)
        #
        
        return stem
    #endMethod
    
    
    def build_ppStem(self, tailtags, Xlist):
        """
        This is a method for building parallel stems, stems that follow
        the rule
        
        stem = [(i, j), (i+1, j+1), (i+2, j+2) .... ]
        
        in short, (i,j) -> (i + ki, j + kj), where ki and kj > 0 and
        generally, ki = kj = 1.
        
        
        The important aspect of this part of the code is that we
        define the criteria (typically rather simple) for what
        constitutes a "parallel stem". This comes from my days with
        RNA where I came to understand that a stem need not be a
        contiguous segment. For example, the following parallel stem
        is surely contiguous.
        
        ....ABCD......abcd....
        
        However, what can we say about the following stem?
        
        ...A..B..C..D......a...b...c...d...
        
        It would surely depend on the type of interactions, which
        should be module sepecific. Nevertheless, perhaps my intuition
        would be that it is weakly bound and perhaps the pairs Aa, Bb,
        Cc and Dd are largely independent.
        
        How about this example?
        
        ...A.B.C.D......a.b.c.d...
        
        Its regularity and closeness would tend to pursuade me that it
        is effectively a single stem, not four separated and isolated
        stems. 
        
        In the very first external call to this general program, the
        scan of the sequence only identifies stems with a contiguous
        set of base pairs. In this method, some simple rules for
        parsing for effective stems is introduced to help resolve
        these issues.
        """
        
        
        debug_build_ppStem = False # True #
        if debug_build_ppStem:
            print "Enter build_ppStem: "
        #
        
        
        stemlist = []
        for sst in tailtags:
            for k in range(0, len(Xlist)):
                ib = Xlist[k].i; jb = Xlist[k].j
                if debug_build_ppStem:
                    print "sst = ", sst, "ijb = ", (ib, jb)
                #
                
                if ib == sst.i and jb == sst.j:
                    stem = self.scanFor_ppStem(Xlist, k)
                    # print stem, len(stem)
                    # sys.exit(0)
                    stemlist += [stem]
                    break
                #
                
            #endfor
            
        #endfor
        
        if debug_build_ppStem:
            for k in range(0, len(stemlist)):
                print k, stemlist[k]
            #
            
        #
        
        kr = 0
        while kr < (len(stemlist)):
            if debug_build_ppStem:
                print stemlist[kr]
            #
            
            tail = 0
            i_tr = stemlist[kr][tail].i; j_tr = stemlist[kr][tail].j
            head = len(stemlist[kr])-1
            i_hr = stemlist[kr][head].i; j_hr = stemlist[kr][head].j
            if debug_build_ppStem:
                print "stem{(%2d,%2d):(%2d,%2d)}" % (i_tr, j_tr, i_hr, j_hr)
            #
            
            kt = kr + 1
            while kt < len(stemlist):
                i_tt = stemlist[kt][tail].i; j_tt = stemlist[kt][tail].j
                if debug_build_ppStem:
                    print "build_ppStem, test (k=%2d)(%2d,%2d)" % (kt, i_tt, j_tt)
                #
                
                if (i_tr <= i_tt and
                    i_tt <= i_hr and 
                    j_hr <= j_tt and
                    i_tt <= j_tr):
                    if debug_build_ppStem:
                        print "del %s" % (stemlist[kt])
                    #
                    del stemlist[kt]
                else:
                    kt += 1
                #
                
            #endwhile
            
            if debug_build_ppStem:
                print "kr(%d) -> %d" % (kr, kr + 1)
            #
            
            kr += 1
        #endwhile
        
        return stemlist
    #endMethod
    
    
    def scanFor_ppStem(self, Xlist, k):
        """
        This is a support method for build_ppStem(). Parallel stems have a
        property that follows the rule
        
        stem = [(i, j), (i+1, j+1), (i+2, j+2) .... ]
        
        in short, (i,j) -> (i + ki, j + kj), where ki and kj > 0 and
        generally, ki = kj = 1.
        """
        
        debug_scanFor_ppStem = False # True #
        
        max_gap = 8
        stem_len = 1
        i_t = Xlist[k].i; j_t = Xlist[k].j
        i_p = i_t; j_p = j_t
        stem = [Xlist[k]]
        kr = k
        flag_cont = True
        gap = 2
        while flag_cont and kr < (len(Xlist)-1):
            kr += 1
            i_nx = Xlist[kr].i; j_nx = Xlist[kr].j
            if debug_scanFor_ppStem:
                print "kr = ", kr, "ij_t = ", (i_t, j_t), ", ij_nx = ", (i_nx, j_nx)
            #
            
            if i_nx - i_p <= gap and j_nx - j_p <= gap:
                if i_t <= i_nx and j_t <= j_nx:
                    stem += [Xlist[kr]]
                    stem_len += 1
                    gap = stem_len / 2 + 1
                    if gap > max_gap:
                        gap = max_gap
                    #
                    
                    if debug_scanFor_ppStem:
                        print "gap  = ", gap
                        print "slen = ", stem_len
                    #
                    
                #
                
                i_p = i_nx; j_p = j_nx
                
            else:
                flag_cont = False
            #
            
        #endwhile
        
        if debug_scanFor_ppStem:
            print "pp stem_len: ", stem_len
            print "pp gap:      ", gap
            print "pp stem:     ", stem
        #
        
        return stem
    #endMethod
    
    
    
    # ###################################^^^^^^^^^^^^^^^^^^^
    # ####  Stem construction tools  ####^^^^^^^^^^^^^^^^^^^
    # ###################################^^^^^^^^^^^^^^^^^^^
    
    
    
    
    # Used to process stem tail junctures that are stored in
    # PKlist. Looks for all the stem tails in PKlist. Presently can
    # only work with antiparallel stems
    def findStemsIn_ppList(self, ppPKtails):
        debug_findStemsIn_ppList  = False # True # 
        debug_findStemsIn_ppListx = False # True # 
        new_way = True
        if debug_findStemsIn_ppList:
            print "Enter findStemsIn_ppList: "
            print "ppPKtails: (before)"
            for rpkk in ppPKtails:
                print rpkk
            #endfor
            
            print "BPlist:"
            for vv in self.vsBPlist:
                print vv.disp_Pair()
            #endfor
            
            print "PKlist:"
            for vv in self.vsPKlist:
                print vv.disp_Pair()
            #endfor
        #
        
        pkppTailList = copyList(ppPKtails)
        pkppTailList = sortPairListWRT_n(pkppTailList, 0)
        ssppTailList = tuple2List(copyList(self.ppBPtails))
        ssppTailList = sortPairListWRT_n(ssppTailList, 0)
        if debug_findStemsIn_ppList:
            print "before doing this additional sorting into ss and pk stems:"
            print "0: ssppTailList: ", ssppTailList
            print "0: pkppTailList: ", pkppTailList
            #print "stop at 0 in findStemsIn_ppList"; sys.exit(0)
        #
        
        
        # ####################################################
        # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
        
        """
        I am not sure that any of this will really happen in folding a
        structure, but whereas you can build a list where normally the
        root can be defined on the same range as the anti-parallel
        secondary stucture, occasionally, this is not so. In this
        first step, I attempt to remove parallel PK roots formed off
        of an anti-parallel region.
        
        Specificially, the following structure leads to troubles if I
        do not carry out this step:
        
                  0         10        20        30        40        50        60    
                  |    .    |    .    |    .    |    .    |    .    |    .    |    .
        ss_seq = ".((((((.A.AAAAA.<BC..))))).((((.>bc....a..aaaaa....))))..)...."
        
        The positions [(16,32),(17,33),(18,34)] represent a parallel
        stem somehow associating with this structure. Though rather
        contorted and unlikely to form by any imaginable process (it's
        even hard to visualize), it is representable in this structure
        diagram. Unfortunately, ssppTailList will have (16,32) in the
        list and pkppTailList will have (17,33). Whereas I generally
        prefer to represent typical parallel stems as a kind of PK
        with a root, that is not possible here, so the final list of
        pkppTailList should only have (16,32) without (17,33).
        
        This section essentially addresses this type of matter, though
        hopefully such crazy things are highly unlikely from a
        practical point of view..
        
        """
        shiftinfo = []
        for kq in range(0, len(self.ssaptails)):
            iaprt = self.ssaptails[kq].i; japrt = self.ssaptails[kq].j
            kr = 0
            while kr < len(ssppTailList):
                ipprt = ssppTailList[kr][0]; jpprt = ssppTailList[kr][1]
                if debug_findStemsIn_ppListx:
                    print "iaprt(%2d) < ipprt(%2d) < japrt(%2d) < jpprt(%2d)" \
                        % (iaprt, ipprt, japrt, jpprt)
                #
                
                if iaprt < ipprt and japrt < jpprt and ipprt < japrt:
                    if debug_findStemsIn_ppListx:
                        print "rearranging to pklist", ssppTailList[kr]
                    #
                    
                    pkppTailList += [(ipprt, jpprt)]
                    shiftinfo += [(ipprt, jpprt)]
                    del ssppTailList[kr]
                    
                else:
                    kr += 1
                #
                
            #endwhile
            
        #endfor
        
        pkppTailList = sortPairListWRT_n(pkppTailList, 0)
        if debug_findStemsIn_ppListx:
            print "pkppTailList(initial): ", pkppTailList
            print "shift info:            ", shiftinfo
        #
        
        if len(shiftinfo) > 0 and new_way:
            for kp in range(0, len(shiftinfo)):
                isf = shiftinfo[kp][0]; jsf = shiftinfo[kp][1]
                kt = 0
                while kt < len(self.vsBPlist):
                    it = self.vsBPlist[kt].i;
                    jt = self.vsBPlist[kt].j
                    ijv = self.vsBPlist[kt].v
                    
                    if (isf == it and jt == jsf):
                        if debug_findStemsIn_ppListx:
                            print "del matching BP(%d,%d) in PK list" % (it, jt)
                        #
                        
                        self.vsPKlist += [self.vsBPlist[kt]]
                        del self.vsBPlist[kt]
                        
                    elif (isf < it and jsf < jt and it < jsf and ijv == 'p'):
                        # ordered list!!!
                        if debug_findStemsIn_ppListx:
                            print "move from BP(%d,%d) to PK list" % (it, jt)
                        #
                        
                        self.vsPKlist += [self.vsBPlist[kt]]
                        del self.vsBPlist[kt]
                        
                    else:
                        kt += 1
                    #
                    
                #endwhile
                
            #endfor
            
            for kr in range(0, len(self.vsPKlist)):
                it = self.vsPKlist[kr].i; jt = self.vsPKlist[kr].j
                ijv= self.vsPKlist[kr].v
                kt = 0
                while kt < len(self.vsBProots):
                    irt = self.vsBProots[kt].i; jrt = self.vsBProots[kt].j
                    
                    if debug_findStemsIn_ppListx:
                        print "ijrt = (%2d,%2d)" % (irt, jrt)
                    #
                    
                    if (irt == it and jrt == jt and ijv == 'p'):
                        # ordered list!!!
                        if debug_findStemsIn_ppListx:
                            print "delete BProots(%d,%d)" % (it, jt)
                        #
                        
                        del self.vsBProots[kt]
                        
                    else:
                        kt += 1
                    #
                    
                #endwhile
                
            #endfor
            
        #
        
        if debug_findStemsIn_ppListx:
            print "updated full lists"
            print "BPlist:"
            for vv in self.vsBPlist:
                print vv.disp_Pair()
            #endfor
            
            print "BProots:"
            for vv in self.vsBProots:
                print vv.disp_Pair()
            #endfor
            
            print "PKlist:"
            for vv in self.vsPKlist:
                print vv.disp_Pair()
            #endfor
            
            print "PKroots:"
            for vv in self.vsPKroots:
                print vv.disp_Pair()
            #endfor
            
            #print "stop at 1 in findStemsIn_ppList"; sys.exit(0)
        #
        
        # Now we address potential overlaps with closely neighoring
        # parallel pairs that are present in the pk list.
        
        # --> keep only the stemtails and remove all the extraneous
        # --> stuff.
        
        kr = 0
        while kr < (len(pkppTailList)-1):
            ikr0 = pkppTailList[kr  ][0]; jkr0 = pkppTailList[kr  ][1]
            ikr1 = pkppTailList[kr+1][0]; jkr1 = pkppTailList[kr+1][1]
            if debug_findStemsIn_ppListx:
                print "ikr0(%2d) < ikr1(%2d) < jkr0(%2d) < jkr1(%2d)" \
                    % (ikr0, ikr1, jkr0, jkr1)
            #
            
            if ikr0 < ikr1 and jkr0 < jkr1 and ikr1 < jkr0:
                if (ikr1 - ikr0) == 1 and (jkr1 - jkr0) == 1:
                    if debug_findStemsIn_ppListx:
                        print "removing from pklist", pkppTailList[kr+1]
                    #
                    
                    del pkppTailList[kr+1]
                    
                else:
                    kr += 1
                #
                
            else:
                kr += 1
            #
            
        #endwhile
        
        # 190328: had to make some major fixes here to resolve a
        # problem with the counter in TreeNode2DotBracket where the
        # same stem was appearing multiple times. I think these hacks
        # fix this, but it seems that issues may remain with parallel
        # stems. At least it remans a possibility.
        
        if debug_findStemsIn_ppList:
            print "pkppTailList(final): ", pkppTailList
            print "PKlist:"
            for vv in self.vsPKlist:
                print vv.disp_Pair()
            #
            
            print "======="
        #
        
        if not new_way:
            self.vsPKlist += matchTupleList2VsList(shiftinfo, self.vsPKlist)
            self.vsPKlist += matchTupleList2VsList(shiftinfo, self.vsBPlist)
        #
        
        self.vsPKlist = self.sortvsList(self.vsPKlist, 'i')
        
        # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        # ####################################################
        
        
        
        # build a vslist (sstails) from self.vsBPlist
        self.sspptails =  matchTupleList2VsList(ssppTailList, self.vsBPlist)
        # build a vslist (pktails) from self.vsPKlist and self.vsBPlist
        self.pkpptails =  matchTupleList2VsList(pkppTailList, self.vsPKlist)
        self.pkpptails += matchTupleList2VsList(pkppTailList, self.vsBPlist)
        # second component (self.vsBPlist and self.vsPKlist) are used as
        # reference to pkppTailList and written to pktails.
        self.pkpptails = self.sortvsList(self.pkpptails, 'i')
        
        if new_way:
            
            kr = 0
            while kr < (len(self.pkpptails)-1):
                ippr = self.pkpptails[kr].i; jppr = self.pkpptails[kr].j;
                
                if debug_findStemsIn_ppListx:
                    print "pkpptails[kr= %2d]: %s" % (kr, self.pkpptails[kr])
                #
                
                kt = kr+1
                while kt < len(self.pkpptails):
                    ippt = self.pkpptails[kt].i; jppt = self.pkpptails[kt].j;
                    
                    if debug_findStemsIn_ppListx:
                        print "pkpptails[kt= %2d]: %s" % (kt, self.pkpptails[kt])
                    #
                    
                    if ippr == ippt and jppr == jppt:
                        del self.pkpptails[kt]
                    else:
                        kt += 1
                    #
                    
                #endwhile
                
                kr += 1
            #endwhile
            
        #
        
        
        
        if debug_findStemsIn_ppList:
            print "======="
            print "setup finished"
            # print "pkppTailList: ", pkppTailList
            # print "ssppTailList: ", ssppTailList
            print "pp sstails: "
            for vv in self.sspptails:
                print vv.disp_Pair()
            #endfor
            
            print "pp pktails"
            for vv in self.pkpptails:
                print vv.disp_Pair()
            #endfor
            
            self.disp_ViennaLists("ViennaList: findStemsIn_ppList")
            print "-------"
            #print "stop at 2 in findStemsIn_ppList"; sys.exit(0)
        #
        
        
        
        # build stems based on the current information:
        ss_ppstemlist  = self.build_ppStem(self.sspptails, self.vsBPlist)
        ss_ppstemlist  = self.ins_sort_stemlist(ss_ppstemlist)
        
        # print "xx"
        if debug_findStemsIn_ppList:
            for vv in self.vsPKlist:
                print vv.disp_Pair()
            #endfor
            
            print self.pkpptails
        #
        
        pk_ppstemlist  = self.build_ppStem(self.pkpptails, self.vsBPlist)
        pk_ppstemlist += self.build_ppStem(self.pkpptails, self.vsPKlist)
        pk_ppstemlist  = self.ins_sort_stemlist(pk_ppstemlist)
        
        if debug_findStemsIn_ppList:
            
            print "pp ssstem (ss_ppstemlist):"
            for k in range(0, len(ss_ppstemlist)):
                print "%2d: %s" % (k, self.disp_stem(ss_ppstemlist[k]))
            #endfor
            
            print "pp pkstem (pk_ppstemlist): "
            for k in range(0, len(pk_ppstemlist)):
                print "%2d: %s" % (k, self.disp_stem(pk_ppstemlist[k]))
            #endfor
            
            #print "stop at 3 (end of) findStemsIn_ppList"; sys.exit(0)
        #
        
        return ss_ppstemlist, pk_ppstemlist
        ##################################
    #endMethod
    
    
    
    
    
    def build_MPstruct(self):
        """
        This step takes the island (MultiPair) data that has been dug out
        from the input structure
        """
        
        debug_build_MPstruct = False # True # 
        if debug_build_MPstruct:
           print "Enter build_MPstruct"
        #
        
        if not self.initialized_bplist:
            print "ERROR: antiparallel secondary structure or PK data not initialized"
            sys.exit(1)
        #
        
        bpTailList = []
        for vv in self.vsBProots:
            bpTailList += [vv]
        #endfor
        
        #print "bpTailList: ", bpTailList
        # sys.exit(0)
        
        """
        Example:
        
        .{...|......}...{..((...))...((...))..|...((...))...}...
        |    '    |    '    |    '    |    '    |    '    |    '
        0         10        20        30        40        50
        
        CTCF                                 branches 
        {  1, 12}[(  1,  5), (  5, 12)]      None
        { 16, 52}[( 16, 38), ( 38, 52)]      [(19,25), (29,35), (42,48)]
        
        output:
        arches: 
        ctcf(  1, 12): 
        island(  1,  5)[-]: branches[-] --> None
        island(  5, 12)[-]: branches[-] --> None
        ctcf( 16, 52): 
        island( 16, 38)[-]: branches[P] --> ( 19, 25), ( 29, 35), 
        island( 38, 52)[-]: branches[S] --> ( 42, 48), 
        
        """
        
        if debug_build_MPstruct:
            print "vsMPlist: ", self.vsMPlist
        #
        
        cldict = {}
        self.MPlist = []
        for kcl in range(0, len(self.vsMPlist)):
            
            cl = self.vsMPlist[kcl]
            MPisland = self.vs.expand_island(cl)
            ww = MultiPair(cl.i, cl.j, cl.v)
            if debug_build_MPstruct:
                print "cl:       ", cl
                print "MPisland: ", MPisland
                print "ww:       ", ww.disp_MultiPair()
                #sys.exit(0)
            #
            
            if len(MPisland) == 1:
                cik = MPisland[0]
                if debug_build_MPstruct:
                    print "1 MP island:", len(MPisland)
                    print "cik = ", cik, cik.i, cik.j
                #
                
                
                mbllist = []
                for bpk in bpTailList:
                    if cik.i < bpk.i and  bpk.j < cik.j:
                        mbllist += [bpk]
                    #
                    
                #endfor
                if debug_build_MPstruct:
                    print "mbllist: ", mbllist
                #
                
                a = ArchMP(cik.i, cik.j, cik.v)
                itype = ''
                if len(mbllist) == 1:
                    itype = 'J'
                elif len(mbllist) > 1:
                    itype = 'P'
                else:
                    itype = '-'
                #
                
                cldict[(cik.i,cik.j)] = itype
                
                a.internal = mbllist
                a.btype    = itype
                ww.arches  = [a]
                
            else:
                if debug_build_MPstruct:
                    print "n MP islands:", len(MPisland)
                #
                
                a = ArchMP(cl.i, cl.j, cl.v)
                a.internal = []
                a.btype    = 'X'
                ww.arches  = [a]
                
                for kclj in range(1, len(MPisland)):
                    cik = MPisland[kclj]
                    if debug_build_MPstruct:
                        print "assign pair: (%2d,%2d) " % (cik.i, cik.j)
                    #
                    
                    
                    mbllist = []
                    for bpk in bpTailList:
                        if cik.i < bpk.i and  bpk.j < cik.j:
                            mbllist += [bpk]
                        #
                        
                    #endfor
                    
                    if debug_build_MPstruct:
                        print "mbllist: ", mbllist
                    #
                    
                    
                    a = ArchMP(cik.i, cik.j, cik.v)
                    itype = ''
                    if len(mbllist) == 1:
                        itype = 'J'
                    elif len(mbllist) > 1:
                        itype = 'P'
                    else:
                        itype = '-'
                    #
                    
                    cldict[(cik.i,cik.j)] = itype
                    a.internal = mbllist
                    a.btype    = itype
                    ww.arches += [a]
                #endfor
                
            #
            
            self.MPlist += [ww]
            
        #
        
        
        
        if debug_build_MPstruct:
            if 0: # produces the same output as the next operation
                print "islands: "
                for ww in self.MPlist:
                    print disp_Motif(ww)
                #endfor
                
            #
            
            print "islands: "
            for ww in self.MPlist:
                print ww.disp_MultiPair()
            #endfor
            
            print "len(self.vsMPlist): ", len(self.vsMPlist)
            for k in range(0, len(self.vsMPlist)):
                print self.vsMPlist[k].disp_Pair()
            #endfor
            
        #
    #endMethod
    
#



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
    # the argument "test0" is just used to give this data set a name.
    chdt.print_ChPairData()
    # no argument in print_ChPairData() means that the output will
    # only be displayed.
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
#

def test0(cl):
    aList = [123, 'xyz', 'zara', 'abc']
    aList.insert( 3, 2009)
    print "Final List : ", aList
#


def test1(cl):
    vs = Vstruct()
    #          0         10        20        30        40        50        60        70        80        90
    #          |         |         |         |         |         |         |         |         |         |
    ss_seq  = ".(((.(((.(((....))).))).)))."
    #ss_seq = ".(((.(((((((.[[...)))))..]].((((...))))..))))).....([)]([)].......(..)..((..)).(..)..........................."
    #ss_seq = "{(((.(((((((.[[...)))))..]].((((...))))..)))))..|..([)]([)]...}.{.(..)..((..)).(..)..}...{.....|...|....|....}"
    #ss_seq = "{((((((.A.AAAAA......))))).BBBB........a..aaaaa....bbbb..).|..([)]([)].ABC.abc.}....{.....}"
    #ss_seq = "{((((((.A.AAAAA......))))).((((........a..aaaaa....))))..).|..([)]([)].ABC.abc.}....{.....}"
    #ss_seq = "{((((((.A.AAAAA......))))).((((........a..aaaaa....))))..).|..([)]([)]...}....{.A.B.C..a.b.c.....}"
    #ss_seq = "{((((((.A.AAAAA......))))).((((........a..aaaaa....))))..).|..([)]([)]...}....{.ABC..DEF..abc..def...}"
    #ss_seq = "{((((((.A.AAAAA.<BC..))))).((((.>bc....a..aaaaa....))))..).|..([)]([)]...}....{.A.B.C..a.b.c.....}"
    #ss_seq = "{((((((.A.AAAAA.<BC..))))).((((.>bc.DE.a..aaaaa..de))))..).|..([)]([)].........}....{.....}"
    #ss_seq  = "{((((((.A.AAAAA.<BC..))))).((((.>bc....a..aaaaa....))))..).|..([)]([)]...}....{.ABC..DEF..abc..def...}"
    
    if len(cl) > 1:
        ss_seq = cl[1]
    vs.parse_fullDotBracketStructure(ss_seq, True)
    v2t = Vienna2TreeNode(vs)
    v2t.vienna2tree()
    print "planned exit"; sys.exit(0);
    
    tf = LThreadBuilder(v2t)
    tf.visit(v2t.genTree)
    print "LThread notation: "
    tf.disp_lt()
    # print "planned exit"; sys.exit(0);
    
    # print dt.makeLThreadDotBracket_VARNA(ds.lt, 0)
    vf = LThread2Vienna()
    vf.lt2vs(tf.lt)
    vf.set_vstr(v2t.vstr)
    vf.set_vseq(v2t.vseq)
    
    print "structure sequence: "
    print vf.vstr
    print "vsBPlist: "
    for bpk in vf.vsBPlist:
        print bpk.disp_Pair()
    #
    print "vsPKlist: "
    for bpk in vf.vsPKlist:
        print bpk.disp_Pair()
    #
    print "vsMPlist: "
    for bpk in vf.vsMPlist:
        print bpk.disp_Pair()
    #
    
#    

def usage_test2():
    print "Usage: %s [structure_sequence]" % PROGRAM
    print "       %s [sequence    structure_sequence]" % PROGRAM
#
def test2(cl):
    
    # Evidently, Vienna() still has a few problems presently because
    # it cannot convert ".ABCD.......abcd...." properly.
    
    rnaseq = ""
    # regular secondary structure
    #          0         10        20        30        40        50        60        70
    #          |    .    |    .    |    .    |    .    |    .    |    .    |    .    |
    #ss_seq  = "(..............)"
    #rnaseq  = "GuuuuuuuuuuuuuuC"
    #ss_seq  = ".(..............).."
    #rnaseq  = "uGuuuuuuuuuuuuuuCuu"
    #ss_seq  = ".((............)).."
    #rnaseq  = "uGGuuuuuuuuuuuuCCuu"
    #ss_seq  = ".((((........)))).."
    #rnaseq  = "uGGGGuuuuuuuuCCCCuu"
    #ss_seq  = ".((((..(...).)))).."
    #rnaseq  = "uGGGGuuGuuuCuCCCCuu"
    #ss_seq   = ".(((....)))."
    #rnaseq   = "uGGGuuuuCCCu"
    #ss_seq   = ".(((.(((.(((....))).))).)))."
    #rnaseq   = "uGGGuGGGuGGGuuuuCCCuCCCuCCCu"
    #ss_seq  = ".((((..(((..........))).)))).."
    #rnaseq  = "uGGGGuuGGGuuuuuuuuuuCCCuCCCCuu"
    #ss_seq   = ".(((.(((......(((....))).......))).)))."
    #rnaseq   = "uGGGuGGGuuuuuuGGGuuuuCCCuuuuuuuCCCuCCCu"
    ss_seq  = ".((((........))))...((((........)))).."
    rnaseq  = "uGGGGuuuuuuuuCCCCuuuGGGGuuuuuuuuCCCCuu"
    #ss_seq  = "((((........))))...((((........))))"
    #rnaseq  = "GGGGuuuuuuuuCCCCuuuGGGGuuuuuuuuCCCC"
    #ss_seq  = "(((.((((........))))...((((........)))).)))"
    #rnaseq  = "GGGuGGGGuuuuuuuuCCCCuuuGGGGuuuuuuuuCCCCuCCC"
    #ss_seq  = ".((((..(((..(...)..))).)))).."
    #rnaseq  = "uGGGGuuGGGuuGuuuCuuCCCuCCCCuu"
    #          0         10        20        30        40        50        60        70
    #          |    .    |    .    |    .    |    .    |    .    |    .    |    .    |
    #ss_seq  = "...(((..(..((((..(...).))))..)...))).."
    #rnaseq  = "uuuGGGuuGuuGGGGuuGuuuCuCCCCuuCuuuCCCuu"
    #ss_seq  = "...(((..((((........))))..(((.........)))..))).."
    #rnaseq  = "uuuGGGuuGGGGuuuuuuuuCCCCuuGGGuuuuuuuuuCCCuuCCCuu"
    #ss_seq  = "...(((..((((..(...).))))..(((.(.....).)))..))).."
    #rnaseq  = "uuuGGGuuGGGGuuGuuuCuCCCCuuGGGuGuuuuuCuCCCuuCCCuu"
    #ss_seq  = ".(.(.(.(..).).).)..(.(.........).).."
    #rnaseq  = "uGuGuGuGuuCuCuCuCuuGuGuuuuuuuuuCuCuu"
    #ss_seq  = ".(.(.(.(.....).).).)..(.(.........).).."
    #rnaseq  = "uGuGuGuGuuuuuCuCuCuCuuGuGuuuuuuuuuCuCuu"
    #ss_seq  = ".((((........))))..(((.........))).."
    #rnaseq  = "uGGGGuuuuuuuuCCCCuuGGGuuuuuuuuuCCCuu"
    #ss_seq  = ".((((..(...).))))..(((.(.....).))).."
    #rnaseq  = "uGGGGuuGuuuCuCCCCuuGGGuGuuuuuCuCCCuu"
    #ss_seq  = "..(((((...(((..((((..(...).))))..(((.(.....).)))..)))..)))))."
    #rnaseq  = "uuGGGGGuuuGGGuuGGGGuuGuuuCuCCCCuuGGGuGuuuuuCuCCCuuCCCuuCCCCCu"
    
    #ss_seq  = ".((((...[[[[.))))..]]]].."
    #rnaseq  = "uGGGGuuuGGGGuCCCCuuCCCCuu"
    #ss_seq  = ".((((...AAAA.......BBBB..))))..aaaa.......bbbb"
    #rnaseq  = "uGGGGuuuGGGGuuuuuuuGGGGuuCCCCuuCCCCuuuuuuuCCCC"
    #ss_seq  = ".((((...AAAA.......BBBB..))))..bbbb.......aaaa"
    #rnaseq  = "uGGGGuuuGGGGuuuuuuuGGGGuuCCCCuuCCCCuuuuuuuCCCC"
    
    
    #          0         10        20        30        40        50        60        70        80        90
    #          |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |
    #ss_seq  = "..(((((((..(((((...(((..((((..(...).))))..(((.(.....).)))..)))..))))).(((....)))..)))))))..((...)).."
    #ss_seq  = "..(((((((..(((((..(((...((((..(...).))))..(((.(.....).)))..)))..))))).(((....)))..)))))))..((...)).."    
    
    # parallel stems and mixtures
    #          0         10        20        30        40        50        60        70        80        90
    #          |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |
    #ss_seq  = ".ABCD...........abcd.."
    #ss_seq  = "..A..B..C..D..a..b..c..d.."
    #ss_seq  = "..A..B..C..D..a..b..c..d....E..F..G....e...f...g..."
    #ss_seq  = ".ABCD..EFG......abcd......efg...."
    #@@@
    #ss_seq  = ".ABCD..EFG..H...abcd......efg..h."
    #@@@
    #ss_seq  = ".ABCD.......abcd....EFG...efg."
    #ss_seq  = ".ABCD.......abcd....((...))."
    #ss_seq  = "..((...))...ABCD......abcd.."
    #ss_seq  = ".ABCD...EF....abcd..ef..........."
    #ss_seq  = ".ABCD...EF....abcd..ef...((...))."
    #ss_seq  = ".ABCD....EFG...efg.....HIJ...hij....abcd.."
    #ss_seq  = ".ABCD....EFG...efg..((.....)).....HIJ...hij....abcd.."
    #@@@
    #ss_seq  = ".ABCD....EFG...efg..((..[[.)).]]..HIJ...hij....abcd.."
    #@@@
    #ss_seq  = ".AB..CD....EFG...efg.....HI..JK...hi..jk....ab..cd.."
    
    #          0         10        20        30        40        50        60        70        80        90
    #          |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    
    #ss_seq  = "ABCD...((((...abcd...))))"
    #ss_seq  = "((((...ABCD...))))...abcd"
    #ss_seq  = "AB...((...ab...))"
    #ss_seq  = "(...)"
    #ss_seq  = "AB...ab"
    #ss_seq   = "([...)]"
    #ss_seq  = "(..AB...ab..)"
    #ss_seq  = "(..(B...)b..)"
    #ss_seq  = "(.....)..AB...ab.."
    #ss_seq  = "((...))..AB...ab.."
    #ss_seq  = "AB...ab..(.....).."
    #ss_seq  = "A..B...C...D....a...b...c...d."
    #ss_seq  = "(..B...C...D....)...b...c...d."
    #ss_seq  = ".A..B..C..D...........a..b..c..d.."
    #ss_seq  = "ABC....abc"
    #ss_seq  = "(((((...ABCD....[[[[...)))))...]]]]...abcd.."
    #ss_seq  = "(((((...AB.....ab...[[...)))))...]]...."
    #ss_seq  = "AB....(((((...ab...[[...)))))...]]...."
    #ss_seq  = "[[....(((((...]]...AB...)))))...ab...."
    
    #          0         10        20        30        40        50        60        70        80        90
    #          |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .
    #ss_seq  = ".AAA...BBB...CCC....aaa...bbb....ccc...."
    #ss_seq  = ".AAA...BBB...aaa...CCC....bbb....ccc...."
    #ss_seq  = ".(((...BBB...)))...(((....bbb....)))...."
    
    # with PK contacts
    #          0         10        20        30        40        50        60        70        80        90
    #          |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |
    #ss_seq  = "..([)].."
    #ss_seq  = "..(.[)].."
    #ss_seq  = "..(.[.)].."
    #ss_seq  = "..(.[.).].."
    #ss_seq  = "..([)]..([)].."
    #ss_seq  = "(.).(.([)])(.)(.)(())(..)(.).(())(..).((.))"
    #@@@
    #ss_seq  = "(.).(.([)])"
    #@@@
    #ss_seq  = "(.).(.ABab)"
    #ss_seq  = "(....(.[..)..]...)..."
    #ss_seq  = "((...(.[..)..]..))..."
    #ss_seq  = "((...(.[<.)..].>))..."
    #ss_seq  = "((...A.BC.a..b.c))..."
    
    #ss_seq  = "((...A..(.).B..(.)..C..(.)....a...b...c..))"
    #ss_seq  = "((...(..(.).B..(.)..C..(.)....)...b...c..))"
    
    
    #          0         10        20        30        40        50        60        70        80        90
    #          |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |
    #ss_seq = ".((....[[......))..]].."
    #ss_seq = ".((....AA......))..((.....BB....))...aa..bb.."
    #ss_seq = ".((....AA......))..((.....BB....))...bb..aa.."
    #ss_seq = ".((....AA..BB..))..((...........))...aa..bb.."
    #ss_seq = ".((....AA..BB..))..((...........))...bb..aa.."
    
    
    #ss_seq = ".((....AA...))..BB....aa....bb." # same
    #ss_seq = ".((....[[...))..((....]]....))." # same
    #ss_seq = ".((....AA...))..((....aa....CC...))....cc.." # structure 1 R/K-type
    #ss_seq = ".((....AA...))..BB....aa....CC...bb....cc.." # structure 1 K-type only
    #ss_seq = ".((....AA...))..((....CC....aa...))....cc.." # structure 1 R/K-type
    #ss_seq  = ".((....AA...))..BB....CC....aa...bb....cc.." # structure 1 K-type only
    
    #ss_seq  = ".((...((..AA.))..aa...((..BB.))..bb..))..((....CC.........))....cc.." # structure 1 R/K-type
    #                -------------   -------------      -------------------------
    
    #ss_seq = ".((....[[...))..((....]]....[[...))...((...]]....)).." # structure 2 R/K-type
    #ss_seq =  ".((....AA...))..BB....aa....CC...bb...DD...cc....dd.." # structure 2 K-type only
    
    
    #ss_seq = ".((....AA...))..((....aa....CC...))...((...cc....)).." # structure 2 R/K-type
    #ss_seq  = ".((....AA...))..BB....aa....CC...bb...DD...cc....dd.." # structure 2 K-type only
    
    #ss_seq = ".((....AA..BB..CC..))..((....aa..bb..cc...DD...EE..FF..)).....dd..ee..ff." # structure 2 R/K-type
    #@@@@@
    #ss_seq = ".((....AA..BB..CC..))..DD....aa..bb..cc...EE...FF..GG..dd.....ee..ff..gg." # structure 2 K-type only
    #ss_seq = ".......AA..BB..CC......((....aa..bb..cc...DD...EE..FF..))..GG...dd..ee..ff...gg.." # structure 2 K-type only
    #ss_seq = ".....(((....AAA....)))...(((....[[[...)))...(((...]]]....)))....aaa...."
    #@@@@@
    
    
    #          0         10        20        30        40        50        60        70        80        90
    #          |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |
    #ss_seq = ".((....AA..BB..CC..))..((....aa..bb..cc...DD...EE..FF..))...((..dd..ee..ff...)).." # structure 2 R/K-type
    #ss_seq = ".((....AA..BB..CC..))..DD....aa..bb..cc...EE...FF..GG..dd...HH..ee..ff..gg...hh.." # structure 2 K-type only
    #ss_seq = ".((....AA..BB..CC..))..DD....aa..bb..cc...EE...FF..GG..dd...((..ee..ff..gg...)).." # structure 2 K-type only
    
    
    
    #          0         10        20        30        40        50        60        70        80        90
    #          |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |
    #ss_seq = ".(((.(((..AA..)))...(((..BB...))).)))..(((.(((..aa..)))...(((..bb...))).)))"
    #ss_seq = ".(((.(((..AA..)))...(((..BB...))).)))..(((.(((..bb..)))...(((..aa...))).)))"
    #ss_seq = ".(((.(((..[[..)))...(((..]]...))).)))..(((.(((..[[..)))...(((..]]...))).)))"
    #ss_seq = ".(((.(((..[[..)))...]]..(((...))).)))..(((.(((..[[..)))...]].(((....))).)))"
    #ss_seq = ".((((....[[[...))))...((((.....))))...((...)).]]]..."
    #ss_seq  = "(((((.....((((....[[[...))))...((((.....))))...((...)).]]].......)))))...."
    #ss_seq = ".(((.(((......)))...(((.......))).)))..(((.(((......)))...(((.......))).)))"
    
    #ss_seq = ".(((.(((..AA..)))...(((..BB...))).)))..((....)).(((.(((..aa..)))...(((..bb...))).)))" # structure 3 R-type only
    #ss_seq = ".(((.(((..AA..)))...(((..CC...))).)))..((....)).EEE.BBB..aa..bbb...DDD..cc...ddd.eee" # structure 3 K-type only
    
    #ss_seq = ".(((.(((..AA..)))...(((..BB...))).)))..((....)).(((.(((..bb..)))...(((..aa...))).)))" # structure 4 R-type only
    #ss_seq = ".(((.(((..CC..)))...(((..AA...))).)))..((....)).EEE.BBB..aa..bbb...DDD..cc...ddd.eee" # structure 4 K-type only
    
    
    
    #          0         10        20        30        40        50        60        70        80        90
    #          |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |
    #ss_seq  = ".(((.(((..AA..)))...(((.......))).)))..((.BB.)).(((.(((..aa..)))...(((..bb...))).)))"
    #ss_seq  = ".(((.(((..AA..)))...(((.......))).)))..((.BB.)).(((.(((..bb..)))...(((..aa...))).)))"
    #ss_seq  = ".(((.(((..AB..)))...(((.......))).)))..((.CD.)).(((.(((..ab..)))...(((..cd...))).)))"
    
    #ss_seq = ".(((...[[.....)))..(((.]]....))).."
    #ss_seq = "...(((..((((..(...[[[..).))))..(((.(..]]].....).)))..))).."
    #ss_seq = ".((....[[......))..]]...((..[[..))..]]."
    
    #          0         10        20        30        40        50        60        70        80        90        100       110       120
    #          |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |
    #ss_seq = "...(((..((((..(...[[[..).))))..(((.(..]]].....).)))..)))...(((..((((..(...[[[..).))))..(((.(..]]]........).)))..))).."
    #ss_seq = "...(((..((((..(...[[[..).))))..(((.(..]]].....).)))..)))...(((..((((..(...AAA..).))))..(((.(..BBB..aaa...).)))..)))..bbb"
    #ss_seq = ".(((.(((((((.[[...)))))..]].((((...))))..))))).....([)]([)].......(..)..((..)).(..)..........................."
    #ss_seq = ".(((.(((((((.[[...)))))..((((.]]...))))..))))).....([)]([)].......(..)..((..)).(..)..........................."
    #ss_seq  = ".((((..((((...((((....))))..((((....))))...[[[[..))))...))))...((((((...]]]]...))))))..."
    #ss_seq  = ".((((..((((...((((....))))..((((...[[[[...))))...))))...))))...((((((...]]]]...))))))..."
    #ss_seq  = ".((((..((((...((((....))))..((((...[[[[...))))...))))...))))...((((((...((((..]]]]...))))...(((....)))...))))))..."
    #ss_seq  = ".((((..((((...((((....))))..((((...[[[[...))))...))))...))))...((((((.]]]]..((((.....))))...(((....)))...))))))..."
    
    # with CTCF contacts  
    #          0         10        20        30        40        50        60        70        80        90
    #          |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    
    #ss_seq = "{....................}"
    #ss_seq = "{.............|......}"
    #ss_seq = "{..((...))....|......}"
    #ss_seq = "{..((...))....|...((....))...}"
    #ss_seq = "{..((...))...((....))...|...((....))...}"
    #ss_seq = "{..((...))....|.....|...}"
    #ss_seq =  ".{...|......}...{..((...))...((...))..|...((...))...}..."
    #ss_seq =  ".{...|......}...{..((.[.))...((.].))..|...((...))...}..."
    #ss_seq = ".{.(...).|.((...))..}...{..((...))...((...))..|...((...))...}..."
    #ss_seq = ".{.((.)).|.((...))..}...{..((...))...((...))..|...((...))...}..."
    #ss_seq = ".{.(...).(....)..|.((...))..}...{..((...))...((...))..|...((...))...}..."
    #ss_seq = "{(((.(((((((.[[...)))))..]].((((...))))..)))))..|..([)]([)]...}.{.(..)..((..)).(..)..}...{.....|...|....|....}"
    #ss_seq  = "{((((((.A.AAAAA.<BC..))))).((((.>bc....a..aaaaa....))))..).|..([)]([)]...}....{.ABC..DEF..abc..def...}"
    #ss_seq  = ".((((((.A.AAAAA.<BC..))))).((((.>bc....a..aaaaa....))))..)."
    vs = Vstruct()
    vs.set_system("RNA")
    #vs.set_system("Chromatin")
    
    if len(cl) == 2:
        ss_seq = cl[1] # only the structure sequence given
        rnaseq = len(ss_seq)*'A' # no information, so just make it poly A
    elif len(cl) == 3:
        rnaseq = cl[1]
        ss_seq = cl[2]
    else:
        if len(cl) > 2:
            print "Too many arguments for this test!"
            print cl
            usage_test2()
            sys.exit(1)
        else:
            print "using defaults for test2"
        #
        
    #
    
    
    try:
        print "main: input structure sequence:"
        if not rnaseq == "": print rnaseq
        print ss_seq
        
    except(UnboundLocalError):
        print "ERROR, ss_seq is not assigned -- you idiot!"
        usage_test2()
        sys.exit(1)
    #
    vs.parse_fullDotBracketStructure(ss_seq, rnaseq, True) # True = print results
    #print "planned exit after running parse_fullDotBracketStructure"; sys.exit(0);
    
    print "next: calculate the tree structure:"
    v2t = Vienna2TreeNode(vs, None, rnaseq)
    v2t.vienna2tree()
    
    print v2t.dispTreeStruct(v2t.ssTree,  "Secondary Structure Tree")
    print v2t.dispTreeStruct(v2t.genTree, "Full Structure Tree")
    #print "planned exit before running TreeNode2Motif"; sys.exit(0);
    
    
    t2m = TreeNode2Motif(v2t)
    t2m.visit(v2t.genTree)
    # print v2t.MPlist
    t2m.post_graftMP()
    print "finished TreeNode2Motif"
    print "%10.2f" % t2m.dGmin
    print "planned exit before running LThreadBuilder"; sys.exit(0);
    
    
    # Test LThread building abilities
    print "LThreadBuilder test"
    ltb = LThreadBuilder(v2t)
    ltb.visit(v2t.genTree)
    print "LThread notation: "
    ltb.disp_lt()
    
    print "review structures in vs"
    vs.disp_allLists()
    print "planned exit after running LThreadBuilder"; sys.exit(0);
    
    
    # Test transformation between LThread and Vienna Pair
    print "now try generating LThread output"
    vf = LThread2Vienna()
    vf.lt2vs(ltb.lt)
    vf.set_vstr(v2t.vstr)
    vf.set_vseq(v2t.vseq)
    
    print "make new Vstruct using LThread2Vienna data (test conversion):"
    vsp = Vstruct()
    vsp.init_vs_from_lt(vf)
    print "review structures in vsp"
    vsp.disp_allLists()
    # sys.exit(0)
    
    
    print "analyze the bplist and make dot-bracket notation"
    v2tp = Vienna2TreeNode(vsp)
    v2tp.vienna2tree()
    tn2db = TreeNode2DotBracket(v2tp)
    tn2db.make_dotbracket(True)
    print "display the thread -> dot-bracket results (vsp)"
    print tn2db.vSeq
    print tn2db.vstr

    print "display the original (vs)" 
    tn2dbo = TreeNode2DotBracket(v2t)
    tn2dbo.make_dotbracket(True)
    print tn2db.vstr
    print tn2db.vSeq
    # print "planned exit after running LThread2Vienna"; sys.exit(0);
    
    
    # Test transformation between LThread and Vienna Pair
    vf = LThread2Vienna()
    vf.lt2vs(ltb.lt)
    vf.set_vstr(v2t.vstr)
    vf.set_vseq(v2t.vseq)
    
    print "structure sequence: "
    print vf.vstr
    print "vsBPlist: "
    for bpk in vf.vsBPlist:
        print bpk.disp_Pair()
    #
    print "vsPKlist: "
    for bpk in vf.vsPKlist:
        print bpk.disp_Pair()
    #
    print "vsMPlist: "
    for bpk in vf.vsMPlist:
        print bpk.disp_Pair()
    #
    
    # print "planned exit before running TreeNode2DotBracket"; sys.exit(0);
    
    tn2db = TreeNode2DotBracket(v2t)
    tn2db.make_dotbracket(True)
    #tn2db.visit(v2t.genTree)
    #print tn2db.makeFinal(tn2db.sssv)
    print tn2db.vSeq
    print tn2db.vstr
    #print "planned exit before running TreeNode2Motif"; sys.exit(0);
    
    
    t2m = TreeNode2Motif(v2t)
    t2m.visit(v2t.genTree)
    # print v2t.MPlist
    t2m.post_graftMP()
    print "finished TreeNode2Motif"
#   




def main(cl):
    if TEST == 0:
        test0(cl)
    elif TEST == 1:
        test1(cl)
        
    elif TEST == 2:
        test2(cl)
        
    elif TEST == 3:
        test3(cl)
#
   
# Main
if __name__ == '__main__':
    main(sys.argv)
#

