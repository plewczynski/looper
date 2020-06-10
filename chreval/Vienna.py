#!/usr/bin/env python3

"""@@@

Main Module:   Vienna.py 

Classes:       Vstruct
               ViennaData

Functions:

Author:        Wayne Dawson
creation date: ~2014/2015

last update:   200603 major revision: vastly improved consistency
               between 1D seq-declared and thread-built structures.
               It is still not clear if it can handle every possible
               input arrangement, but its ability to handle them has
               been vastly improved.

version:       1.0

#####################################################################
VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

note "issue xxxx 191103"

Not sure if the issue is actually a "problem" at present, but it is
clumsy.  In general, Vstruct is always called with no arguments in all
the programs.

200206: I have started addressing this. Actually, I didn't notice the
note here, but anyway.

200412: Unified the consistency between reading an input structure and
reading the LThread data.

AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
#####################################################################

Purpose:

This tool reads the Vienna formatted pairing sequence and converts it
into a set of pair-contacts, or in the case of chromatin, both pair
contacts and ctcf contacts.

Originally, this was intended as a tool for handling Vienna package
data with SimRNA. Now, it is used by some other programs, particularly
chreval.py and its derivatives. It also has been significantly
upgraded to handle and represent a very diverse variety of 1D
structures containing pseudoknot like structures and also multipair
interactions.

The two functions that carry this action out are

parse_SecondaryStructure       -- only dot bracket structures .((...)).
parse_fullDotBracketStructure  -- all types of complex dot bracket structures


This package can be tested in the following way

vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
Python 2.7.6 (default, Mar 22 2014, 22:59:56) 
[GCC 4.8.2] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> from Vienna import Vstruct
>>> vs = Vstruct()
>>> vs.parse_SecondaryStructure("(((((....)))))")
(((((....)))))
(    0,   13)
(    1,   12)
(    2,   11)
(    3,   10)
(    4,    9)
0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

"""

import sys
from copy import deepcopy
from collections import OrderedDict

# used for PKs and parallel stems
# this notation at least works with VARNA
from Constants import num2lpr
from Constants import num2rpr
from Constants import lpr2num
from Constants import rpr2num
from Constants import PKfndx
from Constants import PKrndx
from Constants import pointer
from Constants import stack
from Constants import counter
from Constants import Xlist

from BasicTools import sortPairListWRT_n
from BasicTools import copyList
from BasicTools import tuple2List

# fundamental representation of pairing
from Pair import Pair
from Pair import SortPair
from Pair import vsPair2list
from Pair import rlist2vsPair

# other constants and parameters
from MolSystem import MolSystem
from MolSystem import sysDefLabels # system RNA, Chromatin
from MolSystem import genChrSeq
from MolSystem import genRNASeq



debug = False # True # 

PROGRAM = "Vienna.py"
"""
tests:
  1 = test parse_SecondaryStructure for a simple
      secondary structure

  2 = test parse_fullDotBracketStructure for more
      complex structures that include pseudoknots. 
"""

dtestNames = { "test0"  :  0, 
               "test1"  :  1,
               "test2"  :  2,
               "test3"  :  3 }


def usage():
    print ("USAGE: %s -t testname [-mstr structure]" % PROGRAM)
    print ("testName options:")
    for v in dtestNames.keys():
        print (v)
    #
#

# error handling 
class MyException(Exception):
    def __init__(self, st):
        self.stmnt = st
    #
    def __str__(self):
        return self.stmnt
    #

    """ @
    Example:
      if not self.set_var:
          mssg = "\nERROR, value of var undefined\n"
          ex=MyException(mssg)
          raise ex
      #
    """
#


def LThread2Pair(lt, title = "undefined"):
    """coverts a list of LThread data to a list of Pair data"""
    
    BPlist = []
    for tr in lt.thread:
        btp = tr.btp    #  s, sp, sa, wyspa, bgn, end, etc.
        if not (btp == 'bgn' or btp == 'end' or btp == 'wyspa'):
            v = tr.ij_ndx
            i = v[0]; j = v[1]
            
            bp = Pair()
            # Vienna is specifically designed for secondary structure
            # of single strand (ss) pairs. Therefore, we use put_ssPair
            bp.put_ssPair(i, j, 'bp', btp)
            BPlist += [bp]
        #
    #
    ordering = SortPair()
    # sorting is not necessarily required but it makes it easier to
    # read at some level.
    BPlist = ordering.sortvsList(BPlist)
    
    return BPlist
#

def prune_pairlist(plist):
    
    # hopefully, this doesn't have to be used all that much, but it is
    # a good example
    plist = sortPairListWRT_n(plist, 0)
    k = 0
    
    while k < (len(plist) - 1):
        ik = plist[k  ][0]; jk = plist[k  ][1]
        il = plist[k+1][0]; jl = plist[k+1][1]
        
        # print ("ijk(%2d, %2d), ijl(%2d,%2d)" % (ik, jk, il, jl))
        
        if ik == il and jk == jl:
            del plist[k+1]
        else:
            k += 1
        #
    #

    return plist
#

def get_CPK_jbranches(ir, jr, il, jl, BProots):
    # BProots -> [class Pair]
    fpntr = []
    flag_forward = True
    if not jr < jl:
        flag_forward = False
    #
    
    

    #print ("flag_forward: ", flag_forward)
    for v in BProots:
        iv = v.i; jv = v.j
        #print ("ijv = (%2d,%2d)" % (iv, jv))
        if flag_forward:
            #print ("jr(%d) < iv(%d) && jv(%d) < jl(%d)" % (jr, iv, jv, jl))
            if jr < iv and jv < jl:
                fpntr += [(iv, jv)]
            #
            
        else:
            #print ("il(%d) < iv(%d) && jv(%d) < ir(%d)" % (il, iv, jv, ir))
            if il < iv and jv < ir:
                fpntr += [(iv, jv)]
            #
            
        #
        
    #|endfor
    
    #print (fpntr)
    
    return fpntr
#


def get_XPK_jbranches(ir1, jr1, ir2, jr2, BProots):
    fpntr = []
    
    for v in BProots:
        iv = v.i; jv = v.j
        if jr1 < iv and jv < ir2:
            fpntr += [(iv, jv)]
        #
        
    #|endfor
    
    return fpntr
#


def compress_to_jbranches(plist):
    
    # first remove any matching pairs
    
    #print (plist)
    #sys.exit(0)
    
    # finally, we have to find only those parts that form branches
    # in the structure
    
    k = 0
    while k < (len(plist) - 1):
        
        ik = plist[k][0]; jk = plist[k][1]
        #print ("k(%2d) (%2d,%2d)" % (k, ik, jk))
        l = k + 1
        flag_cont = True
        while l < len(plist) and flag_cont:
            il = plist[l][0]; jl = plist[l][1]
            #print ("l(%2d) (%2d,%2d)" % (l, il, jl))
            if ik < il and jl < jk:
                del plist[l]
            elif jk < il:
                flag_cont = False
            #
            
        #|endwhile
        
        k += 1
        
    #|endwhile
    
    return plist
#


class PKdmn(object):
    def __init__(self, ref_root, linkage):
        self.pktype = "" # C = Core, X = eXtended
        self.ipk = -1
        self.jpk = -1
        self.is_extended = None
        
        self.ref_root = [] # reference root stemtail(s)
        self.linkage = linkage
        
        self.root = [] # the actual root stemtail(s)
        
        # describe what is inside the root stem(s)
        self.r1levels = []
        self.r2levels = []
        
        # if the root stem is itself a PK, save this info
        self.r1pk = None # takes up class PKdmn if not None
        self.r2pk = None # takes up class PKdmn if not None
        
        self.jbranches = [] # branches in the junction region 
        
        self.dpkstack = {}
        if len(ref_root) == 1:
            self.is_extended = False
            #print ("core")
            self.pktype = "C" # core PK
            self.ref_root = ref_root
            vr = ref_root[0]
            self.ipk = vr[0]; self.jpk = vr[1]
        elif len(ref_root) == 2:
            self.is_extended = True
            #print ("extended")
            self.pktype = "X" # core PK
            self.ref_root = ref_root
            vr1 = ref_root[0]
            vr2 = ref_root[1]
            self.ipk = vr1[0]; self.jpk = vr2[1]
        else:
            print ("ERROR: root stem should consist of one or two tails")
            print ("input: ", ref_root)
            sys.exit(1)
        #
        
        
        
        
        
    #
    
    
    def update_PKdmn(self):
        #print ("update_PKdmn:")
        if self.pktype == "X":
            #print ("extended PK")
            #print ("type: r1", self.root[0], type(self.r1pk).__name__)
            if not type(self.r1pk).__name__ == 'NoneType':
                self.r1pk.update_PKdmn()
                #print ("X PK: using r1pk(%d,%d) %s" \
                #       % (self.r1pk.ipk, self.r1pk.jpk, self.r1pk))
                self.ipk = self.r1pk.ipk
            else:
                self.ipk = self.root[0][0]
            #
            
            #print ("type: r2", self.root[1], type(self.r2pk).__name__)
            if not type(self.r2pk).__name__ == 'NoneType':
                self.r2pk.update_PKdmn()
                #print ("X PK: using r2pk(%d,%d) %s" \
                #       % (self.r2pk.ipk, self.r2pk.jpk, self.r2pk))
                self.jpk = self.r2pk.linkage[0][1]
            else:
                self.jpk = self.root[1][1]
            #
            
        else:
            #print ("core PK")
            #print ("type: r1", self.root[0], type(self.r1pk).__name__)
            if not type(self.r1pk).__name__ == 'NoneType':
                #print ("C PK: using r1pk(%d,%d) %s" \
                #       % (self.r1pk.ipk, self.r1pk.jpk, self.r1pk))
                self.ipk = self.r1pk.ipk
                self.jpk = self.r1pk.jpk
            else:
                
                ir = self.root[0][0];    jr = self.root[0][1]
                ipkx = ir; jpkx = jr
                for lk in self.linkage:
                    il = lk[0]; jl = lk[1]
                    
                    if jl > jpkx:
                        # found a linkage that is larger than the
                        # current maximum range for the list
                        # containing the linkage(s)
                        self.jpk = jl
                        
                        # if only one value then becomes (ir, lk[1])
                        
                    #
                    
                    if il < ipkx:
                        # found a linkage that is smaller than the
                        # current minimum range for the list
                        # containing the linkage(s)
                        
                        self.ipk = il
                        
                        # if only one value then becomes (lk[0], jr)
                    #
                    
                #
                
                
                #print ("ijr(%2d,%2d), ijl(%2d,%2d) ijpk(%2d,%2d)" \
                #       % (ir, jr, il, jl, self.ipk, self.jpk))
                #sys.exit(0)
                
            #
            
        #
        
        #print ("ijpk(%3d,%3d)" % (self.ipk, self.jpk))
        #sys.exit(0)
    #
    
    
    def get_roots(self):
        
        act_roots = []
        
        if not type(self.r1pk).__name__ == 'NoneType':
            self.r1pk.update_PKdmn()
            #print ("X PK: using r1pk(%d,%d) %s" \
            #       % (self.r1pk.ipk, self.r1pk.jpk, self.r1pk))
            act_roots += [(self.r1pk.ipk, self.r1pk.jpk)]
        else:
            act_roots += [self.root[0]]
        #
            
        if self.pktype == "X":
            
            #print ("type: r2", self.root[1], type(self.r2pk).__name__)
            if not type(self.r2pk).__name__ == 'NoneType':
                self.r2pk.update_PKdmn()
                #print ("X PK: using r2pk(%d,%d) %s" \
                #       % (self.r2pk.ipk, self.r2pk.jpk, self.r2pk))
                act_roots += [(self.r2pk.ipk, self.r2pk.jpk)]
            else:
                act_roots += [self.root[1]]
            #
            
        #
        
        return act_roots
    #
    
    
    def add_pk_to_root1(self, pkdmn1): # class PKdmn
        #print ("writing r1: ", pkdmn1)
        self.r1pk = pkdmn1
        self.update_PKdmn()
        #print (self.__str__())
        #print ("r1 -> ijpk(%3d,%3d)" % (self.ipk, self.jpk))
        #print ("stop 1"); sys.exit(0)
    #
    
    
    def add_pk_to_root2(self, pkdmn2): # class PKdmn
        #print ("writing r2: ", pkdmn2)
        self.r2pk = pkdmn2
        self.update_PKdmn()
        #print ("r2 -> ijpk(%3d,%3d)" % (self.ipk, self.jpk))
        #print ("stop 2"); sys.exit(0)
        self.is_extended = True
    #


    def disp_dpkstack(self):
        s = ''
        skeys = list(self.dpkstack.keys())
        if len(skeys) > 0:
            for k in range(len(skeys)-1):
                s += "%s\n" % (self.dpkstack[skeys[k]])
            #|endfor
            s += "%s" % (self.dpkstack[skeys[len(skeys)-1]])
        #
        
        return s
    #
    
    def disp_linkage(self):
        n = len(self.linkage)
        s = "linkage["
        for k in range(n-1):
            v = self.linkage[k]
            s += "(%3d,%3d) " % (v[0], v[1])
        #|endfor
        
        v = self.linkage[n - 1]
        s += "(%3d,%3d)]" % (v[0], v[1])
        
        return s
    #
    
    
    def __str__(self):
        s = '(%3d,%3d) ' % (self.ipk, self.jpk)
        if len(self.root) > 1:
            
            if not type(self.r1pk).__name__ == "NoneType":
                s += "root1(%3d,%3d) " % (self.r1pk.ipk, self.r1pk.jpk)
            else:
                s += "root1(%3d,%3d) " % (self.root[0][0], self.root[0][1])
            #
            
            if not type(self.r2pk).__name__ == "NoneType":
                s += "root2(%3d,%3d) " % (self.r2pk.ipk, self.r2pk.jpk)
            else:
                s += "root2(%3d,%3d) " % (self.root[1][0], self.root[1][1])
            #
            
            s += self.disp_linkage()
            
        elif len(self.root) == 1:
            
            if not type(self.r1pk).__name__ == "NoneType":
                s += "root1(%3d,%3d) " % (self.r1pk.ipk, self.r1pk.jpk)
            else:
                s += "root1(%3d,%3d) " % (self.root[0][0], self.root[0][1])
            #
            
            s += self.disp_linkage()
            
        elif len(self.ref_root) > 1:
            # initial structure is one that is not fully resulved
            s += "ref_root1(%3d,%3d) " % (self.ref_root[0][0], self.ref_root[0][1])
            s += "ref_root2(%3d,%3d) " % (self.ref_root[1][0], self.ref_root[1][1])
            s += self.disp_linkage()
            
        elif len(self.ref_root) == 1:
            s += "ref_root(%3d,%3d) " % (self.ref_root[0][0], self.ref_root[0][1])
            s += self.disp_linkage()
            
        else:
            print ("what? ... PKdmn")
            sys.exit(1)
        #
        
        if self.is_extended:
            s += " extended"
        else:
            s += " core"
        #
        
        return s
    #

    def __repr__(self):
        return self.__str__()
    #
#



class Weights(object):
    def __init__(self, i, j): 
        
        self.i = i
        self.j = j
        
        self.nsst = 0
        # nsst: Number of Secondary Structure stem Tails. This
        # basically counts the number of stem tails that would (in
        # principle) correspond to the secondary structure in the
        # region between i and j (excluding i and j). If there are no
        # additional stem tails between i and j, then this value is
        # zero. The count is obtained by evaluating subtract_crossovers()
        
        self.stem = 1
        # It is important to know where indendenent stems are
        # located. This would help serve as a beginning point from
        # BProots, though it would need to be sorted so that the PKs
        # are plucked out.
        
        self.npairs  = 1
        # this is the main measure of weight. The reason the number is
        # one is because a domain containes at least itself as a
        # pair. Weights is only used in conjuction with existing
        # BPlist. If there is no pair, this cannot be called.
        
        self.rpntr = []
        self.fpntr = []
        self.spntr = []
    #
    
    def show_Weights(self, show_header = False):
        if show_header:
            self.show_header()
        #
        
        s = ("(%3d,%3d)   %3d     %3d     %-20s   %-20s" \
             % (self.i, self.j, self.stem, self.npairs, self.rpntr, self.fpntr))
        return s
    #

    def show_header(self):
        s = "(i,  j)     stem    npairs  rpntr                  fpntr"
        return s
    #
    
    def __str__(self):
        s = ("%3d %3d   %3d   %3d  %3d   %s  %s" \
             % (self.i, self.j, self.nsst, self.stem, self.npairs, self.rpntr, self.fpntr))
        return s
    #
    
    def __repr__(self):
        return self.__str__()
    #
#


class Vstruct(SortPair):
    def __init__(self, molsys): # class MolSystem
        ###########################
        self.molsys  = molsys
        self.dinput  = None  # class InputSettings
        ###########################
        """@
        
        I'd say that molsys and dinput are mostly here for the purpose
        of the "cargo" they carry; they provide the basic tools for
        informing other modules what sorts of energy functions are
        required and where the data came from. The most critical
        parameters that should be fully available to Vstruct are vstr,
        vseq, N and wt.
        
        """
        
        # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        self.system       = molsys.system
        self.useTurner    = molsys.useTurner
        self.useViS       = molsys.useViS
        self.usegMatrix   = molsys.usegMatrix
        self.useChromatin = molsys.useChromatin
        self.jobtype      = molsys.jobtype
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        
        self.vstr    = '' # structure
        self.vseq    = '' # sequence (optional)
        
        self.N       = -1 # sequence length
        
        self.wt      = -1 # used with my_generation to add specified weights to heatmaps
        
        # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        
        # Probably should insist that all these things be defined
        
        if len(molsys.mstr) > 0:
            self.vstr = molsys.mstr
        else:
            print ("ERROR: no structure defined")
            sys.exit(0)
        #
        
        if len(molsys.mseq) > 0:
            self.vseq = molsys.mseq
            
        else:
            if self.jobtype == "generator":
                # we don't care about the sequence
                if self.system == "Chromatin":
                    self.vseq = genChrSeq(self.vstr)
                elif self.system == "RNA":
                    self.vseq = genRNASeq(self.vstr)
                else:
                    print ("ERROR: undefined system %s" % self.system)
                    sys.exit(1)
                #
                
            else:
                print ("ERROR: no sequence defined")
                sys.exit(0)
            #
            
        #
        
        if len(self.vseq) == len(self.vstr):
            self.N = len(self.vseq)
        else:
            print ("ERROR: sequence length and structure length don\'t match")
            print ("seq(%d) != str(%d)" % (len(self.vseq), len(self.vstr)))
            sys.exit(1)
        #
        
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        
        self.use_old_build = False
        self.ppgap = 2
        self.apgap = 3
        # presently, I would say ppgap = 2, but a waffle on 2 or 3.
        
        # various dictionaries for statistics
        self.prncpl_wts  = OrderedDict() # PRiNCiPaL Weights (list of class Weights)
        self.pkeylevels  = OrderedDict() # Principal KEY LEVELS
        self.dPKlinkages = {} # place holder to store initial PK information
        self.dQroots     = {} # list of root stems used in PKlinkage development
        
        self.dCPKlinfo = {} # Core PseudoKnot Linkage INFOrmation (raw)
        self.dCPK_r2l  = {} # Core PseudoKnot: root to linkage info
        self.dCPKinfo  = {} # the final Core PK INFOrmation
        
        self.dXPKlinfo = {} # eXtended PseudoKnot Linkage INFOrmation (raw)
        self.dXPK_xp2l = {} # eXtended PK: ex PK to linkages
        self.dXPK_xp2r = {} # eXtended PK: ex PK to root stems
        self.dXPKinfo  = {} # the final extended PK INFOrmation
        # both are dictionaries that contain objects of class PKdmn
        
        self.pxlist_ref = []
        self.flag_show_pxlist = False # True # 
        """@ 
        
        There is a broad proliferation of dictionaries here, so I
        provide a more expanded explanation of what the contents are.
        
        prncpl_wts (PRiNCiPaL WeighTS) 
        
        -- various quanties about domains
        
        pkeylevels (Principal KEY LEVELS)
        
        -- this actually achieves a lot of the goal of TreeNode for
           the secondary structure.
        
        dPKlinkages (Dict of PseudoKnot Linkages) 
        
        -- place holder to store initial PK information. Actually,
           this list is constructed in separate_sheep_from_goats and
           in findMaxStemTails it is split into dCPKlinfo and
           dXPKlinfo. Therefore, I think that I might just eliminate
           it in go directly to these two dictionaries in
           separate_sheep_from_goats. It eventually becomes these
           other two dictionaries anyway.
        
        dQoots (Dict of Qualified ROOTS) 
        
        -- place holder dictionary of root stems that is used in
           constructing dPKlinkage so that root stems claimed by other
           PKs in the structure are not inadvertantly reclassified as
           linkages. It is important when there are multiple linkages
           in a structure
        
        
        dCPKlinfo (dict of Core PseudoKnot Linkage INFOrmation (raw))
        dXPKlinfo (dict of eXtended PseudoKnot Linkage INFOrmation (raw))
        
        -- dCPKlinfo and dXPKlinfo contain an object of class PKdmn. 
        
           I found that I have to split the linkages into core and
           extended PKs from the start and a gradually build up a more
           accurate construction of the PK by successively working
           through the ordered list. In this first step, 
        
           These serve as very important building tools, holders if
           you will, for constructing a full PK construction that can
           be passed from here to Vienna2TreeNode. These dictionaries
           categorize and save the linkage and root information for
           each entry in self.PKroots. In addition, the store related
           critical information from prncpl_wts.
        
        dCPKinfo (dict final Core PK INFOrmation)
        dXPKinfo (dict final eXtended PK INFOrmation)
        
           This is the final result of the process. Where the
           information of the root stem is stored in dCPKlinfo and
           dXPKlinfo, these have the coordinates of the PK itself.
        
        dCPK_r2l       (Core PseudoKnot: root to linkage info)
        self.dXPK_xp2l (eXtended PK: ex PK to linkages)
        self.dXPK_xp2r (eXtended PK: ex PK to root stems)
        
           These are direct hashes. For dCPK_r2l, we use the root key
           to provide a reference to linkage stem 5'/3' ends.
        
           For dXPK_xp2l and dXPK_xp2r we are providing the pseudoknot
           boundaries as the reference and these then point either to
           the beginning of the linage nad root stms.
        
        dstcount  (Dictionary of STructure COUNTed)
        
           in making prncpl_wts, we need to do them just once.
        
        save_branches = {}
        
           Presently (200515), not sure if I will need this. I am
           having some trouble with tandem PKs and thought this might
           be a way to keep track of the correct tandems. However, I
           begin to see that the problem is a bit more complex.
        
        """
        
        self.sskeylist  = []
        self.pkkeylist  = []
        # standard secondary structure (ss)
        self.BPlist  = []
        self.BProots = []
        
        # a both ss and pseudoknot (PK) in same line
        self.PKlist  = []
        self.PKroots = []
        
        # CTCF islands, or triple helices
        self.MPlist  = []
        
        """@
        
        170306wkd: I know the following looks a little strange because
        I import Xlist (etc.) from Constants.py; however, for some
        reason that I cannot explain, the module Xlist is changed when
        these objects are changed, even though self.Xlist should be a
        independent object and Xlist should be just a
        template. Nevertheless, if I use
        
          self.Xlist.update({xl : Xlist[xl]}),
          self.stack.update({sl : stack[sl]}),
          etc.
        
        The module values are updated with previous entries generated
        from prior created Vstruct objects.  To make a long story
        short, finally what I had to do to make sure that self.Xlist
        didn't end up with misinformation was to hard wire write it
        independent of the contents of the module Xlist. So, in the
        end, the keys from the module are used, but not the
        contents. Strange, but anyway.....
        
        """
        
        
        
        # reset the key reference dictionary of pointers 
        self.Xlist = {}
        for xl in Xlist.keys():
            self.Xlist.update({xl : [] })
        #|endfor
        
        self.stack = {}
        for sl in stack.keys():
            self.stack.update({sl : 0 })
        #|endfor
        
        self.counter = {}
        for cl in counter.keys():
            self.counter.update({cl : 0 })
        #|endfor
        
        self.pointer = {}
        for pl in pointer.keys():
            self.pointer.update({pl : 0 })
        #|endfor
        
    #
    
    
    def set_system(self, iSetUp): # class InputSettings
        
        ###########################
        self.dinput = iSetUp
        ###########################
        sname = iSetUp.system
        if not sname in sysDefLabels:
            print ("ERROR: unrecognized system name (%s)" % sname)
            print ("       allowed names ")
            for snm in sysDefLabels.keys():
                print (snm)
            #|endfor
            
            sys.exit(1)
            
        #
        
        ###########################
        self.molsys = iSetUp.molsys
        ###########################
        
        
        self.system = sname
        
        
        if self.system == "RNA":
            self.vseq    = self.dinput.rnaseq  # sequence (optional) 
            self.vstr    = self.dinput.rnastr  # structure
            
        elif self.system == "Chromatin":
            self.vseq    = self.dinput.chrseq  # sequence (optional) 
            self.vstr    = self.dinput.chrstr  # structure
        #
        
        self.N       = len(self.vseq)      # sequence length
        
        #print ("rSeq", self.dinput.rSeq)
        #print (self.system)
        
        
    #
    
    # init_vs_from_Xlf(self, vf):
    def init_vs_from_XlistFile(self, vf):
        # Basically, the only difference between init_vs_from_lt and
        # init_vs_from_XlistFile is the origin of the input object
        # vf. In this one, we requre List2Vienna, in the other we
        # require LThread2Vienna.
        debug_init_vs_from_XlistFile = False # True # 
        if debug_init_vs_from_XlistFile:
            print ("Entered init_vs_from_XlistFile")
        #
        
        object_name = type(vf).__name__
        print (object_name)
        if not object_name == "List2Vienna":
            print ("ERROR(init_vs_from_XlistFile) requires type List2Vienna")
            sys.exit(1)
        #
        
        ###########################
        self.molsys = vf.molsys
        """self.dinput = ?? #  """
        ###########################
        
        self.reset_VstructLists()
        # system dependent parameters
        self.system       = vf.molsys.system
        self.useTurner    = vf.molsys.useTurner
        self.useViS       = vf.molsys.useViS
        self.usegMatrix   = vf.molsys.usegMatrix
        self.useChromatin = vf.molsys.useChromatin
        
        self.N         = vf.N         
        self.vstr      = vf.molsys.mstr      
        self.vseq      = vf.molsys.mseq
        
        
        self.BPlist    = vf.vsBPlist
        self.BProots   = vf.vsBProots 
        self.PKlist    = vf.vsPKlist   # PKs
        self.PKroots   = vf.vsPKroots 
        self.MPlist    = vf.vsMPlist   # e.g., ctcf etc
        
        
        
        if debug_init_vs_from_XlistFile:
            self.disp_allLists("from init_vs_from_XlistFile (before corrections)", False)
        #
        
        if self.use_old_build:
            # post base pair list processing
            self.PKlist = self.sortvsList(self.PKlist, 'i')
            self.PKlist = self.assign_PKdirection(self.PKlist)
            self.compress_to_BPlist()
            if debug_init_vs_from_XlistFile:
                self.disp_allLists("from init_vs_from_XlistFile (after adjustments)", False)
                #print ("stop at end of init_vs_from_XlistFile"); sys.exit(0)
            #
            
        else:
            self.build1DStructureLists(self.BPlist, self.PKlist)
            if debug_init_vs_from_XlistFile:
                print ("Results using new_build")
                self.disp_allLists("from init_vs_from_XlistFile (after adjustments)", False)
                print ("stop at end of init_vs_from_XlistFile"); sys.exit(0)
            #
            
        #
            
        
    #endMethod
    
    def init_vs_from_lt(self, vf):
        # basically, the only difference between init_vs_from_lt and
        # init_vs_from_XlistFile is the origin of the input object
        # vf. In this one, we requre LThread2Vienna, in the other we
        # require List2Vienna.
        debug_init_vs_from_lt = False # True # 
        if debug_init_vs_from_lt:
            print ("Entered init_vs_from_lt")
        #
        
        object_name = type(vf).__name__
        #print (object_name)
        if not object_name == "LThread2Vienna":
            print ("ERROR(init_vs_from_lt) requires type LThread2Vienna")
            sys.exit(1)
        #
        
        ###########################
        self.molsys = vf.molsys
        """self.dinput = ?? #  """
        ###########################
        
        self.reset_VstructLists()
        # system dependent parameters
        self.system       = vf.molsys.system
        self.useTurner    = vf.molsys.useTurner
        self.useViS       = vf.molsys.useViS
        self.usegMatrix   = vf.molsys.usegMatrix
        self.useChromatin = vf.molsys.useChromatin
        
        self.N         = vf.N         
        self.vstr      = vf.molsys.mstr      
        self.vseq      = vf.molsys.mseq
        
        self.BPlist    = vf.vsBPlist
        self.BProots   = vf.vsBProots 
        self.PKlist    = vf.vsPKlist   # PKs
        self.PKroots   = vf.vsPKroots 
        self.MPlist    = vf.vsMPlist   # e.g., ctcf etc
        
        
        
        if debug_init_vs_from_lt:
            self.disp_allLists("from init_vs_from_lt (before corrections)", False)
        #
        
        if self.use_old_build:
            # post base pair list processing
            self.PKlist = self.sortvsList(self.PKlist, 'i')
            self.PKlist = self.assign_PKdirection(self.PKlist)
            self.compress_to_BPlist()
            if debug_init_vs_from_lt:
                self.disp_allLists("from init_vs_from_lt (after adjustments)", False)
                #print ("stop at end of init_vs_from_lt"); sys.exit(0)
            #
            
        else:
            self.build1DStructureLists(self.BPlist, self.PKlist)
            if debug_init_vs_from_lt:
                print ("Results using new_build")
                self.disp_allLists("from init_vs_from_lt (after adjustments)", False)
                print ("stop at end of init_vs_from_lt"); sys.exit(0)
            #
            
        #
        
        """@@@
        
        It turns out, with the structure "(.).(.([)])", the resulting
        BPlist and PKlist from chreval.py is somewhat different from
        when the list is read from the sequence.
        
        ERROR(buildTreeDmns near 10): somehow, core PK is not properly defined
        iz( 0)   < pLmin( 7) < jz( 2)   < qLmax( 9)
        <or>   pLmin( 7) < iz( 0)   < qLmax( 9) < jz( 2)
        buildTreeDmns near 10
        structural sequence
        
        secondary structure
        (    0,    2)[a]
        (    4,   10)[a]
        (    6,    8)[a]   <=== 
        (    7,    9)[a]   <=== these should both be [p]
        secondary structure roots
        (    0,    2)[a]
        (    4,   10)[a]
        
        The output from vs is 
        
        structure sequence: 
        (.).(.([)])
        vsBPlist: 
        (    0,    2)[a]
        (    4,   10)[a]
        (    6,    8)[p]  
        (    6,    8)[p]  <=== evidently double counts but at least both are marked [p]
        vsPKlist: 
        (    7,    9)[p]
        vsMPlist: 
        xcycxcxxyyy
        (.).(.[<]>)
        
        .... So it seems that I will have to run an independent check of
        the structure through LThread, and LThreadBuilder and here.
        
        """
    #endMethod
    
    
    def set_Vstruct(self, s):
        self.vstr = s
        self.N = len(s)
        return 0
    #endMethod
    
    
    def set_Vseq(self, s):
        self.vseq = s
        return 0
    #endMethod
    
    
    def set_VstructSystem(self, molsys): # class MolSystem
        """@
        
        In general, this should not be necessary to you ever, except
        perhaps in testing this program from idiotic things that might
        happen. If you start with a particular potential, it seems
        rather odd that you would suddenly change it.
        
        """
        
        ###########################
        self.molsys  = molsys
        self.dinput  = None
        """# dinput -> class InputSettings"""
        ###########################
        
        # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        self.system       = molsys.system
        self.useTurner    = molsys.useTurner
        self.useViS       = molsys.useViS
        self.usegMatrix   = molsys.usegMatrix
        self.useChromatin = molsys.useChromatin
        self.jobtype      = molsys.jobtype
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        
        self.wt      = -1
        # used qwith my_generation to add specified weights to heatmaps
        
        
        # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        
        # Probably should insist that all these things be defined
        
        if len(molsys.mstr) > 0:
            self.vstr = molsys.mstr
        else:
            print ("ERROR: no structure defined")
            sys.exit(0)
        #
        
        if len(molsys.mseq) > 0:
            self.vseq = molsys.mseq
            
        else:
            if   self.jobtype == "generator":
                # we don't care about the sequence
                if self.system == "Chromatin":
                    self.vseq = genChrSeq(self.vstr)
                elif self.system == "RNA":
                    self.vseq = genRNASeq(self.vstr)
                else:
                    print ("ERROR: undefined system %s" % self.system)
                    sys.exit(1)
                #
                
            elif self.jobtype == "testing":
                # specifics about the sequence will be defined later
                if self.system == "Chromatin":
                    self.vseq = "cccccc"
                elif self.system == "RNA":
                    self.vseq = "UUUUUU"
                else:
                    print ("ERROR: undefined system %s" % self.system)
                    sys.exit(1)
                #
                
            else:
                print ("ERROR: no sequence defined")
                sys.exit(0)
            #
            
        #
        
        if len(self.vseq) == len(self.vstr):
            self.N = len(self.vseq)
        else:
            print ("ERROR: sequence length and structure length don\'t match")
            print ("seq(%d) != str(%d)" % (len(self.vseq), len(self.vstr)))
            sys.exit(1)
        #
        
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        
        self.reset_VstructLists()
    #
    
    def reset_VstructLists(self):
        """@
        
        Previously, I was resetting everything including system, N,
        vseq and vstr. However, this creates problems when working
        with Chromatin because everything will default back to RNA. So
        I separated reset_VstructLists into reset_VstructLists for the
        various contact information and and reset_VstructSystem()
        which resets the entire system and potential functions.
        
        """
        
        # various dictionaries for statistics
        self.prncpl_wts  = OrderedDict() # PRiNCiPaL Weights (list of class Weights)
        self.pkeylevels  = OrderedDict() # Principal KEY LEVELS
        self.dPKlinkages = {} # place holder to store initial PK information
        self.dQroots     = {} # list of root stems used in PKlinkage development
        
        self.dCPKlinfo = {} # Core PseudoKnot Linkage INFOrmation (raw)
        self.dCPK_r2l  = {} # Core PseudoKnot: root to linkage info
        self.dCPKinfo  = {} # the final Core PK INFOrmation
        
        self.dXPKlinfo = {} # eXtended PseudoKnot Linkage INFOrmation (raw)
        self.dXPK_xp2l = {} # eXtended PK: ex PK to linkages
        self.dXPK_xp2r = {} # eXtended PK: ex PK to root stems
        self.dXPKinfo  = {} # the final extended PK INFOrmation
        # both are dictionaries that contain objects of class PKdmn
        self.pxlist_ref = []
        
        #self.flag_show_pxlist = False # It seems better not to reset it

        
        self.sskeylist  = []
        self.pkkeylist  = []
        # standard secondary structure (ss)
        self.BPlist  = []
        self.BProots = []
        # a both ss and pseudoknot (PK) in same line
        self.PKlist  = []
        self.PKroots = []
        
        # CTCF islands
        self.MPlist = []
        
        
        """@
        
        170306wkd: I know the following looks a little strange because
        I import Xlist (etc.) from Constants.py; however, for some
        reason that I cannot explain, the module Xlist is changed when
        these objects are changed, even though self.Xlist should be an
        independent object and Xlist should be just a
        template. Nevertheless, if I use
        
          self.Xlist.update({xl : Xlist[xl]}),
          self.stack.update({sl : stack[sl]}),
          etc.
        
        The module values are updated with previous entries generated
        from prior created Vstruct objects.  To make a long story
        short, finally what I had to do to make sure that self.Xlist
        didn't end up with misinformation was to hard wire write it
        independent of the contents of the module Xlist. So, in the
        end, the keys from the module are used, but not the
        contents. Strange, but anyway.....
        
        """
        
        # reset the key reference dictionary of pointers 
        self.Xlist = {}
        for xl in Xlist.keys():
            self.Xlist.update({xl : []})
        #|endfor
        
        self.stack = {}
        for sl in stack.keys():
            self.stack.update({sl : 0})
        #|endfor
        
        self.counter = {}
        for cl in counter.keys():
            self.counter.update({cl : 0})
        #|endfor
        
        self.pointer = {}
        for pl in pointer.keys():
            self.pointer.update({pl : 0})
        #|endfor
        
        return 0
    #endMethod
    
    
    def scan_vs(self, i):
        
        # very primative secondary structure reader (the most original
        # in the family), can do secondary structure
        if debug:
            print ("scan_vs(%d):" % i)
        #
        
        mm = -1
        if i < self.N:
            s = self.vstr[i]
            # separate '(' and ')' from '[', ']', and '.'
            if not (s == '(' or s == ')' or s == '[' or s == ']' or s == '.'):
                print ('ERROR: improperly defined Fontana formated sequence')
                print ('       offending character \'%c\' located at position %d' \
                       % (s, i+1))
                sys.exit(1)
            #
            
            if s == '(':
                # add to the stack
                b = Pair()
                b.put_ssPair_i(i, "bp")
                self.BPlist += [b]
                if debug:
                    print (b.disp_Pair())
                    # print ('i = %5d, %s' % (i, s))
                #
                
                self.stack[0] += 1 
                self.pointer[0] = self.stack[0]
                self.counter[0] +=1
                mm = i + 1
                self.scan_vs(mm)
            elif s == ')':
                self.pointer[0] = self.find_next_Xpoint_n(self.BPlist, self.pointer[0])
                self.BPlist[self.pointer[0]-1].put_ssPair_j(i, "bp")
                if debug:
                    print (self.BPlist[self.pointer[0]-1].disp_Pair())
                    # print ('j = %5d, %s' % (i, s))
                #
                
                self.pointer[0] -= 1
                self.counter[0] -= 1
                mm = i + 1
                self.scan_vs(mm)
            else:
                mm = i
                # ignore pseudoknot brackets [[[...]]], if present.
                # Converts these to unpairsed strand data.
                if (s == '[') or (s == ']'):
                    s = '.'
                #
                
                # Either look for next closing bracket in sequence or
                # terminate at the end of the sequence if nothing is
                # found.
                while (s == '.') and (mm < self.N):
                    if debug:
                        print ('i = %5d, %s' % (mm, s))
                    #
                    
                    mm += 1
                    if mm == self.N:
                        break
                    #
                    
                    s = self.vstr[mm]
                    
                #|endwhile
                
                self.scan_vs(mm)
            #
            
        else:
            if debug:
                print ('point = %d' % self.pointer[0])
            #
            
            if not self.counter[0] == 0:
                case = 0
                if self.counter[0] > 0:
                    case = 0
                else:
                    case = 1
                #
                
                print ('jcount = %d' %  self.counter[0])
                print ('ERROR!!! Fontana notation is not correct!')
                if case == 0:
                    print ('         %d too many \'(\' brackets!' % self.counter[0])
                else:
                    print ('         %d too many \')\' brackets!' % (-self.counter[0]))
                #
                
                print ('input structure:')
                print (self.vstr)
                sys.exit(1)
            #
            
            if not self.check_Xlist_n(self.BPlist):
                print ('ERROR!!! Fontana notation is not correct!')
                print ('         At least one structure of type \')...(\' was found')
                print ('input structure:')
                print (self.vstr)
                sys.exit(1)
            #
            
        #
        
        return mm
    #endMethod
    
    
    def scan_allTypes(self, i, layer):
        
        # this is far more general 1D structure reader
        
        debug_bp   = False # True # 
        debug_PK   = False # True # 
        debug_ctcf = False # True # 
        if debug_bp or debug_PK or debug_ctcf:
            print ("enter scan_allTypes(%d):" % i, self.pointer[0])
        #
        
        if layer > self.N:
            # if it does this, something is definitely wrong! The
            # variable layer is mainly used to track the recursion
            # level and make sure that things have not gone weird. It
            # seems like this part of the program works fine now. I
            # have not encoutered a stop in ages. However, the
            # debugging code is still here because I don't know when
            # it might be needed.
            print ("ERROR(scan_allTypes): layer(%d) exceeds the sequence length(%d)!" \
                % (layer, self.N))
            print ("                      Something is wrong. ... Last call i = %d"   % (i))
            sys.exit(1)
        #
        
        mm = -1
        if i < self.N:
            s = self.vstr[i]
            if not ((s in lpr2num) or (s in rpr2num) or s == '|' or s == '.'):
                print ('ERROR(scan_allTypes): improperly defined Fontana formated sequence')
                print ('                      offending character \'%c\' located at' % (s))
                print ('position %d' % (i+1))
                sys.exit(1)
            #
            
            if s == '(': # standard bp argument 
                # add to the stack
                BPndx = lpr2num[s]
                if debug_bp:
                    print ("BPndx(%2d) => %s" % (BPndx, s))
                #
                
                b = Pair()
                b.put_ssPair_i(i, "bp")
                self.Xlist[BPndx] += [b]
                if debug_bp:
                    print (b.disp_Pair())
                    print ('i = %5d, %s' % (i, s))
                    print ("element(%2d) => counter(%2d), pointer(%s), stack(%2d)" \
                        % (BPndx, self.counter[0], self.pointer[0], self.stack[0]))
                #
                
                mm = i + 1
                self.stack[BPndx]   += 1
                self.pointer[BPndx]  = self.stack[BPndx]
                self.counter[BPndx] += 1
                self.scan_allTypes(mm, layer + 1)
                
            elif s == ')':
                BPndx = rpr2num[s]
                if debug_bp:
                    print ("BPndx(%2d) => %s" % (BPndx, s))
                    flag_details = False
                    if flag_details:
                        print ("list of current contacts")
                        for bp in self.Xlist[BPndx]:
                            print (bp.disp_Pair())
                        #|endfor
                        
                        print ("------")
                    #
                    
                #
                
                self.pointer[BPndx] \
                    = self.find_next_Xpoint_n(self.Xlist[BPndx], self.pointer[BPndx])
                self.Xlist[BPndx][self.pointer[BPndx]-1].put_ssPair_j(i, "bp")
                if debug_bp:
                    print (self.Xlist[BPndx][self.pointer[BPndx]-1].disp_Pair())
                    print ('j = %5d, %s' % (i, s))
                    print ("element(%2d) => counter(%2d), pointer(%s), stack(%2d)" \
                        % (BPndx, self.counter[0], self.pointer[0], self.stack[0]))
                #
                
                self.pointer[BPndx] -= 1
                self.counter[BPndx] -= 1
                mm = i + 1
                self.scan_allTypes(mm, layer + 1)
                
            elif s == '[': # PK
                PKndx = lpr2num[s]
                if debug_PK:
                    print ("PKndx(%2d) => %s" % (PKndx, s))
                #
                
                mm = i
                b = Pair()
                b.put_ssPair_i(i, "pk")
                self.Xlist[PKndx] += [b]
                if debug_PK:
                    print (b.disp_Pair())
                    print ('i = %5d, %s' % (i, s))
                #
                
                self.stack[PKndx]   += 1 
                self.pointer[PKndx]  = self.stack[PKndx]
                self.counter[PKndx] +=1
                mm = i + 1
                self.scan_allTypes(mm, layer + 1)
                
            elif s == ']':
                PKndx = rpr2num[s]
                if debug_PK:
                    print ("PKndx(%2d) => %s" % (PKndx, s))
                #
                
                self.pointer[PKndx] \
                    = self.find_next_Xpoint_n(self.Xlist[PKndx], self.pointer[PKndx], "PK")
                self.Xlist[PKndx][self.pointer[PKndx]-1].put_ssPair_j(i, "pk")
                if debug_PK:
                    print (self.Xlist[PKndx][self.pointer[PKndx]-1].disp_Pair())
                    print ('j = %5d, %s' % (i, s))
                #
                
                self.pointer[PKndx] -= 1
                self.counter[PKndx] -= 1
                mm = i + 1
                self.scan_allTypes(mm, layer + 1)
                
                
                
                # ########################################
                # CTCF
                # ########################################
                # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
                
            elif s == '{':  # xxxx  CTCF  xxxx
                CTCFndx = lpr2num[s]
                if debug_ctcf:
                    print ("CTCFndx(%2d) => %s" % (CTCFndx, s))
                #
                
                mm = i
                b = Pair()
                b.put_ssPair_i(i, "ctcf")
                self.Xlist[CTCFndx] += [b]
                if debug_ctcf:
                    print (b.disp_Pair())
                    print ('i = %5d, %s' % (i, s))
                #
                
                self.stack[CTCFndx]   += 1 
                self.pointer[CTCFndx]  = self.stack[CTCFndx]
                self.counter[CTCFndx] +=1
                mm = i + 1
                self.scan_allTypes(mm, layer + 1)
                
            elif s == '}':
                CTCFndx = rpr2num[s]
                if debug_ctcf:
                    print ("CTCFndx(%2d) => %s" % (CTCFndx, s))
                #
                
                self.pointer[CTCFndx] = self.find_next_Xpoint_n(self.Xlist[CTCFndx],
                                                                self.pointer[CTCFndx],
                                                                "CTCF")
                self.Xlist[CTCFndx][self.pointer[CTCFndx]-1].put_ssPair_j(i, "ctcf")
                if debug_ctcf:
                    print (self.Xlist[CTCFndx][self.pointer[CTCFndx]-1].disp_Pair())
                    print ('j = %5d, %s' % (i, s))
                #
                
                
                self.pointer[CTCFndx] -= 1
                self.counter[CTCFndx] -= 1
                mm = i + 1
                self.scan_allTypes(mm, layer + 1)
                
            elif s == '|':
                CTCFndx = lpr2num['{']
                if debug_ctcf:
                    print ("CTCFndx(%2d) => |" % (CTCFndx))
                #
                
                self.Xlist[CTCFndx][self.pointer[CTCFndx]-1].put_contacts(i)
                mm = i + 1
                if debug_ctcf:
                    print ('j = %5d, %s' % (i, s))
                #
                
                self.scan_allTypes(mm, layer + 1)
                # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                # ########################################
                # CTCF
                # ########################################
                
                
                
                # ########################################
                # expanded PK notations
                # ########################################
                # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
            elif s in lpr2num:
                l_mate = s
                r_mate = num2rpr[lpr2num[s]]
                PKndx = lpr2num[s]
                if debug_PK:
                    print ("PKndx(%2d) => %s, search for pattern %s and %s" \
                        % (PKndx, s, l_mate, r_mate))
                    
                    print ("left(i = %d) ==> %s" % (i, s))
                    print ("lpr2num[s]:                          ", lpr2num[s])
                    print ("num2lpr[lpr2num[s]]:                 ", num2lpr[lpr2num[s]])
                    print ("PKrndx[lpr2num[s]]:                  ", PKrndx[lpr2num[s]])
                    print ("num2lpr[PKfndx[PKrndx[lpr2num[s]]]]: ", \
                        num2lpr[PKfndx[PKrndx[lpr2num[s]]]])
                    print ("search for connecting element:          '%s'" \
                        % num2rpr[lpr2num[s]])
                    print ("element(%2d) => counter(%2d), pointer(%s), stack(%2d)" \
                        % (PKndx, self.counter[PKndx],
                           self.pointer[PKndx], self.stack[PKndx]))
                #
                
                mm = i
                b = Pair()
                b.put_ssPair_i(i, "pk")
                self.Xlist[PKndx] += [b]
                if debug_PK:
                    print (b.disp_Pair())
                    print ('i = %5d, %s' % (i, s))
                #
                
                self.stack[PKndx]   += 1 
                self.pointer[PKndx]  = self.stack[PKndx]
                self.counter[PKndx] += 1
                mm = i + 1
                self.scan_allTypes(mm, layer + 1)
                
            elif s in rpr2num:
                l_mate = num2lpr[rpr2num[s]]
                r_mate = s
                PKndx = rpr2num[s]
                ndx = PKrndx[rpr2num[s]]
                if debug_PK:
                    print ("right(i = %d) ==> %s" % (i, s))
                    print ("rpr2num[s]:                          ", rpr2num[s])
                    print ("num2rpr[rpr2num[s]]:                 ", num2rpr[rpr2num[s]])
                    print ("PKrndx[rpr2num[s]]:                  ", PKrndx[rpr2num[s]])
                    print ("num2rpr[PKfndx[PKrndx[rpr2num[s]]]]: ", \
                        num2rpr[PKfndx[PKrndx[rpr2num[s]]]])
                    print ("search for connecting element:          '%s'" \
                        % num2rpr[rpr2num[s]])
                    print ("element(%2d) => counter(%2d), pointer(%s), stack(%2d)" \
                        % (PKndx, self.counter[PKndx], \
                           self.pointer[PKndx], self.stack[PKndx]))
                #
                
                self.pointer[PKndx] \
                    = self.find_next_Xpoint_n(self.Xlist[PKndx],
                                              self.pointer[PKndx], "PK", False)
                self.Xlist[PKndx][self.pointer[PKndx]-1].put_ssPair_j(i, "pk")
                
                if debug_PK:
                    print (self.Xlist[PKndx][self.pointer[PKndx]-1].disp_Pair())
                    print ('j = %5d, %s' % (i, s))
                #
                
                self.pointer[PKndx] -= 1
                self.counter[PKndx] -= 1
                mm = i + 1
                self.scan_allTypes(mm, layer + 1)
                
                # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                # ########################################
                # expanded PK notations
                # ########################################
                
            else:
                # exit out the non-interacting point
                mm = i
                while (s == '.') and (mm < self.N):
                    if debug_bp or debug_PK or debug_ctcf:
                        print ('i = %5d, %s' % (mm, s))
                    #
                    
                    mm += 1
                    if mm == self.N:
                        break
                    #
                    
                    s = self.vstr[mm]
                #|endwhile
                
                self.scan_allTypes(mm, layer + 1)
                
            #
            
        else:
            
            # ######################################################
            # ######  When reach this point, means we are finished
            # ######  searching the sequence
            # ######  ##############################################
            
            #print ("came here!!!")
            
            if self.use_old_build:
                # secondary structure contacts 
                self.BPlist = self.Xlist[0]
                self.BPlist = self.sortvsList(self.BPlist, 'i')
                # self.BPlist = sorted(self.Xlist[0], key=self.getKey_i)
                self.BProots = self.findroots(self.BPlist, False)
                
                # post base pair list processing
                
                # pseudoknot contacts 
                self.PKlist = []
                pkndx = PKrndx.keys()
                for x_k in pkndx:
                    self.PKlist += self.Xlist[x_k]
                #|endfor
                
                self.PKlist = self.assign_PKdirection(self.PKlist)
                self.PKlist = self.sortvsList(self.PKlist, 'i')
                self.compress_to_BPlist()
            else:
                
                # build initial input BPlist
                self.BPlist = self.Xlist[0]
                self.BPlist = self.sortvsList(self.BPlist, 'i')
                
                # build intial input PKlist
                # pseudoknot contacts 
                self.PKlist = []
                pkndx = PKrndx.keys()
                for x_k in pkndx:
                    self.PKlist += self.Xlist[x_k]
                #|endfor
                
                self.build1DStructureLists(self.BPlist, self.PKlist)
            #
            
            if debug_bp:
                print ("BPlist: ")
                for nn in self.BPlist:
                    print (nn.disp_Pair())
                #|endfor
                
                print ("BProots: ")
                for nn in self.BProots:
                    print (nn.disp_Pair())
                #|endfor
                
                print ("PKlist: ")
                for nn in self.PKlist:
                    print (nn.disp_Pair())
                #|endfor
                
                print ("PKroots: ")
                for nn in self.PKroots:
                    print (nn.disp_Pair())
                #|endfor
                
                if self.use_old_build:
                    print ("stop at 1a: scan_allTypes (after compress_to_BPlist)")
                    sys.exit(0)
                else:
                    print ("stop at 1b: scan_allTypes (after build1DStructureLists)")
                    sys.exit(0)
                #
                
            #
            
            # MP contacts; e.g., CTCFs
            self.MPlist = self.Xlist[2]
            for cl in range(0, len(self.MPlist)):
                self.MPlist[cl].v = '-'
            #|endfor
            
            # display the layout of results
            if debug_bp or debug_PK or debug_ctcf:
                print ("index   pointer    counter       stack")
                for x_k in self.Xlist.keys():
                    print (" %2d      %3d        %3d        %3d" \
                        % (x_k, self.pointer[x_k], self.counter[x_k], self.stack[x_k]))
                #|endfor
                
            #
            
            for x_k in self.Xlist.keys():
                if not self.counter[x_k] == 0:
                    case = 0
                    if self.counter[x_k] > 0:
                        case = 0
                    else:
                        case = 1
                    #
                    print ('counter[%d] = %d' %  (x_k, self.counter[x_k]))
                    print ('ERROR!!! Fontana notation is not correct!')
                    if case == 0:
                        print ('         %d too many \'%s\' brackets!' \
                            % (self.counter[x_k], num2lpr[x_k], num2rpr[x_k]))
                    else:
                        print ('         %d too many \'%s\' brackets!' \
                            % (-self.counter[x_k], num2lpr[x_k], num2rpr[x_k]))
                    #
                    print ('input structure:')
                    print (self.vstr)
                    sys.exit(1)
                #
                
            #|endfor
            
            for x_k in self.Xlist.keys():
                if not self.check_Xlist_n(self.Xlist[x_k]):
                    print ("ERROR!!! Fontana notation is not correct!")
                    print ("         At least one structure of type \'%s...%s\' was found" \
                        % (num2rpr[x_k], num2lpr[x_k]))
                    print ("         i.e., the order should be %s...%s." \
                        % (num2lpr[x_k], num2rpr[x_k]))
                    print ('input structure:')
                    print (self.vstr)
                    sys.exit(1)
                #
                
            #|endfor
            
            if debug_bp or debug_PK or debug_ctcf:
                print ("mostly finished scan_allTypes(layer=%2d, ndx=%4d, N=%4d" \
                    % (layer, i, self.N))
                # print ("stop at 2 in scan_allTypes"); sys.exit(0)
            #
        #
        
        if debug_bp or debug_PK or debug_ctcf:
            print ("exit scan_allTypes(%d):" % i, self.pointer[0])
            # print ("stop at 3 (exit from) scan_allTypes"); sys.exit(0)
        #
        
        return mm
    #
    
    # pknotes = self.search_outPKstems(pknotes, k+2, jref)
    def checkForDmnRoot(self, k, jpk):
        debug_ad = False # True # 
        
        if debug_ad:
            print ("Enter checkForDmnRoot k=%d, jpk=%d:" % (k, jpk))
        #
        
        n = len(self.BPlist)
        flag_cont = True
        dmn_size = 0
        dmnroot = None
        while (k < n-1) and flag_cont:
            ibp0 = self.BPlist[k  ].i; jbp0 = self.BPlist[k  ].j
            ibp1 = self.BPlist[k+1].i; jbp1 = self.BPlist[k+1].j
            
            if debug_ad:
                print ("ibp0(%2d) < ibp1(%2d) < jbp0(%2d) < jbp1(%2d)" \
                    % (ibp0, ibp1, jbp0, jbp1))
                print ("ibp0(%2d) < jpk(%2d) < jbp0(%2d)" \
                    % (ibp0, jpk, jbp0))
            #
            
            if jpk < ibp0:
                if debug_ad:
                    print ("jpk(%2d) < ibp0(%2d)" \
                           % (jpk, ibp0))
                #
                
                flag_cont = False
            #
            
            if ibp0 < jpk < jbp0:
                
                # In think in an ordered list, the first item
                # encountered should be the largest domain.
                
                test_dmn_size = jbp0 - ibp0 + 1
                # print ("test_dmn_size: ", test_dmn_size)
                if test_dmn_size > dmn_size:
                    dmnroot = self.BPlist[k] # replace with the largest domain
                    dmn_size = test_dmn_size
                    # print ("assign dmnroot: ", dmnroot)
                #
            #
            
            k += 1
            
        #|endfor
        
        if debug_ad:
            
            if not dmnroot == None:
                print ("finished checkForDmnRoot -> dmnroot:", dmnroot.disp_Pair())
            else:
                print ("finished checkForDmnRoot: nothing found")
            #
            
            #print ("stop at end of checkForDmnRoot"), sys.exit(0)
        #
        
        #print (self.PKlist)
        return dmnroot
        
    #
    
    
    def show_pxlist(self, pxlist, title = "pxlist:"):
        # sometimes needed for debugging
        print (title)
        print (" k   ----- pair -----")
        for k in range(len(pxlist)):
            print ("%3d  %s" % (k, pxlist[k].disp_Pair()))
        #
        
    #
    
    
    def build1DStructureLists(self, inBPlist, inPKlist):
        # inBPlist, inPKlist: list of objects of class Pair
        
        """@
        
        build1DStructureLists: Sets up a standard form for all the
        pairs by distributing them into the two lists BPlist(BProots)
        and PKlist(PKroots) and decides if the stem is parallel or
        antiparallel.
        
        Here are a couple of examples of tests runs 
        
           > Vienna.py -mstr <structure>
        or 
           > Vienna2TreeNode.py -t test1 -mstr <structure>
        
        
        example 1
                   0         10        20        30        40 
                   |    .    |    .    |    .    |    .    |  
        ss_seq  = ".A..B..C..D...........a..b..c..d.."
        result:
        (   1,   22)[p]
        (   4,   25)[p]
        (   7,   28)[p]
        (  10,   31)[p]
        
        example 2:
                   0         10        20        30        40 
                   |    .    |    .    |    .    |    .    |  
        ss_seq  = ".AAA...BBB...CCC...aaa....bbb....ccc...."
        input:
        (   1,   21)[a]
        (   2,   20)[a]
        (   3,   19)[a]
        (   7,   28)[a]
        (   8,   27)[a]
        (   9,   26)[a]
        (  13,   35)[a]
        (  14,   34)[a]
        (  15,   33)[a]
        
        Results:
        structural sequence
        cxxxcccxxxcccxxxcccyyycccyyycccyyyc
        .AAA...BBB...CCC...aaa...bbb...ccc.
        secondary structure
        (    1,   21)[a]
        (    2,   20)[a]
        (    3,   19)[a]
        secondary structure roots
        (    1,   21)[a]
        pseudoknot linkages: 
        (    7,   27)[a]
        (    8,   26)[a]
        (    9,   25)[a]
        (   13,   33)[a]
        (   14,   32)[a]
        (   15,   31)[a]
        pseudoknot linkage roots
        (    7,   27)[a]
        (   13,   33)[a]
        
        example 3:
                   0         10        20        30        40 
                   |    .    |    .    |    .    |    .    |  
        ss_seq  = ".AAA...BBB...aaa...CCC....bbb....ccc...." # same
        ss_seq  = ".(((...BBB...)))...(((....bbb....)))...." # same
        Input:
        (   1,   15)[a]
        (   2,   14)[a]
        (   3,   13)[a]
        (   7,   28)[a]
        (   8,   27)[a]
        (   9,   26)[a]
        (  19,   35)[a]
        (  20,   34)[a]
        (  21,   33)[a]
        
        Results:
        structural sequence
        cxxxcccxxxcccyyycccxxxcccyyycccyyyc
        .AAA...BBB...aaa...CCC...bbb...ccc. 
    ==> .(((...[[[...)))...(((...]]]...))). <==
        secondary structure
        (    1,   15)[a]
        (    2,   14)[a]
        (    3,   13)[a]
        (   19,   33)[a]
        (   20,   32)[a]
        (   21,   31)[a]
        secondary structure roots
        (    1,   15)[a]
        (   19,   33)[a]
        pseudoknot linkages: 
        (    7,   27)[a]
        (    8,   26)[a]
        (    9,   25)[a]
        pseudoknot linkage roots
        (    7,   27)[a]
        (   19,   33)[a]
        
        It appears to be able to distinguish consecuitive pairs that
        would satisfy a parallel pattern and they do not have to be
        contiguous either (as shown in example 1). However, if it is
        just an odd sort of pseudoknot where you have antiparallel
        stems but a kind of ladder structure, it will not call this
        ladder arrangement "parallel" (assigning the "p" to the
        brackets [a] or [p].
        
        190219: appears to work ok with the more odd or difficult
        examples I could think of.
        
        """
        debug_b1Dstr = False # True # 
        
        # mash inBPlist and inPKlist together list into plist
        pxlist = deepcopy(inBPlist) + deepcopy(inPKlist)
        pxlist = self.sortvsList(pxlist, 'i')
        self.pxlist_ref = deepcopy(inBPlist) + deepcopy(inPKlist)
        
        
        if self.flag_show_pxlist:
            # sometimes needed for debugging. This leaves the option
            # that the raw input is shown when testing but other
            # debugging is not turned on.
            self.show_pxlist(pxlist, "input pxlist:")
        #
        
        if debug_b1Dstr:
            print ("Enter build1DStructureLists: ")
            print ("inBPlist: (before)")
            for vv in inBPlist:
                print (vv.disp_Pair())
            #|endfor
            
            print ("inBPlist: (before)")
            for vv in inPKlist:
                print (vv.disp_Pair())
            #|endfor
            
            self.show_pxlist(pxlist, "build1DStructureLists pxlist at 0.0:")
            
            #print ("stop at 0.0 in build1DStructureLists"); sys.exit(0)
            
        #
        
        #  prroots = self.findroots(prlist, debug_b1Dstr)
        
        # findPrincipalDMNs tries to organize the table so that it is
        # BProots and PKroots. It is a more sophisticated version of
        # findroots that actually does this subdivision.
        
        self.findPrincipalDMNs(pxlist)
        # builds
        # self.BProots -> self.BPlist,
        # self.PKroots -> self.PKlist
        
        if debug_b1Dstr:
            print ("finished findPrincipalDMNs")
            #print ("stop at 0.1 in build1DStructureLists"); sys.exit(0)
        #
        
        self.BPlist  = self.assign_ppStems(self.BPlist, debug_b1Dstr)
        self.BProots = self.assign_ppRoots(self.BPlist, self.BProots)
        self.PKlist  = self.assign_ppStems(self.PKlist, debug_b1Dstr)
        self.PKroots = self.assign_ppRoots(self.PKlist, self.PKroots)
        
        if debug_b1Dstr:
            
            print (self.disp_allLists("finished assign_ppStems/Roots", False))
            
            print ("inBPlist: (before)")
            for vv in inBPlist:
                print (vv.disp_Pair())
            #
            
            print ("inBPlist: (before)")
            for vv in inPKlist:
                print (vv.disp_Pair())
            #|endfor
            
            
            print ("dXPKinfo:")
            #print (self.dXPKlinfo)
            for vk in list(self.dXPKlinfo):
                pk = self.dXPKlinfo[vk]
                lk = pk.linkage[0]
                rt1 = pk.ref_root[0]
                rt2 = pk.ref_root[1]
                print ("linkage: ", self.prncpl_wts[lk])
                print ("root1  : ", self.prncpl_wts[rt1])
                print ("root2  : ", self.prncpl_wts[rt2])
                
            #
            
            print ("dCPKlinfo:")
            for vk in list(self.dCPKlinfo):
                print (self.dCPKlinfo[vk])
            #
            
            #print ("stop at 0.2 in build1DStructureLists"); sys.exit(0)
        #
        
        self.findPKboundaries()
        #debug_b1Dstr = True
        if debug_b1Dstr:
            
            
            
            print ("BPlist: (after ap/pp assignment)")
            for vv in self.BPlist:
                print (vv.disp_Pair())
            #|endfor
            
            print ("BProots: (after ap/pp assignment)")
            for vv in self.BProots:
                print (vv.disp_Pair())
            #|endfor
            
            print ("PKlist: (after ap/pp assignment)")
            for vv in self.PKlist:
                print (vv.disp_Pair())
            #|endfor
            
            print ("PKroots: (after ap/pp assignment)")
            for vv in self.PKroots:
                print (vv.disp_Pair())
            #|endfor
            
            print ("dXPKinfo:")
            for v in list(self.dXPKinfo.keys()):
                print (self.dXPKinfo[v])
            #
            
            print ("dCPKinfo:")
            for v in list(self.dCPKinfo.keys()):
                print (self.dCPKinfo[v])
            #
            
            print ("stop at end of build1DStructureLists"), sys.exit(0)
        #
        
    #
    
    
    def findPKboundaries(self):
        debug = False # True # 
        
        if debug:
            print ("Enter findPKboundaries()")
            print (self.dXPKlinfo)        
            print (self.dCPKlinfo)
            #print ("stop at 1 in findPKboundaries"); sys.exit(0)
        #
        
        """@ 
        
        First, we have to figure out what the true "root stem"
        actually is.  This step builds the basic pseudoknot
        information based on the linkage stmes. 
        
        builds
        dCPKlinfo (Core PseudoKnot Linkage INFOrmation)
        dXPKlinfo (eXtended PseudoKnot Linkage INFOrmation)
        
        """
        
        
        for vv in list(self.PKroots):
            
            vk = (vv.i, vv.j)
            if vk in self.dXPKlinfo:
                #print ("dXPKinfo:")
                #print (self.dXPKlinfo[vk])
                
                pk = self.dXPKlinfo[vk]
                lk = pk.linkage[0]
                rt1 = pk.ref_root[0]
                rt2 = pk.ref_root[1]
                #print ("prncpl_wts", rt1, self.prncpl_wts[rt1].rpntr)
                for rv in self.prncpl_wts[rt1].rpntr:
                    iv = rv[0]; jv = rv[1]
                    if iv <= rt1[0] and rt1[0] <= jv and jv < lk[1]:
                        
                        pk.r1levels += [rv]
                        """@
                        
                        example:
                        
                          iv3                                           jv3
                           |           |         <----------|            |
                           |           il                  jl            |
                         ..(...(....(..[...)...)....(...(...]..)...).....)...
                              iv2  iv1    jv1 jv2  irt2          jrt2   
                               A               A
                               |               |
                        
                        Hence pk.r1levels = [(iv2, jv2), (iv1, jv1)].  
                        
                        The base (iv3,jv3) encompasses both root1,
                        root2 and the linkage stem, therefore is
                        excluded. As a result, rlevels will contain
                        (iv1, jv1) and (iv2, jv2).
                        
                        """
                                                
                    #
                    
                #|endfor
                
                if debug:
                    print (pk)
                    print (pk.r1levels)
                #
                
                rt1 = pk.r1levels[0] # update rt1
                pk.root += [rt1]
                
                
                
                for rp in self.BProots:
                    ip = rp.i; jp = rp.j
                    if rt1[1] < ip and ip < lk[1] and lk[1] < jp:
                        
                        pk.r2levels += [(ip, jp)]
                        #print ("ijp(%2d,%2d), r2levels = %s" % (ip, jp, pk.r2levels))
                        
                        """@
                        
                        example:
                        
                                                                          
                           |           |<---  linkage  ---->|            |
                           |           il                  jl            |
                         ..(...(....(..[...)...)....(...(...]..)...).....)...
                              rt1i            rt1j ip2 ip1    jp1 jp2   
                                                    A              A
                                                    |              |
                        
                        Hence pk.r2levels = [(ip2, jp2), (ip1, jp1)].  
                        
                        Like the previous example, rlevels will
                        contain (ip1, jp1) and (ip2, jp2).
                        
                        """
                        
                    #
                    
                #|endfor
                
                if debug:
                    print (pk)
                    print (pk.r2levels)
                #
                
                rt2 = pk.r2levels[0] # update rt1
                pk.root += [rt2]
                pk.update_PKdmn()
                
                if debug:
                    print ("linkage: ", self.prncpl_wts[lk])
                    print ("root1  : ", self.prncpl_wts[rt1])
                    print ("root2  : ", self.prncpl_wts[rt2])
                    print (pk.r1levels)
                    print (pk.r2levels)
                    print (pk.root)
                #
                
                self.dXPKlinfo[vk] = pk
                
                #print (self.dXPKlinfo[vk])
                #sys.exit(0)
            elif vk in self.dCPKlinfo:
                if debug:
                    print ("dCPKlinfo:")
                    print (self.dCPKlinfo[vk])
                #
                
                pk = self.dCPKlinfo[vk]
                lk = pk.linkage[0]
                rt1 = pk.ref_root[0]
                #print (self.prncpl_wts[rt1].rpntr)
                for rv in self.prncpl_wts[rt1].rpntr:
                    #print ("rv -> ", self.prncpl_wts[rv])
                    
                    iv = rv[0]; jv = rv[1]
                    if debug:
                        print ("iv(%d) <= rt1(%d) <= jv(%d) < lk.j(%d)?" \
                               % (iv, rt1[0], jv, lk[1]))
                    #
                    
                    if iv <= rt1[0] and rt1[1] <= jv and jv < lk[1]:
                        
                        pk.r1levels += [rv]
                        #print (rv, "pk.r1levels", pk.r1levels)
                        """@
                        
                        example:
                        
                           
                                       |<----------->|      
                                      lk.i          lk.j     
                         ......(....(..[...)...).....].....
                              iv2  iv1    jv1 jv2   
                               A               A
                               |               |
                              iv              jv
                        
                        Hence pk.r1levels = [(iv2, jv2), (iv1, jv1)].  
                        
                        
                        """
                                                
                        # if iv <= rt1[0] and rt1[1] <= jv and jv < lk[1]:
                        
                    elif lk[0] < iv and iv <= rt1[0] and rt1[1] <= jv:
                        
                        pk.r1levels += [rv]
                        #print (rv, "pk.r1levels", pk.r1levels)
                        """@
                        
                        example:
                        
                           
                          |<---------->|      
                          lk.i         lk.j     
                         .[....(....(..]...)...).....].....
                              iv2  iv1    jv1 jv2   
                               A               A
                               |               |
                              iv              jv
                        
                        Hence pk.r1levels = [(iv2, jv2), (iv1, jv1)].  
                        
                        
                        """
                    #
                    
                #|endfor
                
                if debug:
                    print (pk)
                    print (pk.r1levels)
                #
                
                rt1 = pk.r1levels[0] # update rt1
                pk.root += [rt1]
                pk.update_PKdmn()
                
                if debug:
                    print ("linkage: ", self.prncpl_wts[lk])
                    print ("root1  : ", self.prncpl_wts[rt1])
                    print (pk.r1levels)
                    print (pk.root)
                #
                
                self.dCPKlinfo[vk] = pk
                
                if debug: print (self.dCPKlinfo[vk])
                #sys.exit(0)
            else:
                print ("somethings wrong in findPKboundaries")
                sys.exit(1)
                
            #
            
        #
        
        if debug:
            
            print ("results from building dCPKlinfo and dXPKlinfo")
            
            print (self.dXPKlinfo)
            for vk in list(self.dXPKlinfo):
                print (vk, self.dXPKlinfo[vk].r1levels)
                print (vk, self.dXPKlinfo[vk].r2levels)
            #
            
            print (self.dCPKlinfo)
            for vk in list(self.dCPKlinfo):
                print (vk, self.dCPKlinfo[vk].r1levels)
            #
            
            #print ("stop at 2 in findPKboundaries"); sys.exit(0)
            print ("next step: build a dictionary of root stems")
        #
        
        """@
        
        In this step, we need to build a reference system where the
        root stem points to the linkages, the inverse of above where
        PKdmn helps us find the root stems. """
        
        self.dCPK_r2l  = {}
        self.dXPK_xp2l = {}
        self.dXPK_xp2r = {}
        for vv in list(self.PKroots): # scan the linkages
            
            vk = (vv.i, vv.j) # linkage 
            if vk in list(self.dXPKlinfo):
                pk = self.dXPKlinfo[vk]
                rt1 = pk.r1levels[0]
                rt2 = pk.r2levels[0]
                ir1 = rt1[0]; jr2 = rt2[1]
                #print (rt1, rt2, vk)
                
                if (ir1, jr2) in self.dXPK_xp2l:
                    self.dXPK_xp2l[(ir1, jr2)] += [vk]
                else:
                    self.dXPK_xp2l.update({(ir1, jr2) : [vk] })
                #
                
                if (ir1, jr2) in self.dXPK_xp2r:
                    self.dXPK_xp2r[(ir1, jr2)] += [rt1, rt2]
                else:
                    self.dXPK_xp2r.update({(ir1, jr2) : [rt1, rt2] })
                #
                
                
            elif vk in self.dCPKlinfo:
                pk = self.dCPKlinfo[vk]
                rt1 = pk.r1levels[0]
                
                # for root key, store linkage
                if rt1 in self.dCPK_r2l:
                    self.dCPK_r2l[rt1] += [vk]
                else:
                    self.dCPK_r2l.update({rt1 : [vk] })
                #
            
            else:
                print ("somethings wrong in findPKboundaries")
                sys.exit(1)
                
            #
            
        #
        
        if debug:
            print ("core pk linkages wrt root domains")
            for vr in list(self.dCPK_r2l.keys()):
                print (vr, self.dCPK_r2l[vr])
            #
            
            print ("extended pk linkages wrt root domains")
            for vr in list(self.dXPK_xp2l.keys()):
                print (vr, self.dXPK_xp2l[vr])
            #
            
            print ("extended pk roots")
            for vr in list(self.dXPK_xp2r.keys()):
                print (vr, self.dXPK_xp2r[vr])
            #
            
            #print ("stop at 3 in findPKboundaries"); sys.exit(0)
            
        #
        
        
        # ########################################################
        # ####   Cluster the contents around the root stems.  ####
        # ########################################################
        
        # First, we built the final core pseudoknot
        
        self.dCPKinfo = {}
        for vk in list(self.dCPK_r2l):
            ir = vk[0]; jr = vk[1]
            if len(self.dCPK_r2l[vk]) > 1:
                lks = self.dCPK_r2l[vk]
                pk = deepcopy(self.dCPKlinfo[lks[0]])
                fpntr = get_CPK_jbranches(ir, jr, lks[0][0], lks[0][1], self.BProots)
                
                for k in range(1, len(lks)):
                    pk.linkage += [lks[k]]
                    fpntr += get_CPK_jbranches(ir, jr, lks[k][0], lks[k][1], self.BProots)
                    
                #|endfor
                
                fpntr = prune_pairlist(fpntr)
                pk.jbranches = compress_to_jbranches(fpntr)
                
                #print (pk)
                pk.update_PKdmn()
                #print (pk, pk.jbranches)
                self.dCPKinfo.update({(pk.ipk, pk.jpk) : pk })
                
                #sys.exit(0)
            else:
                lks = self.dCPK_r2l[vk]
                pk = deepcopy(self.dCPKlinfo[lks[0]])
                self.dCPKinfo.update({(pk.ipk, pk.jpk) : pk })
                pk = deepcopy(self.dCPKlinfo[lks[0]])
                fpntr = get_CPK_jbranches(ir, jr, lks[0][0], lks[0][1], self.BProots)
                pk.jbranching = compress_to_jbranches(fpntr)
                #print (pk)
                pk.update_PKdmn()
                #print (pk, pk.jbranches)
                self.dCPKinfo.update({(pk.ipk, pk.jpk) : pk })
            #
            
            
        #
        
        if debug:
            print ("4a")
            print ("core pk linkages wrt root domains")
            for cpk in list(self.dCPKinfo.keys()):
                print (self.dCPKinfo[cpk], self.dCPKinfo[cpk].jbranches)
            #
            print ("----------------------")
            
            print ("\nextended pk linkage info:")
            for vk in self.dXPKlinfo:
                print (self.dXPKlinfo[vk])
            #
            
            #print ("stop at 4a in findPKboundaries"); sys.exit(0)
            
        #
        
        
        
        """@
        
        The matter seems to be a lot more complicated for extended
        pseudoknots, requiring two essential parts to get to the end
        of it.
        
        In this first step, we have to look at the root stems. If they
        are themselves core PKs, we should also include this in the in
        the final construction.
        
        I have opted to place this step here so that the root stem is
        built up into the final linkage information based solely on
        the linkage information. This can then be written into a
        tandem PK.  """
        
        # wrap up extended PK info
        for vk in list(self.dXPKlinfo):
            pk = self.dXPKlinfo[vk]
            r1 = pk.root[0]; r2 = pk.root[1]
            
            if r1 in self.dCPK_r2l:
                
                for lk in self.dCPK_r2l[r1]:
                    
                    if debug:
                        print ("r1(%3d,%3d) -> lk(%3d,%3d): " \
                               % (r1[0], r1[1], lk[0], lk[1]))
                        print ("r1.i(%d) < l.i(%d) && r1.j(%d) < l.j(%d) && l.j(%d) < r2.i(%d)" \
                               % (r1[0], lk[0], r1[1], lk[1], lk[1], r2[0]))
                    #
                    
                    if r1[0] < lk[0] and r1[1] < lk[1] and lk[1] < r2[0]:
                        
                        """@
                                  vk
                           |<--------------|
                           |               |
                        (..<..[..)..]..(...>..)
                        ir1   | jr1 | ir2    jr2
                              |<--->|
                                 lk
                        """
                        
                        if lk in self.dCPKlinfo:
                            #print ("add r1 C", self.dCPKlinfo[lk])
                            pk.add_pk_to_root1(self.dCPKlinfo[lk])
                        else:
                            #print ("add r1 X", self.dXPKlinfo[lk])
                            pk.add_pk_to_root1(self.dXPKlinfo[lk])
                        #
                        
                    #
                    
                #
                
            #
            
            if r2 in self.dCPK_r2l:
                
                for lk in self.dCPK_r2l[r2]:

                    if debug:
                        print ("r2(%3d,%3d) -> lk(%3d,%3d): "
                               % (r2[0], r2[1], lk[0], lk[1]))
                    #
                    
                    
                    if r2[0] < lk[0] and r2[1] < lk[1]:
                        
                        
                        """@
                                  vk
                           |<--------------|
                           |               |
                        (..<..[..)..]..(..A>..)....a
                        ir1     jr1   ir2 |   jr2  |
                                          |<------>|
                                              lk
                        """
                        if lk in self.dCPKlinfo:
                            #print ("add r2 C", self.dCPKlinfo[lk])
                            pk.add_pk_to_root2(self.dCPKlinfo[lk])
                        else:
                            #print ("add r2 X", self.dXPKlinfo[lk])
                            pk.add_pk_to_root2(self.dXPKlinfo[lk])
                        #
                        
                        
                    #
                    
                #
                
            #
            
            #print (pk)
            self.dXPKlinfo[vk] = pk
            #print ("stop at vk=(%d,%d) 4a -> 4b" % (vk[0], vk[1])); sys.exit(0)
            
        #
        
        
        if debug:
            print ("4b")
            for vk in list(self.dXPKlinfo):
                print (vk, self.dXPKlinfo[vk])
            #
            
            #print ("stop at 4b in findPKboundaries"); sys.exit(0)
        #
        
        """@
        
        Now we need to take care of tandem extended PKs, We have to
        combine together all the tandem parts together so that the
        final assembly is of the object of class PseudoKnot can be
        constructed.
        
        """
        
        self.dXPKinfo = {}
        pxkeys = sortPairListWRT_n(list(self.dXPK_xp2r.keys()), 0)
        #print (pxkeys)
        
        k = 0;
        ir1 = self.N; jr2 = 0
        
        while k < len(pxkeys):
            rk = self.dXPK_xp2r[pxkeys[k]]
            ir1k = rk[0][0]; jr2k = rk[1][1]
            
            ir1 = ir1k; jr2 = jr2k
            
            #print ("ir1(%d),jr2(%d)" % (ir1, jr2))
            
            roots    = self.dXPK_xp2r[(ir1, jr2)]
            linkages = self.dXPK_xp2l[(ir1, jr2)]
            pk = deepcopy(self.dXPKlinfo[linkages[0]])
            pk.dpkstack.update({linkages[0]: deepcopy(pk) })
            #print ("roots: ", roots, ", linkages: ", linkages)
            if len(linkages) > 1:
                for k in range(1, len(linkages)):
                    pk.linkage += [linkages[k]]
                    pk.dpkstack.update({linkages[k]:
                                        deepcopy(self.dXPKlinfo[linkages[k]])})
                #|endfor
                
            #
            
            # have to use pk(ipk,jpk) because the boundaries may not
            # be just the root stems but pseudoknots comprising the
            # root stem.
            self.dXPKinfo.update({(pk.ipk, pk.jpk) : pk })
            #print (self.dXPKinfo)
            
            if k == len(pxkeys) - 1:
                if debug:
                    print ("end ir1(%d),ir2(%d) -> pk(%d,%d)" \
                           % (ir1, jr2, pk.ipk, pk.jpk))
                #
                
            #
            
            #sys.exit(0)
            # each new cycle of k, we set (ir1,jr2) = (ir1k, jr2k)
            
            l = k + 1; flag_cont = True
            
            while l < len(pxkeys) and flag_cont:
                rl = self.dXPK_xp2r[pxkeys[l]]
                ir1l = rl[0][0]; jr2l = rl[1][1]
                
                if debug:
                    print ("test ir1(%d),ir2(%d) vs jr1l(%d),jr2l(%d)" \
                           % (ir1, jr2, ir1l, jr2l))
                #
                
                if jr2 < ir1l:
                    #print ("1 out of range jr2k(%d) < ir1l(%d)" % (jr2k, ir1l))
                    flag_cont = False
                    
                elif ir1 <= ir1l and jr2 < jr2l:
                    if debug:
                        print ("2 satisfied ir1(%d) <= ir1l(%d) < jr2(%d) <= jr2l(%d)" \
                               % (ir1, ir1l, jr2, jr2l))
                    #
                    
                    
                    ir1_o = ir1; jr2_o = jr2
                    pk_o     = deepcopy(self.dXPKinfo[(ir1_o, jr2_o)])
                    roots_o  = pk_o.root
                    
                    # must expand the size of the PK, update and
                    # delete the old reference.
                    jr2 = jr2l
                    roots    = [( roots_o[0][0], roots_o[0][1] )]
                    roots_n  = self.dXPK_xp2r[(ir1l, jr2l)]
                    roots   += [( roots_n[1][0], roots_n[1][1] )]
                    pk_o.root = roots
                    #print ("assign new pk(%2d,%2d)" % (ir1,   jr2))
                    linkages = self.dXPK_xp2l[(ir1l, jr2l)]
                    #print ("roots: ", roots, ", linkages: ", linkages)
                    pk_o.linkage += [linkages[0]]
                    pk_o.dpkstack.update({linkages[0]:
                                          deepcopy(self.dXPKlinfo[linkages[0]])})
                    if len(linkages) > 1:
                        for k in range(1, len(linkages)):
                            pk_o.linkage += [linkages[k]]
                            pk_o.dpkstack.update({linkages[k]:
                                                  deepcopy(self.dXPKlinfo[linkages[k]])})
                        #|endfor
                        
                    #
                    
                    #if debug: print ("update_PKdnm vvvvv")
                    pk_o.update_PKdmn() # sets to the new boundaries
                    #if debug: print ("update_PKdnm ^^^^^")
                    
                    #self.dXPKinfo.update({(ir1, jr2) : deepcopy(pk_o) })
                    
                    # I think I need to use pk_o(ipk, jpk) here!
                    self.dXPKinfo.update({(pk_o.ipk, pk_o.jpk) : deepcopy(pk_o) })
                    #print ("delete old pk(%2d,%2d)" % (ir1_o, jr2_o))
                    del self.dXPKinfo[(ir1_o, jr2_o)]
                    
                    #print (self.dXPKinfo)
                    l += 1
                    #sys.exit(0)
                    
                    
                else:
                    
                    if debug:
                        print ("2 subset ir1(%d) <= jr2l(%d) < jr1l(%d) <= jr2(%d)" \
                               % (ir1, jr2l, ir1l, jr2))
                        print ("assign new pk(%2d,%2d)" % (ir1,   jr2))
                    #
                    
                    lnks = self.dXPK_xp2l[(ir1l, jr2l)]
                    #print ("roots: ", roots, ", linkages: ", lnks)
                    v = (ir1, jr2)
                    self.dXPKinfo[v].linkage += [lnks[0]]
                    self.dXPKinfo[v].dpkstack.update({lnks[0]:
                                   deepcopy(self.dXPKlinfo[lnks[0]]) })
                    if len(lnks) > 1:
                        for k in range(1, len(lnks)):
                            self.dXPKinfo[v].linkage += [lnks[k]]
                            self.dXPKinfo[v].dpkstack.update({lnks[k]:
                                   deepcopy(self.dXPKlinfo[lnks[k]]) })
                        #|endfor
                        
                    #
                    
                    #print ("update_PKdnm vvvvv")
                    self.dXPKinfo[(ir1, jr2)].update_PKdmn()
                    #print ("update_PKdnm ^^^^^")
                    
                    if debug: print (self.dXPKinfo)
                    
                    l += 1
                    #sys.exit(0)
                #
                
            #|endwhile
            
            k = l
            
        #|endwhile
        
        
        # Now we finalize the junction branches on the self.dXPKinfo.
        
        #print("self.dXPKinfo:", self.dXPKinfo)
        for xpk in list(self.dXPKinfo.keys()):
            r = self.dXPKinfo[xpk].get_roots()
            #print (r, ", vs ", self.dXPKinfo[xpk].root)
            ir1 = r[0][0]; jr1 = r[0][1]
            ir2 = r[1][0]; jr2 = r[1][1]
            #print (self.dXPKinfo[xpk])
            #print (self.dXPKinfo[xpk].disp_dpkstack())
            jbranches = get_XPK_jbranches(ir1, jr1, ir2, jr2, self.BProots)
            self.dXPKinfo[xpk].jbranches = compress_to_jbranches(jbranches)
            #print (self.dXPKinfo[xpk].jbranches)
        #
        
        
        
        if debug:
            print ("4c")
            print ("extended pk wrt root domains")
            for xpk in list(self.dXPKinfo.keys()):
                print (xpk, self.dXPKinfo[xpk])
                print (self.dXPKinfo[xpk].disp_dpkstack())
                print (self.dXPKinfo[xpk].jbranches)
            #
            
            print ("\nextended pk linkage info:")
            for vk in self.dXPKlinfo:
                print (self.dXPKlinfo[vk])
            #
            
            #print ("stop at 4c in findPKboundaries"); sys.exit(0)
            
        #
        
        #debug = True
        if debug:
            print ("\n### final tally:")
            print ("core pk linkages wrt root domains:")
            for vr in list(self.dCPK_r2l.keys()):
                print (vr, self.dCPK_r2l[vr])
            #
            
            print ("\nextended pk linkages wrt root domains:")
            for vr in list(self.dXPK_xp2l.keys()):
                print (vr, self.dXPK_xp2l[vr])
            #
            
            print ("\nfinal core pk linkage info:")
            for vk in self.dCPKlinfo:
                print (self.dCPKlinfo[vk])
            #
            
            print ("\nfinal extended pk linkage info:")
            for vk in self.dXPKlinfo:
                print (self.dXPKlinfo[vk])
            #
            
            print ("\nfinal core pk info:")
            for vk in self.dCPKinfo:
                print (vk, self.dCPKinfo[vk])
            #
            
            print ("\nfinal extended pk info:")
            for vk in self.dXPKinfo:
                print (vk, self.dXPKinfo[vk])
            #
            
            
            print ("Exiting findPKboundaries()")
            #print ("stop at 4 in findPKboundaries"); sys.exit(0)
        #
        
    #
    
    
    
    def assign_ppRoots(self, pxlist, pxroots):
        
        # should be run AFTER assign_ppStems!!!!
        
        # reassigns all Pair.v (pxroots) that are of type 'p' in
        # pxlist.
        
        for pxrk in pxroots:
            irk = pxrk.i; jrk = pxrk.j; vrk = pxrk.v
            for pxlk in pxlist:
                ilk = pxlk.i; jlk = pxlk.j; vlk = pxlk.v
                if irk == ilk and jrk == jlk:
                    if not vlk == vrk:
                        pxrk.v = vlk
                    #
                    
                #
                
            #
            
        #
        
        return pxroots
    #
    
    
    def assign_ppStems(self, pxlist, debug = False):
        
        # reassigns each Pair.v in the stem to 'p' if the stem is
        # clearly parallel
        
        debug_b1Dstr = debug
        #debug_b1Dstr = False # True #
        
        if debug:
            print ("Entered assign_ppStems")
        #
        
        # since the lists are ordered, we can check through
        # self.BPlist and identify parallel stems by the simple rule
        # of parallel-ness: ibp0 < ibp1 < jbp0 < jbp1.
        
        n = len(pxlist)
        pknotes = []
        dmnnotes = []
        zone_lock = False
        i_lock = self.N
        j_lock = 0
        
        for k in range(0, n - 1):
            
            ibp0 = pxlist[k  ].i; jbp0 = pxlist[k  ].j
            ibp1 = pxlist[k+1].i; jbp1 = pxlist[k+1].j
            
            zone_lock = (jbp1 <= j_lock)
            if debug_b1Dstr:
                print ("ibp0(%2d) < ibp1(%2d) < jbp0(%2d) < jbp1(%2d)" \
                    % (ibp0, ibp1, jbp0, jbp1))
                print ("ij_lock(%2d,%2d) = %s" % (i_lock, j_lock, zone_lock))
            #
            
            if ibp0 < ibp1 and ibp1 < jbp0 and jbp0 < jbp1:
                flag_pass = False
                #print ("xxx1")
                test_22_is_pp    = False
                test_p2_is_pp    = False
                test_p2_is_shift = False
                if k <= n-3:
                    ibp2 = pxlist[k+2].i; jbp2 = pxlist[k+2].j
                    # ..(.....[.....<......).....].....>..
                    #   ^     ^     ^      ^     ^     ^
                    # ibp0  ibp1   ibp2  jbp0  jbp1  jbp2
                    test_p2_is_pp = (ibp1 < ibp2 and jbp1 < jbp2 and ibp2 < jbp1)
                    test_p2_is_shift = (jbp1 < ibp2)
                    test_22_is_pp  = (ibp1 - ibp0 < 2 and jbp1 - jbp0 < 2)
                    
                elif k <= n-2:
                    #print ("xxx")
                    test_22_is_pp  = (ibp1 - ibp0 < 2 and jbp1 - jbp0 < 2)
                    #print (test_22_is_pp); sys.exit(0)
                #
                
                test_m1_is_pp    = False
                test_m1_is_shift = False
                if k > 0:
                    # ..(.....[.....<......).....].....>..
                    #   ^     ^     ^      ^     ^     ^
                    # ibm1  ibp0   ibp1  jbm1  jbp0  jbp1
                    ibm1 = pxlist[k-1].i; jbm1 = pxlist[k-1].j
                    test_m1_is_pp = (ibm1 < ibp0 and jbm1 < jbp0 and ibp0 < jbm1)
                    test_m1_is_shift = (jbm1 < ibp0)
                #
                
                if test_m1_is_pp or test_p2_is_pp or test_22_is_pp:
                    # if at least 1 of these is also a parallel stem,
                    # then this passes. In short There must be at
                    # least three of these together to make this work.
                    flag_pass = True
                #
                
                
                # parallel relationship. This is the nearest
                # neighboring Pair and it is in BPlist (the dominant
                # one). Therefore, it should have priority and be
                # assigned properly.
                
                if flag_pass:
                    if debug_b1Dstr:
                        print ("Reassign to p:")
                        print (pxlist[k])
                        print (pxlist[k+1])
                    #
                    
                    pxlist[k  ].v = 'p'
                    pxlist[k+1].v = 'p'
                    
                else:
                    if debug_b1Dstr:
                        print ("zone_lock: ", zone_lock)
                    #
                    
                    if not zone_lock:
                        pknotes += [pxlist[k+1]]
                        iref = pxlist[k+1].i
                        jref = pxlist[k+1].j
                        if debug_b1Dstr:
                            print ("Search for other pks within pk(%d,%d), next k = %d" \
                                   % (iref, jref, k+2))
                        #
                        
                        dmnroot = self.checkForDmnRoot(k + 2, jref)
                        if not dmnroot == None:
                            i_lock = dmnroot.i
                            j_lock = dmnroot.j
                            dmnnotes += [dmnroot]
                        #
                        
                    else:
                        
                        if debug_b1Dstr:
                            print ("zone locked between %d and %d, <- (%d,%d)" \
                                   % (i_lock, j_lock, ibp1, jbp1))
                        #
                        
                    #
                    
                #
                
            #
            
        #|endfor
        
        if debug_b1Dstr:
            self.show_pxlist(pxlist, "assign_ppStems(pxlist) at 1:")
            #print ("stop at 1 in assign_ppStems"); sys.exit(0)
            
            if len(pknotes) > 0:
                
                # I think pknotes is no longer necesary with the newer
                # approach, which is what this rountine was built for.
                
                print ("pknotes:")
                for vv in pknotes:
                    print (vv.disp_Pair())
                #|endfor
                
            #
            
            #print ("stop at 2 in assign_ppStems"); sys.exit(0)
        #
        
        return pxlist
    #
    
    
    def find_rpntr(self, k, ik, jk, vkeys):
        
        
        vk = (ik, jk)
        n = len(vkeys)
        
        l = 0; flag_cont = True
        
        while l < n and flag_cont:
            
            vl = vkeys[l]; il = vl[0]; jl = vl[1]
            if il < ik and jk < jl:
                self.prncpl_wts[vk].rpntr += [vl]
                l += 1
                
            elif jk < il:
                flag_cont = False
                
            else:
                l += 1
            #
            
        #
        
        
        if len(self.prncpl_wts[vk].rpntr) > 1:
            self.prncpl_wts[vk].rpntr = sortPairListWRT_n(self.prncpl_wts[vk].rpntr, 0)
            # sort wrt index 0
        #
            
    #
    
    
    def revise_pointers(self):
        
        # this __should__ be used ___after___ running
        # findMaxStemTails(!!!) 
        
        sskeys = self.sskeylist
        n = len(sskeys)
        
        flag_cont = True
        k = 0
        while k < n and flag_cont:
            vk = sskeys[k]
            ik = vk[0]; jk = vk[1];
            vk = sskeys[k] 
            
            if k == n - 1:
                # an isolated last one means it can only be the stem itself.
                flag_cont = False
            #
            
            l = k + 1; flag_cont_l = True
            
            while l < n and flag_cont_l:
                
                vl = sskeys[l]; il = vl[0]; jl = vl[1]
                if ik < il and jl < jk:
                    self.prncpl_wts[vl].rpntr += [vk]
                    l += 1
                    
                else:
                    flag_cont_l = False
                #
                
            #
            
            k += 1
            
        #
        
        for v in self.sskeylist:
            if len(self.prncpl_wts[v].rpntr) > 1:
                self.prncpl_wts[v].rpntr = sortPairListWRT_n(self.prncpl_wts[v].rpntr, 0)
                # sort wrt index 0
            #
            
        #|endfor
        
        # need to include something to remove the fpntr etc from pk
        for v in self.pkkeylist:
            #self.prncpl_wts[v].fpntr = [] # reset to zero
            self.prncpl_wts[v].npairs = self.prncpl_wts[v].stem
        #|endfor
        
        
    #
    
    def subtract_crossovers(self, fpntr, debug = False):
        n = len(fpntr)
        if debug: print ("n = ", n)
        dmaskPK = {}
        # dmaskPK: used as a place holder to mask potential pk
        # references that have already be marked.
        for k in range(len(fpntr)):
            ik = fpntr[k][0]; jk = fpntr[k][1]
            if (ik, jk) in dmaskPK:
                continue
            #
            
            l = k + 1
            flag_cont = True
            while l < len(fpntr) and flag_cont:
                il = fpntr[l][0]; jl = fpntr[l][1]
                if debug: print ("ijk(%2d,%2d) v ijl(%2d,%2d)" % (ik, jk, il, jl))
                if il > jk:
                    if debug: print ("finished")
                    flag_cont = False
                elif ik < il and jl < jk:
                    if debug: print ("inside")
                    l += 1
                elif ik < il and il < jk and jk < jl:
                    if debug: print ("crossover")
                    dmaskPK.update( { (il, jl) : True } )
                    n -= 1
                    l += 1
                else:
                    #print ("other")
                    l += 1
                #
            #|endwhile
            
        #|endfor
        
        return n
    
    #
    
    # currently the favored tool
    def separate_sheep_from_goats(self, vkeys):
        debug = False # True # 
        
        if debug:
            print ("Enter separate_sheep_from_goats")
            for v in vkeys:
                print ("(%3d, %3d)" % (v[0], v[1]))
            #
            
        #
        
        kpntr = deepcopy(vkeys)
        n = len(kpntr)
        if debug: print ("n = ", n)
        self.dPKlinkages = {}
        for k in range(len(kpntr)):
            ik = kpntr[k][0]; jk = kpntr[k][1]
            if (ik, jk) in self.dPKlinkages:
                continue
            #
            
            
            l = k + 1
            flag_cont = True
            while l < len(kpntr) and flag_cont:
                il = kpntr[l][0]; jl = kpntr[l][1]
                if debug:
                    print ("ijk(%2d,%2d) v ijl(%2d,%2d)" % (ik, jk, il, jl))
                #
                
                if   il > jk:
                    if debug:
                        print ("finished")
                    #
                    
                    flag_cont = False
                elif (il, jl) in self.dPKlinkages or (ik, jk) in self.dPKlinkages:
                    # skip it
                    if debug:
                        print ("skip ijl(%2d,%2d)" % (il, jl))
                    #
                    
                    l += 1
                    #sys.exit(0)
                    
                elif ik < il and jl < jk:
                    
                    if debug:
                        print ("inside")
                    #
                    
                    l += 1
                elif ik < il and il < jk and jk < jl:
                    if debug:
                        print ("crossover")
                    #
                    
                    pkdmn = self.untangle_crossover(ik, jk, k, il, jl, l, vkeys)
                    
                    ilr = pkdmn.linkage[0][0]
                    jlr = pkdmn.linkage[0][1]
                    
                    if debug: print ("pkdmn: ", pkdmn)
                    
                    self.dPKlinkages.update( { (ilr, jlr) : pkdmn } )
                    
                    if debug: print (self.dPKlinkages)
                    #sys.exit(0)
                    flag_stop_at = False # True # 
                    if flag_stop_at:
                        #i_t =  7; j_t = 45
                        #i_t =  5; j_t = 45
                        i_t = 15; j_t = 65
                        #i_t =  5; j_t = 25
                        #i_t =  5; j_t = 85
                        #i_t = 22; j_t = 48
                        if i_t == ilr and j_t == jlr:
                            print ("stop when ijl(%d,%d)" % (il, jl))
                            sys.exit(0)
                        #
                        
                    #
                    
                    
                    l += 1
                else:
                    if debug:
                        print ("other")
                    #
                    
                    l += 1
                #
                
            #|endwhile
            
        #|endfor
        
        if debug:
            print ("exit separate_sheep_from_goats")
            for v in list(self.dPKlinkages.keys()):
                print ("(%3d, %3d)" % (v[0], v[1]), self.dPKlinkages[v])
            #
            
        #
        
        return self.dPKlinkages
    
    #
    
    
    # currently the favored tool
    def findMaxStemTails(self, vkeys):
        debug = False # True # 
        
        if debug:
            print ("Enter findMaxStemTails")
            print ("dXPKlinfo", self.dXPKlinfo)
            print ("dCPKlinfo", self.dCPKlinfo)
            #print ("stop at 0 in findMaxStemTails"); sys.exit(0)
        #
        
        
        for k in range(0, len(vkeys)):
            vk = vkeys[k]
            ik = vk[0]; jk = vk[1]
            
            self.prncpl_wts[vk].npairs = self.prncpl_wts[vk].stem
            self.prncpl_wts[vk].rpntr += [vk]
            self.maximize_stemtails(ik, jk, k, vkeys)
            self.prncpl_wts[vk].nsst \
                = self.subtract_crossovers(self.prncpl_wts[vk].fpntr)
            #self.find_rpntr(k, ik, jk, vkeys) 
            if debug: print (self.prncpl_wts[vk])
        #
        
        self.separate_sheep_from_goats(vkeys)
        
        for vv in list(self.dPKlinkages.keys()):
            pk = deepcopy(self.dPKlinkages[vv])
            if self.dPKlinkages[vv].is_extended:
                self.dXPKlinfo.update({vv : pk })
            else:
                self.dCPKlinfo.update({vv : pk })
            #
        #
        
        if debug:
            print ("### finished part 1 of findMaxStemTails")
            print ("dCPKlinfo: ", self.dCPKlinfo)
            print ("dXPKlinfo: ", self.dXPKlinfo)
            #print ("stop at part 1 findMaxStemTails"); sys.exit(0)
        #
        
        
        bss, bpk = self.maximize_sstails(vkeys)
        
        self.sskeylist = sortPairListWRT_n(bss, 0) # sort wrt index 0
        self.pkkeylist = sortPairListWRT_n(bpk, 0) # sort wrt index 0
        
        
        if debug:
            print ("### finished part 2 of findMaxStemTails")
            print (bss, bpk)
            print ("dCPKlinfo: ", self.dCPKlinfo)
            print ("dXPKlinfo: ", self.dXPKlinfo)
            #print ("stop at part 2 findMaxStemTails"); sys.exit(0)
        #
        
        # now we build the pair list shifting the contests of
        # pkkeylist to PKlist and retaining the remaining material in
        # BPlist.
        
        self.revise_pointers()
        
        
        if debug:
        
            print ("### finished part 3 of findMaxStemTails")
            print (self.display_prncpl_wts("ss tails:", self.sskeylist))
            print (self.display_prncpl_wts("pk tails:", self.pkkeylist))
            
            
            print ("stop at the end of findMaxStemTails"); sys.exit(0)
            
        #
        #print ("stop at the end of findMaxStemTails"); sys.exit(0)
        
                
    #
    
    
    # support for findMaxStemTails
    def maximize_stemtails(self, ib, jb, k, vkeys):
        # vkeys: OrderedDict set of keys
        
        debug = False # True # 
        """@
        
        Within the boundaries of ib,jb, obtains references to the
        structures (stem loops) that are contained inside. At this
        stage, it does not discriminate if it is a PK or not, only
        whether the stemsd exist inside the boundaries.
        
        This data is filled into prncpl_wts in the handles npairs
        and fpntr
        """
        
        if debug:
            print ("maximize_stemtails(ijb(%2d,%2d), k(%2d))" % (ib, jb, k))
            
            vk = vkeys[k]; ik = vk[0]; jk = vk[1]
            print ("reference k = %2d -> ijb(%2d,%2d)" % (k, ik, jk))
            """
            print ("processing vkeys: ")
            for vv in vkeys:
                print ("(%3d, %3d)" % (vv[0], vv[1]))
            #
            """
            
        #
        
        
        fbranching = []
        
        l = k + 1; n = len(vkeys)
        flag_cont = True
        
        while l < n and flag_cont:
            vl = vkeys[l]; il  = vl[0]; jl  = vl[1]
            if debug:
                print ("l = %2d; -> ib(%d) < il(%d) < jl(%d) < jb(%d)" \
                       % (l, ib, il, jl, jb))
            #
            
            if ib < il and jl < jb:
                fbranching += [(il, jl)]
                
                if debug: print (l, fbranching)
                
                self.prncpl_wts[vkeys[k]].npairs += self.prncpl_wts[vkeys[l]].stem
                l += 1
                
            elif jb <= il:
                flag_cont = False
            else:
                if debug: print ("skip ijl(%2d,%2d)" % (il, jl))
                l += 1
            #
            
        #
        
        
        
        
        if len(fbranching) > 0:
            self.prncpl_wts[vkeys[k]].fpntr += fbranching
                
        #
        
        if debug:
            print ("exiting maximize_stemtails: ")
            print ("prncpl_wts:")
            self.prncpl_wts[vkeys[k]].show_Weights(True)
            
            #if k == 3: sys.exit(0)
            
            #print ("stop at end of maximize_stemtails"); sys.exit(0)
            
        #
        
        return l, fbranching
    #
    
    
    
    # support for findMaxStemTails
    def scanFor_extendedPK(self,
                           ik, jk, k,
                           il, jl, l,
                           vkeys, debug_sfe = False):
        debug = debug_sfe
        #debug = True # False #
        
        if debug:
            print ("Enter scanFor_extendedPK(ijk(%2d,%2d), ijl(%2d,%2d))" \
                   % (ik, jk, il, jl))
        #
        
        is_extended = False
        m = l + 1; flag_cont = True
        n = len(vkeys)
        while m < n and flag_cont:
            vm = vkeys[m]; im = im = vm[0]; jm = vm[1]
            
            
            if jl < im:
                
                """@
                
                # (im, jm) OUT OF RANGE of (il, jl)
                
                        il            jl
                        |              |
                 ((((...[[....))))....]]...(((....)))
                 |               |         |        |
                 ik              jk        im       jm
                
                (im,jm) is out of range of the PK.  """
                
                flag_cont = False
                
            elif ik < im and im < jk and jk < jm:
                
                """@
                
                # im leg of (im, jm) is INSIDE (ik, jk) 
                
                        il              jl
                        |                |
                 ((((...[[..<<..))))....]]..>>
                 |          |      |         |
                 ik         |     jk         |
                            im               jm
                
                (ik, jk) and (im, jm) also form a PK, so it (im,jm)
                overlaps with (ik,jk).  Hence, whatever it is, it is
                not an extended pseudoknot. """
                
                if debug:
                    print ("skipping: ik(%d) < im(%d) < jk(%d) < jm(%d)" % (ik, im, jk, jm))
                #
                
                m += 1
                
            elif il < im and im < jl and jl < jm:
                
                """@
                
                # potentially "just right" (but need to check things too)
                
                        il                 jl
                        |                   |
                 ((((...[[....))))....(((..]]..)))
                 |               |    |          |
                 ik              jk   im        jm
                
                The extended PK conditions are met.  """
                
                nk = self.prncpl_wts[(ik, jk)].nsst
                nl = self.prncpl_wts[(il, jl)].nsst
                nm = self.prncpl_wts[(im, jm)].nsst
                
                if not nl > nk:
                    # if nk has lots of contacts, it is probably the
                    # best choice for a doman.
                    
                    #            root1     root2      linkage
                    pk = PKdmn([(ik, jk), (im, jm)], [(il, jl)])
                    
                    print ("k-nsst: ", self.prncpl_wts[(ik, jk)].nsst)
                    print ("l-nsst: ", self.prncpl_wts[(il, jl)].nsst)
                    print ("m-nsst: ", self.prncpl_wts[(im, jm)].nsst)
                    # The handle "rpntr" needs the information about
                    # ss stem tails, so at this point, it is not
                    # sufficiently defined yet
                    
                    self.dXPKlinfo.update({ (il, jl) : pk })
                    
                    
                    if debug:
                        print ("found extended PK:")
                        print ("ik(%3d) < il(%3d) < jk(%3d) < im(%3d) < jl(%3d) < jm(%3d)" \
                               % (ik, il, jk, im, jl, jm))
                        print ("kkkkkkk   -------   kkkkkkk   mmmmmmm   -------   mmmmmmm")
                        #print ("exit in scanFor_extendedPK"); sys.exit(0)
                    #
                    
                    
                    is_extended = True
                    flag_cont = False
                    # we're only looking at this particular linkage
                    # (il,jl) not all linkages associated with
                    # (ik,jk).
                    
                else:
                    tk, tm = self.is_tangled_linkage(ik, jk, il, jl, im, jm)
                    
                    """@
                    
                    should be balanced
                    
                    [
                    """
                    
                    if debug:
                        print ("ijk: ", self.prncpl_wts[(ik, jk)].fpntr)
                        print ("ijl: ", self.prncpl_wts[(il, jl)].fpntr)
                        print ("ijm: ", self.prncpl_wts[(im, jm)].fpntr)
                        print ("nl(%d) > nk(%d) or nl(%d) > nm(%d)" % (nl, nk, nl, nm))
                        sys.exit(0)
                    #
                    
                    m += 1
                #
                
                """@
                
                Both root1 and root2 contain pointers to themselves
                presently. However, it is not necessarily that the
                root stem is located at the base of the domain, it can
                be inside. We must have the full hierarchy available
                at least, in case the PK connects to some stem loop
                inside of a larger domain. In assembling the motif,
                this information needs to be present.
                
                pk.linkage = 
                this can only be the linkage, but it is posible that
                there is more than one linkage stem that is
                connected together.
                
                """
                
            else:
                m += 1
            #
        #
        
        if debug:
            print ("exiting scanFor_extendedPK()")
        #
        
        return is_extended
    #
    
    def is_tangled_linkage(self, ik, jk, il, jl, im, jm):
        lpntr = self.prncpl_wts[(il,jl)].fpntr
        tangled_with_k = False
        tangled_with_m = False
        for lp_k in lpntr:
            ilp = lp_k[0]; jlp = lp_k[1]
            
            if ilp < jk and jk < jlp:
                print ("ilp(%d) < jk(%d) < jlp(%d)" % (ilp, jk, jlp))
                tangled_with_k = True
            #
            
            if ilp < im and im < jlp:
                print ("ilp(%d) < im(%d) < jlp(%d)" % (ilp, im, jlp))
                tangled_with_m = True
            #
                
            
        #
        
        print ("tangled_with_k(%s), tangled_with_m(%s)" % (tangled_with_k, tangled_with_m))
        sys.exit(0)
        return tangled_with_k, tangled_with_m
    #
    
    
    
    def untangle_crossover(self, ikr, jkr, kr, ilr, jlr, lr, vkeys):
        """@
        
        This routine tries to unscramble whether the PK linkage is a
        core PK or an extended PK. It turns out that this was a rather
        difficult problem to solve. I wouldn't have thought so, all
        you have to do is look at the pairs and you can decide, but it
        seems that it is not so simple. This routine is meant to
        generate the initial stem tail arrangement and whether those
        stem tails reflect a core pseudoknow or an extended
        pseudoknot.
        
        """
        
        debug = False # True # 
        if debug:
            print ("untangle_crossover(ijkr(%2d,%2d), ijlr(%2d,%2d))" \
                   % (ikr, jkr, ilr, jlr))
        #
        
        
        is_extended, id_mn, jd_mn = self.check_is_extendedPK(ikr, jkr, kr,
                                                             ilr, jlr, lr,
                                                             vkeys, debug)
        
        if debug:
            print ("is_extended(%s) ijd_mn(%2d, %2d)" % (is_extended, id_mn, jd_mn))
        #
        
        
        l = kr + 1
        d_mx = jd_mn - id_mn + 1
        id_mx = id_mn; jd_mx = jd_mn
        v = {}
        for l in range(l, len(vkeys)):
            idp = vkeys[l][0]; jdp = vkeys[l][1]
            #print ("ijdp(%d, %d)" % (idp, jdp))
            if idp <= id_mn and jd_mn <= jdp:
                if debug:
                    print ("idp(%d) <= id_mn(%d), jd_mn(%d) <= jdp(%d)" \
                           % (idp, id_mn, jd_mn, jdp))
                #
                
                if is_extended:
                    if idp > jkr:
                        v.update({(idp, jdp) : True })
                    #
                    
                else:
                    v.update({(idp, jdp) : True })
                #
                
            #
            
        #
        
        del v[(id_mn, jd_mn)]
        
        vv = list(v.keys())
        
        if debug:
            print (v)
            for vvk in vv:
                print (vvk)
            #|endfor
            
        #
        
        id_mx, jd_mx = self.pruneStemTailList(vv, id_mn, jd_mn, debug)
        d_mx = jd_mx - id_mx + 1
        im = id_mx; jm = jd_mx
        
        nk = self.prncpl_wts[(ikr,jkr)].nsst
        nl = self.prncpl_wts[(ilr,jlr)].nsst
        nm = self.prncpl_wts[(id_mx,jd_mx)].nsst
        
        flag_weight_rhs = False
        flag_weight_lhs = False
        flag_weight_mid = False
        
        stroots = []
        fpntr = []
        
        ir1 = ikr;   jr1 = jkr
        ir2 = -1;    jr2 = -1
        il  = ilr;   jl  = jlr
        debugC = debug # True # False # 
        debugX = debug # True # False # 
        
        
        if jkr < id_mx:
            # effectively some sort of extended result popped out
            
            ir1 = ikr;    jr1 = jkr
            ir2 = im;     jr2 = jm
            il  = ilr;    jl  = jlr
            stroots = rlist2vsPair(vkeys, 'a', "bp")
            fpntr = get_XPK_jbranches(ir1, jr1, ir2, jr2, stroots) 
            #print ("fpntr: ", fpntr)
            
            kmscore, lmscore = self.score_Xedit(ir1, jr2, kr,
                                                ilr, jlr,
                                                im, jm, vkeys, debugX)
            
            #print ("nk/l/m: ", nk, nl, nm)
            #print ("km/lmscore: %s %s" % (kmscore, lmscore))
            
            if kmscore > lmscore:
                
                if debug:
                    print ("X: kmscore(%d) > lmscore(%d), nk(%d)/nl(%d)/nm(%d)" \
                           % (kmscore, lmscore, nk, nl, nm))
                #
                
                flag_weight_lhs = True
                flag_weight_rhs = True
                is_extended     = True
                self.dQroots.update({ (ir1, jr1) : True })
                self.dQroots.update({ (ir2, jr2) : True })
                
                
            elif kmscore <= lmscore:
                kscore, lscore = self.score_Cedit(ikr, jkr, kr,
                                                  ilr, jlr, vkeys, debugX)
                
                if debug:
                    print ("C: kscore(%d) <= lscore(%d), nk(%d)/nl(%d)" \
                           % (kscore, lscore, nk, nl))
                #
                
                if nk == 0 and nl > 0:
                    ir1 = ilr;    jr1 = jlr
                    il  = ikr;    jl  = jkr
                    flag_weight_lhs = False
                    flag_weight_mid = True
                    is_extended     = False
                elif kscore < lscore:
                    ir1 = ilr;    jr1 = jlr
                    il  = ikr;    jl  = jkr
                    flag_weight_lhs = False
                    flag_weight_mid = True
                    is_extended     = False
                elif kscore >= lscore:
                    ir1 = ikr;    jr1 = jkr
                    il  = ilr;    jl  = jlr
                    flag_weight_lhs = True
                    flag_weight_mid = False
                    is_extended     = False
                else:
                    
                    print ("nk/l/m: ", nk, nl, nm)
                    print ("km/lmscore: %s %s" % (kmscore, lmscore))
                    print ("klscore: %s %s" % (kscore, lscore))
                    print ("??? don't know why we're here 1")
                    sys.exit(0)
                    
                #
                
                # final correction, in case one of these roots was
                # already spoken for. ... ... ...
                if (il, jl) in self.dQroots:
                    ir1 = ikr;    jr1 = jkr
                    il  = ilr;    jl  = jlr
                    flag_weight_lhs = True
                    flag_weight_mid = False
                    is_extended     = False
                #
                
            else:
                
                print ("nk/l/m: ", nk, nl, nm)
                print ("km/lmscore: %s %s" % (kmscore, lmscore))
                print ("klscore: %s %s" % (kscore, lscore))
                print ("??? don't know why we're here 2")
                sys.exit(0)
            #
            
            
        else:
            
            kscore, lscore = self.score_Cedit(ikr, jkr, kr, ilr, jlr, vkeys, debugC)
            stroots = rlist2vsPair(vkeys, 'a', "bp")
            
            if kscore > lscore:
                if debug:
                    print ("C: kscore(%d) > lscore(%d), nk(%d)/nl(%d)" \
                           % (kscore, lscore, nk, nl))
                #
                
                ir1 = ikr;    jr1 = jkr
                il  = ilr;    jl  = jlr
                flag_weight_lhs = True
                flag_weight_mid = False
                is_extended     = False
                self.dQroots.update({ (ir1, jr1) : True })
                self.dQroots.update({ (ir2, jr2) : True })
            else:
                if debug:
                    print ("C: kscore(%d) <= lscore(%d), nk(%d)/nl(%d)" \
                           % (kscore, lscore, nk, nl))
                #
                
                ir1 = ilr;    jr1 = jlr
                il  = ikr;    jl  = jkr
                flag_weight_lhs = False
                flag_weight_mid = True
                is_extended     = False
                self.dQroots.update({ (ir1, jr1) : True })
                self.dQroots.update({ (ir2, jr2) : True })
            #
            
            # final correction, in case one of these roots was
            # already spoken for. ... ... ...
            if (il, jl) in self.dQroots:
                tmp = (ir1, jr1)
                ir1 = il;      jr1 = jl
                il  = tmp[0];  jl  = tmp[1]
                if il < ir1:
                    flag_weight_lhs = False
                    flag_weight_mid = True
                else:
                    flag_weight_lhs = True
                    flag_weight_mid = False
                #
                
                is_extended     = False
            #
                
        #
        
        
        
        if debug:
            
            if is_extended:
                print ("ijr1(%2d,%2d), ijr2(%2d,%2d), ijl(%2d, %2d)" \
                       % (ir1, jr1, ir2, jr2, il, jl))
            else:
                print ("ijr1(%2d,%2d), ijl(%2d, %2d)" \
                       % (ir1, jr1, il, jl))
            #
            
            print ("ijd_mx(%2d,%2d), d_mx(%2d), is_extended(%s)" \
                   % (id_mx, jd_mx, d_mx, is_extended))
            print ("nk(%d), nl(%d), nm(%d)" % (nk, nl, nm))
            
            print ("jkr(%d) < id_mx(%d)" % (jkr, id_mx))
            
            print ("flag_weight_lhs: ", flag_weight_lhs)
            print ("flag_weight_mid: ", flag_weight_mid)
            print ("flag_weight_rhs: ", flag_weight_rhs)
            print ("is_extended:     ", is_extended)
            #sys.exit(0)
            
        #
        
        
        roottails = [(ir1, jr1)]
        if is_extended:
            roottails += [(ir2, jr2)]
        #
        
        linkage = [(il, jl)]
        
        pk = PKdmn(roottails, linkage)
        pk.jbranches = fpntr
        
        
        if debug:
            
            print (pk)
            
            if is_extended:
                print ("extended PK")
                print ("ijr1(%3d,%3d), ijl(%3d,%3d), ijr2(%3d,%3d)" \
                       % (ir1, jr1, il, jl, ir2, jr2))
            else:
                print ("core PK")
                print ("ijr1(%3d,%3d), ijl(%3d,%3d)" \
                       % (ir1, jr1, il, jl))
            #
            
        #
        
        #sys.exit(0)
        flag_stop_at = False # True # 
        if flag_stop_at:
            #i_t = 7; j_t = 45
            i_t = 5; j_t = 45
            i_t =15; j_t = 65
            #i_t = 5; j_t = 25
            #i_t = 5; j_t = 85
            if ilr == i_t and jlr == j_t:
                print ("planned stop at ijlr(%d,%d) in untangle_crossover" \
                       % (i_t, j_t))
                sys.exit(0)
            #
            
            
        #
        
        return pk
    #
    
    
    
    def score_Cedit(self, ikr, jkr, kr, ilr, jlr, vkeys, debug = False):
        
        if debug:
            print ("Enter score_Cedit(ijkr(%2d, %2d), ijlr(%2d, %2d))" \
                   % (ikr, ikr, ilr, jlr))
            print ("vkeys: ", vkeys)
            
        #
        
        kshow = {(ikr, jkr) : True}
        lshow = {(ilr, jlr) : True}
        
        
            
        if debug:
            print ("kshow ijkr(%2d, %2d)" % (ikr, jkr))
            print ("kshow ijlr(%2d, %2d)" % (ilr, jlr))
        #
        
        k = kr
        
        kscore = 0 # len(vkeys)
        lscore = 0 # len(vkeys)
        
        ik = ikr; jk = jkr
        
        if debug:         
            print ("build k-score")
        #
        
        l = kr; flag_cont = True
        while l < len(vkeys) and flag_cont:
            iv = vkeys[l][0]; jv = vkeys[l][1]
            kscore += 1
            
            if debug: print ("ijv(%2d,%2d)" % (iv, jv))
            
            
            if   jlr < iv: # scan all the way from ik to jm
                if debug: print ("l finished")
                flag_cont = False
                
            elif jv > jlr:
                if debug: print ("skip iv jlr(%d) < jv(%d)" % (jlr, jv))
                
                
            elif ik <= iv and jv <= jk:
                
                if debug:
                    print ("pass ik(%2d) <= iv(%2d) < jv(%2d) <= jk(%2d)" \
                           % (ik, iv, jv, jk))
                #
                
            elif ik < iv and iv < jk and jk < jv:
                
                if debug:
                    print ("ik(%2d) < iv(%2d) < jk(%2d) < jv(%2d)" \
                           % (ik, iv, jk, jv))
                #
                
                if (ik, jk) in kshow:
                    # when k are removed, this is a penalty for l
                    
                    if debug:
                        print ("k(%2d,%2d) entangled with l(%2d,%2d), delete from kscore" \
                               % (ik, jk, iv, jv))
                    #
                    
                    kscore -= 1
                        
                #
                
                
            elif iv < ik and ik < jv and jv < jk:
                    
                if debug:
                    print ("iv(%2d) < ik(%2d) < jv(%2d) < jk(%2d)" \
                           % (iv, ik, jv, jk))
                #
                
                if (iv, jv) in kshow:
                    # when k are removed, this is a penalty for l
                    
                    if debug:
                        print ("v(%2d,%2d) entangled with k(%2d, %2d), delete from kscore" \
                               % (iv, jv, ik, jk))
                    #
                    
                    kscore -= 1
                    
                #
                
            else:
                
                if debug: print ("ijk(%2d,%2d), ijv(%2d,%2d)" % (ik, jk, iv, jv))
                #sys.exit(0)
                
            #
                
            l += 1
            
        #|endwhile
        
        if debug: 
            print ("klscore: ", kscore, lscore)
            #sys.exit(0)
            print ("build l-score")
        #
        
        lscore = 0 # len(vkeys)
        
        il = ilr; jl = jlr
        
        l = kr; flag_cont = True
        while l < len(vkeys) and flag_cont:
            iv = vkeys[l][0]; jv = vkeys[l][1]
            lscore += 1
            
            if debug: print ("ijv(%2d,%2d)" % (iv, jv))
            
            if   jlr < iv: # scan all the way from il to jm
                if debug: print ("l finished")
                flag_cont = False
            elif jv > jlr:
                if debug: print ("skip iv jlr(%d) < jv(%d)" % (jlr, jv))
                
                
            elif il <= iv and jv <= jl:
                if debug: 
                    print ("pass il(%2d) <= iv(%2d) < jv(%2d) <= jl(%2d)" \
                           % (il, iv, jv, jl))
                #
                
            elif il < iv and iv < jl and jl < jv:
                
                if debug: 
                    print ("il(%2d) < iv(%2d) < jl(%2d) < jv(%2d)" \
                           % (il, iv, jl, jv))
                #
                
                if (il, jl) in lshow:
                    # when l are removed, this is a penalty for k/m
                    
                    if debug:
                        print ("1 l(%2d,%2d) entangled with l(%2d,%2d), delete from lscore" \
                               % (il, jl, iv, jv))
                    #
                    
                    lscore -= 1
                #
                
                
            elif iv < il and il < jv and jv < jl:
                
                if debug: 
                    print ("iv(%2d) < il(%2d) < jv(%2d) < jl(%2d)" \
                           % (iv, il, jv, jl))
                #
                
                
                if (il, jl) in lshow:
                    # when l are removed, this is a penalty for k/m
                    
                    if debug:
                        print ("2 l(%2d,%2d) entangled with l(%2d,%2d), delete from lscore" \
                               % (il, jl, iv, jv))
                    #
                    
                    lscore -= 1
                #
                
            else: 
                if debug:
                    print ("ijl(%2d,%2d), ijv(%2d,%2d)" \
                           % (il, jl, iv, jv))
                #
                
            #
                
            l += 1
            
        #|endwhile
            
        
        
        if debug:
            
            print ("Exit score_Cedit()")
            print ("klscore: ", kscore, lscore)
            #sys.exit(0)
            
        #
        
        return kscore, lscore
    #
    
    
    def score_Xedit(self, ikr, jkr, kr, ilr, jlr, imr, jmr, vkeys, debug = False):
        
        if debug:
            print ("Enter score_Xedit(ijkr(%2d, %2d), ijlr(%2d, %2d), ijmr(%2d,%2d))" \
                   % (ikr, ikr, ilr, jlr, imr, jmr))
            print ("vkeys: ", vkeys)
            
        #
        
        kshow = {(ikr, jkr) : True}
        lshow = {(ilr, jlr) : True}
        mshow = {(imr, jmr) : True}
        
        if debug:
            print ("kshow ijkr(%2d, %2d)" % (ikr, jkr))
            print ("lshow ijlr(%2d, %2d)" % (ilr, jlr))
            print ("mshow ijmr(%2d, %2d)" % (imr, jmr))
            # kshow scean from (ikr to jkr)
            
        #
        
        kmscore = 0
        lmscore  = 0 
        
        k = kr
        ik = ikr; jk = jkr
        
        if debug: print ("build k part of score")
        l = kr; flag_cont = True
        while l < len(vkeys) and flag_cont:
            iv = vkeys[l][0]; jv = vkeys[l][1]
            kmscore += 1
            if debug: print ("ijv(%2d,%2d)" % (iv, jv))
            
            
            if   jmr < iv: # scan all the way from ik to jm
                if debug: print ("l finished")
                flag_cont = False
                
            elif jv > jmr:
                if debug: print ("skip iv jmr(%d) < jv(%d)" % (jmr, jv))
                
                
            elif ik <= iv and jv <= jk:
                
                if debug:
                    print ("pass ik(%2d) <= iv(%2d) < jv(%2d) <= jk(%2d)" \
                           % (ik, iv, jv, jk))
                #
                
            elif ik < iv and iv < jk and jk < jv:
                
                if debug:
                    print ("ik(%2d) < iv(%2d) < jk(%2d) < jv(%2d)" \
                           % (ik, iv, jk, jv))
                #
                
                include_linkage = False
                if iv < jmr and jmr < jv:
                    include_linkage = True
                #
                
                if (ik, jk) in kshow and include_linkage:
                    # when k are removed, this is a penalty for l
                    
                    if debug:
                        print ("k(%2d,%2d) entangled with l(%2d,%2d), delete from kmscore" \
                               % (ik, jk, iv, jv))
                    #
                    
                    kmscore -= 1
                        
                elif not include_linkage:
                    
                    if debug: 
                        print ("NOTE: linkage l(%2d,%2d) does not satisfy" \
                               % (iv, jv))
                        print ("        im(%d) < jl(%d) < jm(%d), skipped" \
                               % (imr, jv, jmr))
                    #
                    
                #
                
            elif iv < ik and ik < jv and jv < jk:
                # this is probably not necessary
                
                if debug:
                    print ("iv(%2d) < ik(%2d) < jv(%2d) < jk(%2d)" \
                           % (iv, ik, jv, jk))
                #
                
                include_linkage = False
                if iv < jmr and jmr < jv:
                    include_linkage = True
                #
                
                if (iv, jv) in kshow:
                    # when k are removed, this is a penalty for l
                    
                    if debug:
                        print ("v(%2d,%2d) entangled with k(%2d, %2d), delete from kmscore" \
                               % (iv, jv, ik, jk))
                    #
                    
                    kmscore -= 1
                    
                #
                
            else:
                
                if debug:
                    print ("ijk(%2d,%2d), ijv(%2d,%2d)" \
                           % (ik, jk, iv, jv))
                #
                
                #sys.exit(0)
                
            #
                
            l += 1
            
        #|endwhile
        
        
        if debug:
            print ("kmscore (after calculating k part): ", kmscore)
            #sys.exit(0)
            print ("build m part of score")
            
        #
        
        
        k = kr
        im = imr; jm = jmr
        
        l = kr; flag_cont = True
        while l < len(vkeys) and flag_cont:
            iv = vkeys[l][0]; jv = vkeys[l][1]
            if debug: print ("ijv(%2d,%2d)" % (iv, jv))
            
            
            if   jmr < iv: # scan all the way from im to jm
                if debug: print ("l finished")
                flag_cont = False
                
            elif jv > jmr:
                
                if debug:  print ("skip iv jmr(%d) < jv(%d)" % (jmr, jv))
                
                
            elif im <= iv and jv <= jm:
                
                if debug:
                    print ("pass im(%2d) <= iv(%2d) < jv(%2d) <= jm(%2d)" \
                           % (im, iv, jv, jm))
                #
                
            elif im < iv and iv < jm and jm < jv:
                
                if debug:
                    print ("im(%2d) < iv(%2d) < jm(%2d) < jv(%2d)" \
                           % (im, iv, jm, jv))
                #
                
                if (im, jm) in mshow:
                    # when k are removed, this is a penalty for l
                    
                    if debug:
                        print ("k(%2d,%2d) entangled with l(%2d,%2d), delete from kmscore" \
                               % (im, jm, iv, jv))
                    #
                    
                    kmscore -= 1
                    
                #
                
                
            elif iv < im and im < jv and jv < jm:
                    
                if debug:
                    print ("iv(%2d) < im(%2d) < jv(%2d) < jm(%2d)" \
                           % (iv, im, jv, jm))
                #
                
                if (im, jm) in mshow:
                    # when k are removed, this is a penalty for l
                    
                    if debug:
                        print ("v(%2d,%2d) entangled with k(%2d, %2d), delete from kmscore" \
                               % (iv, jv, im, jm))
                    #
                    
                    kmscore -= 1
                    
                #
                
            else:
                
                if debug:
                    print ("ijm(%2d,%2d), ijv(%2d,%2d)" \
                           % (im, jm, iv, jv))
                    #sys.exit(0)
                #
                
            #
                
            l += 1
            
        #|endwhile
        
        if debug: 
            print ("kmscore (after calculating m part): ", kmscore)
            #sys.exit(0)
            print ("build lm part of score")
        #
        
        il = ilr; jl = jlr
        
        l = kr; flag_cont = True
        while l < len(vkeys) and flag_cont:
            iv = vkeys[l][0]; jv = vkeys[l][1]
            lmscore += 1
            
            if debug: print ("ijv(%2d,%2d)" % (iv, jv))
            
            if   jmr < iv: # scan all the way from il to jm
                if debug: print ("l finished")
                flag_cont = False
            elif jv > jmr:
                if debug: print ("skip iv jmr(%d) < jv(%d)" % (jmr, jv))
                
                
            elif il <= iv and jv <= jl:
                if debug: 
                    print ("pass il(%2d) <= iv(%2d) < jv(%2d) <= jl(%2d)" \
                           % (il, iv, jv, jl))
                #
                
            elif il < iv and iv < jl and jl < jv:
                if debug: 
                    print ("il(%2d) < iv(%2d) < jl(%2d) < jv(%2d)" \
                           % (il, iv, jl, jv))
                #
                
                if (il, jl) in lshow:
                    # when l are removed, this is a penalty for k/m
                    
                    if debug:
                        print ("1 l(%2d,%2d) entangled with l(%2d,%2d), del from lmscore" \
                               % (il, jl, iv, jv))
                    #
                    
                    lmscore -= 1
                #
                
                
            elif iv < il and il < jv and jv < jl:
                if debug: 
                    print ("iv(%2d) < il(%2d) < jv(%2d) < jl(%2d)" \
                           % (iv, il, jv, jl))
                #
                
                
                if (il, jl) in lshow:
                    # when l are removed, this is a penalty for k/m
                    
                    if debug:
                        print ("2 l(%2d,%2d) entangled with l(%2d,%2d), del from lmscore" \
                               % (il, jl, iv, jv))
                    #
                    
                    lmscore -= 1
                #
                
            else:
                if debug: 
                    print ("ijl(%2d,%2d), ijv(%2d,%2d)" % (il, jl, iv, jv))
                #
                
            #
                
            l += 1
            
        #|endwhile
            
        
        
        if debug:
            
            print ("Exit score_Xedit()")
            print ("km/lmscore: ", kmscore, lmscore)
            #sys.exit(0)
            
        #
        
        
        return kmscore, lmscore
    #
    
    
    
    def pruneStemTailList(self, stlist, im, jm, debug = False):
        if debug:
            print ("Enter pruneStemTailList(ijm(%d,%d), %s" % (im, jm, stlist))
        #
        
        k = 0
        id_mx = im; jd_mx = jm
        flag_cont_k = True
        while k < len(stlist) and flag_cont_k:
            vk = stlist[k]
            ik = vk[0]; jk = vk[1]
            l = k + 1; flag_cont_l = True
            while l < len(stlist):
                vl = stlist[l]
                il = vl[0]; jl = vl[1]
                if ik < il and il < jk and jk < jl:
                    if debug: print ("ik(%d) < il(%d) < jk(%d) < jl(%d)" % (ik, il, jk, jl))
                    ddk = jk - ik
                    ddl = jl - il
                    # keep the shortest one
                    if ddk < ddl:
                        del stlist[l]
                    else:
                        del stlist[k]
                    #
                    
                else:
                    l += 1
                #
                
            #|emndwhile
            
            k += 1
            
        #|endwhile

        
        for stk in stlist:
            ik = stk[0]; jk = stk[1]
            if ik < id_mx and jd_mx < jk:
                id_mx = ik; jd_mx = jk
            #
            
        #|endfor
        
        if debug:
            print (stlist)
            print ("exit pruneStemTailList(ijm(%d,%d) -> ijd_mx(%d,%d)" \
                   % (im, jm, id_mx, jd_mx))
            #print ("stop at exit"); sys.exit(0)
        #
        
        return id_mx, jd_mx
    #
            
    
    def findMinStemLoop(self, im, jm, jlr, lr, vkeys, debug = False):
        if debug:
            print ("findMinStemLoop(im(%d), jlr(%d), jm(%d)" % (im, jlr, jm))
        #
        
        if im < jlr and jlr < jm:
            if debug: print ("pass")
        else:
            
            print ("findMinStemLoop failed: not im(%d) < jlr(%d) < jm(%d)"
                   % (im, jlr, jm))
            sys.exit(1)
        #
        
        
        l = lr
        flag_cont = True
        d_mn = 1e9
        id_mn = im; jd_mn = jm; l_mn = lr
        flag_cont = True
        while l < len(vkeys) and flag_cont:
            idp = vkeys[l][0]; jdp = vkeys[l][1]
            dd = jdp - idp + 1
            if debug: print ("ijdp(%2d,%2d)" % (idp, jdp))
            
            if id_mn < idp and  idp < jlr and jlr < jdp and jdp < jd_mn and d_mn > dd:
                
                if debug:
                    print ("id_mn(%d) < idp(%d) < jlr(%d) < jdp(%d) < jd_mn(%d)" \
                           % (id_mn, idp, jlr, jdp, jd_mn))
                #
                
                id_mn = idp; jd_mn = jdp; d_mn = dd; l_mn = l
                l += 1
            
            else:
                l += 1
                
            #
            
        #
        
        if debug:
            print ("findMinStemLoop: ijd_mn(%d,%d)" % (id_mn, jd_mn))
        #
        
        return id_mn, jd_mn
    
    #
    
        
    
    # support for findMaxStemTails (adapted from the an earlier and
    # now obsolete version scanFor_extendedPK)
    def check_is_extendedPK(self,
                            ik, jk, k,
                            il, jl, l,
                            vkeys, debug_sfe = False):
        debug = debug_sfe
        #debug = True # False #
        
        if debug:
            print ("Enter check_is_extendedPK(ijk(%2d,%2d), ijl(%2d,%2d))" \
                   % (ik, jk, il, jl))
        #
        
        is_extended = False
        m = l + 1; flag_cont = True
        n = len(vkeys)
        imv = il; jmv = jl
        while m < n and flag_cont:
            vm = vkeys[m]; im = im = vm[0]; jm = vm[1]
            
            
            if jl < im:
                
                """@
                
                # (im, jm) OUT OF RANGE of (il, jl)
                
                        il            jl
                        |              |
                 ((((...[[....))))....]]...(((....)))
                 |               |         |        |
                 ik              jk        im       jm
                
                (im,jm) is out of range of the PK.  """
                
                flag_cont = False
                
            elif ik < im and im < jk and jk < jm:
                
                """@
                
                # im leg of (im, jm) is INSIDE (ik, jk) 
                
                        il              jl
                        |                |
                 ((((...[[..<<..))))....]]..>>
                 |          |      |         |
                 ik         |     jk         |
                            im               jm
                
                (ik, jk) and (im, jm) also form a PK, so it (im,jm)
                overlaps with (ik,jk).  Hence, whatever it is, it is
                not an extended pseudoknot. """
                
                if debug:
                    print ("skipping: ik(%d) < im(%d) < jk(%d) < jm(%d)" % (ik, im, jk, jm))
                #
                
                m += 1
                
            elif il < im and im < jl and jl < jm:
                
                """@
                
                # potentially "just right" (but need to check things too)
                
                        il                 jl
                        |                   |
                 ((((...[[....))))....(((..]]..)))
                 |               |    |          |
                 ik              jk   im        jm
                
                The extended PK conditions are met.  """
                
                if debug:
                    print (" ..(.......A...........).........(..........a.........)")
                    print ("ik(%3d) < il(%3d) < jk(%3d) < im(%3d) < jl(%3d) < jm(%3d)" \
                           % (ik, il, jk, im, jl, jm))
                #
                
                imv, jmv = self.findMinStemLoop(im, jm, jl, k, vkeys, debug)
                
                is_extended = True
                flag_cont = False
                
            else:
                m += 1
            #
        #
        
        if debug:
            print ("ijmv(%2d,%2d)" % (imv, jmv))
            print ("exiting check_is_extendedPK()")
            #print ("stop at end of check_is_extendedPK"); sys.exit(0)
        #
        
        return is_extended, imv, jmv 
    #
    
    
    
    
    # support for findMaxStemTails
    def maximize_sstails(self, fbranching):
        debug = False # True # 
        
        if debug:
            print ("maximize_sstails", fbranching)
        #
        
        n = len(fbranching)
        bss = deepcopy(fbranching)
        bpk = []
        dpk_account = {} 
        
        k = 0
        while k < len(bss)-1:
            ik = bss[k][0]; jk = bss[k][1]
            
            # print ("k: ", k, bss[k])
            
            l = k + 1
            flag_cont = True
            while l < len(bss):
                
                # print ("l: ", l, bss[l])
                il = bss[l][0]; jl = bss[l][1]
                
                if debug:
                    print ("k(%2d), l(%2d): ik(%2d) < il(%2d) < jk(%2d) < jl(%2d)?" \
                           % (k, l, ik, il, jk, jl))
                #
                
                if (il, jl) in self.dXPKlinfo:
                    # remove extended PKs
                    if not (il, jl) in dpk_account:
                        bpk += [(il, jl)]
                        del bss[l]
                        dpk_account.update({ (il, jl) : True } )
                    else:
                        l += 1
                    #
                    
                elif ik < il and il < jk and jk < jl:
                    
                    # crossover found
                    
                    if (il, jl) in self.dCPKlinfo:
                        if debug:
                            print ("maximize_sstails:")
                            print ("root    stem(k = %2d): (%3d,%3d)" % (k, ik, jk))
                            print ("linkage stem(l = %2d): (%3d,%3d)" % (l, il, jl))
                        #

                        if not (il, jl) in dpk_account:
                            bpk += [(il, jl)]
                            del bss[l]
                            dpk_account.update({ (il, jl) : True } ) 
                        else:
                            l += 1
                        #
                        
                    else:
                        if debug:
                            print ("maximize_sstails:")
                            print ("root    stem(l = %2d): (%3d,%3d)" % (l, il, jl))
                            print ("linkage stem(k = %2d): (%3d,%3d)" % (k, ik, jk))
                        #
                        
                        if not (ik, jk) in dpk_account:
                            bpk += [(ik, jk)]
                            del bss[k]
                            dpk_account.update({ (ik, jk) : True } )
                        else:
                            l += 1
                        #
                        
                    #
                    
                    
                    if debug:
                        print ("updated in maximize_sstails:")
                        print ("ss: ", bss)
                        print ("pk: ", bpk)
                        #sys.exit(0)
                    #
                    
                else:
                    l += 1
                #
                
            #|endwhile
            
            k += 1
            #print ("new k: ", k)
            
        #
        
        if debug:
            print (bss, bpk)
            #print ("stop at end of maximize_sstails:"); sys.exit(0)
        #
        
        return bss, bpk
    #
    
    
    def obsolete_optimizeStemTailCount(self, pxlist, debug):
        
        # ###################################################################
        # ####  Now we find the weight on each of these principal stems  ####
        # ###################################################################
        
        vkeys = list(self.prncpl_wts.keys())
        n = len(vkeys)
        
        for k in range(0, n):
            vk = vkeys[k]
            ik = vk[0]; jk = vk[1];
            self.prncpl_wts[vk].npairs = self.prncpl_wts[vk].stem
            if k == n - 1:
                # an isolated last one means it can only be the stem itself.
                break
            #
            
            l = k + 1
            flag_cont_l = True
            
            
            while l < n and flag_cont_l:
                vl = vkeys[l]
                il  = vl[0]; jl  = vl[1]
                
                if debug:
                    print ("enter: k(%2d) l(%2d), (lijk(%2d,%2d), ijl(%2d,%2d)" \
                           % (k, l, ik, jk, il, jl))
                #
                
                if ik < il  and jl < jk:
                    # ap stem
                    self.prncpl_wts[(ik, jk)].npairs += self.prncpl_wts[(il, jl)].stem
                    l += 1
                    
                    if debug:
                        print ("ap stem: ijl(%2d,%2d)" % (il, jl))
                    #
                    
                elif  il - ik < self.ppgap and jl - jk < self.ppgap and il < jk:
                    
                    # pp stem (have to be more restrictive
                    
                    # ...([<ABC..DEFG...)]>abc..defg...
                    #    ^    ^  ^      ^    ^  ^
                    #    ik  ilx il    jk   jlx jl
                    
                    self.prncpl_wts[(ik, jk)].npairs += self.prncpl_wts[(il, jl)].stem
                    l += 1
                    
                    if debug:
                        print ("pp stem: ijl(%2d,%2d)" % (il, jl))
                    #
                    
                else:
                    if debug:
                        print ("shift: ijk(%2d,%2d), ijl(%2d,%2d)" % (ik, jk, il, jl))
                    #
                    
                    flag_cont_l = False
                    
                #
                
            #|endwhile
            
            k += 1
            print ("end: k = ", k)
        #|endfor
        
        
        if debug:
            print ("finished part 1 obsolete_optimizeStemTailCount <- findPrincipalDMNs")
            for v in list(self.prncpl_wts.keys()):
                print ("(%2d, %2d): %s" % (v[0], v[1], self.prncpl_wts[v]))
            #|endfor
            
            print (self.prncpl_wts.keys())
            #print ("stop here at part 2"); sys.exit(0)
            
        #
        
        pkeys = list(self.prncpl_wts.keys())
        self.make_apLevelList(-1, self.N+1, pkeys, 0, 0, debug)
        
        if debug:
            print ("finished part 2 obsolete_optimizeStemTailCount <- findPrincipalDMNs")
            print ("prncpl_wts:   ")
            for v in pkeys:
                print ("(%2d, %2d): %s" % (v[0], v[1], self.prncpl_wts[v]))
            #|endfor
            
            print ("pkeylevels: ")
            for v in list(self.pkeylevels.keys()):
                if len(self.pkeylevels[v]) > 0:
                    print ("level %2d: %s" % (v, self.pkeylevels[v]))
                #
                
            #|endfor
            
            #print ("stop here at part 3"); sys.exit(0)
        #
        
        # ########################################################
        # ####  decompose pkeylevels into a 1D list of pairs  ####
        # ########################################################
        
        pBPkeys = []
        for k in list(self.pkeylevels.keys()):
            
            for vr in self.pkeylevels[k]:
                iv = vr[0]; jv = vr[1]
                pBPkeys += [(iv, jv)]
            #|endfor
            
        #|endfor
        
        # sort the list of pairs
        pBPkeys = sortPairListWRT_n(pBPkeys, 0) # sort wrt index 0
        if debug:
            print ("pBPkeys:")
            print (pBPkeys)
        #
        
        
        # ###############################################
        # ####  make newBProots and newPKroots list  ####
        # ###############################################
        
        nb = len(pxlist)
        flag_cont = True
        ku = 0
        kp = 0
        newBProots = []
        newPKroots = []
        
        iu = 0; ju = 0
        for v in pkeys:
            iv = v[0]; jv = v[1]
            
            #print ("v: ", v)
            
            is_pk = False
            if ku < len(pBPkeys): 
                iu = pBPkeys[ku][0]; ju = pBPkeys[ku][1]
                
                if debug:
                    print ("pBPkeys[%2d]: (%2d,%2d)" % (ku, iu, ju))
                #
                
                if iu == iv and ju == jv:
                    ku += 1
                else:
                    is_pk = True
                #
                
            else:
                is_pk = True
            #
                
            if debug:
                print ("pkeys        (%2d,%2d), is_pk = %s" % (iv, jv, is_pk))
            #
            
            flag_cont = True
            while kp < nb and flag_cont:
                ib = pxlist[kp].i; jb = pxlist[kp].j
                
                if debug:
                    print ("BPlist(%2d): (%2d, %2d)" % (kp, ib, jb))
                #
                
                if is_pk:
                    if debug:
                        print ("pk <- (%2d,%2d)" % (ib, jb))
                    #
                    
                    if ib == iv and jb == jv:
                        newPKroots += [pxlist[kp]]
                        flag_cont = False
                    #
                    
                else:
                    if debug:
                        print ("ss <- (%2d,%2d)" % (ib, jb))
                    #
                    
                    if ib == iu and jb == ju:
                        newBProots += [pxlist[kp]]
                        flag_cont = False
                    #
                    
                #
                
                kp += 1
                #print ("next: -> ", pxlist[kp])
                
            #|endwhile
            
        #|endfor
        
        
        if debug:
            print ("newBProots:")
            for v in newBProots:
                print (v)
            #|endfor
            
            print ("newPKroots:")
            for v in newPKroots:
                print (v)
            #|endfor
            
            print ("stop at 4x in obsolete_optimizeStemTailCount"); sys.exit(0)
        #
        
        self.BProots = deepcopy(newBProots)
        self.PKroots = deepcopy(newPKroots)
        
        # ###########################################################
        # #### Finally, we move the pseudoknot pairs into newPKlist
        # #### and remover them from pxlist
        # ###########################################################
        
        npkr = len(newPKroots)
        kb   = 0
        newPKlist = []
        for vpk in newPKroots:
            i_tpk = vpk.i; j_tpk = vpk.j
            if debug:
                print ("ij_tpk(%2d,%2d)" % (i_tpk, j_tpk))
            #
            
            flag_cont = True
            kb = 0
            while kb < len(pxlist) and flag_cont:
                
                iv = pxlist[kb].i; jv = pxlist[kb].j
                if debug:
                    print ("ijv(%2d, %2d)[%2d]" % (iv, jv, kb))
                #
                
                if i_tpk == iv and jv == j_tpk:
                    for ks in range(0, self.prncpl_wts[(iv,jv)].stem):
                        if debug:
                            print ("write to newPKlist ",
                                   pxlist[kb],
                                   self.prncpl_wts[(i_tpk, j_tpk)].stem)
                        #
                        
                        newPKlist += [deepcopy(pxlist[kb])]
                        del pxlist[kb]
                    #|endfor
                    
                elif iv > j_tpk:
                    # BPlist is out of range
                    flag_cont = False
                    
                else:
                    kb += 1
                    
                #
                
            #|endwhile
            
        #|endfor
        
        self.BPlist = deepcopy(pxlist)
        self.PKlist = deepcopy(newPKlist)
        
        if debug:
            print ("Exit obsolete_optimizeStemTailCount")
            #print ("stop here at part 5"); sys.exit(0)
        #
        
    #
    
    
    def new_optimizeStemTailCount(self, pxlist, debug = False):
        
        debug = False # True # 
        
        if debug:
            print ("entered new_optimizeStemTailCount:")
            print (list(self.prncpl_wts.keys()))
            
        #
        
        pkeys = list(self.prncpl_wts.keys())
        
        self.findMaxStemTails(pkeys)
        
        
        if debug:
            print ("finished part 1 new_optimizeStemTailCount <- findPrincipalDMNs")
            print (self.display_prncpl_wts("prncpl_wts:", pkeys))
            
            #print ("stop here at part 2"); sys.exit(0)
            
        #
        
        
        pBPkeys = self.sskeylist
        pPKkeys = self.pkkeylist
        self.make_apLevelList(-1, self.N+1, pBPkeys, 0, 0, debug)
        
        if debug:
            print ("finished part 2 new_optimizeStemTailCount <- findPrincipalDMNs")
            print (self.display_prncpl_wts("prncpl_wts (ss):   ", pBPkeys))
            print (self.display_prncpl_wts("prncpl_wts (pk):   ", pPKkeys))
            
            print ("pkeylevels: ")
            for v in list(self.pkeylevels.keys()):
                if len(self.pkeylevels[v]) > 0:
                    print ("level %2d: %s" % (v, self.pkeylevels[v]))
                #
                
            #|endfor
            
            #print ("stop here at part 3"); sys.exit(0)
        #
        
        
        
        # ###############################################
        # ####  make newBProots and newPKroots list  ####
        # ###############################################
        
        newBProots = self.make_PairRoots(pBPkeys, pxlist)
        newPKroots = self.make_PairRoots(pPKkeys, pxlist)
        
        
        if debug:
            print ("newBProots:")
            for v in newBProots:
                print (v.disp_Pair())
            #|endfor
            
            print ("newPKroots:")
            for v in newPKroots:
                print (v.disp_Pair())
            #|endfor
            
            #print ("stop here at part 4"); sys.exit(0)
            
        #
        
        
        self.BProots = deepcopy(newBProots)
        self.PKroots = deepcopy(newPKroots)
        
        # ###########################################################
        # #### Finally, we move the pseudoknot pairs into newPKlist
        # #### and remover them from pxlist
        # ###########################################################
        
        newBPlist, newPKlist = self.make_finalPairList(newPKroots, pxlist)
        
        self.BPlist = deepcopy(newBPlist)
        self.PKlist = deepcopy(newPKlist)
        
        if debug:
            
            print (self.disp_allLists("final tally from new_optimizeStemTailCount", False))
            print ("Exit new_optimizeStemTailCount")
            #print ("stop here at part 5"); sys.exit(0)
        #
        
    #
    
    
    def display_prncpl_wts(self, title, pkeys):
        s = "%s\n" % (title)
        n = len(pkeys)
        for k in range(0, n):
            vk = pkeys[k]
            if k == 0:
                s += (self.prncpl_wts[vk].show_header() + '\n')
            #
            
            if k < n -1:
                s +=  (self.prncpl_wts[vk].show_Weights() + '\n')
            else:
                s +=  (self.prncpl_wts[vk].show_Weights())
            #
            
        #|endfor
        return s
    #
    
    
    def findPrincipalDMNs(self, pxlist):
        
        # I think this version does what I really want it to do over
        # all in psyching out parallel stems as well as the
        # antiparallel stems.
        
        """@
        
        This is used to identify the domains and help establish how to
        assign linkages. For example, a large domain should be
        weighted more heavily than a could of pairs with a large
        span. This helps say which pairs have priority.
        
        sets the following internal variables
        
        self.BProots ==> self.BPlist
        self.PKroots ==> self.PKlist
        
        """
        
        debug = False # True # 
        
        if debug:
            print ("Enter findPrincipalDMNs()")
        #
        
        flag_show_pxlist = False # True # 
        if flag_show_pxlist:
            self.show_pxlist(pxlist, "findPrincipalDMNs (input pxlist):")
        #
        
        # ###########################################################
        # ####  first find principal stems tails to latch onto   ####
        # ###########################################################
        
        self.prncpl_wts = OrderedDict()
        n = len(pxlist)
        k = 0
        flag_cont_k = True
        while k < n and flag_cont_k:
            
            ik = pxlist[k].i; jk = pxlist[k].j
            if debug:
                print ("ijk (%2d, %2d) vvvvvvvv" % (ik, jk))
                print ("k while: k = %d, %s" % (k, self.prncpl_wts))
            #
            
            self.prncpl_wts.update({(ik,jk) : Weights(ik, jk) })
            flag_cont_l = True
            l = k+1
            
            while l < n and flag_cont_l:
                
                # gather together contiguous and semi-contiguous parts
                # of the immediate stem. This is not following the
                # stricter rules of Kuhn length etc. but it is trying
                # to more or less get the most obvious of the types of
                # contiguous and effective stems.
                
                il  = pxlist[l].i;    jl  = pxlist[l].j
                ilx = pxlist[l-1].i;  jlx = pxlist[l-1].j
                
                if debug:
                    print ("ijk (%2d, %2d)" % (ik, jk))
                    print ("ijlx(%2d, %2d)" % (ilx, jlx))
                    print ("ijl (%2d, %2d)" % (il,  jl))
                    print ("l while: l = %d, %s" % (l, self.prncpl_wts))
                #
                
                if (ilx < il  and
                    jl < jlx and
                    il - ilx < self.apgap and
                    jlx - jl < self.apgap):
                    
                    # satsfies the gap requirements for antiparallel stems
                    
                    if debug: print ("x1 ap")
                    # ap stem
                    
                    # ...((((((..(((( . . . ))))..))))))....
                    #    ^    ^  ^             ^  ^    ^
                    #    ik  ilx il           jl jlx   jk
                    
                    self.prncpl_wts[(ik, jk)].stem += 1
                    l += 1
                    
                elif  (ik < il and il < jk and jk < jl):
                    
                    if debug:
                        sx  = "x2 pp: il(%d)-ilx(%d) = %d, " % (il, ilx, il-ilx)
                        sx += "jl(%d)-jlx(%d) = %d, " % (jl, jlx, jl - jlx)
                        sx += "il(%d) < jlx(%d)" % (il, jlx)
                        print (sx)
                    #
                    
                    if  il - ilx < self.ppgap and jl - jlx < self.ppgap and il < jlx:
                        
                        # satisfies the gap requrements for parallel stems
                        
                        if debug:
                            print ("contiguous: ")
                        #
                        
                        # pp stem
                        
                        # ...([<ABC..DEFG...)]>abc..defg...
                        #    ^    ^  ^      ^    ^  ^
                        #    ik  ilx il    jk   jlx jl
                        
                        if ik <= ilx and jk <= jlx:
                            self.prncpl_wts[(ik, jk)].stem += 1
                            l += 1
                            
                            if debug:
                                print ("x2a:")
                                print (self.prncpl_wts)
                            #
                            
                        elif (il, jl) in self.prncpl_wts:
                            # I think this shouldn't happen anymore
                            self.prncpl_wts[(ik, jk)].stem += 1
                            l += 1
                            
                            if True: # debug:
                                print ("x2b:")
                                print (self.prncpl_wts)
                                print ("thought this would never happen!")
                                print ("stop at x2b"); sys.exit(0)
                            #
                            
                        else:
                            
                            self.prncpl_wts.update({(il, jl) : Weights(il, jl) })
                            flag_cont_l = False
                            
                            if debug:
                                print (self.prncpl_wts)
                                print ("x2c: make ijl(%d,%d)" % (il, jl))
                                #print ("stop at x2c"); sys.exit(0)
                            #
                            
                        #
                        
                        
                    else:
                        
                        # if the gap is too big, the condition is not
                        # satisfied and we stop here. This only looks
                        # for the primative types of stems. Additional
                        # tools in Vienna2TreeNode are required to
                        # include Kuhn length rules. Here, we just use
                        # the very basic primative type of stem rules
                        # only.
                        
                        
                        if debug:
                            print ("make ijl(%d,%d)" % (il, jl))
                        #
                        
                        self.prncpl_wts.update({(il, jl) : Weights(il, jl) })
                        flag_cont_l = False
                        
                    #
                    
                    if debug:
                        print ("l = %d" % l)
                    #
                    
                else:
                    
                    # nothing is satisfied, so we must quit this cycle
                    # and look for another.
                    
                    if debug:
                        print ("x3 none")
                    #
                    
                    flag_cont_l = False
                #
                
            #|endwuile
            
            k = l
            
            if debug:
                print ("k = %d" % k)
            #
            
        #|endwhile
        
        # 
        
        if debug:
            print ("finished part 1 findPrincipalDMNs")
            for v in list(self.prncpl_wts.keys()):
                print ("(%2d, %2d): %s" % (v[0], v[1], self.prncpl_wts[v]))
            #
            
            #print ("stop here at part 1"); sys.exit(0)
        #
        
        #debug = True
        
        # ###################################################################
        # ####  Now we find the weight on each of these principal stems  ####
        # ###################################################################
        
        flag_use_old_method = False # True # 
        if flag_use_old_method:
            self.obsolete_optimizeStemTailCount(pxlist, debug)
        else:
            self.new_optimizeStemTailCount(pxlist, debug)
        #
        
        #debug = True
        if debug:
            print ("new BPlist:")
            for v in self.BPlist:
                print (v.disp_Pair())
            #|endfor
            
            print ("new BProots:")
            for v in self.BProots:
                print (v.disp_Pair())
            #|endfor
            
            print ("new PKlist:")
            for v in self.PKlist:
                print (v.disp_Pair())
            #|endfor
            
            print ("new PKroots:")
            for v in self.PKroots:
                print (v.disp_Pair())
            #|endfor
            
            print ("Exiting findPrincipalDMNs()")
            #print ("stop at exit of findPrincipalDMNs"); sys.exit(0)
            
        #
    
    #
    
    
    
    def make_finalPairList(self, newPKroots, pxlist, debug = False):
        
        # newPKroots: list of objects of class Pair
        # pxlist:     list of objects of class Pair
        
        # note the PKroots here. This operation transfers the data in
        # pxlist and drops into the list newPKlist. In the course of
        # this operation, all references to PKs are removed from
        # pxlist. Therefore, this should not be done stupidly or you
        # will end up with garbage!
        
        debug = False # True #
        
        npkr = len(newPKroots)
        kb   = 0
        newPKlist = []
        for vpk in newPKroots:
            i_tpk = vpk.i; j_tpk = vpk.j
            if debug:
                print ("ij_tpk(%2d,%2d)" % (i_tpk, j_tpk))
            #
            
            flag_cont = True
            kb = 0
            while kb < len(pxlist) and flag_cont:
                
                iv = pxlist[kb].i; jv = pxlist[kb].j
                if debug:
                    print ("ijv(%2d, %2d)[%2d]" % (iv, jv, kb))
                #
                
                if i_tpk == iv and jv == j_tpk:
                    for ks in range(0, self.prncpl_wts[(iv,jv)].stem):
                        if debug:
                            print ("write to newPKlist ",
                                   pxlist[kb],
                                   self.prncpl_wts[(i_tpk, j_tpk)].stem)
                        #
                        
                        newPKlist += [deepcopy(pxlist[kb])]
                        del pxlist[kb]
                    #|endfor
                    
                elif iv > j_tpk:
                    # BPlist is out of range
                    flag_cont = False
                    
                else:
                    kb += 1
                    
                #
                
            #|endwhile
            
        #|endfor
        
        return pxlist, newPKlist
    #
    
    def make_PairRoots(self, pkeys, pxlist, debug = False):
        #debug = False # True # 
        
        newroots = []
        
        nb = len(pxlist)
        kb = 0
        flag_cont = True
        
        kv = 0
        for v in pkeys:
            iv = v[0]; jv = v[1]
            
            flag_cont = True
            while kb < nb and flag_cont:
                ib = pxlist[kb].i; jb = pxlist[kb].j
                
                if debug:
                    print ("BPlist[kb=%2d](%2d,%2d) vs ijv[kv=%2d](%2d,%2d)" \
                           % (kb, ib, jb, kv, iv, jv))
                #
                
                if ib == iv and jb == jv:
                    
                    newroots += [pxlist[kb]]
                    flag_cont = False
                    
                #
                
                kb += 1
                #print ("next: -> ", pxlist[kb])
                
            #|endwhile
            kv += 1
            
        #|endfor
        
        if debug:
            for v in newroots:
                print (v.disp_Pair())
            #
            
        #
        
        return newroots
    #
    
    def remove_this_key(self, i, j, x):
        k = 0
        while k < len(x):
            xi = x[k][0]; xj = x[k][1]
            if i == xi and j == xj:
                del x[k]
                break
            else:
                k += 1
            #
            
        #|endwhile
        
        return x
    #
    
    
    def key_exists(self, i, j, x):
        flag_exists = False
        k = 0
        while k < len(x):
            xi = x[k][0]; xj = x[k][1]
            if i == xi and j == xj:
                flag_exists = True
                break
            else:
                k += 1
            #
            
        #|endwhile
        
        return flag_exists
    #
    
    
    
    def hasPKconflict(self, i, j, level):
        debug = False # True # 
        
        if debug:
            print ("enter noPKconflict((%2d,%2d, level=%d)" % (i, j, level))
        #
        
        flag_hasPK = False
        
        for v in self.pkeylevels[level]:
            
            ix = v[0]; jx = v[1]
            if debug:
                print ("ijx(%2d,%2d)" % (ix, jx))
            #
            
            if ix < i and i < jx and jx < j:
                if debug:
                    print ("ix(%2d) < i(%2d) < jx(%2d) < j(%2d)" \
                           % (ix, i, jx, j))
                #
                
                flag_hasPK = True
                
            elif i < ix and ix < j and j < jx:
                if debug:
                    print ("i(%2d) < ix(%2d) < j(%2d) < jx(%2d)" \
                           % (i, ix, j, jx))
                #
                
                flag_hasPK = True
            #
            
        #
        
        return flag_hasPK
    #
    
    
    def make_apLevelList(self, i, j, pkeys, kstart, level, debug_x = False):
        
        debug = debug_x
        
        #debug = True # False # 
        
        if level > 20:
            print ("Something is wrong! level = %n" % level)
            sys.exit(1)
        #
        
        if debug:
            print ("Enter make_apLevelList(ij(%2d, %2d), kstart(%2d), level(%2d)" \
                   % (i, j, kstart, level))
            #print ("stop at entrance"); sys.exit(0)
        #
        
        flag_cont = True
        k = kstart
        if not level in self.pkeylevels: 
            self.pkeylevels.update({level : []})
            if debug:
                print ("make: level(%d) key" % level)
            #
            
        #
        
        npkeys = len(pkeys)
        while (k < npkeys) and flag_cont:
            i1 = pkeys[k][0]; j1 = pkeys[k][1]
            #print ("ijl", (i1, j1))
            if self.hasPKconflict(i1, j1, level):
                #print ("skip", (i1, j1))
                k += 1
                continue
            #
            
            
            if debug: 
                print ("0: k  (%2d) i(%2d) < i1(%2d) < j1(%2d) < j(%2d)" \
                       % (k, i, i1, j1, j))
                #if k > 3: print ("stop at 0xx"); sys.exit(0)
            #
            
            if i < i1 and j1 < j:
                
                # check if this is the last item in the sequence of
                # pairs. It has to be treated different.
                
                if k == npkeys-1:
                    if not self.key_exists(i1, j1, self.pkeylevels[level]):
                        if debug: print ("1z, write: ", (i1, j1))
                        self.pkeylevels[level] += [(i1, j1)]
                        if debug:
                            print ("1z pkeylevels(%2d): %s" \
                                   % (level, self.pkeylevels[level]))
                            print ("1z + after ij2(%2d,%2d), k -> %d" \
                                   % (i1, j1, k))
                        #
                        
                    else:
                        if debug:
                            print ("1z (%2d,%2d) already present in list" \
                                   % (i1, j1))
                            print ("1z pkeylevels(%2d): %s" \
                                   % (level, self.pkeylevels[level]))
                        #
                        
                    #
                    
                    k += 1
                    break
                #
                
                if debug: print ("k = %d, npkeys = %d" % (k, npkeys))
                i2 = pkeys[k+1][0]; j2 = pkeys[k+1][1]
                if debug: 
                    print ("1: k+1(%2d) i(%2d) < i2(%2d) < j2(%2d) < j(%2d)" \
                           % (k + 1, i, i2, j2, j))
                #
                
                if j < i2:
                    
                    if debug: print ("1a: j(%d) < i2(%d)" % (j, i2))
                    if not self.key_exists(i1, j1, self.pkeylevels[level]):
                        if debug: print ("1a k + 0, write: ", (i1, j1))
                        self.pkeylevels[level] += [(i1, j1)]
                        k = self.make_apLevelList(i1, j1, pkeys, k+1, level+1)
                        if debug:
                            print ("1a pkeylevels(%2d): %s" \
                                   % (level, self.pkeylevels[level]))
                            print ("1a after ij1(%2d,%2d), k -> %d, %s" \
                                   % (i1, j1, k, pkeys[k]))
                        #
                        
                    else:
                        if debug:                        
                            print ("1a (%2d,%2d) already present in list" \
                                   % (i1, j1))
                            print ("1a current key k -> %d, %s" \
                                   % (k, pkeys[k])) 
                            print ("1a pkeylevels(%2d): %s" \
                                   % (level, self.pkeylevels[level]))
                        #
                        
                    #
                    
                    flag_cont = False
                    
                elif i1 < i2 and j2 < j1:
                    
                    # antiparallel case
                    if debug: 
                        print ("1c: i1(%2d) < i2(%2d) < j2(%2d) < j1(%2d)" \
                               % (i1, i2, j2, j1))
                    #
                    
                    if not self.key_exists(i1, j1, self.pkeylevels[level]):
                        if debug: print ("1c k + 0, write: ", (i1, j1))
                        self.pkeylevels[level] += [(i1, j1)]
                        #print (self.pkeylevels[level])
                        k = self.make_apLevelList(i1, j1, pkeys, k+1, level+1)
                        if debug: 
                            print ("1c pkeylevels(%2d): %s" \
                                   % (level, self.pkeylevels[level]))
                            if k < len(pkeys):
                                print ("1c after ij1(%2d,%2d), k -> %d, %s" \
                                       % (i1, j1, k, pkeys[k]))
                            else:
                                print ("1c after ij1(%2d,%2d), k -> %d last" \
                                       % (i1, j1, k))
                            #
                            
                        #
                        
                    else:
                        if debug: 
                            print ("1c (%2d,%2d) already present in list" \
                                   % (i1, j1))
                            if k < len(pkeys):
                                print ("1c after ij1(%2d,%2d), k -> %d, %s" \
                                       % (i1, j1, k, pkeys[k]))
                            else:
                                print ("1c after ij1(%2d,%2d), k -> %d last" \
                                       % (i1, j1, k))
                            #
                            
                            print ("1c pkeylevels(%2d): %s" \
                                   % (level, self.pkeylevels[level]))
                        #
                        
                    #
                    
                elif i1 < i2 and i2 < j1 and j1 < j2:
                    
                    # pk like case
                    dmn1 = self.prncpl_wts[pkeys[k]].npairs
                    dmn2 = self.prncpl_wts[pkeys[k+1]].npairs
                    
                    if debug: 
                        print ("1d: i1(%2d) < i2(%2d) < j1(%2d) < j2(%2d)" \
                               % (i1, i2, j1, j2))
                    #
                    
                    if dmn1 < dmn2: # (dmn2 >= dmn1)
                        if debug: print ("1d k + 1, write: ", (i2, j2))
                        self.pkeylevels[level] \
                            = self.remove_this_key(i1, j1, \
                                                   self.pkeylevels[level])
                        self.pkeylevels[level] += [(i2, j2)]
                        
                        k = self.make_apLevelList(i2, j2,
                                                  pkeys, k + 2, level + 1)
                        
                        if debug: 
                            print ("1d pkeylevels(%2d): %s" \
                                   % (level, self.pkeylevels[level]))
                            if k < len(pkeys):
                                print ("1d after ij1(%2d,%2d), k -> %d, %s" \
                                       % (i1, j1, k, pkeys[k]))
                            else:
                                print ("1d after ij1(%2d,%2d), k -> %d last" \
                                       % (i1, j1, k))
                            #
                            
                            print ("1d pkeylevels(%2d): %s" \
                                   % (level, self.pkeylevels[level]))
                        #
                        
                    else:
                        if debug: print ("1d k + 0, write: ", (i1, j1))
                        if not self.key_exists(i1, j1, self.pkeylevels[level]):
                            self.pkeylevels[level] += [(i1, j1)]
                            k = self.make_apLevelList(i1, j1,
                                                      pkeys, k+1, level+1)
                        else:
                            if debug: 
                                print ("1d (%2d,%2d) already present in list" \
                                       % (i1, j1))
                                if k < len(pkeys):
                                    print ("1d after ij1(%2d,%2d), k -> %d, %s" \
                                           % (i1, j1, k, pkeys[k]))
                                else:
                                    print ("1d after ij1(%2d,%2d), k -> %d last" \
                                           % (i1, j1, k))
                                #
                                print ("1d pkeylevels(%2d): %s" \
                                       % (level, self.pkeylevels[level]))
                            #
                            
                        #
                        
                        if debug: 
                            print ("1d pkeylevels(%2d): %s" \
                                   % (level, self.pkeylevels[level]))
                            print ("1d after ij2(%2d,%2d), k -> %d, %s" \
                                   % (i2, j2, k, pkeys[k]))
                        #
                        
                    #
                    
                    
                elif j1 < i2:
                    
                    if debug: print ("1b; j1(%d) < i2(%d)" % (j1, i2))
                    if not self.key_exists(i1, j1, self.pkeylevels[level]):
                        if debug: print ("1b k + 0, write: ", (i1, j1))
                        self.pkeylevels[level] += [(i1, j1)]
                        k = self.make_apLevelList(i1, j1, pkeys, k+1, level+1)
                        # k + 1 because we're looking at the next level
                        if debug: 
                            print ("1b pkeylevels(%2d): %s" \
                                   % (level, self.pkeylevels[level]))
                            print ("1b after ij1(%2d,%2d), k -> %d, %s" \
                                   % (i1, j1, k, pkeys[k]))
                        #
                        
                    else:
                        if debug: 
                            print ("1b (%2d,%2d) already present in list" \
                                   % (i1, j1))
                            print ("1b current key k -> %d, %s" \
                                   % (k, pkeys[k])) 
                            print ("1b pkeylevels(%2d): %s" \
                                   % (level, self.pkeylevels[level]))
                        #
                        
                    #
                    
                    
                    if i < i2 and j2 < j:
                        if debug:
                            print ("1b+ i(%2d) < i2(%2d) < j2(%2d) < j(%2d)"
                                   % (i, i2, j2, j))
                        #
                        
                        if not self.key_exists(i2, j2, self.pkeylevels[level]):
                            if debug: print ("1b + k + 1, write: ", (i2, j2))
                            self.pkeylevels[level] += [(i2, j2)]
                            k = self.make_apLevelList(i2, j2,
                                                      pkeys, k+1, level+1)
                            if debug: 
                                print ("1b pkeylevels(%2d): %s" \
                                       % (level, self.pkeylevels[level]))
                                print ("1b + after ij2(%2d,%2d), k -> %d" \
                                       % (i2, j2, k))
                            #
                            
                        else:
                            if debug: 
                                print ("1b (%2d,%2d) already present in list" \
                                       % (i2, j2))
                                if k < len(pkeys):
                                    print ("1b current key k -> %d, %s" \
                                           % (k, pkeys[k]))
                                else:
                                    print ("1b current key k -> %d exit" \
                                           % (k))
                                #
                                
                                print ("1b pkeylevels(%2d): %s" \
                                       % (level, self.pkeylevels[level]))
                            #
                            
                        #
                        
                        #k += 1
                    #
                    
                    
                else:
                    print ("full pxlist_ref: (input to Vienna.Vstruct())")
                    for k in range(len(self.pxlist_ref)):
                        print ("%2d  %s" % (k, self.pxlist_ref[k].disp_Pair()))
                    #
                    
                    print ("pkeys: (input to make_apLevelList)")
                    for k in range(len(pkeys)):
                        print ("%2d, (%3d,%3d)" % (k, pkeys[k][0], pkeys[k][1]))
                    #
                    
                    print ("1e")
                    print ("stop at 1e"); sys.exit(0)
                    k += 1
                #
                
            else:
                if debug: 
                    print ("2: level %d: ij1(%2d,%2d) finished" \
                           % (level, i1, j1))
                #
                
                flag_cont = False
                
            #
            
            if debug: 
                print ("end of loop, k = %d, n = %d" % (k, len(pkeys)))
            #
            
        #
        
        return k
        
    #
    
    
    def assign_PKdirection(self, plist): # plist <- self.PKlist
        """#########
        !!!!!
        
        Decides if the stem is parallel or antiparallel. Here are a
        couple of examples of tests that I ran from Threads.py.
        
        example 1
                   0         10        20        30        40 
                   |    .    |    .    |    .    |    .    |  
        ss_seq  = ".A..B..C..D...........a..b..c..d.."
        result:
        (   1,   22)[p]
        (   4,   25)[p]
        (   7,   28)[p]
        (  10,   31)[p]
        
        example 2:
                   0         10        20        30        40 
                   |    .    |    .    |    .    |    .    |  
        ss_seq  = ".AAA...BBB...CCC...aaa....bbb....ccc...."
        input:
        (   1,   21)[a]
        (   2,   20)[a]
        (   3,   19)[a]
        (   7,   28)[a]
        (   8,   27)[a]
        (   9,   26)[a]
        (  13,   35)[a]
        (  14,   34)[a]
        (  15,   33)[a]
        
        Results:
        structural sequence
        cxxxcccxxxcccxxxcccyyycccyyycccyyyc
        .AAA...BBB...CCC...aaa...bbb...ccc.
        secondary structure
        (    1,   21)[a]
        (    2,   20)[a]
        (    3,   19)[a]
        secondary structure roots
        (    1,   21)[a]
        pseudoknot linkages: 
        (    7,   27)[a]
        (    8,   26)[a]
        (    9,   25)[a]
        (   13,   33)[a]
        (   14,   32)[a]
        (   15,   31)[a]
        pseudoknot linkage roots
        (    7,   27)[a]
        (   13,   33)[a]
        
        example 3:
                   0         10        20        30        40 
                   |    .    |    .    |    .    |    .    |  
        ss_seq  = ".AAA...BBB...aaa...CCC....bbb....ccc...." # same
        ss_seq  = ".(((...BBB...)))...(((....bbb....)))...." # same
        Input:
        (   1,   15)[a]
        (   2,   14)[a]
        (   3,   13)[a]
        (   7,   28)[a]
        (   8,   27)[a]
        (   9,   26)[a]
        (  19,   35)[a]
        (  20,   34)[a]
        (  21,   33)[a]
        
        Results:
        structural sequence
        cxxxcccxxxcccyyycccxxxcccyyycccyyyc
        .AAA...BBB...aaa...CCC...bbb...ccc.
        secondary structure
        (    1,   15)[a]
        (    2,   14)[a]
        (    3,   13)[a]
        (   19,   33)[a]
        (   20,   32)[a]
        (   21,   31)[a]
        secondary structure roots
        (    1,   15)[a]
        (   19,   33)[a]
        pseudoknot linkages: 
        (    7,   27)[a]
        (    8,   26)[a]
        (    9,   25)[a]
        pseudoknot linkage roots
        (    7,   27)[a]
        (   19,   33)[a]
        
        It appears to be able to distinguish consequitive pairs that
        would satisfy a parallel pattern and they do not have to be
        contiguous either (as shown in example 1). However, if it is
        just an odd sort of pseudoknot where you have antiparallel
        stems but a kind of ladder structure, it will not call this
        ladder arrangement "parallel" (assigning the "p" to the
        brackets [a] or [p].
        
        190219: appears to work ok with the more odd or difficult
        examples I could think of.
        
        """
        debug_ad = False # True # 
        if debug_ad:
            print ("Enter assign_PKdirection: ")
            print ("plist: (before)")
            for vv in plist:
                print (vv.disp_Pair())
            #|endfor
            
            #print ("stop at 0.0 in assign_PKdirection"), sys.exit(0)
        #
        
        # since the lists are ordered, we can check through
        # self.BPlist and identify parallel stems by the simple rule
        # of parallel-ness: ibp0 < ibp1 < jbp0 < jbp1.
        
        n = len(self.BPlist)
        pknotes = []
        dmnnotes = []
        zone_lock = False
        i_lock = self.N
        j_lock = 0
        
        
        for k in range(0, n - 1):
            
            ibp0 = self.BPlist[k  ].i; jbp0 = self.BPlist[k  ].j
            ibp1 = self.BPlist[k+1].i; jbp1 = self.BPlist[k+1].j
            
            zone_lock = (jbp1 <= j_lock)
            if debug_ad:
                print ("ibp0(%2d) < ibp1(%2d) < jbp0(%2d) < jbp1(%2d)" \
                    % (ibp0, ibp1, jbp0, jbp1))
                print ("ij_lock(%2d,%2d) = %s" % (i_lock, j_lock, zone_lock))
            #
            
            if ibp0 < ibp1 and ibp1 < jbp0 and jbp0 < jbp1:
                flag_pass = False
                
                test_p2_is_pp    = False
                test_p2_is_shift = False
                if k <= n-3:
                    ibp2 = self.BPlist[k+2].i; jbp2 = self.BPlist[k+2].j
                    # ..(.....[.....<......).....].....>..
                    #   ^     ^     ^      ^     ^     ^
                    # ibp0  ibp1   ibp2  jbp0  jbp1  jbp2
                    test_p2_is_pp = (ibp1 < ibp2 and jbp1 < jbp2 and ibp2 < jbp1)
                    test_p2_is_shift = (jbp1 < ibp2)
                #
                
                test_m1_is_pp    = False
                test_m1_is_shift = False
                if k > 0:
                    # ..(.....[.....<......).....].....>..
                    #   ^     ^     ^      ^     ^     ^
                    # ibm1  ibp0   ibp1  jbm1  jbp0  jbp1
                    ibm1 = self.BPlist[k-1].i; jbm1 = self.BPlist[k-1].j
                    test_m1_is_pp = (ibm1 < ibp0 and jbm1 < jbp0 and ibp0 < jbm1)
                    test_m1_is_shift = (jbm1 < ibp0)
                #
                
                if test_m1_is_pp or test_p2_is_pp:
                    # if at least 1 of these is also a parallel stem,
                    # then this passes. In short There must be at
                    # least three of these together to make this work.
                    flag_pass = True
                #
                
                
                # parallel relationship. This is the nearest
                # neighboring Pair and it is in BPlist (the dominant
                # one). Therefore, it should have priority and be
                # assigned properly.
                
                if flag_pass:
                    if debug_ad:
                        print ("Reassign to p:")
                        print (self.BPlist[k])
                        print (self.BPlist[k+1])
                    #
                    
                    self.BPlist[k  ].v = 'p'
                    self.BPlist[k+1].v = 'p'
                    
                else:
                    print ("zone_lock: ", zone_lock)
                    if not zone_lock:
                        pknotes += [self.BPlist[k+1]]
                        iref = self.BPlist[k+1].i
                        jref = self.BPlist[k+1].j
                        print ("Search for other pks within pk(%d,%d), next k = %d" \
                               % (iref, jref, k+2))
                        dmnroot = self.checkForDmnRoot(k + 2, jref)
                        if not dmnroot == None:
                            i_lock = dmnroot.i
                            j_lock = dmnroot.j
                            dmnnotes += [dmnroot]
                        #
                        
                    else:
                        
                        print ("zone locked between %d and %d, <- (%d,%d)" \
                               % (i_lock, j_lock, ibp1, jbp1))
                    #
                    
                #
                
            #
            
        #|endfor
        
        if debug_ad:
            print ("BPlist: (after checking)")
            for vv in self.BPlist:
                print (vv.disp_Pair())
            #|endfor
            
            print ("PKlist: (before checking)")
            for vv in self.PKlist:
                print (vv.disp_Pair())
            #|endfor
            
            print ("pknotes:")
            for vv in pknotes:
                print (vv.disp_Pair())
            #|endfor
            
            print ("stop at 0.1 in assign_PKdirection"), sys.exit(0)
        #
        
        # pre-search for parallel stems before 
        ee = {}
        n = len(pknotes) - 1
        kr = 0
        pktemp = []
        while kr < len(pknotes):
            ipk  = pknotes[kr  ].i; jpk  = pknotes[kr  ].j
            if debug_ad:
                print ("ipk(%2d,%2d), kr = %d" % (ipk, jpk, kr))
            #
            
            for kt in range(0, len(self.BPlist)):
                ibp = self.BPlist[kt].i; jbp = self.BPlist[kt].j
                if debug_ad:
                    print ("ipk(%2d) < ibp(%2d) < jpk(%2d) < jbp(%2d)" \
                        % (ipk, ibp, jpk, jbp))
                    print ("ee = ", ee)
                #
                
                #   ..[....(.....]....)...
                #     ^    ^     ^    ^
                #    ipk  ibp   jpk  jbp
                
                if ipk < ibp and ibp < jpk and jpk < jbp:
                    if (ipk, jpk) in ee:
                        ee[(ipk, jpk)] += [(ibp, jbp)]
                    else:
                        ee.update({(ipk, jpk) : [(ibp, jbp)]})
                        pktemp += [(ipk, jpk)]
                        
                    #
                    
                elif jpk < ibp:
                    break
                #
                
            #|endfor
            
            kr += 1
            
        #|endwhile
        
        
        #print (self.PKlist)
        
        # pre-search for parallel stems before 
        dd = {}
        n = len(self.BPlist) - 1
        kr = 0
        bptemp = []
        while kr < len(self.BPlist):
            ibp  = self.BPlist[kr  ].i; jbp  = self.BPlist[kr  ].j
            if debug_ad:
                print ("ibp(%2d,%2d), kr = %d" % (ibp, jbp, kr))
            #
            
            for kt in range(0, len(self.PKlist)):
                ipk = self.PKlist[kt].i; jpk = self.PKlist[kt].j
                if debug_ad:
                    print ("ibp(%2d) < ipk(%2d) < jbp(%2d) < jpk(%2d)" \
                        % (ibp, ipk, jbp, jpk))
                    print ("dd = ", dd)
                #
                
                #   ..[....(.....]....)...
                #     ^    ^     ^    ^
                #    ipk  ibp   jpk  jbp
                
                if ibp < ipk and jbp < jpk and ipk < jbp:
                    if (ibp, jbp) in dd:
                        dd[(ibp, jbp)] += [(ipk, jpk)]
                    else:
                        dd.update({(ibp, jbp) : [(ipk, jpk)]})
                        bptemp += [(ibp, jbp)]
                        
                    #
                    
                elif jbp < ipk:
                    break
                #
                
            #|endfor
            
            kr += 1
            
        #|endwhile
        
        if debug_ad:
            print ("dd results(before)")
            for ddk in bptemp:
                print (ddk, dd[ddk])
            #|endfor
            
            print ("stop at 0.2 in assign_PKdirection"), sys.exit(0)
            
            print ("ee results(before)")
            for eek in pktemp:
                print (eek, ee[eek])
            #|endfor
            
            print ("stop at 0.2 in assign_PKdirection"), sys.exit(0)
        #
        
        # Now we have to filter the initial data so that only parallel
        # stems are examine.
        
        for ddk in bptemp:
            pkpptest = dd[ddk]
            
            for kt in range(0, len(pkpptest)-1):
                ipp1 = pkpptest[kt  ][0]; jpp1 = pkpptest[kt  ][1]
                ipp2 = pkpptest[kt+1][0]; jpp2 = pkpptest[kt+1][1]
                
                if not (ipp1 < ipp2 and jpp1 < jpp2 and ipp2 < jpp1):
                    # if it is not parallel (ipp1 < ipp2 < jpp1 <
                    # jpp2), then we are not interested at this
                    # point. We are only interested in cases where it
                    # is definitely a parallel stem.
                    dd[ddk] = []
                    break
                #
                
            #|endfor
            
        #|endwhile
        
        if debug_ad:
            self.disp_allLists("initial list:", False)
            #print ("stop at 0.3 in assign_PKdirection"), sys.exit(0)
        #
        
        # delete results that don't satisfy the requirement
        
        for k in range(0, len(bptemp)):
            if bptemp[k] in dd:
                #print ("bptemp: ", bptemp, dd[bptemp[k]])
                if len(dd[bptemp[k]]) == 0: 
                    del dd[bptemp[k]]
                #
                
            #
            
        #|endfor
        
        if debug_ad:
            self.disp_allLists("partiall refined list:", False)
            #print ("stop at 0.4 in assign_PKdirection"), sys.exit(0)
        #
        
        #sys.exit(0)
        for ddk in dd.keys():
            i = ddk[0]; j = ddk[1]
            
            # find the corresponding Pair class object in BPlist
            for k in range(0, len(self.BPlist)):
                if i == self.BPlist[k].i and j == self.BPlist[k].j:
                    self.BPlist[k].v = 'p'
                #
                
            #|endfor
            
            for pkk in dd[ddk]:
                # print (pkk)
                ipk = pkk[0]; jpk = pkk[1]
                # find the corresponding Pair class object in PKlist
                for k in range(0, len(self.PKlist)):
                    if ipk == self.PKlist[k].i and jpk == self.PKlist[k].j:
                        self.PKlist[k].v = 'p'
                    #
                    
                #|endfor
                
            #|endfor
            
        #|endfor
        
        
        if debug_ad:
            print ("dd results(after)")
            for ddk in dd.keys():
                print (ddk, dd[ddk])
            #|endfor
            
            print ("BPlist:")
            for vv in self.BPlist:
                print (vv.disp_Pair())
            #|endfor
            
            print ("PKlist:")
            for vv in self.PKlist:
                print (vv.disp_Pair())
            #|endfor
            
            #print ("stop at 0.5 in assign_PKdirection"), sys.exit(0)
        #
        
        
        # Note: "plist" is the current PK list, not all bps
        
        plist = self.sortvsList(plist, 'i')
        # plist.sort(key=self.getKey_i)
        for k in range(0, len(plist)-1):
            i_kp0 = plist[k  ].i; j_kp0 = plist[k  ].j
            i_kp1 = plist[k+1].i; j_kp1 = plist[k+1].j
            if debug_ad:
                print ("a) ij_kp0[=%2d](%2d,%2d) ... ij_kp1[=%2d](%2d,%2d)" % \
                    (k, i_kp0, j_kp0, k+1, i_kp1, j_kp1));
            #
            
            if i_kp0 <  i_kp1 and i_kp1 < j_kp0 and j_kp0 < j_kp1:
                """
                First, if it is a parallel stem, then it must be that it looks
                like a PK with the linkage to the right side.
                
                "...A.....B.........a........B..." 
                    ^     ^         ^        ^
                 i_kp0  i_kp1     i_kp1    j_kp0
                
                """
                    
                
                test_km1 = False; test_kp2 = False
                if k == 0 and len(plist) == 2:
                    # we only have these two pairs so there is no k-1
                    # or k+2 case to test this against.
                    test_km1 = True
                    test_kp2 = True
                #
                
                """
                We are _given_ the following
                
                i_kp0 < i_kp1 < j_kp0 < j_kp1
                
                ============================================
                this case should be accepted
                    k=0
                       i_kp1  i_kp2      j_kp2      j_kp1   
                         v      v          v          v            
                      ..AB......[[........]].........ab...
                        A                            A             
                      i_kp0                        j_kp0
                
                    i_kp0 < i_kp1 < j_kp0 < j_kp1 (given)
                    i_kp1 < i_kp2 < j_kp2 < j_kp1 **
                
                
                this case should be rejected
                    k=1
                       i_km1  i_kp1      j_kp1     j_km1   
                        v       v          v         v            
                      ..AB......[[........]].........ab...
                         A       ^        ^           A             
                       i_kp0   i_kp2    j_kp2       j_kp0
                    
                    i_km1 < i_kp0 < j_km1 < j_kp0
                    i_kp0 < i_kp1 < j_kp1 < j_kp0 fails here
                    i_kp1 < i_kp2 < j_kp2 < j_kp1 **
                
                ============================================
                this case should be rejected
                    k=1
                      i_km1   i_kp1     j_kp1       j_km1   
                        v       v         v           v            
                      ..[[......AB........ab.........]]...
                         A       ^         ^         A             
                       i_kp0   i_kp2     j_kp2     j_kp0
                
                    
                    i_km1 < i_kp0 < j_kp0 < j_kp1 **
                    i_kp0 < i_kp1 < j_kp1 < j_kp0 fails here
                    i_kp1 < i_kp2 < j_kp1 < j_kp2 
                
                this case should be accepted
                    k=2
                       i_km1   i_kp1     j_kp1     j_km1   
                         v       v         v         v            
                      ..[[......AB........ab.........]]...
                                A         A             
                              i_kp0     j_kp0
                
                    i_km1 < i_kp0 < j_kp0 < j_kp1 **
                    i_kp0 < i_kp1 < j_kp0 < j_kp1
                
                
                ============================================
                this case should be accepted
                    k=0
                       i_kp1       j_kp1 i_kp2      j_kp2   
                         v           v     v          v            
                      ..AB..........ab.....[[........]]...
                        A           A             
                      i_kp0       j_kp0
                    
                    i_kp0 < i_kp1 < j_kp0 < j_kp1 (given)
                    i_kp1 < j_kp1 < i_kp2 < j_kp2 **
                            -------------
                
                    
                ============================================
                this case should be accepted
                    k=2
                       i_km1     j_km1   i_kp1        j_kp1   
                         v        v         v           v            
                      ..[[........]].......AB..........ab.
                                           A           A             
                                         i_kp0       j_kp0
                    
                    i_km1 < j_km1 < i_kp0 < j_kp0 **
                    i_kp0 < i_kp1 < j_kp0 < j_kp1 (given)
                            -------------
                
                ============================================
                this case should be rejected
                    k=0
                       i_kp1     i_kp2   j_kp1        j_kp2   
                         v        v        v            v            
                      ..AA........BB.......aa..........bb.
                        A                   A             
                      i_kp0               j_kp0
                    
                    i_kp0 < i_kp1 < j_kp1 < j_kp0 fails here
                            -------------
                    i_kp1 < i_kp2 < j_kp1 < j_kp2 
                
                this case should be rejected
                    k=1
                       i_km1   i_kp1,2    j_km1      j_kp1,2   
                        v         vx        v          xv            
                      ..AA........BB.......aa..........bb.
                         A                 A             
                       i_kp0             j_kp0
                    
                    i_km1 < j_kp0 < i_kp0 < j_km1 **
                            -------------
                    i_kp0 < i_kp1 < j_kp0 < j_kp1 (given)
                    i_kp1 < i_kp2 < j_kp2 < j_kp1 **
                            -------------
                
                this case should be rejected
                    k=2
                       i_km1     i_kp1   j_km1      j_kp1   
                         v         v       v           v            
                      ..AA........BB.......aa..........bb.
                                  A                     A             
                                i_kp0                 j_kp0
                    
                    i_km1 < j_kp0 < i_km1 < j_kp0 
                    i_kp0 < i_kp1 < j_kp1 < j_kp0 fails here
                            -------------
                
                """
                if k > 0:
                    i_km1 = plist[k-1].i; j_km1 = plist[k-1].j
                    """
                    We are _given_ i_kp0 < i_kp1 < j_kp0 < j_kp1
                    
                    k=1,2 
                    
                    Both cases fail the above condition so they don't
                    end up here
                    
                    k=3 (first case satisfied)
                         i_kp0  i_kp1        j_kp0       j_kp1
                           V      v           V            v
                      ..[[[[......ABCD........]]]].........abcd...  (i_km1 < i_kp0 and j_kp0 < j_km1) <===
                          ^                    ^                    (i_kp0 < i_kp1  <  jkp0  < j_kp1)[*]
                        i_km1                j_km1                  should be rejected
                    
                    k=4
                               i_kp0,1,2                j_kp0,1,2         
                                  Vvv                      Vvv           
                      ..[[[[......ABCD........]]]].........abcd...  (i_km1 < i_kp0  < j_km1 < j_kp0)
                           ^                  ^                     (i_kp0 < i_kp1  < j_kp0 < j_kp1)[*]
                         i_km1              j_km1                   (i_kp1 < i_kp2  < j_kp1 < j_kp2)
                                                                    should be accepted
                    
                    k=5
                                i_kp0,1,2                j_kp0,1,2         
                                   Vvv                      Vvv            
                      ..[[[[......ABCD........]]]].........abcd...  (i_km1 < i_kp0  < j_km1 < j_kp0)
                                  ^                        ^        (i_kp0 < i_kp1  < j_kp0 < j_kp1)[*]
                                i_km1                    j_km1      (i_kp1 < i_kp2  < j_kp1 < j_kp2)
                                                                    should be accepted
                    
                    k=6
                                 i_kp0,1                   j_kp0,1         
                                    Vv                       Vv            
                      ..[[[[......ABCD........]]]].........abcd...  (i_km1 < i_kp0  < j_km1 < j_kp0)
                                   ^                        ^       (i_kp0 < i_kp1  < j_kp0 < j_kp1)[*]
                                 i_km1                    j_km1     should be accepted
                    
                    """
                    if debug_ad:
                        print ("b) ij_kp0[=%2d](%2d,%2d) ... ij_km1[=%2d](%2d,%2d)" % \
                            (k, i_kp0, j_kp0, k-1, i_km1, j_km1));
                    #
                    """
                    This is the former method, I think it is not so secure
                    
                    k       i_km1 < i_kp0    i_kp0 < j_kp1     j_kp0 < j_km1    j_km1 < j_kp1   ijkp01
                    3             T                T                 T                 T       --> a
                    4             T                T                 F                 T       --> p
                    5             T                T                 F                 T       --> p
                    6             T                T                 F                 T       --> p
                    #if not (i_km1 < i_kp0 and i_kp0 < j_kp1 and j_kp0 < j_km1 and j_km1 < j_kp1):
                    #    test_km1 = True
                    ##
                    """
                    
                    if (not (i_km1 < i_kp0 and j_kp0 < j_km1) or # <=== (i_km1 < i_kp0 and j_kp0 < j_km1) 
                        (j_km1 < i_kp0) or 
                        (i_km1 < i_kp0  and j_km1 < j_kp0 and i_km1 < j_kp0)):
                        test_km1 = True
                    #
                    
                    
                #
                
                if debug_ad:
                    print ("k= %d, len(plist) = %d" % (k, len(plist)))
                #
                
                if k + 2 < len(plist):
                    i_kp2 = plist[k+2].i; j_kp2 = plist[k+2].j
                    if debug_ad:
                        print ("c) ij_kp0[=%2d](%2d,%2d) ... ij_kp1[=%2d](%2d,%2d) ... ij_kp2[=%2d](%2d,%2d)" % \
                            (k, i_kp0, j_kp0, k+1, i_kp1, j_kp1, k+2, i_kp2, j_kp2));
                    #
                    
                    """
                    This is the former method, I think it is not so secure
                    
                    k       i_kp0 < i_kp1    i_kp1 < i_kp2     j_kp0 < j_kp2     j_kp2 < j_kp1  ijkp01
                    3             T                T                 T                 F       --> p
                    4             T                T                 T                 F       --> p
                    5             T                T                 T                 F       --> p
                    #if not (i_kp0 < i_kp1 and i_kp1 < i_kp2 and j_kp0 < j_kp2 and j_kp2 < j_kp1):
                    #    test_kp2 = True
                    ##
                    """
                    if (not(i_kp1 < i_kp2 and j_kp2 < j_kp1) or # <=== (i_kp1 < i_kp2 and j_kp2 < j_kp1) 
                        (j_kp1 < i_kp2) or 
                        (i_kp1 < i_kp2  and j_kp1 < j_kp2 and i_kp2 < j_kp1)):
                        test_kp2 = True
                    #
                    
                #
                
                if debug_ad:
                    print ("test_km1(%s), test_kp2(%s)" % (test_km1, test_kp2))
                #
                
                if test_km1 or test_kp2:
                    plist[k  ].v = 'p'
                    plist[k+1].v = 'p'
                #
                
            #
            
        #|endfor
        
        if debug_ad:
            print ("PKlist: (after)")
            for vv in plist:
                print (vv.disp_Pair())
            #|endfor
            
            print ("BPlist: (after)")
            for vv in self.BPlist:
                print (vv.disp_Pair())
            #|endfor
            
            print ("Exit assign_PKdirection: ")
            #print ("stop at end of assign_PKdirection"); sys.exit(0)
        #
        
        return plist
    # #########
    
    
    def findroots(self, prlist, debug = False):
        debug_findroots = debug
        if debug_findroots:
            print ("Enter findroots:")
            print ("prlist:")
            for vv in prlist:
                print (vv.disp_Pair())
            #|endfor
            
        #
        
        nprlist = []
        for pr_k in prlist:
            # make sure we don't just make a pointer
            nprlist += [pr_k]
        #|endfor
        
        k = 0
        while k < (len(nprlist) - 1):
            i_ref = nprlist[k].i; j_ref = nprlist[k].j; tp_ref = nprlist[k].v
            if debug_findroots:
                print ("ij_ref: ", nprlist[k].disp_Pair())
            #
            
            ll = 0
            while ll < len(nprlist):
                i_t = nprlist[ll].i; j_t = nprlist[ll].j; tp_t = nprlist[ll].v
                if i_t == i_ref and j_t == j_ref:
                    ll += 1
                    continue
                #
                
                if 0: #debug_findroots:
                    print ("k = %2d, ll = %2d: i_ref(%2d) < i_t(%2d) and j_t(%2d) < j_ref(%2d)" \
                        % (k, ll, i_ref,  i_t, j_t, j_ref))
                #
                
                delete_t = False
                if not (tp_t == 'p' and tp_ref == 'p'):
                    # still not sure if this is the only condition or
                    # if this uniquely defines the problem of a
                    # parallel stem, but this is what I can currently
                    # think of to make it so the PK roots contain all
                    # these elements.
                    if i_ref < i_t and j_t < j_ref:
                        if (tp_t == 'a' and tp_ref == 'a'):
                            delete_t = True
                        #
                        
                    #
                    
                #
                
                if delete_t:    
                    if debug_findroots:
                        print ("deleting ", nprlist[ll].disp_Pair())
                    #
                    del nprlist[ll]
                else:
                    ll += 1
                #
                
            #|endwhile
            
            k+= 1
            
        #|endwhile
        
        if debug_findroots:
            print ("end of findroots: vvvvvvvvvvvvv")
            print ("nprlist: ")
            for bproot in nprlist:
                print (bproot.disp_Pair())
            #|endfor
            
            print ("prlist: ")
            for bproot in prlist:
                print (bproot.disp_Pair())
            #|endfor
            
            print ("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
            #print ("stop at the exit point from findroots"); sys.exit(0)
        #
        
        return nprlist
    #
    
    def test_entwined(self, pkroots, ssroots, debug_test_entwined = False):
        """
        Consider this example (a structure in Threads.py)
        
           ".((....AA...))..BB....aa....CC...bb....cc.."
        
        In this routine, ssroots and pkroots leave us with the
        following paired down structure
        
           ".((....AA...))..BB....aa....CC...bb....cc.."
        
             |                  |                   |
             V                  V                   V
        
           ".(.....A.....)..B......a....C.....b.....c.."
        
        and this routine isolates the PK roots
        
           "................B...........C.....b.....c.."
        
        So the output from this routine reference to the PK roots
        corresponding to Bb and Cc.
        
        Likewise, following the same reasoning, this next structure
        is processed in the following way.
        
           ".((....AA...))..BB....CC....aa...bb....cc.." 
        
             |                  |                   |
             V                  V                   V
        
           ".(.....A.....)..B.....C......a....b.....c.."
        
             |                  |                   |
             V                  V                   V
        
           "................B.....C...........b.....c.."
        """
        
        debug_test_entwined = False # True # 
        stop_at_end         = False # True # 
        if debug_test_entwined:
            print ("Enter test_entwined: (before)")
            print ("ssroots: ")
            for vv in ssroots:
                print (vv.disp_Pair())
            #|endfor
            
            print ("pkroots: ")
            for vv in pkroots:
                print (vv.disp_Pair())
            #|endfor
            
            print ("------")
        #
        
        ss_pssbl = []  # possible secondary structure roots
        for pklk in pkroots:
            ss_pssbl += [pklk]
        #|endfor
        
        kp = 0
        while kp < len(ss_pssbl):
            ss_entwined = False
            i_p = ss_pssbl[kp].i; j_p = ss_pssbl[kp].j; v_p = ss_pssbl[kp].v
            kr = 0
            for kr in range(0, len(ssroots)):
                # compare the secondary structure roots with the 
                i_r = ssroots[kr].i; j_r = ssroots[kr].j; v_r = ssroots[kr].v
                
                """
                 ...[.....(.......].......)..     ...(.....[.......).......]..
                    i_p   i_r     j_p     j_r        i_r   i_p     j_r     j_p
                """
                if debug_test_entwined:
                    if (i_p < i_r and j_p < j_r and i_r < j_p):
                        # i_p < i_r < j_p < j_r
                        print ("1 i_p(%2d) < i_r(%2d) & j_p(%2d) < j_r(%2d) & i_r(%2d) < j_p(%2d), v_p(%s) v_r(%s)" \
                            % (i_p, i_r, j_p, j_r, i_r, j_p, v_p, v_r))
                    elif (i_r < i_p and j_r < j_p and i_p < j_r):
                        # i_r < i_p < j_r < j_p
                        print ("2 i_r(%2d) < i_p(%2d) & j_r(%2d) < j_p(%2d) & i_p(%2d) < j_r(%2d), v_p(%s) v_r(%s)" \
                            % (i_r, i_p, j_r, j_p, i_p, j_r, v_p, v_r))
                    elif (i_r < i_p and j_p < j_r): 
                        print ("3 i_r(%2d) < i_p(%2d) < j_p(%2d) < j_r(%2d), v_p(%s) v_r(%s)" \
                            % (i_r, i_p, j_p, j_r, v_p, v_r))
                    elif (i_p < i_r and j_r < j_p): 
                        print ("4 i_p(%2d) < i_r(%2d) < j_r(%2d) < j_p(%2d), v_p(%s) v_r(%s)" \
                            % (i_p, i_r, j_r, j_p, v_p, v_r))
                    elif (j_r < i_p): 
                        print ("5 i_r(%2d) < j_r(%2d) < i_p(%2d) < j_p(%2d), v_p(%s) v_r(%s)" \
                            % (i_r, j_r, i_p, j_p, v_p, v_r))
                    elif (j_p < i_r):
                        print ("6 i_p(%2d) < j_p(%2d) < i_r(%2d) < j_r(%2d), v_p(%s) v_r(%s)" \
                            % (i_p, j_p, i_r, j_r, v_p, v_r))
                    elif (i_r == i_p and j_r == j_p):
                        print ("7 ij_p(%2d,%2d) <=> ij_r(%2d,%2d), v_p(%s) v_r(%s)" \
                            % (i_p, j_p, i_r, j_r, v_p, v_r))
                    else:
                        print ("don't know what is going on")
                        print ("8 ij_p(%2d,%2d) <??> ij_r(%2d,%2d), v_p(%s) v_r(%s)" \
                            % (i_p, j_p, i_r, j_r, v_p, v_r))
                        sys.exit(0)
                    #
                    
                #
                
                if ((i_p < i_r and j_p < j_r and i_r < j_p) or
                    (i_r < i_p and j_r < j_p and i_p < j_r)):
                    if debug_test_entwined:
                        print ("deleting ss_pssbl: ", ss_pssbl[kp].disp_Pair())
                    #
                    
                    del ss_pssbl[kp] # pk meshed with ssroots
                    ss_entwined = True
                    break
                    
                #
                
            #|endfor
            
            if not ss_entwined:
                kp += 1
            #
            
        #|endwhile
        
        if debug_test_entwined:
            print ("check overlaps with BPlist:  vvvvvvvvvvvvvv")
        #
        
        kv = 0
        while kv < len(ss_pssbl):
            ic = ss_pssbl[kv].i; jc = ss_pssbl[kv].j
            if debug_test_entwined:
                print ("ijc(%2d,%2d)" % (ic, jc))
            #
            
            flag_del = False
            for bpk in self.BPlist:
                iv = bpk.i; jv = bpk.j
                if    iv < ic and jv < jc and ic < jv:
                    if debug_test_entwined:
                        print ("overlap with region iv(%2d) < ic(%2d) < jv(%2d) < jc(%2d)" \
                            % (iv, ic, jv, jc))
                    #
                    del ss_pssbl[kv]
                    flag_del = True
                    break
                elif ic < iv and jc < jv and iv < jc:
                    if debug_test_entwined:
                        print ("overlap with region ic(%2d) < iv(%2d) < jc(%2d) < jv(%2d)" \
                            % (ic, iv, jc, jv))
                    #
                    del ss_pssbl[kv]
                    flag_del = True
                    break
                #
                
            #|endfor
            
            if not flag_del:
                kv += 1
            #
            
        #|endwhile
        
        if debug_test_entwined:
            print ("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
            print ("Exit test_entwined: (after)")
            print ("ssroots: ")
            for vv in ssroots:
                print (vv.disp_Pair())
            #|endfor
            
            print ("pkroots: ")
            for vv in pkroots:
                print (vv.disp_Pair())
            #|endfor
            
            print ("ss_pssbl: ")
            for nn in ss_pssbl:
                print (nn.disp_Pair())
            #|endfor
            
            print ("------")
            if stop_at_end:
                print ("requested to stop at the end of test_entwined()")
                sys.exit(0)
            #
            
        #
        
        return ss_pssbl
    #
    
    
    def compress_to_BPlist(self):
        """
        Compress_to_BPlist: organizes BPlist, and PKlist so that it fits a
        kind of standard format of secondary structure and PK
        notation. It also generates two new important pieces of
        information: BProots and PKroots. These latter lists can be
        used to identify the base of each secondary structure domain
        and the arrangement of the pseudoknots that hangs over them.
    
    
    Example 1:
                  0         10        20        30        40 
                  |    .    |    .    |    .    |    .    |  
       ss_seq  = ".A..B..C..D...........a..b..c..d.."
       
       result:
       new BPlist:
       (   1,   22)[p]
       new BProots:
       (   1,   22)[p]
       new PKlist:
       (   4,   25)[p]
       (   7,   28)[p]
       (  10,   31)[p]
       new PKroots:
       (   4,   25)[p]
       (   7,   28)[p]
       (  10,   31)[p]
    
        [Note that for parallel stems, this operation does not really
        do all that much to reduce the size of the set of roots of
        parallel stems. It seems like this could be handled, as it is
        just difference of 1 between (i1,j1) and (i2,j2)
        
        i.e., (i2 - i1) = 1 and (j2 - j1) = 1                           
    
        However, presently, the output will always include except the
        root stem.]
    
    Example 2:
                  0         10        20        30        40 
                  |    .    |    .    |    .    |    .    |  
       ss_seq  = ".AAA...BBB...CCC....aaa...bbb....ccc...."
       
       result:
       new BPlist:
       (   1,   22)[a]
       (   2,   21)[a]
       (   3,   20)[a]
       new BProots:
       (   1,   22)[p]
       new PKlist:
       (   7,   28)[a]
       (   8,   27)[a]
       (   9,   26)[a]
       (  13,   35)[a]
       (  14,   34)[a]
       (  15,   33)[a]
       new PKroots:
       (   7,   28)[a]
       (  13,   35)[a]
    
    Example 3:
                  0         10        20        30        40 
                  |    .    |    .    |    .    |    .    |  
       ss_seq  = ".AAA...BBB...aaa...CCC....bbb....ccc...." # same
       ss_seq  = ".(((...BBB...)))...(((....bbb....)))...." # same
   
       result:
       new BPlist:
       (   1,   15)[a]
       (   2,   14)[a]
       (   3,   13)[a]
       (  19,   35)[a]
       (  20,   34)[a]
       (  21,   33)[a]
       new BProots:
       (   1,   15)[a]
       (  19,   35)[a]
       new PKlist:
       (   7,   28)[a]
       (   8,   27)[a]
       (   9,   26)[a]
       new PKroots:
       (   7,   28)[a]

        """
        debug_cbpl  = False # True # 
        stop_at_end = False # True # 
        if debug_cbpl:
            print ("enter compress_to_BPlist:")
        #  
        
        # First, we need to build direct copies of self.PKlist and
        # self.BProots and sort the lists with respect to the first
        # element.
        ssroots = copyList(self.BProots)
        ssroots = self.sortvsList(ssroots, 'i')
        #ssroots = sorted(ssroots, key=self.getKey_i)
        
        pklist = copyList(self.PKlist)
        # pklist is already sorted. 
        pkroots = self.findroots(pklist, debug_cbpl) 
        pkroots = self.sortvsList(pkroots, 'i')
        #self.PKroots = sorted(pkroots, key=self.getKey_i)
        
        # print out the current list of ss-domains and PK regions
        if debug_cbpl:
            print ("Before: vvvvvvvvvvvvvvvvvv")
            print ("base pair lists: ")
            print ("BPlist:")
            for sspr_k in self.BPlist:
                print (sspr_k.disp_Pair())
            #|endfor
            
            print ("root lists: ")
            print ("ssroots:")
            for pkpr_k in self.BProots:
                print (pkpr_k.disp_Pair())
            #|endfor
            
            print ("pseudoknot lists")
            print ("PKlist:")
            for pkpr_k in self.PKlist:
                print (pkpr_k.disp_Pair())
            #|endfor
            
            print ("pkroots:")
            for pkpr_k in self.PKroots:
                print (pkpr_k.disp_Pair())
            #|endfor
            
            print ("^^^^^^^^^^^^^^^^^^ :Before")
            #print ("stop at 0 in compress_to_BPlist"); sys.exit(0)
        #
        
        """
        Second, we need to isolate the regions where the various PK root
        entries are free of any secondary structure root
        entanglement. The resulting PK root entries may still be an
        entangled list of PKs, but after this operation, the remaining
        root PKs are independent of any secondary structure roots. In
        the first operation, we prune out all PK roots that are
        entangled with the root secondary structure.
        
        In this first step, we scan the list of secondary structure
        against the PK list and look for regions where there are
        clearly additional possible roots that may also be
        equivalently represented as secondary structure.
        
        """
        
        ss_pssbl = self.test_entwined(pkroots, ssroots, debug_cbpl)
        
        """
           ".((....AA...))..BB....aa....CC...bb....cc.."
        
        In this routine, ssroots and pkroots leave us with the
        following paired down structure
        
           ".((....AA...))..BB....aa....CC...bb....cc.."
           
             |                  |                   |
             V                  V                   V
           
           ".(.....A.....)..B......a....C.....b.....c.."
           
             |                  |                   |
             V                  V                   V
           
           "................B...........C.....b.....c.."
           
        In this next step, Bb and Cc are examined and, in this
        example, this algorithm will dig out Bb.
           
             |                  |                   |
             V                  V                   V
           
           "................B.................b........"
           
           ss_pssbl:
           (  16,  34)[a]
           (  28,  40)[a]
           revised ss_psbbl:
           (  16,  34)[a]
        """
       
        kp = 0
        while kp < (len(ss_pssbl)-1):
            i_p = ss_pssbl[kp].i; j_p = ss_pssbl[kp].j
            # remember, this is an ordered list!!!!
            kn = 0
            while kn < len(ss_pssbl):
                # compare the secondary structure roots with the 
                i_n = ss_pssbl[kn].i; j_n = ss_pssbl[kn].j
                """
                ...[.....(.......].......)..     ...(.....[.......).......]..
                   i_p   i_n     j_p     j_n        i_n   i_p     j_n     j_p
                """
                if ((i_p < i_n and i_n < j_p and j_p < j_n) or
                    (i_n < i_p and i_p < j_n and j_n < j_p)):
                    if debug_cbpl:
                        if   (i_p < i_n and i_n < j_p and j_p < j_n):
                            print ("i_p(%2d) < i_n(%2d) < j_p(%2d) < j_n(%2d)" \
                                % (i_p, i_n, j_p, j_n))
                        elif (i_n < i_p and i_p < j_n and j_n < j_p):
                            print ("i_n(%2d) < i_p(%2d) < j_n(%2d) < j_p(%2d)" \
                                % (i_n, i_p, j_n, j_p))
                        else:
                            print ("problems! cbpl")
                            print ("ij_n(%2d,%2d), ij_p(%2d,%2d)" % (i_n, j_n, i_p, j_p))
                        #
                        
                    #
                    
                    del ss_pssbl[kn] # pk meshed with ssroots
                else:
                    kn += 1
                #
                
            #|endwhile
            
            kp += 1
        #|endwhile
        
        if debug_cbpl:
            print ("revised ss_pssbl: ")
            for nn in ss_pssbl:
                print (nn.disp_Pair())
            #|endfor
            
            #print ("stop at 1 in compress_to_BPlist"); sys.exit(0)
        #
        
        """
        What have we done? Suppose that the following PK roots were
        extracted in the first step
        
           A.....B.......a...C.....b.......c
        
        This is the first result we would find in ss_pssbl.
        
        Because Bb is entangled between both Aa and Cc, this is 
        restructured to the following
        
           A.....B.......a...C.....b.......c -> (.....B.......)...(.....b.......)
        """
        pk2ssroots = copyList(ss_pssbl)
        pk2ssroots = self.sortvsList(pk2ssroots, 'i')
        #pk2ssroots = sorted(pk2ssroots, key=self.getKey_i)
        
        if debug_cbpl:
            print ("PKroots:")
            for vv in pkroots:
                print (vv.disp_Pair())
            #|endfor
            
            print ("pk2ssroots:")
            for vv in pk2ssroots:
                print (vv.disp_Pair())
            #|endfor
            
        #
        
        kv = 0
        for pk2ss in pk2ssroots:
            it = pk2ss.i; jt = pk2ss.j
            while kv < len(pkroots):
                ipkr = pkroots[kv].i; jpkr = pkroots[kv].j
                if debug_cbpl:
                    print ("ijt(%2d,%2d) vs ijpkr(%2d,%2d)" % (it, jt, ipkr, jpkr))
                #
                if it == ipkr and jt == jpkr:
                    if debug_cbpl:
                        print ("delete: ", pkroots[kv])
                    #
                    del pkroots[kv]
                else:
                    kv += 1
                #
                
            #|endwhile
            
        #|endfor
        
        del self.PKroots
        self.PKroots  = copyList(pkroots)
        
        if debug_cbpl:
            print ("------")
            print ("pkroots: ")
            for nn in self.PKroots:
                print (nn.disp_Pair())
            #|endfor
            
            print ("revised ss_pssbl: ")
            for nn in ss_pssbl:
                print (nn.disp_Pair())
            #|endfor
            
            print ("pk2ssroots: (should be the same as ss_pssbl)")
            for nn in pk2ssroots:
                print (nn.disp_Pair())
            #|endfor
            
            print ("------")
            
            #print ("stop at 2 in compress_to_BPlist"); sys.exit(0)
        #
        
        # now we have to rebuild the pk and ss structures from the
        # root information we have acquired and transfer this information. 
        
        
        BPgroup = copyList(self.BPlist)
        # We don't bother to sort this list because we may end up
        # adding to it from pklist and (moreover) the order of the
        # BPgroup list doesn't matter in the operations being done
        # here.
        for k in range(0,len(pk2ssroots)):
            i_pkr = pk2ssroots[k].i; j_pkr = pk2ssroots[k].j; v_pkr = pk2ssroots[k].v 
            if debug_cbpl:
                print ("k  = %2d, ij_pkr(%2d,%2d)[%s]" % (k, i_pkr, j_pkr, v_pkr))
            #
            
            kl = 0
            while kl < len(pklist):
                i_pkt = pklist[kl].i; j_pkt = pklist[kl].j; v_pkt = pklist[kl].v
                if debug_cbpl:
                    print ("kl = %2d, ij_pkt(%2d,%2d)[%s]" % (kl, i_pkt, j_pkt, v_pkt))
                    #sys.exit(0)
                #
                
                if v_pkr == 'a' and v_pkt == 'a':
                    if debug_cbpl:
                        print ("v_pkr = %s, v_pkt = %s" % (v_pkr, v_pkt))
                        if i_pkr <= i_pkt and j_pkt <= j_pkr:
                            print ("add to BPgroup: ", pklist[kl])
                        #
                        
                    #
                    
                    if i_pkr <= i_pkt and j_pkt <= j_pkr:
                        BPgroup += [pklist[kl]]
                        del pklist[kl]
                    else:
                        kl += 1
                    #
                    
                elif i_pkr == i_pkt and j_pkt == j_pkr:
                    if debug_cbpl:
                        print ("ij_pkr(%2d,%2d) <=> ij_pkt(%2d,%2d)" \
                            % (i_pkr, j_pkr, i_pkt, j_pkt))
                    #
                    
                    BPgroup += [pklist[kl]]
                    del pklist[kl]
                else:
                  kl += 1
                #
                
            #|endwhile
            
        #|endfor
        
        # Now we sort this list (after possibly updating it)
        
        BPgroup = self.sortvsList(BPgroup, 'i')
        #BPgroup.sort(key=self.getKey_i)
        
        # update the PKlist and rest the initial self.BPlist and
        # self.PKlist, finally reconstruct these two lists.
        PKgroup = copyList(pklist)
        PKgroup = self.sortvsList(PKgroup, 'i')
        #PKgroup.sort(key=self.getKey_i)
        del self.BPlist
        del self.PKlist
        
        self.BPlist = copyList(BPgroup)
        self.PKlist = copyList(PKgroup)
        self.BProots = self.sortvsList(self.findroots(self.BPlist), 'i')
        #self.BProots = sorted(self.findroots(self.BPlist, False), key=self.getKey_i)
        
        
        """
        This first sequence does exactly what I want it to do....
        
        ss_seq  = ".ABCD....EFG...efg..((..[[.)).]]..HIJ...hij....abcd.."
        > python Threads.py
              ...
        updated base pair lists: 
        new BPlist:
        (    1,   47)[p]
        (    9,   15)[p]
        (   20,   28)[a]
        (   21,   27)[a]
        (   34,   40)[p]
        new BProots:
        (    1,   47)[p]
        (    9,   15)[p]
        (   20,   28)[a]
        (   34,   40)[p]
        updated pseudoknot lists: 
        new PKlist:
        (    2,   48)[p]
        (    3,   49)[p]
        (    4,   50)[p]
        (   10,   16)[p]
        (   11,   17)[p]
        (   24,   31)[a]
        (   25,   30)[a]
        (   35,   41)[p]
        (   36,   42)[p]
        new PKroots:
        (    2,   48)[p]
        (    3,   49)[p]
        (    4,   50)[p]
        (    9,   15)[p]
        (   10,   16)[p]
        (   11,   17)[p]
        (   24,   31)[a]
        (   34,   40)[p]
        (   35,   41)[p]
        (   36,   42)[p]
        ^^^^^^^^^^^^^^^^^^ :After
        
        The following sequence is a bit more nasty and it required
        some fixes to test_entwined to get it to work a little
        better. The current output follows.
        
        ss_seq  = ".((((((.A.AAAAA.<BC..))))).((((.>bc....a..aaaaa....))))..)."
        
        After: vvvvvvvvvvvvvvvvvvv
        updated base pair lists: 
        new BPlist:
        (    1,   57)[a]
        (    2,   25)[a]
        (    3,   24)[a]
        (    4,   23)[a]
        (    5,   22)[a]
        (    6,   21)[a]
        (   27,   54)[a]
        (   28,   53)[a]
        (   29,   52)[a]
        (   30,   51)[a]
        new BProots:
        (    1,   57)[a]
        updated pseudoknot lists: 
        new PKlist:
        (    8,   46)[a]
        (   10,   45)[a]
        (   11,   44)[a]
        (   12,   43)[a]
        (   13,   42)[a]
        (   14,   39)[a]
        (   16,   32)[p]
        (   17,   33)[p]
        (   18,   34)[p]
        new PKroots:
        (    8,   46)[a]
        (   16,   32)[p]
        (   17,   33)[p]
        (   18,   34)[p]
        ^^^^^^^^^^^^^^^^^^ :After
        
        """
        if debug_cbpl:
            print ("After: vvvvvvvvvvvvvvvvvvv")
            print ("updated base pair lists: ")
            print ("new BPlist:")
            for bpg in self.BPlist:
                print (bpg.disp_Pair())
            #|endfor
            
            print ("new BProots:")
            for bpr in self.BProots:
                print (bpr.disp_Pair())
            #|endfor
            
            print ("updated pseudoknot lists: ")
            print ("new PKlist:")
            for pkg in self.PKlist:
                print (pkg.disp_Pair())
            #|endfor
            
            print ("new PKroots:")
            for pkr in self.PKroots:
                print (pkr.disp_Pair())
            #|endfor
            
            print ("^^^^^^^^^^^^^^^^^^ :After")
            print ("Exit compress_to_BPlist:")
            #print ("stop at 3 (exit from) compress_to_BPlist"); sys.exit(0)
            if stop_at_end:
                sys.exit(0)
            #
            
        #
        
        return 0
    #
    
    
    
    def expand_island(self, CTCFisland):
        MPlist = [CTCFisland]
        n_k = len(CTCFisland.contacts)
        # print (n_k)
        if n_k == 1:
            for l in CTCFisland.contacts:
                b1 = Pair()
                b1.put_ssPair(CTCFisland.i, l, "ctcf", '-')
                b2 = Pair()
                b2.put_ssPair(l, CTCFisland.j, "ctcf", '-')
                MPlist += [b1,b2]
            #|endfor
            
        elif n_k > 1:
            b = Pair()
            i = CTCFisland.i
            j = CTCFisland.contacts[0]
            b.put_ssPair(i, j, "ctcf", '-')
            MPlist += [b]
            for l in range(1,len(CTCFisland.contacts)):
                for k in range(0,l):
                    if l > k + 1:
                        continue
                    #
                    
                    b = Pair()
                    i = CTCFisland.contacts[k]
                    j = CTCFisland.contacts[l]
                    # print (i, j)
                    b.put_ssPair(i, j, "ctcf", '-')
                    MPlist += [b]
                #|endfor
                
            #|endfor
            
            b = Pair()
            n = len(CTCFisland.contacts)
            i = CTCFisland.contacts[n-1]
            j = CTCFisland.j
            b.put_ssPair(i, j, "ctcf", '-')
            MPlist += [b]
        #
        
        return MPlist
    #
    
    
    """
    Since the pair elements are derived from the order in which they
    appear on the stack, after you encounter a "mate", you search back
    through the stack until you find one unassigned element. So
    essentially, you first put a bunch of elements (e.g., 'X') on the
    stack (of course, assuming you have a string of them), now, you
    have just encountered the partner of this example ('x') and you
    look at its position on the stack using p (i.e., the
    pointer[rpr2num['x']] position last incremented). From here of
    course, you know that the position on the stack must decrease by
    exactly the number of 'X' items you put on the stack (more to the
    point, it BETTER match!!!)
    
    For a more concrete example, suppose you have the structure
    
       ".X.XXXXX......x..xxxxx.."
    
    where element X has the index 27
    
    we reach the top of the stack with counter[27] = pointer[27] = stack[27] = 6
    
    self.Xlist[27] should now look like this:
    
      ( 1, -1)
      ( 3, -1)
      ( 4, -1)
      ( 5, -1)
      ( 6, -1)
      ( 7, -1)
    
    Now we go through this list looking for the -1. In this case, it
    is at the top of the list, so it is the first element. So we
    assign pp = 5 (because the list starts from 0 so the last element
    is 5).
    
     ( 1, -1)
     ( 3, -1)
     ( 4, -1)
     ( 5, -1)
     ( 6, -1)
     ( 7, 14)
    
    As we find the remaining elements, pointer[27] decreases. So after
    6 iterations, we obtain the correct assignment of elements on the
    stack.
    
     ( 1, 21)
     ( 3, 20)
     ( 4, 19)
     ( 5, 18)
     ( 6, 17)
     ( 7, 14)
     ".X.XXXXX......x..xxxxx.."
    
    Personally, I think the looping should not occur, because that
    means that the stack is out of order, but maybe this search can be
    employed to detect problems in either building the stack or
    closing it. So I have left this somewhat incongruous feature
    within the tool for the moment (since it is still under
    development). Better safe than sorry as they say.

    """
    
    def find_next_Xpoint_n(self, plist, p, nm = "BP", debug = False):
        if debug:
            print ("find_next_Xpoint_n, p = %d, name(%s)" % (p, nm))
            print ("plist:")
            for pk in plist:
                print (pk.disp_Pair())
            #|endfor
            
        #
        
        pp = p
        j = plist[p-1].j
        if debug:
            print ("pp(%3d), j(%3d)" % (pp, j))
        #
        
        while j > -1 and pp >= 0:
            pp -= 1
            j = plist[pp-1].j
            if debug:
                print ("pp(%3d), j(%3d)" % (pp, j))
            #
            
        #|endfor
        
        if debug:
            print ("result: pp = %d" % pp)
        #
        
        return pp
    #
    
    
    def check_Xlist_n(self, plist):
        pass_X = True
        for b in plist:
            if b.j < 0:
                pass_X = False
            #
            
        #|endfor
        
        return pass_X
    #
    
    
    def print_vstr(self):
        print (self.vstr)
    #
    
    
    def print_vseq(self):
        
        if not self.vseq == "":
            print (self.vseq)
        #
        
        #print (self.vseq)
        
    #
    
    
    def print_Xlist_n(self, plist):
        for b in plist:
            print (b.disp_Pair())
        #|endfor
        
    #
    
    
    # only used for secondary structure
    def parse_SecondaryStructure(self, r_ss, rseq = ""):
        self.reset_VstructLists()
        self.set_Vstruct(r_ss)
        if not rseq == "":
            self.set_Vseq(rseq)
        #
        
        self.scan_vs(0)
        
        self.print_vseq()
        self.print_vstr()
        self.print_Xlist_n(self.BPlist)
        return 0
    #
    
    
    def disp_allLists(self, sttmnt = "Results:", show_seqs = True):
        print (sttmnt)
        if show_seqs:
            print ("structural sequence")
            self.print_vseq()
            self.print_vstr()
        #
        
        print ("secondary structure")
        self.print_Xlist_n(self.BPlist)
        print ("secondary structure roots")
        self.print_Xlist_n(self.BProots)
        
        if len(self.PKlist) > 0:
            print ("pseudoknot linkages: ")
            self.print_Xlist_n(self.PKlist)
            print ("pseudoknot linkage roots")
            self.print_Xlist_n(self.PKroots)
        #
        
        print ("dXPKinfo:")
        for v in list(self.dXPKinfo.keys()):
            print (self.dXPKinfo[v])
        #
        
        
        print ("dCPKinfo:")
        for v in list(self.dCPKinfo.keys()):
            print (self.dCPKinfo[v])
        #
        
        
        if len(self.MPlist) > 0:
            print ("CTCF connects")
            self.print_Xlist_n(self.MPlist)
        #
        
        if len(self.MPlist) > 0:
            islands = []
            for cl in self.MPlist:
                islands += [self.expand_island(cl)]
            #|endfor
            
            print ("CTCF breakdown:")
            kk = 1
            for island_k in islands:
                print ("island(%d):" % kk)
                for ii in island_k:
                    print (ii.disp_Pair())
                #|endfor
                
                kk += 1
                
            #|endfor
            
        #
        
    #
    
    
    # can do RNA or chromatin with all types of structures
    def parse_fullDotBracketStructure(self, r_ss, rseq = "", print_result=False):
        pass_check = True
        if  not type(r_ss) == str:
            print ("Vstruct: input structure not a string:")
            print ("input type:    ", type(r_ss))
            print ("input content: ", r_ss)
            pass_check = False
        #
        
        if not type(rseq) == str:
            print ("Vstruct: input sequence not a string:")
            print ("input type:    ", type(rseq))
            print ("input content: ", rseq)
            pass_check = False
        #
        
        if not type(print_result) == bool:
            print ("Vstruct: print results can only be True or False:")
            print ("input type:    ", type(print_result))
            print ("input content: ", print_result)
            pass_check = False
        #
        
        if not pass_check:
            sys.exit(1)
        #
        
        
        
        debug = False
        self.reset_VstructLists()
        self.set_Vstruct(r_ss)
        if not rseq == "":
            self.set_Vseq(rseq)
        #
        
        self.scan_allTypes(0, 0)
        if print_result:
            if debug:
                self.disp_allLists("Results at finishing Vienna")
                print ("Exit parse_fullDotBracketStructure")
            else:
                self.disp_allLists()
            #
            
        #
        
        return 0
    #
#



# reads in Vienna package data created by the SimRNA_trafl2pdbs
# processing function
class ViennaData(object):
    
    def __init__(self):
        # secondary structure created from SimRNA_trafl2pdbs processing
        self.ss_flnm = ''
        self.ss_data = []
        self.n_structs = -1  
    #
    
    ################################
    ########   Functions    ########  
    ################################
    
    def reset_ViennaData(self):
        self.ss_flnm = ''
        self.ss_data = []
        self.n_structs = -1  
        return 0
    #
    
    # reading in trajectory data (*.trafl) from SimRNA run
    def get_ViennaData(self, fl):
        self.ss_flnm = fl
        try: 
            input = open(self.ss_flnm, 'r')
        except:
            print ('ERROR: %s does not exist!' % self.ss_flnm)
            sys.exit(1)
        #
        
        while True:
            string = input.readline().strip()        # ((((.....)))) data
            if not string: break
            self.ss_data += [string]
        #|endwhile
        
        input.close()
        self.n_structs   = len(self.ss_data)
        if self.n_structs == 0:
            mssg = '\nERROR: no data read\n'
            raise MyException(mssg)
        #
        
        return 0
    #
    
    # reading in trajectory data (*.trafl) from SimRNA run
    def read_SimRNA_ss(self, data):
        
        # print ("read_SimRNA_ss(): \n%s" % data)
        
        self.ss_data = data
        self.n_structs   = len(self.ss_data)
        if self.n_structs == 0:
            mssg = '\nERROR: no data read\n'
            raise MyException(mssg)
        #
        
        # print (self.n_structs, self.ss_data)
        return 0
    #
    
#

def test0():
    ms = MolSystem()
    ms.set_system("RNA")
    ms.set_JobType("evaluation")
    ms.set_program("Vienna::test0")
    ms.set_mstr("((((....))))")
    ms.set_mseq("gggguuuucccc")
    ms.set_ParamType("Turner") # no parameter file will be specified
    
    vs = Vstruct(ms)
    fpntr = [(0, 23), (3, 5), (7, 17), (9, 11), (10, 14), (13, 15), (19, 21)]
    print (vs.subtract_crossovers(fpntr))
    
    
def test1(mstr):
    ss_seq = "(((((....)))))"
    if len(mstr) > 1:
        ss_seq = mstr
    #
    
    rnaseq = genRNASeq(ss_seq)
    
    ms = MolSystem()
    ms.set_system("RNA")
    ms.set_JobType("evaluation")
    ms.set_program("Vienna::test1")
    ms.set_mstr(ss_seq)
    ms.set_mseq(rnaseq)
    ms.set_ParamType("Turner") # no parameter file will be specified
    # other possible arguments in to use set_ParamType(xxx)
    # "ViS"
    # "ParFile", "/home/dawson/python/RNA/ViSparams.par"
    # "gMatrix", "/home/dawson/python/RNA/test3s_1mm+1mmp1+2mm+3mm_3nt-w5x5.gMtrx"
    
    """@
    
    Admittedly, compared to the previous way where I just entered "vs
    = struct()", selected the structure (often without even bothering
    to specify the sequence), and that was it, this comes out as
    bureaucratic as an application for an HP part. However, maybe now,
    I finally understand why HP component and part applications were
    so complicated. Over the years, they too had expanded into a very
    diverse set of products and needed to be sure that the application
    was specified correctly. Either they have to have someone do that,
    or they need the customer to do so. So maybe this is a reality
    that I have not considered before. On the other hand, this way, we
    are sure that the data is properly set up for use and there is no
    ambiguity as the was with simply Vstruct() and no other
    information.  """
    
    vs = Vstruct(ms)
    vs.parse_SecondaryStructure(ss_seq, rnaseq)
    # output 
    # (((((....)))))
    # (    0,   13)
    # (    1,   12)
    # (    2,   11)
    # (    3,   10)
    # (    4,    9)
    # 0
#   

def test2(mstr):
    
    #ss_seq = "((((.....))))...((....))."
    #ss_seq = "{(((((.[[...)))))..]]...|..([)]([)]...}"
    #ss_seq =  "{(((((.A.AAAAA.....)))))......a..aaaaa...|..([)]([)]...}....{.....}"
    #ss_seq =  "{(((((.A.AAAAA.CB..)))))..bc..a..aaaaa...|..([)]([)]...}....{.....}"
    #ss_seq =  "{((((((.A.AAAAA......))))).BBBB........a..aaaaa....bbbb..).|..([)]([)].ABC.abc.}....{.....}"
    #ss_seq =  "{((((((.A.AAAAA......))))).((((........a..aaaaa....))))..).|..([)]([)].ABC.abc.}....{.....}"
    ss_seq  =  "{((((((.A.AAAAA.<BC..))))).((((.>bc....a..aaaaa....))))..).|..([)]([)]...}....{.A.B.C..a.b.c.....}"
    #ss_seq =  "{((((((.A.AAAAA.<BC..))))).((((.>bc.DE.a..aaaaa..de))))..).|..([)]([)]...}....{.A.B.C..a.b.c.....}"
    #ss_seq  =  "((((...AAAA...))))...BBBB....aaaa...((((...bbbb..))))..."
    if len(mstr) > 1:
        ss_seq = mstr
    #
    
    chrseq = genChrSeq(ss_seq)
    
    ms = MolSystem()
    ms.set_system("Chromatin")
    ms.set_JobType("evaluation")
    ms.set_program("Vienna::test2")
    ms.set_mstr(ss_seq)
    ms.set_mseq(chrseq)
    ms.set_ParamType("genheat") # no parameter file will be specified
    # other possible arguments in to use set_ParamType(xxx)
    # "ViS"
    # "ParFile", "/home/dawson/python/RNA/ViSparams.par"
    # "gMatrix", "/home/dawson/python/RNA/test3s_1mm+1mmp1+2mm+3mm_3nt-w5x5.gMtrx"
    vs = Vstruct(ms)
    vs.set_Vseq(chrseq)
    vs.parse_fullDotBracketStructure(ss_seq, chrseq, True)
    # output 
    # |(((((.[[...)))))..]]...|..([)]([)]...|
    # secondary structure
    # (    1,   16)
    # (    2,   15)
    # (    3,   14)
    # (    4,   13)
    # (    5,   12)
    # (   27,   29)
    # (   31,   33)
    # pk connects
    # (    7,   20)
    # (    8,   19)
    # (   28,   30)
    # (   32,   34)
    # CTCF connects
    # (    0,   24)
    # (    0,   38)
    # (   24,   38)
#

def test3():
    # answer should be [(1, 10), (12, 20)]
    plist = [(1, 10), (2, 8), (4, 7), (12, 20), (14, 18), (30, 40), (32, 38)]
    plist += deepcopy(plist)
    
    print (plist)
    plist = prune_pairlist(plist)
    print (plist)
    
    
    pref  = [(1, 10), (2, 8), (4, 7), (12, 20), (14, 18), (30, 40), (32, 38)]
    BProots = []
    
    for prefk in pref:
        bp = Pair()
        # Vienna is specifically designed for secondary structure
        # of single strand (ss) pairs. Therefore, we use put_ssPair
        bp.put_ssPair(prefk[0], prefk[1], 'bp', 'a')
        BProots += [bp]
    #
    
    ir = 1; jr = 10; il = 5; jl = 25
    fpntr = get_CPK_jbranches(ir, jr, il, jl, BProots)
    plist = compress_to_jbranches(fpntr)
    print (plist)
    # answer should be [(12, 20)]
    
    ir = 12; jr = 20; il = 0; jl = 15
    fpntr = get_CPK_jbranches(ir, jr, il, jl, BProots)
    plist = compress_to_jbranches(fpntr)
    print (plist)
    # answer should be [(1, 10)]
    
    
    ir1 = 1; jr1 = 10; ir2 = 30; jr2 = 40
    fpntr = get_XPK_jbranches(ir1, jr1, ir2, jr2, BProots)
    plist = compress_to_jbranches(fpntr)
    print (plist)
    # answer should be [(12, 20)]
#



def testOptions():
    print ("test0:    test various functions associated with class Vstruct")
    print ("test1:    test only 1D secondary structure sequences")
    print ("test2:    test all types of 1D structure sequences")
    print ("test3:    test various simple functions")
#


def main(cl):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', default="test3", type=str,
                        dest='testType',
                        help='select a particular test.')
    
    parser.add_argument('-mstr', default="", type=str,
                        dest='mstr',
                        help='select a particular 1D structure.')
    
    parser.add_argument('-h.showTests', action='store_true', default=False,
                        dest='show_tests',
                        help='shows a list to the current tests available.')
    
    args = parser.parse_args()
    
    test = args.testType
    mstr = args.mstr
    
    if args.show_tests:
        testOptions()
        sys.exit(0)
    #

    v2t_test = -1
    if test in dtestNames:
        v2t_test = dtestNames[test]
    else:
        print ("ERROR: unrecognized test (%s)" % test)
        usage()
        sys.exit(1)
    #

    
    
    if   v2t_test == 0: # test simple functions in Vstruct
        test0()
        
    elif v2t_test == 1: # test only 1D secondary structure sequences
        test1(mstr)
        
    elif v2t_test == 2: # test all types of 1D structure sequences
        test2(mstr)
        
    elif v2t_test == 3: # test simple functions in Vienna
        test3()
        
    else:
        print ("should not be here! ")
        print ("command line", cl)
        sys.exit(1)
    #
#

# Main
if __name__ == '__main__':
    main(sys.argv)
#
