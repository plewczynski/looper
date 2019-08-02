#!/usr/bin/env python

"""@@@

Main Module:   LThreadBuilder.py 

Objects:       LThreadBuilder
               NodeVisitor


Author:        Wayne Dawson
creation date: 170725 (originally part of Threads.py & RThreads.py and later TreeNode) 
last update:   190725
version:       0.1

Purpose:

This module is currently used only by "sarabande.py" to analyze RNA
structure. Because there are some differences in the analysis routines
of the chromatin, RNA and protein based program, I have had to make
these separate.

Comments:

190705:

   This started out as Threads, moved to RThreads (for RNA
   problems). That became a huge unmanageable kluge, so I started
   process of breaking things down, first by moving it to TreeNode and
   finally reducing it into the submodule LThreadBuilder. Mostly, it
   is easier to manage TreeNode without fighting with
   LThreadBuilder. In the future, I may eventually remove LThread from
   the program, but presently, it is still too integral a part of this
   package.


   The LThread data is rather independent of the data sets. It is very
   similar to Pair data and I am thinking that eventually, I may merge
   these two data sets, as it is an inadequate representation to be
   highly useful and I don't much like the proliferation of very
   similar looking data sets.


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

"""for Vienna Pair object representation"""
from Pair import find_this_ij
from Pair import is_ap_stem
from Pair import is_pp_stem
from Pair import Pair
from Pair import SortPair
from Pair import vsPair2list

from Vienna import Vstruct

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

"""for Motif data representation"""
from Motif      import Motif # a core object of these building programs
from Motif      import Link
from Motif      import LGroup
from Motif      import Map # main map for all the FE and structures

"""contains free energy parameters"""
#from RFreeEnergy import FreeEnergy

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


# ################################################################
# ############  Local-global functions and constants  ############
# ################################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# 1. control parameters tests etc in the program

# Test   function
#  0     a _very basic_ test of this module
#  1     insert into a list
#  2     parser for simple structures

TEST = 2


# 2. special local global function

PROGRAM      = "LThreadBuilder.py"  # name of the program

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
For an explanation why!?! we need object, please see
https://stackoverflow.com/questions/4015417/python-class-inherits-object

The short answer is that it is not necessary for Python 3, but it is
necessary for backward compatibility in Python 2.x.

it seems according that source, it seems the best strategy is to use
this style because it gives you access to mnay different python
features. In particular, descriptors. See

https://docs.python.org/3/howto/descriptor.html

The write up in the above refence adds ... One of the downsides of
new-style classes is that the class itself is more memory
demanding. Unless you're creating many class objects, though, I
doubt this would be an issue and it's a negative sinking in a sea of
positives.

Note, instead of the system "NodeVisitor", we could have written
visit(node) in the following way. Here is one way to do it...


def visit(node):
    node_type = type(node).__name__
    if node_type == 'BinOp':
        return self.visit_BinOp(node)
    elif node_type == 'Num':
        return self.visit_Num(node)
    elif ...
    # ...


Another alternative way of doing this is as follows:


def visit(node):
    if isinstance(node, BinOp):
        return self.visit_BinOp(node)
    elif isinstance(node, Num):
        return self.visit_Num(node)
    elif ...


The "visit(node)" is explained in some detail in part 7;
https://ruslanspivak.com/lsbasi-part7/

The "Visit" concept is actually quite a sophisticated approach. There
is considerable discussion and examples at the following link.
https://en.wikipedia.org/wiki/Visitor_pattern

In the wiki reference, they use yet a different way to achieve this
same approach with module "abc".  Apparently, Visit is basically a way
to permit flexible code so you can use various objects and reference
them from a rather simple procedure. ... "my brain is full". 

In this program, NodeVisitor is used by LThreadBuilder and TreeNode2Motif

"""

class NodeVisitor(object):
    def visit(self, node):  # HERE is where visit goes first
        """
        This is really cool! I could use this with class Motif() to permit
        objects to be built for stem, PK, etc without fighting the
        type issues. The fact that Parse in this example actually
        remembers the tree after builting it is an encouraging point.

        """
        debug_visit = False
        method_name = 'visit_' + type(node.name).__name__
        if debug_visit:
            print method_name
        #
        
        # method_name takes the form of the following strings:
        # visit_BinOp, visit_UnaryOp, visit_Num
        
        visitor = getattr(self, method_name, self.generic_visit)
        """
        If method_name is valid, this generates a number
        
           getattr(...)
           getattr(object, name[, default]) -> value
        
        Get a named atrribute from an object; getattr(x, 'y') is
        equivalent to x.y. When a default argument is given, it is
        returned when the attribute doesn't exist; without it, an
        exception is raise in that case.
        
        In this case
        
        getattr(self, method_name ... ) -> self.<method_name>
        
        """
        
        return visitor(node)
    #
    
    def generic_visit(self, node):
        raise Exception('No visit_{} method'.format(type(node).__name__))
    #
#


"""
Originally, I though that I would use LThread as an intermediate
for moving between the Vienna format to the Motif format. However, I
found LThread to be a very inadequete appoach to this problem and
eventually discovered that building trees was far more effective.

Nevertheless, although LThread is not all that much better than the
Vienna notation when it comes to easily parsing complex pseudoknot
structures, it carries information on free energy wich is useful for
building a structure computational engine. Likewise, since LThread is
used in the structure prediction part to help paste different
suboptimal structures together, it is still a convenient intermediate
tool. Finally, one of the current problems with the structure
prediction part is that rather complicated twisted up structures can
be generated during the folding process. Whereas some such structures
may still be hard to understand, the general tools developed here to
convert Vienna to TreeNode representation can be extended further to
convert LThread to TreeNode representations -- which means that the
calculated results can be put in some standard form.

Therefore, whereas LThread (though systematic) is not a really good
notation system for describing the structure of RNA, chromatin or
other biomolecules -- it still takes considerable effort (which is
highly error prone) to read and understand it -- LThread is useful for
building the engine and converting predicted structures into some
standard representation. In short, LThread is useful for computing the
free energy of the structure with the intermediary of applying
Vienna2TreeNode.

"""

class LThreadBuilder(NodeVisitor):
    def __init__(self, v2t):
        inputObject = type(v2t).__name__
        if not inputObject == "Vienna2TreeNode":
            print "Error: LThreadBuilder requires objects of class Vienna2TreeNode."
            print "       object entered: %s" % inputObject
            sys.exit(1)
        #
        self.ssStems  = v2t.ssStems # helps with searching out branching
        self.MPlist   = v2t.MPlist
        self.N        = v2t.N
        self.lt       = LThread(self.N)
        self.fe       = v2t.fe # basic FE/entropy function for chromatin

        self.write_MultiPairs()
        
        if len(self.ssStems) == 0:
            print "WARNING: structure list is empty"
        #
        self.initalized = True
        """
        LNode
        def __init__(self, ij_ndx, dGij_B, ctp, btp):
           self.ij_ndx = ij_ndx # (i, j) This is very important!
           # probably should have ij_ndx.i, ij_ndx.j calling it a pair.
           
           self.dGij_B = dGij_B # the FE of the specific bond at ij
           self.ctp    = ctp    # connection type
           self.btp    = btp    # bond type
        #
        """
    #
    
    
    
    def write_MultiPairs(self):
        if len(self.MPlist) > 0:
            for island in self.MPlist:
                 
                i_W = island.it; j_W = island.jt
                self.lt.thread += [LNode((i_W,j_W), 0.0, 'W', 'bgn')]
                dG = -100/float(len(island.arches))
                for k in range(0, len(island.arches)):
                    ctp = 'w'
                    if k == 0:
                        ctp = 'W'
                    #
                    
                    ww =  island.arches[k]
                    i_w = ww.it; j_w = ww.jt
                    self.lt.thread += [LNode((i_w,j_w), dG, ctp, 'wyspa')]
                    # print "ww.internal: ij_w(%2d,%2d)" % (i_w, j_w), ww.internal
                    if len(ww.internal) > 0:
                        for vv in ww.internal:
                            self.lt.thread += [LNode((vv.i,vv.j), 0.0, ww.btype, 'wyspa')]
                        #
                        
                    #
                    
                #endfor            
                self.lt.thread += [LNode((i_W,j_W), 0.0, 'W', 'end')]
            #endfor
        #
    #
    
    
    def visit_Stem(self, node):
        vstem = node.name
        if not vstem.vtype == "root":
            self.write_StemThread(vstem)
        #
        
        if not node.children == None:
            if len(node.children) > 0:
                for child in node.children:
                    self.visit(child)
                #endfor
                
            #
            
        #
        
    #
    
    def visit_PseudoKnot(self, node):
        vpk = node.name
        ipk = vpk.it; jpk = vpk.jt
        
        
        debug_visit_PseudoKnot = False # True # 
        if debug_visit_PseudoKnot:
            print "Enter TreeBuilder.visit_PseudoKnot ijpk(%2d,%2d)" % (ipk, jpk)
            print vpk
        #
        
        
        pktp = vpk.pktype
        self.lt.thread += [LNode((ipk,jpk), 0.0, pktp, 'bgn')]
        
        if len(vpk.rootstems) == 1:
            self.lt.thread += [LNode((ipk,jpk), 0.0, 'J', 's')]
            
        elif len(vpk.rootstems) > 1:
            self.lt.thread += [LNode((ipk,jpk), 0.0, 'P', 's')]
            
        else:
            print "ERROR(visit_PseudoKnot): (%d,%d)" % (ipk, jpk)
            print vpk.rootstems
            print vpk.linkages
            sys.exit(1)
        #
            
        for stem in vpk.rootstems:
            if debug_visit_PseudoKnot:
                print "rootstem: ", stem
            #
            
            self.write_StemThread(stem)
        #endfor
        
        self.lt.thread += [LNode((ipk,jpk), 0.0, 'l', 'bgn')]
        for stem in vpk.linkages:
            if debug_visit_PseudoKnot:
                print "linkages: ", stem
            #
            
            self.write_StemThread(stem, 'l')
        #endfor
        
        self.lt.thread += [LNode((ipk,jpk), 0.0, 'l', 'end')]
        
        self.lt.thread += [LNode((ipk,jpk), 0.0, pktp, 'end')]
        if not node.children == None:
            if len(node.children) > 0:
                for child in node.children:
                    self.visit(child)
                #endfor
                
            #
            
        #
        
    #
    
    
    def write_StemThread(self, vstem, stype = 's'):
        """vstem -> class Stem"""
        
        # stype = 's' (strandard) or 'l' (linkage)
        debug_write_StemThread = False # True # 
        i_t = vstem.it; j_t = vstem.jt
        i_h = vstem.ih; j_h = vstem.jh
        
        branches = self.get_branches(i_h, j_h)
        nbranches = len(branches)
        
        if debug_write_StemThread:
            print "Enter write_StemThread{(%2d,%2d):(%2d,%2d), stype(%s)" \
                % (i_t, j_t, i_h, j_h, stype)
            print "branches: ", branches
        #
        
        # presently, branches is not needed directly, but in the
        # future with the RNA engine, it is more important to have
        # this kind of information for determining the free energy.
        
        ctp = 'S'
        self.lt.thread += [LNode((i_t, j_t),  0.0, ctp, 'bgn')]
        dG = 0.0
        n = vstem.slen - 1
        for k in range(0, n): # bpk in vstem.stem:
            bpk = vstem.stem[k]
            i = bpk.i; j = bpk.j
            # xx!!!! 190707 !!!!xx revise potential ... also, the stem
            # needs some changes in definition
            dG = self.fe.calc_dG(i, j, 5.0, self.fe.T)
            btp = "%s%s" % (stype, bpk.v)
            if not (bpk.v == 'a' or bpk.v == 'p'):
                btp = bpk.name
            #
            
            nth = LNode((i,j), dG, ctp, btp)
            self.lt.thread += [nth]
            
            if debug_write_StemThread:
                print "(%2d,%2d)[ctp=%s][btp=%s]" % (i,j, ctp, btp)
            #
            
        #endfor
        
        bpk = vstem.stem[n]
        # print vstem.stem
        i = bpk.i; j = bpk.j; v = bpk.v
        dG = self.fe.calc_dG(i, j, 5.0, self.fe.T)
        btp = "%s%s" % (stype, bpk.v)
        if not (bpk.v == 'a' or bpk.v == 'p'):
            btp = bpk.name
        #
        
        # write the terminal part of the stem
        
        ctp_term = 'x'
        if nbranches == 1:
            # i-loop
            # print "iloop"
            ctp_term = 'I'
        elif nbranches > 1:
            # print "mloop"
            # m-loop
            ctp_term = 'M'
        else:
            # no children = leaf
            # print "hloop"
            ctp_term = 'B'
        #
        
        # xx!!!! 190707 !!!!xx Stem is self contained
        self.lt.thread += [LNode((i_h, j_h), dG, ctp_term, btp)]
        
        if debug_write_StemThread:
            print "(%2d,%2d)[ctp=%s][btp=%s]" % (i,j, ctp, btp)
        #
        
        # write the structure end note
        self.lt.thread += [LNode((i_t, j_t),  0.0, ctp, 'end')]
    #
    
    def get_branches(self, ib, jb):
        debug_count_branches = False # True # 
        
        branches = []
        for vstem in self.ssStems:
            gstemtype = type(vstem).__name__
            # print gstemtype
            if not gstemtype == "Stem":
                print "ERROR: ssStem contains objects other than type Stem"
                print "       %s: ijss = (%2d,%2d)" % (gstemtype, vstem.it, vstem.jt)
                sys.exit(1)
            #
            
            iss = vstem.it; jss = vstem.jt
            if ib < iss and jss < jb:
                branches += [(iss, jss)]
            #
            
        #endfor
        
        return branches
    #
    
    def disp_lt(self):
        for ltk in self.lt.thread:
            print ltk.disp_lnode()
        #
        
    #
    
#





def test0(cl):
    # this is not really a test that involves anything related to
    # TreeNode. I am not sure that it is possible to test this
    # independent of Vienna2TreeNode because many routines dependent
    # on Vienna2TreeNode as an input. The modules are highly
    # intertwined as Vienna2TreeNode depends on the classes in this
    # module to generate structures, but these classes also depended
    # on the existence of Vienna2TreeNode to function.
    
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

def test1(cl):
    aList = [123, 'xyz', 'zara', 'abc']
    aList.insert( 3, 2009)
    print "Final List : ", aList
#


def test2(cl):
    from Vienna2TreeNode import Vienna2TreeNode
    # The only way to get this to work is to load this module _within_
    # this specific method. When I do that, it will run without
    # causing circular definitions of TreeNode and Vienna2TreeNode. So
    # it accomplishes a certain amound of independence as though this
    # were a separate program nowhere located within TreeNode.py. 
    
    vs = Vstruct()
    vs = Vstruct()
    vs.set_system("RNA")
    #vs.set_system("Chromatin")

    
    
    #          0         10        20        30        40        50        60        70        80        90
    #          |         |         |         |         |         |         |         |         |         |
#    ss_seq = ".(((.(((((((.[[...)))))..]].((((...))))..))))).....([)]([)].......(..)..((..)).(..)..........................."
    #ss_seq = "{(((.(((((((.[[...)))))..]].((((...))))..)))))..|..([)]([)]...}.{.(..)..((..)).(..)..}...{.....|...|....|....}"
    #ss_seq = "{((((((.A.AAAAA......))))).BBBB........a..aaaaa....bbbb..).|..([)]([)].ABC.abc.}....{.....}"
    #ss_seq = "{((((((.A.AAAAA......))))).((((........a..aaaaa....))))..).|..([)]([)].ABC.abc.}....{.....}"
    #ss_seq = "{((((((.A.AAAAA......))))).((((........a..aaaaa....))))..).|..([)]([)]...}....{.A.B.C..a.b.c.....}"
    #ss_seq = "{((((((.A.AAAAA......))))).((((........a..aaaaa....))))..).|..([)]([)]...}....{.ABC..DEF..abc..def...}"
    #ss_seq = "{((((((.A.AAAAA.<BC..))))).((((.>bc....a..aaaaa....))))..).|..([)]([)]...}....{.A.B.C..a.b.c.....}"
    #ss_seq = "{((((((.A.AAAAA.<BC..))))).((((.>bc.DE.a..aaaaa..de))))..).|..([)]([)].........}....{.....}"
    #ss_seq  = "{((((((.A.AAAAA.<BC..))))).((((.>bc....a..aaaaa....))))..).|..([)]([)]...}....{.ABC..DEF..abc..def...}"
    
    rnaseq = ""
    # regular secondary structure
    #          0         10        20        30        40        50        60        70        80        90
    #          |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |
    #rnaseq  = "uGGuuuuuuuuuuuuCCuu"
    #ss_seq  = ".((............)).."
    #rnaseq  = "uGGGGuuuuuuuuCCCCuu"
    #ss_seq  = ".((((........)))).."
    #rnaseq = "uGGGGuuGuuuCuCCCCuu"
    #ss_seq  = ".((((..(...).)))).."
    #rnaseq   = "uGGGuGGGuGGGuuuuCCCuCCCuCCCu"
    #ss_seq   = ".(((.(((.(((....))).))).)))."
    rnaseq   = "uGGGuGGGuuuuuuGGGuuuuCCCuuuuuuuCCCuCCCu"
    ss_seq   = ".(((.(((......(((....))).......))).)))."
    
    if len(cl) > 1:
        ss_seq = cl[1]
    vs.parse_fullDotBracketStructure(ss_seq, True)
    v2t = Vienna2TreeNode(vs, None, rnaseq)
    v2t.vienna2tree()
    
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
    #endfor
    
    print "vsPKlist: "
    for bpk in vf.vsPKlist:
        print bpk.disp_Pair()
    #endfor
    
    print "vsMPlist: "
    for bpk in vf.vsMPlist:
        print bpk.disp_Pair()
    #endfor
    
#    



def main(cl):
    if TEST == 0:
        test0(cl)
        
    elif TEST == 1:
        test1(cl)
    
    elif TEST == 2:
        test2(cl)
    #
    
#
   
# Main
if __name__ == '__main__':
    main(sys.argv)
#

