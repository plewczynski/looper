#!/usr/bin/env python3

"""@@@

Main Module:   TreeNode.py 

Objects:       NodeVisitor
               TreeNode2DotBracket
               TreeNode2Motif


Author:        Wayne Dawson
creation date: 170126 (originally part of Threads.py & RThreads.py) 
last update:   200211 (upgraded to python3), 191225
version:       0.1

Purpose:

This module is currently used only by "sarabande.py" to analyze RNA
structure. Because there are some differences in the analysis routines
of the chromatin, RNA and protein based program, I have had to make
these separate.

Comments:

190726:

   TreeNode has gradually become the true structural intermediate
   between Map and dot-bracket representations. The goals of this
   module are two objectives.

   1. TreeNode serves both to build up a 1D structural sequence
      representation (including at least some derive or generated
      sequence of beads) into an object of class Map (a colletion of
      objects of class Motif) or to reduce an object of class Map down
      to comprehendable dot-bracket (1D) structural sequence
      representation.

   2. The other goal is that the TreeNode programs can create
      intermediate test structures so that I can verify whether
      modules are working. So the design should work both from bottom
      up with the prediction programs and top down with the structure
      input sequences.

   Originally, LThread was intended to serve as the true intermediate,
   and it still serves as the intermediate in breaking down the
   suboptimal structure in objects of class Map into dot-bracket
   (1D) structural data. 

   The TreeNode data structures are more domain specific, so it still
   seems difficult to merge different leafs (in an NaryTree) together
   under a single processing strategy, but maybe. If this can be done,
   then we can eliminate this rather silly LThread representation that
   currently hosts the suboptimal structures between Map and
   dot-bracket representations.

   The remaining advantage of LThread is still that it is greatly
   simplified so it doesn't eat up quite as much memory as TreeNode
   would. The idea situation is that both the LThread and TreeNode
   modules can be operated rather universally. I'm not sure this can
   be achieved, but I am still thinking the matter over.


earlier comments (updated 190705)

   The _intention_ was that this code would be shared among a variety of
   programs and modules.

   The current objective is to always produce human readable output,
   regardless of whether it comes from Calculate (an actual
   computation of a structure) or from an input structural sequence
   (dot bracket format or otherwise) along with a sequence
   representation.

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
from Motif import ins_sort_StemList

"""for Motif data representation"""
from Motif      import Motif # a core object of these building programs
from Motif      import Link
from Motif      import LGroup
from Motif      import Map # main map for all the FE and structures

from LoopRecords import MBLptr

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

from LThreadBuilder import LThreadBuilder

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

PROGRAM      = "TreeNode.py"  # name of the program

# This is not used so much now that I have introduced GetOpts.py,
# but it may still be useful at the very beginning of the program and
# possibly in other parts.

def usage():
    print ("USAGE: %s" % PROGRAM)
#

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ################################################################
# ################################################################



"""For an explanation why!?! we need object, please see
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

In this module, NodeVisitor is used by TreeNode2DotBracket and
TreeNode2Motif

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
            print (method_name)
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




class TreeNode2DotBracket(NodeVisitor):
    def __init__(self, v2t, is_chromatin = True):
        inputObject = type(v2t).__name__
        if not inputObject == "Vienna2TreeNode":
            print ("Error: TreeNode2DotBracket requires objects of class Vienna2TreeNode.")
            print ("       object entered: %s" % inputObject)
            sys.exit(1)
        #
        
        self.debug_TreeNode2DotBracket = False # True # 
        
        self.genTree = v2t.genTree # the main tree we need
        self.ssStems = v2t.ssStems # helps with searching out branching
        self.MPlist  = v2t.MPlist
        self.N       = v2t.N
        
        # sequence stuff 
        self.seqv    = []
        self.sssv    = []
        
        self.counter = 1
        self.counter = 1
        self.reset_counter()
        
        
        self.chromatin = is_chromatin # specifies that the sequence is chromatin
        
        
        if self.debug_TreeNode2DotBracket:
            print ("here, is_chromatin = ", is_chromatin)
            print ("      v2t.vstr: ", v2t.vstr)
            print ("      v2t.vseq: ", v2t.vseq)
        #
        
        
        if len(v2t.vseq) == len(v2t.vstr) and len(v2t.vseq) == self.N:
            self.seqv = list(v2t.vseq)
        else:
            if self.chromatin:
                # no sequence provide, so we just make a blank one
                self.seqv = self.N*['c']
            else:
                print ("ERROR(TreeNode2DotBracket): RNA lacks a defined sequence!")
                sys.exit(1)
            #
            
        #
        
        self.sssv = self.N*['.']
        
        
        self.vstr = self.makeFinal(self.sssv)
        self.vSeq = ''
        #print (self.seqv)
        #print (self.sssv)
        #sys.exit(0)
        
        # presently, MP structures are automatically written. I am
        # thinking that these should be written on a separate line,
        # but it is convenient presently to have the { | | | }
        # notation present.
        self.write_MP()
        
        
        self.initialize_t2m   = True
        #print ("initialized TreeNode2DotBracket(), counter = ", self.counter)
    #
    
    def make_dotbracket(self, is_chromatin = True):
        self.chromatin = is_chromatin
        
        self.reset_counter()
        # chromatin also generates the sequence
        self.visit(self.genTree)
        self.vstr = self.makeFinal(self.sssv)
        
        self.vSeq = self.makeFinal(self.seqv)
        #print (self.vSeq)
        
        return self.vstr
    #        
    
    def set_seqv(self, seqv):
        self.seqv = list(seqv)
    #
    
    
    
    def makeFinal(self, seqList):
        seqx = ''.join(seqList) 
        # also: python> ''.join(sssv); python> ''.join(seqList) 
        return seqx
    #
    
    def reset_counter(self):
        self.counter = 1
    #
    
    def inc_count(self):
        self.counter += 1
        if self.counter == 2:
            """
            avoids current MP notations {...}; i.e., "..{...|...|...}.."
            
            I think in the future, I would prefer to use WXYZwxyz but
            for representation, VARNA (which recognizes this general
            format) only goes up to Nn. So, since the appearance is
            nice, I am still using this rather personal definition for
            islands. At any rate, this problem will emerge whether I
            like it not, even if I go to Ww/Xx/Yy/Zz.
            
            """
            self.counter += 1
        #
        
        if self.counter >= len(num2lpr):
            self.reset_counter() # self.counter = 1
        #
        
        #print ("counter: ", self.counter)
    #
    
    def visit_Stem(self, node):
        debug_visit_Stem = False # True # 
        if debug_visit_Stem:
            print ("Enter visit_Stem: ", node.name)
        #endif
        
        vstem = node.name
        if not vstem.vtype == "root":
            self.write_Stem(vstem)
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
        debug_visit_PseudoKnot = False # True # 
        
        vpk = node.name
        ipk = vpk.it; jpk = vpk.jt
        if debug_visit_PseudoKnot:
            print ("Enter visit_PseudoKnot", vpk, vpk.pktype)
            print (vpk)
        #
        
        for stem in vpk.rootstems:
            if debug_visit_PseudoKnot:
                print ("vpk: rootstem: ", stem, stem.vtype)
            #
            
            if type(stem).__name__ == "PseudoKnot":
                if debug_visit_PseudoKnot:
                    print ("begin pk in pk")
                #
                
                self.write_PseudoKnot(stem)
                #print ("finished pk in pk"); sys.exit(0)
                
            else:
                self.write_Stem(stem)
            #
        #endfor
        
        if debug_visit_PseudoKnot:
            print ("all linkages: ", vpk.linkages)
        #
        
        for stem in vpk.linkages:
            if debug_visit_PseudoKnot:
                print ("vpk: linkage: ", stem, stem.vtype)
            #
            
            self.write_Linkage(stem)
            self.inc_count() # self.counter += 1
        #endfor
        
        if not node.children == None:
            if len(node.children) > 0:
                for child in node.children:
                    self.visit(child)
                #endfor
                
            #
            
        #
        
    #
    
    
    def write_PseudoKnot(self, vpk):
        debug_visit_PseudoKnot = False # True # 
        for stem in vpk.rootstems:
            if debug_visit_PseudoKnot:
                print ("vpk: rootstem: ", stem, stem.vtype)
            #
            
            if type(stem).__name__ == "PseudoKnot":
                if debug_visit_PseudoKnot:
                    print ("begin pk in pk inside pk")
                #
                
                self.write_PseudoKnot(stem)
                #print ("finished pk in pk"); sys.exit(0)
                
            else:
                self.write_Stem(stem)
            #
        #endfor
        
        if debug_visit_PseudoKnot:
            print ("all linkages: ", vpk.linkages)
        #
        
        for stem in vpk.linkages:
            if debug_visit_PseudoKnot:
                print ("vpk: linkage: ", stem, stem.vtype)
            #
            
            self.write_Linkage(stem)
            self.inc_count() # self.counter += 1
        #endfor
    #
    
    
    def write_Linkage(self, vstem):
        # stype = 's' (strandard) or 'l' (linkage)
        
        debug_write_Linkage = False # True # 
        
        i_t = vstem.it; j_t = vstem.jt
        i_h = vstem.ih; j_h = vstem.jh
        
        vtype = vstem.vtype
        if vtype == 'p':
            for stmk in vstem.stem:
                self.sssv[stmk.i] = num2lpr[self.counter]
                self.sssv[stmk.j] = num2rpr[self.counter]
                self.inc_count() # self.counter += 1
                if self.chromatin:
                    self.seqv[stmk.i] = 'x'
                    self.seqv[stmk.j] = 'y'
                #
                
            #endfor
            
        else:
            for stmk in vstem.stem:
                self.sssv[stmk.i] = num2lpr[self.counter]
                self.sssv[stmk.j] = num2rpr[self.counter]
                if self.chromatin:
                    self.seqv[stmk.i] = 'x'
                    self.seqv[stmk.j] = 'y'
                #
                
            #endfor
            
        #
    #
    
    
    def write_Stem(self, vstem):
        # stype = 's' (strandard) or 'l' (linkage)
        
        debug_write_StemThread = False # True # 
        if debug_write_StemThread:
            print ("write_Stem: ", vstem)
        #
        
        i_t = vstem.it; j_t = vstem.jt
        i_h = vstem.ih; j_h = vstem.jh
        
        vtype = vstem.vtype
        if vtype == 'p':
            for stmk in vstem.stem:
                self.sssv[stmk.i] = num2lpr[self.counter]
                self.sssv[stmk.j] = num2rpr[self.counter]
                self.inc_count()  # self.counter += 1
                if self.chromatin:
                    self.seqv[stmk.i] = 'x'
                    self.seqv[stmk.j] = 'y'
                #
                
            #endfor
            
        else:
            for stmk in vstem.stem:
                self.sssv[stmk.i] = '('
                self.sssv[stmk.j] = ')'
                if self.chromatin:
                    
                    self.seqv[stmk.i] = 'x'
                    self.seqv[stmk.j] = 'y'
                #
                
            #endfor
            
            
        #
        
    #
    
    
    def write_MP(self):
        if len(self.MPlist) > 0:
            for wyspa in self.MPlist:
                i_W = wyspa.it; j_W = wyspa.jt
                
                if self.chromatin:
                    self.seqv[i_W] = 'W'; self.seqv[j_W] = 'Z'
                #
                
                for k in range(0, len(wyspa.arches)):
                    iw = wyspa.arches[k].it
                    jw = wyspa.arches[k].jt
                    
                    if k == 0:
                        self.sssv[i_W] = '{'
                        self.sssv[j_W] = '}'
                        if self.chromatin:
                            self.seqv[i_W] = 'W'
                            self.seqv[j_W] = 'Z'
                        #
                        
                    else:
                        if i_W <= iw and jw < j_W:
                            self.sssv[jw] = '|'
                            if self.chromatin:
                                self.seqv[jw] = 'I'
                            #
                            
                        #
                        
                    #
                    
                #endfor
            #endfor
            
        # end [if len(self.MPlist) > 0:]
    #
    
#





class TreeNode2Motif(NodeVisitor):
    def __init__(self, v2t):
        flag_debug = False
        inputObject = type(v2t).__name__
        if not inputObject == "Vienna2TreeNode":
            print ("Error: TreeNode2Motif requires objects of class Vienna2TreeNode.")
            print ("       object entered: %s" % inputObject)
            sys.exit(1)
        #
        
        self.ssStems  = v2t.ssStems  # List of Stem from secondary structure
        self.genStems = v2t.genStems # List of Stem and PseudoKnot types
        # I suspect that we should use genStems, but I'll leave it this way
        self.MPlist   = v2t.MPlist   # List of Islands
        self.N        = v2t.N
        self.fe       = v2t.fe       # currently RNAModules or ChromatinModules
        self.dGmin    = 1000.0
        self.dHmin    = 1000.0
        self.initialize_t2m = True
        
        if flag_debug:
            print (self.fe.rnaseq)
            
            for j in range(0, len(self.fe.btype)):
                for i in range(0, j):
                    print ("%3d, %3d: " % (i, j), self.fe.btype[i][j])
                #
            #
            
            # print ("initialized TreeNode2Motif()")
            #sys.exit(0)
        #
        
    #
    
    def build_wyspa(self, mpstruct):
        debug_build_wyspa = False # True #
        if debug_build_wyspa:
            print (mpstruct)
            print (type(mpstruct).__name__)
        #
        
        if len(mpstruct.arches) > 0:
            wyspa = []
            iW = mpstruct.arches[0].it; jW = mpstruct.arches[0].jt
            for k in range(1, len(mpstruct.arches)):
                ww = mpstruct.arches[k]
                iw = ww.it; jw = ww.jt
                
                if debug_build_wyspa:
                    print ("ijW(%2d,%2d) <- ijw(%2d,%2d)" % (iW, jW, iw, jw))
                #
                
                if jw < jW:
                    wyspa += [jw]
                #
                
            #endfor
        #
        
        return wyspa
    #
    
    
    def post_graftMP(self):
        
        if len(self.MPlist) > 0:
            debug_post_graftMP = False # True #
            
            dGMP = 0.0
            MProots = []
            for island in self.MPlist:
                
                i_W = island.it; j_W = island.jt
                MProots += [(i_W, j_W)]
                join = [(i_W, j_W)]
                wyspa = self.build_wyspa(island)
                # xx!!!! 190707 !!!!xx need a real potential
                dGW = -10.0 # just an arbitrary value presently
                for k in range(0, len(island.arches)):
                    ww =  island.arches[k]
                    i_w = ww.it; j_w = ww.jt
                    wbr = []
                    if k > 0 or len(island.arches) == 1:
                        # get substructure inside the MultiPair region 
                        wbr = self.get_genbranches(i_w, j_w)
                    #
                    
                    for vv in wbr:
                        iv = vv[0]; jv = vv[1]
                        dGW += self.fe.smap.glink[iv][jv].lg[0].get_Vij()
                    #
                    
                #endfor
                
                """@ 
                
                The more I look at this, the more I realize that it
                doesn't make sense the way it is written. This whole
                thing needs to be fixed before it can be integrated
                with the rest of the program.
                
                """
                
                mpmotif = MultiPair(i_W, j_W, 'wyspa')
                mpmotif.Vij = dG
                
                # somehow, need to write in the list "wyspa" into this thing
                
                # original # link_W = Link(i_W, j_W, dGW, 'W', 'wyspa', join, [], wyspa)
                self.fe.smap.glink[i_W][j_W].add_link(Link(mpmotif))
                dGMP += dGW
            #endfor
            
            i_t = 0; j_t = self.N-1
            self.assign_P(i_t, j_t, MProots, 'P', '-')
            
            if debug_post_graftMP:
                print ("result")
                for mm in self.fe.smap.glink[i_t][j_t].lg:
                    print (mm.motif.show_Motif())
                #endfor
                
            #
            
        #
    #
    
    
    
    def visit_Stem(self, node):
        debug_visit_Stem = False # True # 
        if debug_visit_Stem:
            it = node.name.it; jt = node.name.jt
            ih = node.name.ih; jh = node.name.jh
            nm = node.name.name
            print ("Enter visit_Stem{ijt(%2d,%2d)::ijh(%2d,%2d) %s}" % (it, jt, ih, jh, nm))
        #
        
        # depth first
        if not node.children == None:
            if len(node.children) > 0:
                for child in node.children:
                    self.visit(child)
                #endfor
                
            #
            
            """@
            
            190707: I not completely sure how the new embedded
            approach factors into this, but I think the treeNode
            already would be aware of this. It is built from Vienna
            dot-bracket notation, so it is the structure.
            
            """
        #
        
        # << postorder action >>
        
        vstem = node.name
        """@
        
        objects of class "Stem" or "Pseudoknot" are stored in class
        "Node" (which is handled by class "NaryTree"; i.e., N-ary
        instead of bi-nary) under the handle "name". Therefore, 
        
        name -> Stem or PseudoKnot
        
        """
        
        if vstem.vtype == "root":
            if debug_visit_Stem:
                print ("Establishing the root of the TreeNode")
            #
            
            i_t = vstem.ih; j_t = vstem.jh
            """@
            
            The root of Tree only contains referemce to beginning
            of the polymer sequence and the end of the sequence; i.e.,     
            for the root of the tree, vstem (i_t,j_t) <-> (i_h,j_h).
            
            From here, we must expand the actual connections. To do
            that, we first must scan between (i_t,j_t) -- bounds
            inclusive -- for any structure. The minimum free energy
            structure (or the user-set structure) can include anything
            including the point (0, N-1); i.e., the boundaries.
            
            VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
            190731(NOTE): does not compute the 3'/5' tail free
            energy!!!!!  This needs a dangling bond assignment in the
            method get_branchingFE().
            AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
            
            """
            gbranches = []
            if not node.children == None:
                if len(node.children) > 0:
                    for child in node.children:
                        gbranches += [(child.name.it, child.name.jt)]
                    #
                    
                #
            #
            
            ngbranches = len(gbranches)
            
            """@
            
            The root position (0, N-1) should be assigned either a
            X/J/P if the meta-structures don't fall precisely on
            (i_t,j_t) or B/I/M/S/K/R/W if the meta-structure land
            exactly on the boundaries (i_t, j_t).
            
            """
            
            Vpqh, Hpqh = self.get_branchingFE(gbranches)
            ib = 0; jb = self.N-1
            
            if debug_visit_Stem:
                print ("base: (%3d,%3d): gbranches -> %s" % (i_t, j_t, gbranches))
                print ("      Vij = %8.2f, Hij = %8.2f" % (Vpqh, Hpqh))
                #print ("planned exit at 1 in visit_Stem"); sys.exit(0)
            #
            
            
            if ngbranches == 0: # root and no branches
                """@
                
                The whole thing is either completely empty of
                significant meta-structure (X-loop) or contains only
                one isolated B-loop.
                
                """
                
                if self.fe.btype[i_t][j_t].pair > 0:
                    # there is meta-structure consisting of one
                    # isolated pair (B-loop)
                    
                    # add the head dangle and corrections
                    dG_H =  self.fe.cle_HloopE(i_t,  j_t)
                    # add the tail dangle and corrections
                    dG_H += self.fe.dG_stem_dangle(ib,  jb,   # boundaries 
                                                   i_t, j_t,  # 5'3' tail of the stem
                                                   True,      # T -> tail dangle
                                                   False)     # invoke debug option
                    
                    rmotif = XLoop(i_t, j_t, dG_H, 'B', 'sa')
                    
                    # enthalpy
                    rmotif.dH_dVlp = self.fe.get_dHhloop_mm1(i_t, j_t)
                    rmotif.dHpqh   = 0.0
                    rmotif.Hij     = rmotif.dH_dVlp
                    # btp is irrelevant whether 'sa' or 'sp'
                    
                    """@
                    
                    I think this needs some further consideration for
                    parallel and antiparallel cases, so, for example,
                    looking at i_t - 1, j_t + 1 vs i_t -1, j_t -
                    1. One problem is the prior bp is not passed in
                    the current design.
                    
                    """
                    
                    self.fe.smap.glink[i_t][j_t].add_link(Link(rmotif))
                    
                    if debug_visit_Stem:
                        print ("B-loop closing final sequence (%d,%d)" % (i_t, j_t))
                        print (rmotif)
                    #
                    
                else: 
                    """@
                    
                    There is absolutely no meta-structure at all. This
                    can happen for a poorly chosen structure or
                    sequence; however, it is rare.
                    
                    In general, with turner rules, this shouldn't
                    happen and if it does, unless this is evaluated by
                    the gMatrix, it will surely produce an error.  """
                    
                    a = self.fe.aSeq[i_t]
                    b = self.fe.aSeq[j_t]
                    c = self.fe.aSeq[i_t+1]
                    d = self.fe.aSeq[j_t-1]
                    
                    # It would make sense to say something about this one
                    print ("WARNING/ERROR: terminal H-loop closed by nWC bp.")
                    print ("               i(%d, %s), j(%d, %s)" % (i_t, a, j_t, b))
                    print ("               %s%s/%s%s" % (a, c, b, d))
                    
                    # it would generally be strange to want the FE of
                    # a structure like this, in general, so about the
                    # only thing to do is calculate the ends.
                    dG_X = self.fe.cle_HloopE(i_t,  j_t)
                    
                    rmotif = XLoop(i_t, j_t, dG_X, 'X', '-')
                    
                    # enthalpy
                    rmotif.Hij = self.fe.get_dHhloop_mm1(i_t, j_t)
                    # (no entropy contributions in enthalpy)
                    
                    self.fe.smap.glink[i_t][j_t].add_link(Link(rmotif))
                    
                    
                    if debug_visit_Stem:
                        print ("no structure in final sequence (%d,%d)" % (i_t, j_t))
                        print (rmotif)
                    #
                #
                
            elif ngbranches == 1: # root and one branch
                i_h = gbranches[0][0]; j_h = gbranches[0][1]
                if (i_h == i_t and j_h == j_t):
                    """@
                    
                    The ends of the sequence are a stucture: I, M, K,
                    S, R, or W. This would already have been
                    calculated. Therefore, we skip it.
                    
                    """
                    
                    if debug_visit_Stem:
                        print ("ngbranches == 1 and the single structure of class %s" \
                            % (vstem.name))
                        print ("closes the far ends of the sequence: ijt = (%2d,%2d)" \
                            % (i_t, j_t))
                    #
                    
                else:
                    
                    """The ends of the sequence are a J-loop."""
                    
                    if debug_visit_Stem:
                        print ("ngbranches == 1 and the single structure of class %s" \
                            % (vstem.name))
                        print ("forms a domain at ijh = (%2d,%2d)" \
                            % (i_h, j_h))
                    #
                    
                    stmtail = True
                    dG_td = self.fe.dG_stem_dangle(ib,  jb,   # boundaries 
                                                   i_h, j_h,  # 5'3' tail of the stem
                                                   stmtail,   # T -> tail dangle
                                                   False)     # invoke debug option
                    dG_J  = Vpqh + dG_td
                    # Vpqh is just the FE without the dangle 
                    
                    jmotif = MBL(i_t, j_t, dG_J, 'J', '-', gbranches)
                    jmotif.vtype   = vstem.vtype
                    
                    # here, the I-loop is already complete, so we just
                    # use Vpqh and ignore dG_dVlp.
                    
                    jmotif.dG_dVlp = dG_td # branch dangle + 5'3'-head + cle closing
                    jmotif.dGVpq   = Vpqh  # branch FE (uncorrected)
                    
                    # enthalpy
                    dH_td = self.fe.dH_stem_dangle(ib,  jb,   # boundaries 
                                                   i_h, j_h,  # 5'3' tail of stem
                                                   stmtail,   # T -> tail dangle
                                                   False)     # invoke debug option
                    
                    jmotif.dH_dVlp = dH_td # branch dangle + 5'3'-head enthalpy
                    jmotif.dHVpq   = Hpqh  # branch enthalpy (uncorrected)
                    jmotif.Hij     = Hpqh + dH_td
                    # Hpqh is just the enthalpy without the dangle 
                    
                    self.fe.smap.glink[i_t][j_t].add_link(Link(jmotif))
                    
                    if debug_visit_Stem:
                        print ("J-loop finishes final sequence (%d,%d)" % (i_t, j_t))
                        print (jmotif)
                    #
                    
                #
                
            elif ngbranches > 1: # root and more than one branch
                
                # Well ... if there is more than one branch, we cannot
                # have an M-loop, so it must be that all branches are
                # considered, there is no skipping here, and we can
                # assume without further ado that it is a P-loop.
                
                if debug_visit_Stem:
                    print ("ngbranches > 1 and the single structure of class %s" \
                        % (vstem.name))
                    print ("closes the far ends of the sequence: ijh = (%2d,%2d)" \
                        % (i_t, j_t))
                #
                
                new_pMBL = MBLptr(ib, jb)
                new_pMBL.nm = 'P'
                new_pMBL.dr = '-'
                new_pMBL.V = Vpqh
                for gbrk in gbranches:
                    iv = gbrk[0]; jv = gbrk[1]
                    new_pMBL.addBranch(iv, jv)
                #
                
                new_pMBL.n = len(gbranches)
                
                
                """The ends of the sequence are a P-loop."""
                dG_td = self.fe.dG_trim_53pStems(ib, jb,   # 5' to 3' [i/p]MBL bounds 
                                                 new_pMBL, # mbl branch info
                                                 False)    # iMBL or pMBL? B.Cs
                
                dG = Vpqh + dG_td
                pmotif = MBL(i_t, j_t, dG, 'P', '-', gbranches)
                pmotif.vtype   = vstem.vtype
                pmotif.dG_dVlp = dG_td  # branch dangle + 5'3'-head + cle closing
                pmotif.dGVpq   = Vpqh   # branches FE *uncorrected)
                
                # !!!!!
                """The ends of the sequence are a P-loop."""
                dH_td = self.fe.dH_trim_53pStems(ib, jb,   # 5' to 3' [i/p]MBL bounds 
                                                 new_pMBL, # mbl branch info
                                                 False)    # iMBL or pMBL? B.Cs
                
                pmotif.dH_dVlp = dH_td  # branch dangle + 5'3'-head enthalpy
                pmotif.dHVpq   = Hpqh   # branch enthalpy (uncorrected)
                pmotif.dHd53p  = 0.0    # open so there can be no tail 5'3' dangles
                pmotif.Hij     = Hpqh + dH_td
                
                self.fe.smap.glink[i_t][j_t].add_link(Link(pmotif))
                
                if debug_visit_Stem:
                    print ("P-loop finishes final sequence (%d,%d)" % (i_t, j_t))
                    print (pmotif)
                #
                
                
            #
            
            self.dGmin = self.fe.smap.glink[i_t][j_t].lg[0].motif.get_Vij()
            self.dHmin = self.fe.smap.glink[i_t][j_t].lg[0].motif.get_Hij()
            
            if debug_visit_Stem:
                print ("vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv")
                print ("result: dGmin = %8.2f [kcal/mol]" % self.dGmin)
                print ("        dHmin = %8.2f [kcal/mol]" % self.dHmin)
                for mm in self.fe.smap.glink[i_t][j_t].lg:
                    print (mm.motif.show_Motif())
                #endfor
                
                print ("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")

                #print ("planned exit at end of root"); sys.exit(0)
                
            #
            
            
            
        else:
            # in general, we will be writing here
            
            # print ("went in here ")
            
            # 190707: maybe the trick here is to make this part of
            # RNAModule or ChromatinModule, then almost everything
            # else should be the same. If it is accessed through
            # Calculate or RCalculate, then the details on the method
            # are completely invisible to this module at this level.
            self.write_StemMotif(vstem)
            #sys.exit(0)
        #
        
    #
    
    
    
    
    def visit_PseudoKnot(self, node):
        debug_visit_PseudoKnot = True # False # 
        if debug_visit_PseudoKnot:
            print ("enter visit_PseudoKnot")
        #
        # depth first
        if not node.children == None:
            if len(node.children) > 0:
                for child in node.children:
                    self.visit(child)
                #endfor
                
            #
            
        #
        
        # << postorder action >>
        
        vpk = node.name
        ipk = vpk.it; jpk = vpk.jt
        #print (vpk)
        #print (vpk.disp_PseudoKnot())
        
        vtype = 'a'  
        pktype = vpk.pktype
        # have to consider fixing this point depending on the type of
        # linkage
        
        id1 = vpk.rootstems[0].it; jd1 = vpk.rootstems[0].jt;
        ih1 = vpk.rootstems[0].ih; jh1 = vpk.rootstems[0].jh;
        id2 = -1; jd2 = -1
        ctp = 'K'
        if len(vpk.rootstems) == 2:
            # extended PK
            ctp = 'R'
            id2 = vpk.rootstems[1].it; jd2 = vpk.rootstems[1].jt
            
        elif len(vpk.rootstems) > 2:
            print ("ERROR! too many branches")
            print (vpk.rootstems)
            
            sys.exit(1)
        #
        
        #print (vpk.rootstems)
        # #################################
        # ###   process the root stems  ###
        # #################################
        
        rt_stm_brnch = [] # root stem branch
        kpk = 0
        dG_root = 0.0
        dH_root = 0.0
        all_rootStems = []
        for rstem in vpk.rootstems:
            ir = rstem.it; jr = rstem.jt
            # print ("rootstem: ", stem)
            rt_stm_brnch += [(ir, jr)]
            self.write_StemMotif(rstem, 's', [])
            all_rootStems += [self.fe.smap.glink[ir][jr].lg[0].motif.Ob]
            
            # xx!!!! 190707 !!!!xx have to develop further.
            dG_root += self.fe.smap.glink[ir][jr].lg[0].motif.get_Vij()
            dH_root += self.fe.smap.glink[ir][jr].lg[0].motif.get_Hij()
            if debug_visit_PseudoKnot:
                print ("root stem %d:  ijr(%3d,%3d)  dG = %8.2f  dH = %8.2f" \
                    % (kpk + 1, ir, jr,
                       self.fe.smap.glink[ir][jr].lg[0].motif.get_Vij(),
                       self.fe.smap.glink[ir][jr].lg[0].motif.get_Hij()))
            #
            
            kpk += 1
            
        #endfor
        
        if debug_visit_PseudoKnot:
            print ("rt_stm_brnch: ", rt_stm_brnch)
        #
        
        nrt_stm_brnch = len(rt_stm_brnch)
        if nrt_stm_brnch <= 0:
            print ("ERROR: ??? don't understand what to make of this case!")
            print ("       nssbranches = %d, ij_h(%2d,%2d)[%s]" % (i_h, j_h, vtype))
            sys.exit(1)
        #
        
        
        # #####################################
        # ###   process the linkage stems   ###
        # #####################################
        
        fr_zones_l = []
        if len(vpk.linkages) == 1:
            if   len(vpk.rootstems) == 1:
                
                # basically should be only one linkage (not linkageS)
                pL = vpk.linkages[0].it; qL = vpk.linkages[0].jt
                pl = vpk.linkages[0].ih; ql = vpk.linkages[0].jh
                
                # core PK or combination of core pk plus other structure
                if vpk.rootstems[0].name == "Stem":
                    fr_zones_l = [(jd1, ql)]
                elif vpk.rootstems[0].name == "PseudoKnot":
                    fr_zones_l = [(jd1, ql)]
                else:
                    print ("ERROR: undefined object (%s) at (%2d,%2d)" \
                        % (vpk.rootstems[0].name, id1, jd1))
                    sys.exit(1)
                #
                
            elif len(vpk.rootstems) == 2:
                # extended PK
                ctp = 'R'
                fr_zones_l += [(jd1, id2)]
                
            else:
                print ("ERROR: 1 Too few linkages")
                print (vpk.linkages)
                sys.exit(1)
            #
            
        elif len(vpk.linkages) > 1:
            """@
            
            190731: I think it is correct that we can have more than
            one linkage on the same PK. I have to think about how to
            arrange that. In building the tree, the linkages can
            definitely be grouped, but it is not clear how I will
            design the fr_zones for that. It is an ordered list, so it
            seems like there would be a way, but I will have to think
            it over a bit more carefully. In the meantime, I will only
            test single linkages, because these are sure to be the
            most common anyway. It should be remembered that this may
            need some upgrading when testing much harder problems.
            
            """
            vpk.linkages = ins_sort_StemList(vpk.linkages, "jt")
            # this orders the linkages according to the largest index jt.
            
            # have to add code here
            # basically should be only one linkage (not linkageS)
            
            if   len(vpk.rootstems) == 1:
                # core pseudoknot
                """@
                
                To set this up for make sure that the order of the
                list is from smallest to largest with respect to the
                index jt; i.e.,
                
                jd1 < ql0 < qL0 <  ql1 < qL1, etc.
                
                So we build n free strand zones such that
                
                [(jd1, pl0), (pL0, pl1), ... (pL[n-1], pln)]
                
                """
                
                pL = vpk.linkages[0].it; qL = vpk.linkages[0].jt
                pl = vpk.linkages[0].ih; ql = vpk.linkages[0].jh
                fr_zones_l = [(jd1, ql)]
                for k in range(1, len(vpk.linkages)):
                    qL_km1 = vpk.linkages[k-1].jt
                    ql_k   = vpk.linkages[k  ].jh
                    fr_zones_l += [(qL_km1, ql_k)]
                #
                
            elif len(vpk.rootstems) == 2:
                
                # extended pseudoknot
                """@
                
                To set this up for make sure that the order of the
                list is from smallest to largest with respect to the
                index jt; i.e.,
                
                jd1 < ql0 < qL0 <  ql1 < qL1, etc. up to id2
                
                So we build n free strand zones such that
                
                [(jd1, pl0), (pL0, pl1), ... (pL[k-1], id2)]
                
                """
                ctp = 'R'
                id2 = vpk.rootstems[1].it; jd2 = vpk.rootstems[1].jt;
                
                ql_0 = vpk.linkages[0].jh
                if ql_0 < id2:
                    fr_zones_l = [(jd1, ql_0)]
                    
                    for k in range(1, len(vpk.linkages)):
                        qL_km1 = vpk.linkages[k-1].jt
                        ql_k   = vpk.linkages[k  ].jh
                        if ql_k < id2:
                            fr_zones_l += [(qL_km1, ql_k)]
                        else:
                            fr_zones_l += [(qL_km1, id2)]
                            break
                        #
                    #
                    
                else:
                    fr_zones_l = [(jd1, id2)]
                #
                
            else:
                print ("ERROR: 2 Too few linkages")
                print (vpk.linkages)
                sys.exit(1)
            #
            
        #
        
        dG_link = 0.0
        dH_link = 0.0
        linktail = []
        all_linkStems = []
        for lstem in vpk.linkages:
            pL = lstem.it; qL = lstem.jt
            linktail += [(pL, qL)]
            # print ("linkage: ", lstem.stem)
            self.write_StemMotif(lstem, 'l', fr_zones_l)
            all_linkStems += [self.fe.smap.glink[pL][qL].lg[0].motif.Ob]
            dG_link += self.fe.smap.glink[pL][qL].lg[0].motif.get_Vij()
            dH_link += self.fe.smap.glink[pL][qL].lg[0].motif.get_Hij()
            
            if debug_visit_PseudoKnot:
                print ("linkage:       pqL(%3d,%3d)  dG = %8.2f  dH = %8.2f" \
                    % (pL, qL,
                       self.fe.smap.glink[pL][qL].lg[0].motif.get_Vij(),
                       self.fe.smap.glink[pL][qL].lg[0].motif.get_Hij()))
            #
            
        #
        
        if debug_visit_PseudoKnot:
            print ("linktail: ", linktail)
        #
        
        if debug_visit_PseudoKnot:
            print ("dG linkages ijpk(%3d,%3d)[dG=%8.2f)" \
                % (ipk, jpk, dG_link))
        #
        
        dGpk = dG_root + dG_link
        dHpk = dH_root + dH_link
        # print ("pk(%3d,%3d), dG = %8.2f" % (ipk, jpk, dGK))
        pkmotif = PseudoKnot(ipk, jpk, ctp)
        pkmotif.roottmp   = deepcopy(rt_stm_brnch)
        pkmotif.rootstems = deepcopy(all_rootStems)
        pkmotif.linkages  = deepcopy(all_linkStems)
        pkmotif.Vij = dGpk
        pkmotif.Hij = dHpk
        """
        pkmotif.branching = []  # now solved inside of rootstems and linkages
        pkmotif.Vbr = 0.0       # now solved inside of rootstems and linkages

        # see comments in Motif within class PseudoKnot.
        
        """
        
        
        self.fe.smap.glink[ipk][jpk].add_link(Link(pkmotif))
        if debug_visit_PseudoKnot:
            for mm in self.fe.smap.glink[ipk][jpk].lg:
                print (mm.motif.show_Motif())
            #endfor
            
            #sys.exit(0)
        #
        
    #
    
    def assign_B(self, i_h, j_h, stype, vtype):
        debug_assign_B = True # False # 
        plan_stop_on_pass = False
        if debug_assign_B:
            branches = []
            print ("Enter assign_B{ij_h(%d,%d), branches = %s" % (i_h, j_h, branches))
        #
        
        # has pairing at i_h, j_h
        opt_iMBL = True # T -> closed at (i_h,j_h), F -> open
        
        # this is the case of a leaf, there is no other
        # interaction but this one.
        
        ctp = 'B'
        btp = "%s%s" % (stype, vtype)
        if not (vtype == 'a' or vtype == 'p'):
            btp = vtype
        #
        
        
        dG = self.fe.cle_HloopE(i_h,  j_h)
        dH = self.fe.dH_HloopE(i_h,  j_h)
        
        # add any corrections for tetraloop, local loop
        # entropy or whatever here
        
        bmotif = XLoop(i_h, j_h, dG, ctp, btp)
        bmotif.Hij = dH
        self.fe.smap.glink[i_h][j_h].add_link(Link(bmotif))
        
        if debug_assign_B:
            for mm in self.fe.smap.glink[i_h][j_h].lg:
                print (mm.motif.show_Motif())
            #endfor

            if plan_stop_on_pass:
                print ("planned exit in TreeNode.assign_B");
                sys.exit(0)
            #
            
        #
        
        #print ("exit at end of assign_B"); sys.exit(0)
        
    #
    
    def assign_X(self, i_h, j_h, stype, vtype):
        debug_assign_X = True # False #
        plan_stop_on_pass = False
        if debug_assign_X:
            branches = []
            print ("Enter assign_X{ij_h(%d,%d), branches = %s" % (i_h, j_h, branches))
        #
                
        # has no pairing at i_h, j_h
        opt_iMBL = False # T -> closed at (i_h,j_h), F -> open
        
        # this is the case of a leaf with empty pairing
        # interaction. Not only is there is no other interactions,
        # there is no pairing interaction at all.
        
        ctp = 'X'
        btp = '-'
        
        dG = self.fe.find_branch_dGwt(ih, jh, []) # 
        dH = 0.0
        # add any corrections for tetraloop, local loop
        # entropy or whatever here
        
        bmotif = XLoop(i_h, j_h, dG, ctp, btp)
        bmotif.Hij = dH
        self.fe.smap.glink[i_h][j_h].add_link(Link(bmotif))
        if debug_assign_X:
            for mm in self.fe.smap.glink[i_h][j_h].lg:
                print (mm.motif.show_Motif())
            #endfor
            
            if plan_stop_on_pass:
                print ("planned exit in TreeNode.assign_X");
                sys.exit(0)
            #
        #
        
    #
    
    
    def assign_I(self, i_h, j_h, branches, stype, vtype):
        debug_assign_I = False # True # 
        plan_stop_on_pass = False
        
        if debug_assign_I:
            print ("Enter assign_I{ij_h(%d,%d), branches = %s" % (i_h, j_h, branches))
        #
        
        # has pairing at i_h, j_h
        opt_iMBL = True # T -> closed at (i_h,j_h), F -> open
        
        ctp = 'I'
        btp = "%s%s" % (stype, vtype)
        if not (vtype == 'a' or vtype == 'p'):
            btp = vtype
        #
        
        # pairing at i_h, j_h        
        Vpqh, Hpqh = self.get_branchingFE(branches)
        # Vpqh = FE contribution from I/M branches or 0.0 for B 
        Vclose = self.fe.find_branch_dGwt(i_h, j_h, branches, opt_iMBL)
        
        
        # add any corrections for tetraloop, local loop
        # entropy or whatever here
        
        dG = Vpqh + Vclose
        
        imotif = MBL(i_h, j_h, dG, ctp, btp, branches)
        imotif.vtype   = vtype
        imotif.dG_dVlp = Vclose # branch dangle + 5'3'-head + cle closing
        imotif.dGVpq   = Vpqh   # FE from branches within (ph, qh)
        imotif.dG_fs   = 0.0    # free strand correction
        
        # enthalpy
        Hclose = self.fe.find_branch_dHwt(i_h, j_h, branches, opt_iMBL)
        imotif.dH_dVlp = Hclose
        imotif.dHVpq   = Hpqh
        imotif.Hij     = Hpqh + Hclose
        
        self.fe.smap.glink[i_h][j_h].add_link(Link(imotif)) 
        
        if debug_assign_I:
            dG = Vpqh + Vclose
            dH = Hpqh + Hclose
            pI = branches[0][0]; qI = branches[0][1]
            print ("assign_I: ij_h(%2d,%2d)[%s][%s], pqI(%2d,%2d), dG(%8.2f), dH(%8.2f)" \
                % (i_h, j_h, ctp, btp, pI, qI, dG, dH))
            print (self.fe.smap.glink[i_h][j_h].lg[0].motif.show_Motif())
            
            if plan_stop_on_pass:
                print ("planned exit at end of assign_I");
                sys.exit(0)
            #
            
        #
        
    #
    
    
    def assign_J(self, i_h, j_h, branches, stype, vtype):
        
        # I think this one is not call (at least normally). In
        # particular, I noticed that it had a simple error but it was
        # basically ignored.
        debug_assign_J = True # False # 
        plan_stop_on_pass = False
        
        if debug_assign_J:
            print ("Enter assign_J{ij_h(%d,%d), branches = %s" % (i_h, j_h, branches))
        #
        
        # has no pairing at i_h, j_h
        opt_iMBL = False # T -> closed at (i_h,j_h), F -> open
        
        ctp = 'J'
        btp = "%s%s" % (stype, vtype)
        if not (vtype == 'a' or vtype == 'p'):
            btp = vtype
        #
        
        """@
        
        The main difference between an I-loop and a J-loop is that we
        don't have the cap closing the loop at (i_h, j_h). Hence, the
        function
        
            self.fe.find_branch_dGwt(i_h, j_h, branches)
        
        is not called and the dangle function is most certainly not
        used.
        
        """
        
        Vpqh, Hpqh = self.get_branchingFE(branches)
        # Vpqh = FE contribution from I/M branches or 0.0 for B 
        Vbranch_wt = self.fe.find_branch_dGwt(i_h, j_h, branches, opt_iMBL)
        
        
        # add any corrections for tetraloop, local loop
        # entropy or whatever here
        
        dG = Vpqh + Vbranch_wt + dGt
        
        jmotif = MBL(i_h, j_h, dG, ctp, btp, branches)
        jmotif.vtype   = vstem.vtype
        jmotif.dG_dVlp = Vbranch_wt # branch dangle + 5'3'-head + cle closing
        jmotif.dGVpq   = Vpqh       # FE from branches within (ph, qh)
        jmotif.dG_fs   = 0.0        # free strand correction
        
        # enthalpy
        jmotif.dHVpq   = Hpqh       # branch dH contributions without adjustments
        jmotif.dG_dVlp = self.fe.find_branch_dHwt(i_h, j_h, branches, opt_iMBL)
        jmotif.Hij     = Hpqh + jmotif.dG_dVlp  # free strand correction
        
        
        
        self.fe.smap.glink[i_h][j_h].add_link(Link(jmotif)) 
        
        
        if debug_assign_J:
            # pairing at i_h, j_h
            pJ = branches[0][0]; qJ = branches[0][1]
            print ("assign_J: ij_h(%2d,%2d)[%s][%s], pqJ(%2d,%2d), dG(%8.2f)" \
                % (i_h, j_h, ctp, btp, pJ, qJ, dG))
            print (self.fe.smap.glink[i_h][j_h].lg[0].motif.show_Motif())
            
            if plan_stop_on_pass:
                print ("planned exit at end of assign_J");
                sys.exit(0)
            #
            
        #
        
    #
    
    
    def assign_M(self, i_h, j_h, branches, stype, vtype):
        debug_assign_M = False # True # 
        plan_stop_on_pass = False
        if debug_assign_M:
            print ("Enter assign_M{ij_h(%d,%d), branches = %s" % (i_h, j_h, branches))
        #
        
        # has pairing at i_h, j_h
        opt_iMBL = True # T -> closed at (i_h,j_h), F -> open
        ctp = 'M'
        btp = "%s%s" % (stype, vtype)
        if not (vtype == 'a' or vtype == 'p'):
            btp = vtype
        #
        
        # xx!!!! 190707 !!!!xx needs some updating
        #dGh = self.fe.calc_dG(i_h, j_h, 5.0, self.fe.T)
        
        Vpqh, Hpqh = self.get_branchingFE(branches)
        # Vpqh = FE contribution from I/M branches or 0.0 for B 
        Vclose = self.fe.find_branch_dGwt(i_h, j_h, branches, opt_iMBL)
        if debug_assign_M:
            print ("assign_M: Vpqh(%8.2f), Vclose(%8.2f), branches:  %s" \
                % (Vpqh, Vclose, branches))
        #
        
        # add any corrections for tetraloop, local loop
        # entropy or whatever here
        
        dG = Vpqh + Vclose
        
        mmotif = MBL(i_h, j_h, dG, ctp, btp, branches)
        mmotif.vtype   = vtype
        mmotif.dG_dVlp = Vclose # branch dangle + 5'3'-head + cle closing
        mmotif.dGVpq   = Vpqh   # FE from branches within (ph, qh)
        mmotif.dG_fs   = 0.0    # free strand correction
        
        # enthalpy
        Hclose = self.fe.find_branch_dHwt(i_h, j_h, branches, opt_iMBL)
        if debug_assign_M:
            print ("assign_M: Hpqh(%8.2f), Hclose(%8.2f), branches:  %s" \
                % (Hpqh, Hclose, branches))
        #
        
        mmotif.dHVpq   = Hpqh 
        mmotif.dH_dVlp = Hclose
        mmotif.Hij     = Hpqh + Hclose
        
        self.fe.smap.glink[i_h][j_h].add_link(Link(mmotif))                
        if debug_assign_M:
            print ("free energy: ========")
            print ("dGVpq:      %8.2f" % (mmotif.dGVpq))
            print ("dG_dVlp:    %8.2f" % (mmotif.dG_dVlp))
            print ("--------    ----------")
            print ("Vij:        %8.2f" % (mmotif.Vij))
            print ("enthalpy:    ========")
            print ("dHVpq:      %8.2f" % (mmotif.dHVpq))
            print ("dH_dVlp:    %8.2f" % (mmotif.dH_dVlp))
            print ("--------    ----------")
            print ("Hij:        %8.2f" % (mmotif.Hij))
            print ("assign_M: ij_h(%2d,%2d)[%s][%s] -> dG(%8.2f)" \
                % (i_h, j_h, ctp, btp, dG), branches)
            
            for Vk in self.fe.smap.glink[i_h][j_h].lg[0].motif:
                print (Vk.show_Motif())
            #endfor
            
            if plan_stop_on_pass:
                print ("planned exit at end of assign_M");
                sys.exit(0)
            #

        #
        
    #
    
    
    
    def assign_P(self, i_h, j_h, branches, stype, vtype):
        debug_assign_P = True # False # 
        plan_stop_on_pass = False
        if debug_assign_P:
            print ("Enter assign_P{ij_h(%d,%d), branches = %s" % (i_h, j_h, branches))
        #
        
        # has no pairing at i_h, j_h
        opt_iMBL = False # T -> closed at (i_h,j_h), F -> open
        ctp = stype
        btp = "s%s" % (vtype)
        if not (vtype == 'a' or vtype == 'p'):
            btp = vtype
        #
        
        """@
        
        The main difference between iMBL and pMBL is that we don't
        have the cap closing the loop at (i_h, j_h). Hence, the
        function. This means opt_iMBL = False in this case. The result
        is that the function
        
            self.fe.find_branch_dGwt(i_h, j_h, branches, opt_iMBL)
        
        does not call cle_MloopEV() and only reports results from
        dG_trim_53pStems(). 
        
        """
        
        Vpqh, Hpqh = self.get_branchingFE(branches)
        # Vpqh = FE contribution from I/M branches or 0.0 for B
        
        Vbranch_wt = self.fe.find_branch_dGwt(i_h, j_h, branches, opt_iMBL)
        if debug_assign_P:
            print ("assign_P: Vpqh(%8.2f), Vbranch_wt(%8.2f), gbranches:  %s" \
                % (Vpqh, Vbranch_wt, branches))
        #
        
        # add any corrections for tetraloop, local loop
        # entropy or whatever here
        
        dG = Vpqh + Vbranch_wt + dGt
        
        pmotif = MBL(i_h, j_h, dG, ctp, btp, branches)
        pmotif.vtype   = vstem.vtype
        pmotif.dG_dVlp = Vbranch_wt # branch dangle + 5'3'-head + cle closing
        pmotif.dGVpq   = Vpqh       # FE from branches within (ph, qh)
        pmotif.dG_fs   = 0.0        # free strand correction
        
        # enthalpy
        pmotif.dHVpq    = Hpqh
        pmotif.dH_dVlp  = self.fe.dH_trim_53pStems(i, j,      ## 5' to 3' of [i/p]MBL
                                                  branches,  ## (class MBPptr)
                                                  opt_iMBL)  ## T -> iMBL, F -> pMBL
        pmotif.Hij     = Hpqh + mmotif.dH_dVlp
        
        
        self.fe.smap.glink[i_h][j_h].add_link(Link(pmotif))                
        
        if debug_assign_P:
            print ("assign_P: ij_h(%2d,%2d)[%s][%s] -> dG(%8.2f)" \
                % (i_h, j_h, ctp, btp, dG), branches)
            
            for Vk in self.fe.smap.glink[i_h][j_h].lg[0].motif:
                print (Vk.show_Motif())
            #endfor
            
            if plan_stop_on_pass:
                print ("planned exit at end of assign_P");
                sys.exit(0)
            #
        #
        
    #
    
    
    def assign_Stem(self, i_t, j_t, i_h, j_h, vstem, stype, vtype, fr_zones = []):
        """@
        
        Constructs the newer versions of the stem design and handles
        construction of pseudoknots with consideration about the
        fr_zonesing of the relevant regions -- something that was missed
        before.
        
        In this step, it is assumed that the free energy calculations
        already have been done in the depth first process that is
        used. 
        
        "fr_zones" is the regions for extracting branching for
        pseudoknots; the so-called "free zones"
        
        * core pseudoknot free zone definitions for ap stems
        (ir1,jr1)::(ph1,qh1) and (pL,qL)::(pl,ql) 
        fr_zones -> for root1   = [(ph1, pL), (pl, qh1)]
                    for linkage =  (jr1, ql)
        
        * extended pseudoknots free zone definitions for ap stems
        (ir1,jr1)::(ph,qh), (ir2,jr2)::(ph2,qh2) and (pL,qL)::(pl,ql) 
        fr_zones -> for root1   = [(ph1, pL), (pl, qh1)]
                    for root2   = [(ph2, qL), (ql, qh2)]
                    for linkage =  (jr1, ir2)
        
        """
        
        debug_assign_Stem = False # True # 
        plan_stop_on_pass = False
        if debug_assign_Stem:
            print ("Enter assign_Stem{(%2d,%2d):(%2d,%2d)" % (i_t, j_t, i_h, j_h))
            print ("                   vstem: %s" % (vstem))
            print ("                   stype[%s], vtype[%s]" % (stype, vtype))
            print ("                fr_zones: %s}" % (fr_zones))
        #
        
        # First, we search the head of the stem for branching
        
        gbranches  = []
        ngbranches = 0
        close_at_ijh = True
        
        
        if len(fr_zones) > 0:
            # if the stem is associated with a pseudoknot, then it
            # will have fr_zones regions where the linkage stem
            # attaches. These need to be specified to ensure that the
            # proper structure is retained.  If it is just a regular
            # stem, it will not have a fr_zones.
            for fr_zones_k in fr_zones:
                im = fr_zones_k[0]; jm = fr_zones_k[1]
                gbranches += self.get_genbranches(im, jm, close_at_ijh)
            #
            
            ngbanches = len(gbranches)
            
        else:
            # stem is just part of a regular generic stem.
            gbranches = self.get_genbranches(i_h, j_h, close_at_ijh)
            ngbranches = len(gbranches)
        #
        
        # Now find the free energy of the loop region or the 
        
        Vpqh, Hpqh = self.get_branchingFE(gbranches)
        # Vpqh = FE contribution from I/M branches or 0.0 for B 
        Vclose = self.fe.find_branch_dGwt(i_h, j_h, gbranches)
        Hclose = self.fe.find_branch_dHwt(i_h, j_h, gbranches)
        if debug_assign_Stem:
            print ("Vpqh(%8.2f), Vclose(%8.2f), gbranches:  %s" \
                % (Vpqh, Vclose, gbranches))
            print ("Hpqh(%8.2f), Hclose(%8.2f), gbranches:  %s" \
                % (Hpqh, Hclose, gbranches))
        #endif
        
        
        
        ctp = 'S'
        btp = "%s%s" % (stype, vstem.vtype)
        
        dG = Vpqh + Vclose
        dH = Hpqh + Hclose
        
        # total free energy and enthalpy at the head including the
        # closing point at (i_h, j_h)
        
        n = vstem.slen - 2
        if debug_assign_Stem:
            print ("ij_h(%2d,%2d), slen = %2d, dG = %8.2f" % (i_h, j_h, n+2, dG))
            print ("                        dH = %8.2f" % (dH))
            # sys.exit(0)
        #
        
        if ngbranches == 0:
            hmotif = XLoop(i_h, j_h, Vclose, 'B', btp)
            hmotif.vtype = vstem.vtype
            
            hmotif.dH_dVlp = Hclose
            hmotif.dHVpq   = Hpqh
            hmotif.Hij     = dH # Hpqh + Hclose
            self.fe.smap.glink[i_h][j_h].add_link(Link(hmotif))
        elif ngbranches == 1:
            imotif = MBL(i_h, j_h, dG, 'I', btp, gbranches)
            imotif.vtype   = vstem.vtype
            
            imotif.dG_dVlp = Vclose # branch dangle + 5'3'-head + cle closing
            imotif.dGVpq   = Vpqh   # FE from branches within (ph, qh)
            
            imotif.dH_dVlp = Hclose
            imotif.dHVpq   = Hpqh
            imotif.Hij     = dH # Hpqh + Hclose
            
            self.fe.smap.glink[i_h][j_h].add_link(Link(imotif))
            
        elif ngbranches > 1:
            mmotif = MBL(i_h, j_h, dG, 'M', btp, gbranches)
            mmotif.vtype   = vstem.vtype
            
            mmotif.dG_dVlp = Vclose # branch dangle + 5'3'-head + cle closing
            mmotif.dGVpq   = Vpqh   # FE from branches within (ph, qh)
            
            mmotif.dH_dVlp = Hclose
            mmotif.dHVpq   = Hpqh
            mmotif.Hij     = dH # Hpqh + Hclose
            
            self.fe.smap.glink[i_h][j_h].add_link(Link(mmotif))
        else:
            print ("ERROR(assign_Stem): something seriously wrong ngbranches = %d" \
                % ngbranches)
            sys.exit(1)
        #
        
        # now that we have finished the head, go to the rest of the
        # stem till the tail (if the stem is longer than 1 Pair).
        
        
        # build a map entry for every Pair in the stem
        stm = [vstem.stem[n + 1]]
        dGp = [Vclose]; dHp = [Hclose]
        for k in range(n, -1, -1): # bpk in vstem.stem:
            
            bpk   = vstem.stem[k]
            bpkp1 = vstem.stem[k+1]
            iss  = bpk.i;   jss  = bpk.j
            issf = bpkp1.i; jssf = bpkp1.j
            
            stm += [bpk]
            num_bp = len(stm)
            
            # set the local free energy
            xi_stm, L_eff, dGloc_stm, dGloc_fs \
                = self.fe.calc_ddGlocal_Stem(iss, jss, # stem tail (needed for checking)
                                             i_h, j_h, # stem head (needed for checking)
                                             num_bp,   # the number of bp in the stem
                                             False)    # option to debug
            
            ddGl = (dGloc_stm - dGloc_fs) # local free energy Gstm - Gfs
            
            # 190726: generalized description of pairing FE  vvvvvvvvvvvvvvvvvvvv
            if not (issf - iss == 1 and jss - jssf == 1):
                # add connected I-loop between the effective stems
                
                dGp_k = self.fe.cle_IloopE(iss,  jss,   # tail of I-loop 
                                           issf, jssf,  # head of I-loop 
                                           xi_stm)      # [nt] Kuhn length
                
                dHp_k = self.fe.dH_IloopE(iss,  jss,   # tail of I-loop 
                                          issf, jssf,  # head of I-loop 
                                          xi_stm)      # [nt] Kuhn length
                
                if debug_assign_Stem:
                    s = "Iloop(k=%2d): ijss(%2d,%2d) -> ijssf(%2d,%2d)[%s]" \
                        % (k, iss, jss, issf, jssf, self.fe.btype[iss][jss].f2nt)
                    s += "xi_stm(%8.2f), L_eff(%8.2f), dGloc_stm(%8.2f), dGloc_fs(%8.2f)\n"\
                         % (xi_stm, L_eff, dGloc_stm, dGloc_fs)
                    s += "ddGl(%8.2f), dGp_k(%8.2f), dHp_k(%8.2f)" % (ddGl, dGp_k, dHp_k)
                    print (s)
                    #sys.exit(0)
                #
                
            else:
                
                dGp_k  = self.fe.calc_dGbp(iss, jss, xi_stm)
                dHp_k  = self.fe.calc_dHbp(iss, jss, xi_stm)
                
                if debug_assign_Stem:
                    s = "aStem(k=%2d): ijss(%2d,%2d)::ijssf(%2d,%2d)[%s]" \
                        % (k, iss, jss, issf, jssf, self.fe.btype[iss][jss].f2nt)
                    s += "xi_stm(%8.2f), L_eff(%8.2f), dGloc_stm(%8.2f), dGloc_fs(%8.2f)\n"\
                         % (xi_stm, L_eff, dGloc_stm, dGloc_fs)
                    s += "ddGl(%8.2f), dGp_k(%8.2f), dHp_k(%8.2f)" % (ddGl, dGp_k, dHp_k)
                    print (s)
                    #sys.exit(0)
                #endif
                
            #
            
            # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            dG  += dGp_k;   dH  += dHp_k
            dGp += [dGp_k]; dHp += [dHp_k]
            
            
            smotif = Stem(list(reversed(stm)))
            smotif.xi         = xi_stm    # Kuhn length of stem
            smotif.dG_l       = dGloc_stm # local FE from stem formation
            smotif.dG_fs      = dGloc_fs  # equivalent free strand FE
            # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
            smotif.Vij        = dG + ddGl # total FE (including local FE)
            # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            smotif.dGp \
                = list(reversed(dGp))     # individual pair FE contibutions
            smotif.dGloop     = Vclose    # H/I/M loop closing FE weight at (i_h,j_h)
            smotif.Vpqh       = Vpqh      # FE contribution from adjoining X/J/P
            smotif.branching  \
                = deepcopy(gbranches)     # branches closing the stem
            
            
            smotif.dHp \
                = list(reversed(dHp))     # individual pair FE contibutions
            smotif.Hpqh       = Hpqh      # FE contribution from adjoining X/J/P
            smotif.Hij        = dH
            
            #print (len(self.fe.smap.glink))
            #print (iss, jss)
            
            self.fe.smap.glink[iss][jss].add_link(Link(smotif))
            
            if debug_assign_Stem:
                print ("a%s stem: ijss(%3d,%3d):ijh(%3d,%3d), stemlen = %2d: " \
                    % (vtype, iss, jss, i_h, j_h, len(stm)))
                print ("dG(%8.2f) + ddGl(%8.2f) = dGtot(%8.2f) vs  Vij(%8.2f)" \
                    % (dG,
                       ddGl,
                       dG + ddGl,
                       self.fe.smap.glink[iss][jss].lg[0].motif.get_Vij()))
                print ("                                   dH(%8.2f) vs  Hij(%8.2f)" \
                    % (dH, self.fe.smap.glink[iss][jss].lg[0].motif.get_Hij()))
                #print ("         %s" % (list(reversed(stm))))
                #sys.exit(0)
            #
            
        #endfor k
        
        if debug_assign_Stem:
            print ("summary(assign_Stem){(%2d,%2d):(%2d,%2d)}:" % (i_t, j_t, i_h, j_h))
            for bpk in vstem.stem:
                iss = bpk.i; jss = bpk.j
                
                for ijSt in self.fe.smap.glink[iss][jss].lg[0].motif:
                    print (ijSt.show_Motif())
                #endfor
                
            #endfor
             
            if plan_stop_on_pass:
                print ("planned exit at end of assign_S");
                sys.exit(0)
            #
            
        #endif
    #
    
    
    
    def write_StemMotif(self, vstem, stype = 's', fr_zones = []):
        """@
        
        vstem:      object of class Stem
        
                    label     meaning
        stype:       's'    (strandard) 
                     'l'    (linkage)
        
        fr_zones:   free strand zones used on PseudoKnots
        
        """
        debug_write_StemMotif = False # True # 
        
        i_t = vstem.it; j_t = vstem.jt
        i_h = vstem.ih; j_h = vstem.jh
        vtype = vstem.vtype
        if debug_write_StemMotif:
            print ("Enter write_StemMotif{(%2d,%2d):(%2d,%2d)}[%s]" \
                % (i_t, j_t, i_h, j_h, stype))
            print ("      stemlen:     %d" % vstem.slen)
        #
        
        """@
        
        fr_zones = free strand zones
        
        Presently, ssbranches is not needed directly, but in the
        future with the RNA engine, it is more important to have this
        kind of information for determining the free energy.
        
        """
        
        
        if vstem.slen > 1:
            if debug_write_StemMotif:
                print ("write_StemMotif(%2d, %2d) -> build Stem" % (i_t, j_t))
            #
            
            # build the stem
            self.assign_Stem(i_t, j_t, i_h, j_h, vstem, stype, vtype, fr_zones = [])
            
        else:
            
            """!!!! may happen for chromatin !!!!"""
            
            # We start at the head of the stem with an isolated pair
            
            close_at_ijh = True
            gbranches = self.get_genbranches(i_h, j_h, close_at_ijh)
            ngbranches = len(gbranches)
            
            if debug_write_StemMotif:
                print ("ngbranches(%d): %s" % (ngbranches, gbranches))
            #
            
            if ngbranches > 0:
                
                if debug_write_StemMotif:
                    print ("found isolated bp at (%2d, %2d),  gbranches:  %s" \
                        % (i_h, j_h, gbranches))
                #
                
                if ngbranches == 1:
                    # j-type 
                    self.assign_I(i_h, j_h, gbranches, stype, vtype)
                else:
                    # p-type
                    self.assign_M(i_h, j_h, gbranches, stype, vtype)
                #
                
                # self.fe.smap.mergeSortLinks(self.fe.smap.glink[i_h][j_h].lg)
                
                """@
                
                190728: Note, this comment is what I claimed in a
                previous structuring of this program. I find the
                reasoning rather strange, I think the data should be
                sorted. It also looks like the general assembly is
                correct in examining the outputs. However, maybe I am
                just confused.
                
                In general, we should sort this FE so that the optimal
                structure ends up on top. However, presently, the
                final structure is what we ask for, even if it is not
                the optimal on at this position.
                
                """     
            
                
            elif ngbranches == 0:
                if debug_write_StemMotif:
                    print ("write_StemMotif(%2d, %2d) -> build B-loop" % (i_t, j_t))
                #
                
                self.assign_B(i_h, j_h, stype, vtype)
                
            else:
                print ("ERROR: ??? don't understand what to make of this case!")
                print ("       ij_h(%2d,%2d)[%s]" % (i_h, j_h, vtype))
                sys.exit(1)
            #
            
        #
        
    #
    
    # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
    # VVVVVVVVVVVVVVVVVV   Branch building tools   VVVVVVVVVVVVVVVVVV
    # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
    """@
    
    On get_branchingFE obtaining the free energy from the branchings
    comprising M-, P-, I-, J-loops
    
    Certainly for chromatin, there is _no_ difference in treatment
    between I-/M-loops or J-/P-loops. I think this should also be
    viewed the same way for RNA, but momentarily, I may still
    entertain the idea that there is maybe a higher ordering cost for
    a single branch because the internal loop can form a
    semi-contiguous connection if the I-loop is relatively short and
    symmetric
    
    """
    
    
    def get_branchingFE(self, branches): # both type I and M
        dH_V = 0.0
        dG_V = 0.0
        dG_local = 0.0 # presently, this is always zero
        for Vk in branches:
            i_Vk = Vk[0]; j_Vk = Vk[1]
            dG_V += self.fe.smap.glink[i_Vk][j_Vk].lg[0].Vij
            dH_V += self.fe.smap.glink[i_Vk][j_Vk].lg[0].motif.get_Hij()
            
            # add any corrections for the cost of forming the branch
            # (or branches) between (i_h,j_h) and [(pV0,qV0),
            # (pV1,qV1) ... ]; e.g., purine corrections, loop size,
            # etc.
        #
        
        
        # add the FE contributions together 
        dGloop = dG_local + dG_V
        dHloop = dH_V
        return dGloop, dHloop
    #
    
    
    
    
    def get_ssbranches(self, ib, jb, closed = True):
        """I think this is probably less used or maybe not used anymore"""
        debug_get_ssbranches = False # True # 
        if debug_get_ssbranches:
            print ("Enter get_ssbranches(), closed = %s" % closed)
        #
        
        branches = []
        for vstem in self.ssStems:
            gstemtype = type(vstem).__name__
            # print (gstemtype)
            if not gstemtype == "Stem": # this __must__ be class Stem
                print ("ERROR: ssStem contains objects other than type Stem")
                print ("       %s: ijss = (%2d,%2d)" % (gstemtype, vstem.it, vstem.jt))
                sys.exit(1)
            #
            
            iV = vstem.it; jV = vstem.jt
            if debug_get_ssbranches:
                print ("ib(%2d) <  iV(%2d) < jV(%2d) <  jb(%d)" % (ib, iV, jV, jb))
            #
            
            if closed:
                # BC: the branches are ib < iV < jV < jb, they exclude (ib,jb)
                if ib < iV and jV < jb:
                    branches += [(iV, jV)]
                #
            else:
                # BC: the branches are ib <= iV < jV <= jb, they include (ib,jb)
                if ib <= iV and jV <= jb:
                    branches += [(iV, jV)]
                #
            #
            
        #endfor
        
        return branches
    #
    
    
    def get_genbranches(self, ib, jb, closed = True):
        debug_get_genbranches = False # True # 
        if debug_get_genbranches:
            print ("Enter get_genbranches(), closed = %s" % closed)
        #
        
        branches = []
        for vstem in self.genStems:
            gstemtype = type(vstem).__name__
            
            if vstem.name == "Base":
                continue
            #
            
            # print (gstemtype)
            if not (gstemtype == "Stem" or gstemtype == "PseudoKnot"):
                print ("ERROR: genStem contains objects other than type Stem and PseudoKnot")
                print ("       %s: ijV = (%2d,%2d)" % (gstemtype, vstem.it, vstem.jt))
                sys.exit(1)
            #
            
            iV = vstem.it; jV = vstem.jt
            if debug_get_genbranches:
                print ("ib(%2d) <= iV(%2d) < jV(%2d) <= jb(%d)" % (ib, iV, jV, jb))
            #
            
            if closed:
                # BC: the branches are ib < iV < jV < jb, they exclude (ib,jb)
                if ib < iV and jV < jb:
                    branches += [(iV, jV)]
                #
            else:
                # BC: the branches are ib <= iV < jV <= jb, they include (ib,jb)
                if ib <= iV and jV <= jb:
                    branches += [(iV, jV)]
                #
            #

            if len(branches) > 1:
                if debug_get_genbranches:
                    print ("get_genbranches: branches = %s" % branches)
                #
                
                k = 0
                while k < len(branches) - 1:
                    ik1 = branches[k  ][0]; jk1 = branches[k  ][1]
                    ik2 = branches[k+1][0]; jk2 = branches[k+1][1]
                    
                    if ik1 <= ik2 and jk2 <= jk1:
                        if debug_get_genbranches:
                            print ("get_genbranches: branches = %s" % branches)
                            print ("ik1(%2d) <= ik2(%2d) < jk2(%2d) <= jk1(%2d)" \
                                % (ik1, ik2, jk2, jk1))
                            print ("deleting ijk2(%2d,%2d)" % (ik2, jk2))
                        #
                        
                        del branches[k+1]
                    else:
                        k += 1
                    #
                    
                #
                
            #
            
        #
        return branches
    #
    
    # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    # AAAAAAAAAAAAAAAAAA   Branch building tools   AAAAAAAAAAAAAAAAAA
    # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    
    
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
        print (thr.disp_lnode())
    #
    print (dt.makeLThreadDotBracket_VARNA(lt[ndx], 0))
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
    print ("Final List : ", aList)
#


def test2(cl):
    from SettingsPacket  import InputSettings
    from Vienna2TreeNode import Vienna2TreeNode
    # The only way to get this to work is to load this module _within_
    # this specific method. When I do that, it will run without
    # causing circular definitions of TreeNode and Vienna2TreeNode. So
    # it accomplishes a certain amound of independence as though this
    # were a separate program nowhere located within TreeNode.py. 

    iSetUp = InputSettings("RNA")
    iSetUp.set_source("test2 *TreeNode)")
    iSetUp.set_program("TreeNode.py")
    
    #iSetUp.useTurner  = True
    #iSetUp.useViS     = False
    #iSetUp.usegMatrix = args.usegMatrix
    #iSetUp.gMflnm     = "none"
    iSetUp.set_FEparamData() # default Turner
    
    # set up the actual FE Data according to the settings
    iSetUp.set_ViennaParams()
    
    vs = Vstruct()
    vs.set_system(iSetUp)
    
    
    
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
    print ("LThread notation: ")
    tf.disp_lt()
    # print ("planned exit"); sys.exit(0);
    
    # print (dt.makeLThreadDotBracket_VARNA(ds.lt, 0))
    vf = LThread2Vienna()
    vf.lt2vs(tf.lt)
    vf.set_vstr(v2t.vstr)
    vf.set_vseq(v2t.vseq)
    
    print ("structure sequence: ")
    print (vf.vstr)
    print ("vsBPlist: ")
    for bpk in vf.vsBPlist:
        print (bpk.disp_Pair())
    #endfor
    
    print ("vsPKlist: ")
    for bpk in vf.vsPKlist:
        print (bpk.disp_Pair())
    #endfor
    
    print ("vsMPlist: ")
    for bpk in vf.vsMPlist:
        print (bpk.disp_Pair())
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

