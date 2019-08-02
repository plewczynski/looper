#!/usr/bin/env python

"""@@@

Main Module:     Motif.py (adapted to chromatin starting from RNA/Motif-bkp160630.py)

Objects:         Motif 
                 Link 
                 LGroup 
                 ndx 
                 Map
                 --- motif class objects
                 AST 
                 APT 
                 Stem 
                 PseudoKnot 
                 ArchMP 
                 MultiPair

Author:          Wayne Dawson
creation date:   170126
last update:     190718 
version:         0.1

Purpose:

This module is used primarily by "chreval.py". However, as with many
packages that are building up, the package has become complicated
enough, that the objects in this module are used by more than one
program.

The remainder of this module deals with building tree diagrams and
will eventually become a unified format for Motif. 


Comments:

190725

   #################################################################
   Be aware that Motif is under major revision. As a result, it is
   presently in a somewhat odd state of disarray. I am presently in
   the process of changing Motif to an object representation (with the
   classes Stem, PseudoKnot, MBL, XLoop and MultiPair), adding handles
   to access the objects and migrating the evalutation and
   construction into the rest of the program.  Presently, not all of
   the routines have been upgraded to handle and build this
   representation. Hence, most of the processes are still using the
   old methods. I expect this to change in the next few months, but
   unfortunately at this time of release, that step has not been
   completed yet.

   The program is also being tested with RNA structure prediction in a
   similar form as chromatin, and these checks will be used to help in
   the migration process.
   #################################################################

older comments

This was set up as a separate object for testing. I'm gradually
pushing chromatin in the same direction as my RNA program with real
objects, but presently, the are not being called yet; mostly because I
had to finalize the structure of thet class PseudoKnot in my RNA
directory.

For Vienna formatted (base) pair data, all that (Pascal-like)
information of thread is stripped off and the final structure computed
from the (base) pairs only. Though maybe there is something important
in all this, with the exception of the free energy (and the
Pascal-like notation), basically Pair and LNode are identical. I am
still thinking that it would be better to somehow organize this so
that Vienna and Thread are merged (so we don't have a feeble
half-notation). I haven't made up my mind and I still hesitate because
there is a little bit of information that thread carries that Pair
does not. Anyway, both are too incomplete in of themselves to be used
too extensively.

190403 found at least a reasonable strategy that appears to
work. Originally, I thought that the visit approach was best, but I
could not get that to work. On the other than, the form of NewMotif
appears to be working.

Oddly, the result of these objects in NewMotif were an outgrowth of my
effort to take Vienna data and convert it to a representation
recognized by LThread. However, I eventually recognized that LThread
is far too simple to be useful for building a general transformation
code that is recognizable by Vienna and Motif simultaneously. I guess
it is not impossible to make such a script, but the result is is
cumbersome and would still require the sophisticated tree making tools
in Threads to make it work.

"""

import sys

from ChrConstants import xi # [beads, bps, aa, nt]
from ChrConstants import xi_fs # useless but required by some classes anyway

from Pair import Pair            # vital
from Pair import SortPair        # vital
from Pair import vsPair2list     # not used
from Pair import BranchList2list # not used
from Pair import rlist2vsPair    # not used
from Pair import find_this_ij    # not used
from Pair import is_ap_stem      # not used
from Pair import is_pp_stem      # not used


#####################################################################
######################  BASIC LINKAGE OBJECTS  ######################
#####################################################################


# #################################################
# #############  LINKAGE DEFINITIONS  #############
# #################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

"""@

CBL = chromatin binding linkage
ctp = connection type, established at position (i,j). Let

(i,j) current position under analysis (may or may not have a CBL)
(i',j') a prior CBL

###############  connection definitions (ctp and btp):  ###############

######### CONNECTION TYPE (ctp):

'B', (i,j) -> [(i,j)] => connection type is a bond ('B') at (i,j) or
the terminal connection at the end of an antiparallel loop.

'b', (i,j) -> [(i,j)] => connection type is the terminal connection at
the end of an antiparallel loop at (i,j).

'S', (i,j) -> [(i+1,j-1), (i+2,j-2) ... ] or [(i,j),(i+1,j+1),
(i+2,j+2) ... ] or [(i,j),(i-1,j-1), (i-2,j-2) ... ]: "Stem".  this is
similar to the RNA stem I used before; however, in the case of
chromatin, it is important to remember that both antiparallel and
parallel directions are allowed.

=========

'I', (i,j) -> [(i',j')] => connection type is a "I-loop" ('I'): closed
at (i,j) with the next adjacent cross link at (i',j'), where i' > i +
1 and/or j' < j -1. Typically, this structure should be more favorable
than 'B' if it exists, but it is not necessarily so.

'J', -> [(i',j')] => connection type is a "J-loop" ('J'): _open_ at
(i,j) with the next adjacent cross link at (i',j'), where i' > i + 1
and/or j' < j -1. A J-loop can be more favorable than the I-loop when
the binding free energy at (i,j) is positive.

------------------------------------------------------------------------
Note that both 'I' and 'J' could be defined as a single-branch loop
'V', where 'I' has the same base as (i,j), and 'J' has a different
base from (i,j).
------------------------------------------------------------------------

=========

'K', (i, j) -> join = [(ib,jb)], pk = [(ik1,jk1), (ik2,jk2)...]  =>
connection type pseudoKnot (presently all combinations of H
type). There is a bit of a problem with chromatin that I did not
encounter in RNA because the chain can be either parallel or
antiparallel. I will only allow one solution presently. This is
because I ASSUME that the chain will favor one or the other direction
but not both for a good solid pseudoknot. This could be problematical,
but, for the moment, I'll take the risk.

'l': a pseudoknot in the "linkage" part of the pseudoknot. There is
nothing in principle that is different between the linkage and the
main structure. Both are essentially singletons. However, it is
convenient to notate it differently for human readability and due to
the fact that linkages are generally more restricted on their
interaction with the rest of the structure.

=========

'M', (i,j) -> [(i1,j1),(i2,j2)] => connection type is a "M-loop":
closed at (i,j) with branches located at [(i1,j1), (i2,j2)
... (in,jn)], subject to the condition that i < i1, j1 < i2, i2 < j3,
.... jn < j.  M-loop is for "multibranch loop" and this form is
typically more favorable than an 'I' or a 'B', but it depends on the
situation.

['P', [(i1,j1),(i2,j2)]] => connection type is a "P-loop": _open_ at
(i,j) with branches located at [(i1,j1), (i2,j2) ... (in,jn)],
subject to the condition that i < i1, j1 < i2, i2 < j3, .... jn < j.
P-loop is for "principal loop". Typically, the M-loop should be more
favorable, but it could be favored if the binding free energy 'B' at
(i,j) is positive.

------------------------------------------------------------------------
Note that both 'M' and 'P' could be defined as a multi-branch loop
'V', where 'M' has the same base as (i,j), and 'P' has a different
base from (i,j).
------------------------------------------------------------------------


'S': stem.  It is important to specify if the
collection of stems is parallel or antiparallel oriented. Therefore,
we must add a qualifier to this notation; i.e., for structures
containing the key 'S', we should specify 'sa' (for singleton
antiparallel) and 'sp' (for singleton parallel).


=========

['W', [(i1, j1, i2, j2 ....)]] => connection type CTCF-island. This
is unique to CTCF islands. In this type of interaction, all the
items on the list are in effect associated. This requires a
different notation from the singleton that presently I assume can be
written as

    ....{......|.........|..............}.....
        i1     j1        i2             j2

where the meaning is that all these regions interact with each
other.



# ######### BOND TYPE (btp):

'-': none, used in the case of pointers because it represents neither
a singleton nor a PET cluster.

'bgn'/'end': markers for the beginning and ending of stems ('S'),
pseudoknots ('K') and CTCF islands ('W')

'c': This is a CTCF cluster that we define as "convergent"; i.e.,
having the "... > .... < ..." type of interaction that _probably_
results in an _antiparallel_ chain configuration. This configuration
results in a very strongly binding interaction and is most commently
encountered.

'l': This is a CTCF cluster oriented in the tandem left configuration;
; i.e., having the "... < .... < ..." type of interaction that
_probably_ results in an _parallel_ chain configuration. parallel
configuration. This configuration results in a far less strongly
binding interaction (about 1/3) but is also often observed in the
experimental data.

'r': This is a CTCF cluster oriented in the tandem right
configuration; ; i.e., having the "... > .... > ..." type of
interaction that _probably_ results in an _parallel_ chain
configuration. This configuration results in a far less strongly
binding interaction (about 1/3) but is also often observed in the
experimental data.

'rp2':   RNApolII

's': This is the "singleton" bond type. Previously (before 160721), I
called this type 'n', for "normal". What is normal and is singleton
any better? Presently, I don't know, but we read in chromatin data of
weak structural ininteractions between local regions of the structure
and we predict some loops structures resulting from that. These loop
structures are assumed to have the chromatin chains oriented in an
unspecifice orientation (parallel loop or antiparallel loop direction)
and the interactions are presumed to be weak. Note that 's' could also
refer to divergent CTCF interactions; i.e., having the "< ... >" type
of interaction that results in weakly interacting antiparallel chain
configurations similar to the convergent CTCF structures only with far
less binding free energy.

'sa', 'sp': Stem direction specifier.  It is important to specify if the
collection of pairs is parallel or antiparallel oriented. Therefore,
we must add a qualifier to this notation; i.e., for structures
containing the key 'S', we should specify 'sa' (for singleton
antiparallel) and 'sp' (for singleton parallel).


't' (i.e., "tandem" or "parallel"): This refers to tandem
configurations wherein the particular orientation of the CTCF binding
sites is unknown. In such cases, the program will automatically assign
this label, since the binding energy of 'l' and 'r' are basically
identical at the current level of understanding.

'wyspa': CTCF island.


These additional "bond types" permit us to distinguish between
direction and type of interaction: i.e., RNApolII, CTCF "convergent",
CTCF "tandem".

# #########

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#################################################
#############  LINKAGE DEFINITIONS  #############
#################################################

160715wkd: I added this additional object because the program will
have to be expanded to handle convergent and tandem structures that
are coupled, and also because some structures will turn out to be
degenerate, so the program should be able to classify multiple
structures with the same energy within the same class. However, I
need to process the structures properly. Therefore, Link contains a
list of objects that have the same free energy Vij and Motif()
represents the objects contained in Link().

==================================================================
180815: I think we can use dictionary to identify objects in
motif. I is true that the handshakes need to be identical, but they
can be accessed with the proper interface.

motifDict = { 'B' : bond(), 'P' : MBL(), 'S' : Stem() ... }

lg[0].motif[0].get_branches()

get_branches(self):
  return motifDict['M'].get_branches()
# ==================================================================

"""
class ParseMotif(object):
    # not sure if this will be useful, but maybe this is a transition tool
    def __init__(self):
        self.debug = True
    #
    def get_base(self, mtf):
        bonds = []
        if mtf.ctp == 'P' or mtf.ctp == 'J':
            # J and P do not have a bond at the base
            bonds = mtf.branch
        elif mtf.ctp == 'K' or mtf.ctp == 'R':
            bonds = [(mtf.i, mtf.j)]
        elif mtf.ctp == 'S':
            bonds = [(mtf.i, mtf.j)]
            
        else:
            # B, I, K, and M have a bond at the base, and W is
            # effectively a cluster of bonds, on of which outlines the
            # entire structure
            bonds = [mtf.base]
        #
        return bonds
    #
    def get_branches(self, mtf):
        jj = []
        for jjk in mtf.branching:
            jj += [jjk]
        return jj
    #

#   


class Motif(object):
    def __init__(self, i, j, Vij, ctp, btp, branching, pk = [], wyspa = []):
        self.Ob        = None    # <<<=== transition to objects Stem, MBL, XLoop
        self.Vij       = Vij     # the free energy
        self.ctp       = ctp     # connection type 
        self.btp       = btp     # bond type
        self.branching = branching    # what linkage is pointed to
        self.tlen      = 0.0     # for Stem
        self.xi        = xi      # default Kuhn length
        if len(self.branching) == 0:
            if ctp == 'B':
                self.branching = [(i,j)]
            else:
                print "ERROR(Motif(%d,%d)): unrecognized linkage for type %s." \
                    % (i, j, ctp)
                sys.exit(1)
            #
        #
        self.base = (i,j)   # the base reference point
        self.pk    = pk
        # contains a list of ordered pairs for the linkage stem(s),
        # e.g., [( 7,20), ( 8,19), ( 9,18), (10,16)]
        self.wyspa = wyspa  # contains a list; e.g., [(38, 43), (40, 43)]
        #
    #
    
    def get_base(self):
        tp = self.Ob
        if tp == None:
            bonds = []
            if self.ctp == 'P' or self.ctp == 'J':
                # J and P do not have a bond at the base
                bonds = self.branching
            else:
                # B, I, K, and M have a bond at the base, and W is
                # effectively a cluster of bonds, on of which outlines the
                # entire structure
                bonds = [self.base]
            return bonds
        #
        elif tp.name == "Stem" or tp.name == "MBL" or tp.name == "XLoop":
            return self.Ob.get_base()
        else:
            print "Motif error(get base)", self.base
            print type(self.Ob)
            print self.Ob
            sys.exit(0)
        #
    
    #
    
    def get_branches(self):
        tp = self.Ob
        if tp == None:
            jj = []
            if not (self.ctp == 'B' or self.ctp == 'X'): 
                for jjk in self.branching:
                    jj += [jjk]
                #
                
            #
            
            return jj
        
        elif tp.name == "Stem" or tp.name == "MBL" or tp.name == "XLoop":
            return self.Ob.get_branches()
        
        else:
            print "Motif error(get_branches): ", self.base
            print type(self.Ob)
            print self.Ob
            sys.exit(0)
        #
        
    #

    def get_stem(self):
        if self.Ob == None:
            print "ERROR: Stem ijt(%d,%d) not defined!" % (self.base[0], self.base[1])
            print self.show_Motif()
            sys.exit(1)
        else:
            return self.Ob.get_stem()
        #
    #
    
    def get_Vij(self):
        tp = self.Ob
        if tp == None:
            """
            print "Vij is not set for %s, ctp(%s), btp(%s)" \
                % (self.base, self.ctp, self.btp)
            sys.exit(1)
            """
            return self.Vij
        elif (tp.name == "Stem" or tp.name == "MBL" or tp.name == "XLoop"):
            return self.Ob.Vij
        else:
            print "Motif error(get_Vij)", self.base
            print type(self.Ob)
            print self.Ob
            sys.exit(0)
        #
    #
    
    def get_ctp(self):
        tp = self.Ob
        if tp == None:
            """
            print "ctp is not set for %s, ctp(%s), btp(%s)" \
                % (self.base, self.ctp, self.btp)
            sys.exit(1)
            """
            return self.ctp
        elif (tp.name == "Stem" or tp.name == "MBL" or tp.name == "XLoop"):
            return self.Ob.ctp
        else:
            print "Motif error(get_ctp)", self.base
            print type(self.Ob)
            print self.Ob
            sys.exit(0)
        #
    #
    def get_btp(self):
        tp = self.Ob
        if tp == None:
            """
            print "btp is not set for %s, ctp(%s), btp(%s)" \
                % (self.base, self.ctp, self.btp)
            sys.exit(1)
            """
            return self.btp
        elif (tp.name == "Stem" or tp.name == "MBL" or tp.name == "XLoop"):
            return self.Ob.btp
        else:
            print "Motif error(get_btp)", self.base
            print type(self.Ob)
            sys.exit(0)
        #
    #
    
    
    def get_pks(self):
        return self.pk # pk contains a list of links (pairs)
    #
    
    def get_wyspa(self):
        return self.wyspa # wyspa contains a list of interaction points
    #
    
    def show_Motif(self):
        if self.Ob == None:
            return self.base, self.ctp, self.btp, self.Vij, self.branching, self.pk, self.wyspa
        
        elif self.Ob.name == "Stem":
            return self.Ob.disp_Stem()
        elif self.Ob.name == "MBL":
            return self.Ob.disp_MBL()
        elif self.Ob.name == "XLoop":
            return self.Ob.disp_XLoop()
        else:
            print "Motif error(show_Motif)", self.base
            print type(self.Ob)
            print self.Ob
            sys.exit(0)
        #
    #
#

# This carries the fundamental information that is needed to construct
# the various "secondary structure" elements in this approach.
class Link:

    # #############################################################
    # 190618; this should now take on the functions that we have
    # defined in Motif. Right now, we are in transition, but I think
    # this is where the operations should occur, not in the current
    # "Motif". The classes that actually should be "motif" are Stem,
    # MBL, XLoop, PseudoKnot, MultiPair, etc. and then we use the
    # tools currently in Motif to process these functions.
    # #############################################################
    
    def __init__(self, i, j, Vij, ctp, btp, branching, pk = [], wyspa = []):
        # loop can contain a set of objects of type Motif, however,
        # initially, it must be assigned one such object.
        
        self.motif = [Motif(i, j, Vij, ctp, btp, branching, pk, wyspa)]

        """@@@

        In general, most of the time, self.motif will contain only one
        Motif object -- at least with the present design of an exact
        energy (Vij). However, there is the occasional possibility
        that a region will contain more than one unique structure
        possessing an idential energy. Such structures are called
        degenerate and cannot be distinguished using free energy.
        Further, if this approach is expanded to allow a range of
        structures within a tolerance of ${\Delta\Delta G}, then this
        object could contain many elements within this range of
        tolerance.

        """
        self.Vij = Vij
    #
    
    # to add a Motif with the specified parameters. 
    def add_Motif(self, i, j, Vij, ctp, btp, branching, pk = [], wyspa = []):
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
        if not self.Vij == Vij:
            newmotif = Motif(i, j, Vij, ctp, btp, branching, pk)
            print "ERROR(Link) energy existing motif and current motif not in tolerance"
            print "new structure:"
            print newmotif.show_Motif()
            print "existing structure(s): "
            for vv in self.motif:
                print vv.show_Motif()
            #
            print "Vij(current) = %8.2f, Vij(new) = %8.2f" % (self.Vij, float(Vij))
            sys.exit(1)
        #
        # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        self.motif += [Motif(i, j, Vij, ctp, btp, branching, pk, wyspa)]
    #
    
    # add an object that is already assigned Motif (presently not used)
    def add_Motif_ob(self, newmotif):
        if not self.Vij == newmotif.Vij:
            print "ERROR(Link) energy existing motif and current motif not in tolerance"
            print "new structure:"
            print newmotif.show_Motif()
            print "existing structure(s): "
            for vv in self.motif:
                print vv.show_Motif()
            #
            sys.exit(1)
        #
        self.motif += [newmotif]
    #
    
#

# This is used as a push down stack for links at (i,j). It is only
# used in Map().
class LGroup:
    def __init__(self):
        self.lg = []
        self.nlg = 0
        # I think there should also be a maximum number of links on
        # the stack (at least in principle). However, presently, it
        # doesn't seem to be so necessary.
    #
    def add_link(self, link):
        self.lg = [link] + self.lg
        self.nlg += 1
    #
    
#

"""@

160530wkd: This ndx (index) is somewhat unnecessary. Originally, I had
the idea to make Map a 1D array that only manipulates variables in the
upper triangle. However, I had some problems getting it to work well,
so I just opted to solve it with a regular NxN matrix to save time
with messing with it. I also was somewhat resistant to applying the
Vienna package solution because it means (to some extent) involving
their software. So this is also a kind of workaround to that.

"""
class ndx:
    def __init__(self, i,j, ij):
        self.ij = ij
        self.v = (i,j)
    #
    def gij(self):
        return self.ij
    #
    
    def gv(self):
        return self.v
    #
#

# This is the main tool for keeping track of free energy data in this
# program. It stores a set of free energy values at _each_ position
# (i,j).
class Map:
    def __init__(self, N):
        self.N  = N
        if N <= 0:
            print "ERROR: size of requested matrix (%d) makes no sense" % N
            sys.exit(1)
        #
        
        self.ij = []
        self.glink = []
        
        k = 0
        # build NxN matrix
        for i in range(0,N):
            ij_row = []
            glink_row = []
            for j in range(0,N):
                ij_row   += [ndx(i,j,k)]
                glink_row += [LGroup()]
                k += 1
            #
            
            # construct row ij
            self.ij += [ij_row]
            self.glink += [glink_row]
        #
        
    #
    
    # this is used to sort entries after obtaining the free energy for
    # various types of interactions ('B', 'I', 'M', 'J', or 'P'). It
    # is also used to sort the distribution of free energies when
    # obtaining the suboptimal structures
    def mergeSortLinks(self, alist):
        # print("Splitting ",alist)
        if len(alist)>1:
            mid = len(alist)//2 # division for integers in python
            lefthalf = alist[:mid]
            righthalf = alist[mid:]
            
            self.mergeSortLinks(lefthalf)
            self.mergeSortLinks(righthalf)
            
            i=0
            j=0
            k=0
            while i < len(lefthalf) and j < len(righthalf):
                if lefthalf[i].Vij < righthalf[j].Vij:
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
    
#

# ####################################################################
# ####################################################################
# ####################################################################
# ####################################################################
# #####  Postprocessing approach: Tree/structure building tools  #####
# ####################################################################
# ####################################################################
# ####################################################################
# ####################################################################

# #######################  and  ################################## 

# ################################################################
# ##################   New Motif Object Format  ##################
# ################################################################

"""@

I have moved all the items below to Motif (from both Threads and
Vienna) because I intend to upgrade this part of Motif.

The new module (NewMotif) is intended to become the new approach on
building structural objects in Calculate, and the other programs in
the package. I think I have finally found what will work for
me. NewMotif functions something like a "Container" for various types
of objects and will use the same objects as Threads, so it makes the
process of transforming between bp notations of Vienna, or free energy
scripts with bps from LThread, or ChPair a little easier.

The main thing that needs to be added to these objections below are
information on branches and the free energy terms. Otherwise, this is
almost at the point were it can work. It might make it easier to do
the more difficult part of transforming free energy to Motif notation
that was currently a rather awkward transformation in Thread.

Additionally, presently, I have thread as the individual (base) pair
free energy tool with LThread (in Threads.py) and Pair as the
fundamental tools for the Vienna class. LNode does contain some
additional information that is not available from Pair -- the length
of stems involved in the structure. It is kind of the meta-structural
notation (Pascal-like) between Map/Motif and Pair; however, as a
notation, it is rather poor, as I came to realize eventually.


"""
def stemType2PairList(stem):
    stem_vv = []
    vtp = '-'
    if   type(stem[0]) == Pair:
        if len(stem) > 1:
            # this is mainly for checking data consistency
            if   stem[0].i < stem[1].i and stem[1].j < stem[0].j:
                vtp = 'a'
            elif stem[0].i < stem[1].i and stem[0].j < stem[1].j:
                vtp = 'p'
            
            else:
                print "ERROR(stemType2PairList): stem makes no sense"
                print "stem = ", stem
                sys.exit(1)
            #
            
        else:
            vtp = 'a'
        #
        
        for sv in stem:
            if not (sv.v == vtp or sv.v == 'base'):
                print "WARNING(stemType2PairList): inconsistent stem types"
                print "       proposed \'%s\', analyzed \'%s\'" % (sv.v, vtp)
                print "       stem: ", stem
                if len(stem) > 1:
                    print "stem[0]=(%d,%d), stem[1]=(%d,%d)" \
                        % (stem[0].i, stem[0].j, stem[1].i, stem[1].j)
                #
                
                
                #sys.exit(1)
                # probably can force a change here, but presently, it stops
            #
            
            stem_vv += [sv] # already properly formated
        #
        
        # could also use "stem_vv = stem"
        
    elif type(stem[0]) == Branch:
        # almost there, but not quite
        if len(stem) > 1:
            if   stem[0].i < stem[1].i and stem[1].j < stem[0].j:
                vtp = 'a'
            elif stem[0].i < stem[1].i and stem[0].j < stem[1].j:
                vtp = 'p'
            
            else:
                print "ERROR(stemType2PairList): stem makes no sense"
                print "stem = ", stem
                sys.exit(1)
            #
            
        else:
            vtp = 'a'
        #
        
        # stem should always be of same type
        
        for sv in stem:
            a = Pair()
            a.put_ssPair(sv.i, sv.j, "bp", vtp)
            stem_vv += [a]
        #
        
    elif type(stem[0]) == tuple or type(stem[0]) == list:
        # it might appear as a list rather than a tuple
        if len(stem) > 1:
            if   stem[0][0] < stem[1][0] and stem[1][1] < stem[0][1]:
                vtp = 'a'
            elif stem[0][0] < stem[1][0] and stem[0][1] < stem[1][1]:
                vtp = 'p'
            
            else:
                print "ERROR(stemType2PairList): stem makes no sense"
                print "stem = ", stem
                sys.exit(1)
            #
            
        else:
            vtp = 'a'
        #
        
        # stem should always be of same type
        
        for sv in stem:
            a = Pair()
            a.put_ssPair(sv[0], sv[1], "bp", vtp)
            stem_vv += [a]
        #
        
    else:
        print "ERROR(stemType2PairList): unrecognized type"
        print "type = ", type(stem[0])
        print "stem = ", stem
        sys.exit(1)
    #
    
    return stem_vv
#


class NewMotif(object):
    def __init__(self, motif):
        self.motif = motif
    #
    
    def gmotif(self):
        return self.motif
    #
    
    def disp(self):
        return disp_Motif(self.motif)
    #
    
#    


class AST(object):
    """
    AST (abstract structure type): inspired by ideas on building a
    compiler that were mentioned in the series
    
    https://ruslanspivak.com/lsbasi-part1/ 
    to
    https://ruslanspivak.com/lsbasi-part14/. 
    
    Otherwise, since it is immediately passed to APT, it may be just
    as valid to start with APT(object). I am still thinking about that
    point, as I didn't understand everything in Ruslan's discussion at
    this point. It seems like he was planning to add other stuff but
    he stopped at part 14 and didn't get to part 17 (and it has been
    some 3 years now, so I guess he will not get to it).
    """
    pass
#


class APT(AST):
    """
    APT (abstract pair type): This should function something
    resembling a token. It could potential contain type
    information. Probably more information like 'a' or 'p' could
    easily be added. That might even be good because we could define
    things like 'b' (base or root) and 'f' (fragment), which currently
    don't have a definition and might in some ways be useful.
    """
    
    def __init__(self, i, j):
        self.it = i
        self.jt = j
    #
#

class XLoop(APT):
    def __init__(self, i, j, V, ctp, btp = 'sa'):
        APT.__init__(self, i, j)   # Assign inherited APT
        self.name = "XLoop" # B,H,X
        self.vtype = 'a' # always true for RNA!
        # vtype is for the variable v in class Pair
        self.xi   = xi_fs # [nt/bp/aa/ctcfs]
        
        # ---------------------------------------------------
        
        # transition to a new motif format
        
        self.Vij     = V      # the total free energy
        self.ctp     = ctp    # connection type
        self.btp     = btp   # bond type (always 'sa')
        
        
        # ---------------------------------------------------
    #
    
    def get_base(self):
        return [(self.it, self.jt)]
    #
    
    def get_branches(self):
        return [(self.it, self.jt)]
    #
    
    def disp_XLoop(self):
        s = "%s-loop(%3d,%3d), dG = %8.2f" % (self.ctp, self.it, self.jt, self.Vij)
        return s
    #
    
    def __str__(self):
        return self.disp_XLoop()
    #
    
    def __repr__(self):
        return self.__str__()
    #
#    
    

class MBL(APT):
    def __init__(self, i, j, Vij, ctp, btp, branches):
        APT.__init__(self, i, j)   # Assign inherited APT
        self.name = "MBL" # I,J,M,P loops
        self.vtype = 'a' # always true for RNA!
        # I prefer stype, but vtype is easier to understand when
        # converting from vienna formation to other formats. It has a
        # lot of historical baggage that comes with it.
        self.xi   = xi_fs # [nt/bp/aa/ctcfs] initially default xi_fs
        
        # ---------------------------------------------------
        
        # transition to a new motif format
        
        self.MLclose   = 0.0    # global corrected FE
        self.dGd53p    = 0.0    # local corrected FE
        self.dG_fs     = 0.0    # free strand correction
        self.dGVpq     = 0.0    # free energy from branches
        self.Vij       = Vij    # the total free energy
        self.ctp       = ctp    # connection type
        self.btp       = btp    # bond type (always 'sa')
        self.branching = []     # branching points
        for k in range(0, len(branches)):
            self.branching += [branches[k]]
        #
        self.base      = (i, j) # reference point or base
        
        
        # ---------------------------------------------------
    #
    
    def get_base(self):
        base = []
        if self.ctp == 'P' or self.ctp == 'J':
            # J and P do not have a bond at the base, the base is just
            # a reference point.
            base = self.branching
        else:
            # B, I, K, and M have a bond at the base, and W is
            # effectively a cluster of bases, on of which outlines the
            # entire structure
            base = [self.base]
        return base
    #
    
    def get_branches(self):
        jj = []
        for jjk in self.branching:
            jj += [jjk]
        return jj
    #
    
    
    
    def disp_MBL(self):
        s = "%s-loop(%3d,%3d), dG = %8.2f\n" % (self.ctp, self.it, self.jt, self.Vij)
        for k in range(0,  len(self.branching)):
            vv = self.branching[k]
            s += "%2d  (%2d,%2d)\n" % (k, vv[0], vv[1])
        #
        return s
    #
    
    def __str__(self):
        return self.disp_MBL()
    #
    
    def __repr__(self):
        return self.__str__()
    #
#    
    

class Branch(object):
    """
    This makes it possible to store base pairing vectors in the form
    of a structure with tags i and j. The value in this is that when
    you manipulate the contents of an array branch, rather than
    writing it is branch[k][0] and branch[k][1], you express it as
    branch[k].i and branch[k].j
    """
    def __init__(self, i, j):
        self.i = i
        self.j = j
    #
    def __str__(self):
        return "%d  %d " % (self.i, self.j)
    #
    def __str__(self):
        s = "(%3d, %3d)" % (self.i, self.j)
        return s
    #
    def __repr__(self):
        return self.__str__()
    #
    
         
#




class Stem(APT): # inherits abstract pair type (APT)
    """
    Stem: builds a stem that can be passed to Motif or used to build a
    structure map
    """
    
    def __init__(self, stem):
        self.stem = stemType2PairList(stem) # converts all entries to class Pair()
        APT.__init__(self, self.stem[0].i, self.stem[0].j)   # Assign inherited APT
        # access using Stem.it, Stem.jt; "t" indicates "tail" of object.
        self.mark = False
        self.name = "Stem"
        self.vtype = self.stem[0].v
        # I prefer stype, but vtype is easier to understand when
        # converting from vienna formation to other formats. It has a
        # lot of historical baggage that comes with it.
        self.xi   = xi # [nt/bp/aa/ctcfs]
        self.slen = len(self.stem)
        n = self.slen - 1
        self.ih = self.stem[n].i; self.jh = self.stem[n].j
        # average stem length
        self.tlen = float((self.stem[n].i - self.stem[0].i) \
                          + (self.stem[0].j - self.stem[n].j) + 2.0)/2.0
        
        # ---------------------------------------------------
        
        # transition to a new motif format
        
        self.dG_l      = 0.0        # local stem corrected FE
        self.dG_fs     = 0.0        # local free strand corrected FE
        self.Vij       = 0.0        # the total free energy at (it, jt)
        self.dGp       = []         # individual base pair free energy contributions
        self.dGloop    = 0.0        # H/I/M loop correction at (ph,qh)
        self.Vpqh      = 0.0        # the linking free energy at (ph, qh): dGcap = dGloop + Vpqh
        self.ctp       = 'S'        # connection type at (it, jt)
        self.btp       = self.vtype # bond type (vtype but with Stem needs addendum)
        self.branching = []         # branching points (if any)
        self.base      = ()
        
        # the base reference point
        if self.vtype == 'p':
            self.base = (self.it, self.jh)   
        else:
            self.base = (self.it, self.jt)
        #
        
        """@ 
        
        I don't find it so easy to add further details of self.btp
        because this aspect is highly context dependent. The stem is
        simply a stem to the program, but when we create the stem, we
        don't know if that stem is part of a pseudoknot or just a
        regular stem, so all we can say is that the stem is parallel
        or antiparallel, nothing more. So I don't feel that there is
        the information within the structure of stem to allow this to
        be determined. 

        """
        
        # ---------------------------------------------------
    #
    
    def get_base(self):
        base = [self.base]
        return base
    #
    
    def get_stem(self):
        stm = []
        for stmk in self.stem: # class Pair
            p = stmk.i; q = stmk.j
            stm += [(p, q)] # convert to list of tuples
        return stm
    #
    
    
    def get_branches(self):
        jj = []
        for jjk in self.branching:
            jj += [jjk]
        return jj
    #
    
    
    def disp_Stem(self):
        if self.vtype == 'p':
            s = "%s Stem(%3d,%3d) -> " % (self.vtype, self.it, self.jh)
        else:
            s = "%s Stem(%3d,%3d) -> " % (self.vtype, self.it, self.jt)
        #
        
        for vv in self.stem:
            s += "(%2d,%2d), " % (vv.i, vv.j)
        #
        s += '\n'
        s += "xi = %8.2f, slen = %2d, tlen = %8.2f" % (self.xi, self.slen, self.tlen)
        return s
    #
    
    """
    The following two methods (get_tterm and get_hterm) provide
    information about the end points of the stem.
    
               parallel stem            anti-parallel stem
           ....ABCD......abcd....     ....((((......))))....
               it ih     jt jh            it ih    jh  jt
    
    t-term       (it, jh)                   (it, jt)
    h-term       (ih, jt)                   (ih, jh)
    
    """
    
    def get_tterm(self):
        if self.vtype == 'p':
            return self.it, self.jh
        else:
            return self.it, self.jt
        #
    #
    
    def get_hterm(self):
        if self.vtype == 'p':
            return self.ih, self.jt
        else:
            return self.ih, self.jh
        #
    #
    
    def __str__(self):
        s = ''
        if self.vtype == 'p':
            s += "%s%s(%3d, %3d)" % (self.vtype, self.name, self.it, self.jh)
        else:
            s += "%s%s(%3d, %3d)" % (self.vtype, self.name, self.it, self.jt)
        #
        return s
    #
    
    def __repr__(self):
        return self.__str__()
    #
    
        
#

# maybe not so necessary.
def is_pktype(tp):
    match = True
    if not (tp == 'K' or tp == 'R'):
        print "ERROR: unrecognized pseudoknot type (%s)" % tp
        match = False
    #
    return match
#


class PseudoKnot(APT):
    """
    PseudoKnot: builds a stem that can be passed to Motif or used to
    build a structure map
    """ 
    
    def __init__(self, i, j, pktype):
        APT.__init__(self, i, j) # Assign inherited APT
        # access using PK.it, PK.jt; "t" indicates "tail" of object.
        self.name      = "PseudoKnot"
        self.pktype    = pktype # must be defined
        """
        K type: core pseudoknot (anti-parallel)
        R type: extended pseudoknot (anti-parallel)
        """
        if not is_pktype(pktype):
            print "Error: unrecognized pseudoknot type (%2d,%2d)[pktype = %s]" \
                % (i, j, pktype)
            sys.exit(1)
        #
        
        self.roottmp   = [] # placeholder saves stemtail information
        self.rootstems = [] # type Stem
        self.linkages  = [] # type Stem
        
        
        
        self.xi = xi # [nt] default value
        
        # ---------------------------------------------------
        
        # transition to a new motif format
        
        self.vtype = 'pk'
        
        
        self.Vij  = 0.0     # the free energy
        self.ctp  = pktype  # connection type 
        self.btp  = 'pk'    # bond type (redundant, but basic label)
        self.branching = []
        # reference to the pMBL comprising the root stem and the
        # remaining structures arranged between ipk and jpk
        self.base = (i,j)   # tail of the structure (ipk,jpk)
        # ---------------------------------------------------
        
    #
    
    def __str__(self):
        return "%s pk struct(%3d, %3d)" % (self.pktype, self.it, self.jt)
    #
    
    
    def disp_PseudoKnot(self):
        s = ''
        if len(self.rootstems) > 0:
            s  += "%s- pk(%3d,%3d): rootstems -> " % (self.pktype, self.it, self.jt)
            for vv in self.rootstems:
                s += "(%2d,%2d)[%s], " % (vv.it, vv.jt, vv.vtype)
            #
        else:
            s  += "%s- pk(%3d,%3d): roottmp   -> " % (self.pktype, self.it, self.jt)
            for vv in self.roottmp:
                s += "(%2d,%2d), " % (vv[0], vv[1])
            #
        #
        s += "\n"
        s += "                linkages  -> "
        for vv in self.linkages:
            s += "(%2d,%2d)[%s], " % (vv.it, vv.jt, vv.vtype)
        #
        return s
    #
    
    
    def __repr__(self):
        return self.__str__()
    #
#


class ArchMP(APT):
    """Support object for MultiPair"""
    def __init__(self, i, j, v):
        APT.__init__(self, i, j) # Assign inherited APT
        self.internal = []
        self.itype    = v
        self.btype    = '-'
    #
    
    def __str__(self):
        return "island(%3d, %3d)" % (self.it, self.jt)
    #

    def disp_AMP(self):
        s = "arch(%3d,%3d)[%s]: branches[%s] --> " \
            % (self.it, self.jt, self.itype, self.btype)
        if len(self.internal) == 0:
            s += "None"
        else:
            for vv in self.internal:
                s += "(%3d,%3d), " % (vv.i, vv.j)
            #
        #
        return s
    #
    
    def __repr__(self):
        return self.__str__()
    #
#


class MultiPair(APT):
    """
    MultiPair: builds a stem that can be passed to Motif or used to
    build a structure map
    """ 
    def __init__(self, i, j, v):
        APT.__init__(self, i, j) # Assign inherited APT
        
        # access using ww.it, ww.jt; "t" indicates "tail" of object.
        a = ArchMP(i, j, v)
        self.arches = [a]
        self.name = 'ctcf'
        
        # ---------------------------------------------------
        
        # transition to a new motif format
        
        self.vtype = "wyspa"
        self.Vij   = 0.0     # the free energy
        self.ctp   = 'W'     # connection type 
        self.btp   = "wyspa" # bond type (redundant, but basic label)
        self.branching  = []      # what linkage is pointed to
        # reference to the pMBLs comprising the root stem extending from and the
        # remaining structures arranged between ipk and jpk
        self.base = (i,j)   # tail of the structure (i_W,j_W)
        # ---------------------------------------------------
    #
    
    def __str__(self):
        return "island(%3d, %3d)" % (self.it, self.jt)
    #
    
    def disp_MP(self):
        s  = "%s(%3d,%3d): " % (self.name, self.it, self.jt)
        for k in range(1, len(self.arches)):
            s += '\n'
            s += self.arches[k].disp_AMP()
        #
        return s
    #
    
    def __repr__(self):
        return self.__str__()
    #
#


def disp_Motif(otype):
    node_type = type(otype).__name__
    s = "%s(%2d,%2d) -> " % (node_type, otype.it, otype.jt)
    if node_type == "Stem":
        s = otype.disp_Stem()
    elif node_type == "PseudoKnot":
        s = otype.disp_PseudoKnot()
    elif node_type == "MultiPair":
        s = otype.disp_MultiPair()
    else:
        s += "\n!!! undefined structure !!!"
    #
    return s
#


def copyStem(stem):
    bps = []
    for bpk in stem:
        i = bpk.i; j = bpk.j; v = bpk.v; nm = bpk.name
        pr = Pair()
        pr.put_ssPair(i, j, nm, v)
        bps += [pr]
    #
    newStem = Stem(bps)
    return newStem
#


# sort stemlist according to i.
def ins_sort_StemList(StemTypes, order = 'it'):
    
    """@

    This insertion sort method has a worst case time complexity of
    (N-1)N/2. The worst case example is one where we have the list in
    exactly the opposite order of what we desire. In this list, we
    want it ordered from smallest to largest. For a simple case of
    integers, it would as follows:
    
    -------------------------
    Exibit 1 (the worst case scenario): [5, 4, 3, 2, 1] 
    
    start with i = 1
    [5, 4, 3, 2, 1] -> [4, 5, 3, 2, 1]
    increment i to 2
    [4, 5, 3, 2, 1] -> [4, 3, 5, 2, 1] -> [3, 4, 5, 2, 1]
    increment i to 3
    [3, 4, 5, 2, 1] -> [3, 4, 2, 5, 1] -> [3, 2, 4, 5, 1] -> [2, 3, 4, 5, 1]  
    increment i to 4
    [2, 3, 4, 5, 1] -> [2, 3, 4, 1, 5] -> [2, 3, 1, 4, 5] -> [2, 1, 3, 4, 5]  
                                                                -> [1, 2, 3, 4, 5]
    -------------------------
    
    Hence, 4 + 3 + 2 + 1, which is exactly 4*5/2 = 10 steps, or
    N*(N-1)/2 steps (maximum cost), where N = 5 in this case.
    
    If the list is already ordered, then the evaluation is as follows:
    
    -------------------------
    Exhibit 2 (the best case scenario): [1, 2, 3, 4, 5] 
    
    start with i = 1
    no change and j = 1
    increment i to 2
    no change and j = 2
    increment i to 3
    no change and j = 3
    increment i to 4
    no change and j = 4
    i finishes, 
    
    -------------------------
    
    so this takes only N steps, where N = 5 in this case.
    
    """
    
    if order == "it":
        for i in range(1,len(StemTypes)):    
            j = i                    
            while j > 0 and StemTypes[j].it < StemTypes[j-1].it: 
                StemTypes[j], StemTypes[j-1] = StemTypes[j-1], StemTypes[j]
                # syntactic sugar: swap the items
                j=j-1 
            #
            
        #
        
    elif order == "jt":
        for i in range(1,len(StemTypes)):    
            j = i                    
            while j > 0 and StemTypes[j].jt < StemTypes[j-1].jt: 
                StemTypes[j], StemTypes[j-1] = StemTypes[j-1], StemTypes[j]
                # syntactic sugar: swap the items
                j=j-1 
            #
            
        #
    else:
        print "ERROR(ins_sort_StemList): unrecognized option (%s)" % order
        print "                          allowed options: \"it\" or \"jt\"."
        sys.exit(1)
    #
        
                
    return StemTypes
#



    

# ##################
# #####  MAIN  #####
# ##################

def test0():
    """
    initial tests of NewMotif in preparation of this upgrade

    """
    
    a1 = Pair()
    a1.put_ssPair(0, 10, 'bp', 'a')
    a2 = Pair()
    a2.put_ssPair(1, 9, 'bp', 'a')
    root_stem = [a1,a2]
    rss = Stem(root_stem)
    
    print rss.disp_Stem()
    nm1 = NewMotif(rss)
    print nm1.motif.disp_Stem()
    print nm1.motif.name
    print nm1.motif.it, nm1.motif.jt
    
    l1 = Pair()
    l1.put_ssPair(5, 18, 'bp', 'a')
    l2 = Pair()
    l2.put_ssPair(6, 17, 'bp', 'a')
    linkage_stem = [l1, l2]
    lss = Stem(linkage_stem)
    
    print lss.disp_Stem()
    nm2 = NewMotif(lss)
    print nm2.motif.disp_Stem()
    print nm2.motif.name
    print nm2.motif.it, nm2.motif.jt
    print "-----"
    print nm2.gmotif().disp_Stem()
    print "-----"
    
    pk = PseudoKnot(0, 18, 'K')
    pk.rootstems = [rss]
    pk.linkages  = [lss]
    print pk.disp_PseudoKnot()
    print pk.pktype
    nm3 = NewMotif(pk)
    print nm3.gmotif().name, nm3.gmotif().it, nm3.gmotif().jt
    print nm3.gmotif().disp_PseudoKnot()
    print "xxx"
    print nm3.disp()
#

def test1():
    stem_list = []
    for k in range(0, 5):
        stem_list += [(k, 12-k)]
    #
    stem_Pair = stemType2PairList(stem_list)
    print "stem_list:   ", stem_list
    print "stem_Pair 1: ", stem_Pair

    stem_Brnx = []
    for pv in stem_list:
        stem_Brnx += [Branch(pv[0], pv[1])]
    #
    stem_Pair2 = stemType2PairList(stem_Brnx)
    print "stem_Pair 2: ", stem_Pair2
    
    stem_Pair3 = stemType2PairList(stem_Pair)
    print "stem_Pair 3: ", stem_Pair3
#


if __name__ == "__main__":
    #test0()
    test1()
#
