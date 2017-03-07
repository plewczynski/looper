#!/usr/bin/env python

# package name: Motif
# creater:      Wayne K Dawson
# creation date: 170126
# last update:   170126


# This module is used primarily by "chreval.py". However, as with many
# packages that are building up, the package has become complicated
# enough, that the objects in this module are used by more than one program.

# Therefore, since this is part of a common set of objects, it is
# better to share the code among all the components that need
# it. Therefore, I have separated this package from chreval to allow
# independent testing and so that can be used by other programs for
# various things (for example, building heatmaps from input data).

# Moreover, I would rather be able to work with it as an independent
# module so that I can test various things without having to run the
# whole program.

import sys


#####################################################################
######################  BASIC LINKAGE OBJECTS  ######################
#####################################################################


# #################################################
# #############  LINKAGE DEFINITIONS  #############
# #################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# CBL = chromatin binding linkage
# ctp = connection type, established at position (i,j). Let

# (i,j) current position under analysis (may or may not have a CBL)
# (i',j') a prior CBL 

# ###############  connection definitions (ctp and btp):  ###############

# ######### CONNECTION TYPE (ctp):

# 'B', (i,j) -> [(i,j)] => connection type is a bond ('B') at (i,j) or
# the terminal connection at the end of an antiparallel loop.

# 'b', (i,j) -> [(i,j)] => connection type is the terminal connection
# at the end of an antiparallel loop at (i,j).

# 'S', (i,j) -> [(i+1,j-1), (i+2,j-2) ... ] or [(i,j),(i+1,j+1),
# (i+2,j+2) ... ] or [(i,j),(i-1,j-1), (i-2,j-2) ... ]: "Stem".  this
# is similar to the RNA stem I used before; however, in the case of
# chromatin, it is important to remember that both antiparallel and
# parallel directions are allowed.

# 'I', (i,j) -> [(i',j')] => connection type is a "I-loop" ('I'):
# closed at (i,j) with the next adjacent cross link at (i',j'), where
# i' > i + 1 and/or j' < j -1. Typically, this structure should be
# more favorable than 'B' if it exists, but it is not necessarily so.

# 'J', -> [(i',j')] => connection type is a "J-loop" ('J'): _open_ at
# (i,j) with the next adjacent cross link at (i',j'), where i' > i + 1
# and/or j' < j -1. A J-loop can be more favorable than the I-loop
# when the binding free energy at (i,j) is positive.

# 'K', (i, j) -> join = [(ib,jb)], pk = [(ik1,jk1), (ik2,jk2)...]  =>
# connection type pseudoKnot (presently all combinations of H
# type). There is a bit of a problem with chromatin that I did not
# encounter in RNA because the chain can be either parallel or
# antiparallel. I will only allow one solution presently. This is
# because I ASSUME that the chain will favor one or the other
# direction but not both for a good solid pseudoknot. This could be
# problematical, but, for the moment, I'll take the risk.

# 'l': a pseudoknot in the "linkage" part of the pseudoknot. There is
# nothing in principle that is different between the linkage and the
# main structure. Both are essentially singletons. However, it is
# convenient to notate it differently for human readability and due to
# the fact that linkages are generally more restricted on their
# interaction with the rest of the structure.

# 'M', (i,j) -> [(i1,j1),(i2,j2)] => connection type is a "M-loop":
# closed at (i,j) with branches located at [(i1,j1), (i2,j2)
# ... (in,jn)], subject to the condition that i < i1, j1 < i2, i2 <
# j3, .... jn < j.  M-loop is for "multibranch loop" and this form is
# typically more favorable than an 'I' or a 'B', but it depends on the
# situation.

# ['P', [(i1,j1),(i2,j2)]] => connection type is a "P-loop": _open_ at
# (i,j) with branches located at [(i1,j1), (i2,j2) ... (in,jn)],
# subject to the condition that i < i1, j1 < i2, i2 < j3, .... jn < j.
# P-loop is for "principal loop". Typically, the M-loop should be more
# favorable, but it could be favored if the binding free energy 'B' at
# (i,j) is positive.

# ['W', [(i1, j1, i2, j2 ....)]] => connection type CTCF-island. This
# is unique to CTCF islands. In this type of interaction, all the
# items on the list are in effect associated. This requires a
# different notation from the singleton that presently I assume can be
# written as

# ....|......|.........|..............|.....
#     i1     j1        i2             j2

# where the meaning is that all these regions interact with each
# other.



# ######### BOND TYPE (btp):

# '-': none, used in the case of pointers because it represents
# neither a singleton nor a PET cluster.

# 'bgn'/'end': markers for the beginning and ending of stems ('S'),
# pseudoknots ('K') and CTCF islands ('W')

# 'c': This is a CTCF cluster that we define as "convergent"; i.e.,
# having the "... > .... < ..." type of interaction that _probably_
# results in an _antiparallel_ chain configuration. This configuration
# results in a very strongly binding interaction and is most commently
# encountered.

# 'l': This is a CTCF cluster oriented in the tandem left
# configuration; ; i.e., having the "... < .... < ..." type of
# interaction that _probably_ results in an _parallel_ chain
# configuration. parallel configuration. This configuration results in
# a far less strongly binding interaction (about 1/3) but is also
# often observed in the experimental data.

# 'r': This is a CTCF cluster oriented in the tandem right
# configuration; ; i.e., having the "... > .... > ..." type of
# interaction that _probably_ results in an _parallel_ chain
# configuration. This configuration results in a far less strongly
# binding interaction (about 1/3) but is also often observed in the
# experimental data.

# 'rp2':   RNApolII

# 's': This is the "singleton" bond type. Previously (before 160721),
# I called this type 'n', for "normal". What is normal and is
# singleton any better? Presently, I don't know, but we read in
# chromatin data of weak structural ininteractions between local
# regions of the structure and we predict some loops structures
# resulting from that. These loop structures are assumed to have the
# chromatin chains oriented in an unspecifice orientation (parallel
# loop or antiparallel loop direction) and the interactions are
# presumed to be weak. Note that 's' could also refer to divergent
# CTCF interactions; i.e., having the "< ... >" type of interaction
# that results in weakly interacting antiparallel chain configurations
# similar to the convergent CTCF structures only with far less binding
# free energy.

# 160921wkd: Now that I also include the stem ('S') elements, along
# with the singleton ('s') notation, it is important to specify if the
# collection of stems is parallel or antiparallel oriented. Therefore,
# we must add a qualifier to this notation; i.e., for structures
# containing the key 'S', we should specify 'sa' (for singleton
# antiparallel) and 'sp' (for singleton parallel).

# 't' (i.e., "tandem" or "parallel"): This refers to tandem
# configurations wherein the particular orientation of the CTCF
# binding sites is unknown. In such cases, the program will
# automatically assign this label, since the binding energy of 'l' and
# 'r' are basically identical at the current level of understanding.

# 'wyspa': CTCF island.


# These additional "bond types" permit us to distinguish between
# direction and type of interaction: i.e., RNApolII, CTCF
# "convergent", CTCF "tandem".

# #########

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# #################################################
# #############  LINKAGE DEFINITIONS  #############
# #################################################

# 160715wkd: I added this additional object because the program will
# have to be expanded to handle convergent and tandem structures that
# are coupled, and also because some structures will turn out to be
# degenerate, so the program should be able to classify multiple
# structures with the same energy within the same class. However, I
# need to process the structures properly. Therefore, Link contains a
# list of objects that have the same free energy Vij and Motif()
# represents the objects contained in Link().



class Motif:
    def __init__(self, i, j, Vij, ctp, btp, join, pk = [], wyspa = []):
        self.Vij  = Vij     # the free energy
        self.ctp  = ctp     # connection type 
        self.btp  = btp     # bond type
        self.join = join    # what linkage is pointed to
        if len(self.join) == 0:
            if ctp == 'B':
                self.join = [(i,j)]
            else:
                print "ERROR(Motif(%d,%d)): unrecognized linkage for type %s." \
                    % (i, j, ctp)
                sys.exit(1)
            #
        #
        self.base = (i,j)   # the base reference point
        self.pk    = pk
        self.wyspa = wyspa
        #
    #
    
    def get_base(self):
        bonds = []
        if self.ctp == 'P' or self.ctp == 'J':
            # J and P do not have a bond at the base
            bonds = self.join
        else:
            # B, I, K, and M have a bond at the base, and W is
            # effectively a cluster of bonds, on of which outlines the
            # entire structure
            bonds = [self.base]
        return bonds
    #
    
    def get_branches(self):
        jj = []
        for jjk in self.join:
            jj += [jjk]
        return jj
    #
    
    def get_pks(self):
        return self.pk
    #
    
    def get_wyspa(self):
        return self.wyspa
    #
    
    def show_Motif(self):
        return self.base, self.ctp, self.btp, self.Vij, self.join, self.pk, self.wyspa
    #
#

# This carries the fundamental information that is needed to construct
# the various "secondary structure" elements in this approach.
class Link:
    def __init__(self, i, j, Vij, ctp, btp, join, pk = [], wyspa = []):
        # loop can contain a set of objects of type Motif, however,
        # initially, it must be assigned one such object.
        
        self.motif = [Motif(i, j, Vij, ctp, btp, join, pk, wyspa)]
        # In general, most of the time, self.motif will contain only
        # one Motif object -- at least with the present design of an
        # exact energy (Vij). However, there is the occasional
        # possibility that a region will contain more than one unique
        # structure possessing an idential energy. Such structures are
        # called degenerate and cannot be distinguished using free
        # energy.  Further, if this approach is expanded to allow a
        # range of structures within a tolerance of ${\Delta\Delta G},
        # then this object could contain many elements within this
        # range of tolerance.
        self.Vij = Vij
    #
    
    # to add a Motif with the specified parameters. 
    def add_Motif(self, i, j, Vij, ctp, btp, join, pk = [], wyspa = []):
        if not self.Vij == Vij:
            newmotif = Motif(i, j, Vij, ctp, btp, join, pk)
            print "ERROR(Link) energy existing motif and current motif not in tolerance"
            print "new structure:"
            print newmotif.show_Motif()
            print "existing structure(s): "
            for vv in self.motif:
                print vv.show_Motif()
            #
            sys.exit(1)
        #
        self.motif += [Motif(i, j, Vij, ctp, btp, join, pk, wyspa)]
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

# This is used as a push down stack for links at (i,j). It is mainly
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

# 160530wkd: This ndx (index) is somewhat unnecessary. Originally, I
# had the idea to make Map a 1D array that only manipulates variables
# in the upper triangle. However, I had some problems getting it to
# work well, so I just opted to solve it with a regular NxN matrix to
# save time with messing with it. I also was somewhat resistant to
# applying the Vienna package solution because it means (to some
# extent) involving their software. So this is also a kind of
# workaround to that.
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
        self.ij = []
        self.glink = []
        
        k = 0
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
            mid = len(alist)//2
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
