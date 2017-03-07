#!/usr/bin/env python

# chreval.py (CHRomatin EVALuation program for heatmaps)
#
#    Chreval is a dynamic programming method for determining the most
#    probably chromatin structure arrangement as well as the
#    distribution of chromatin structure arrangements as a function of
#    free energy inside of loops
#
#    command line example:
#    > chreval.py -f chr17_55674644_55751596_res5kb.heat -add_PET_wt 
#
#    Let the file chr17_55674644_55751596_res5kb.heat be called
#    "chrN_x_y_res5kb.heat" for short. The program requires that
#    "chrN_x_y_res5kb.heat" consist of a matrix of integer values that
#    indicate the instances of an interaction. The input file
#    "chrN_x_y_res5kb.heat" must have the extension "heat". The output
#    consists of a directory with the same name as the input file
#    minus the extension ".heat"; i.e., the directory name
#    "chrN_x_y_res5kb". This directory contains (in separate files)
#    all the calculated loop structures. These files have the
#    extension "DBN" and can be read by the 3rd party program
#    varna. Additionally, Chevral deposits two additional files:
#    chrN_x_y_res5kb_BDwt.clust contains a matrix with the Boltzmann
#    probabilities for different interactions and
#    chrN_x_y_res5kb_summary.txt contains a shorthand list of the
#    secondary structures.



# Presently, I think the name is rather unoriginal (i.e., "CHRomatin
# EVALuate"). Anyway, what is important is that it works, and maybe we
# can come up with a better name later.


from math import exp, log, ceil, floor
import sys
import os
import string

# main tool objects
from Functions import KahanSumExp
from GetOpts import GetOpts

# Free energy parameters
from FreeEnergy import FreeEnergy

# Motif object representation
from Motif import Motif # a core object of these building programs
from Motif import Link
from Motif import LGroup
from Motif import Map # main map for all the FE and structures

# LThread object representation
from Threads import DispThreads
from Threads import LThread


# ################################################################
# ######################  Global constants  ######################
# ################################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# for entropy calculations
from Constants import kB # [kcal/molK]
from Constants import xi # [bps]
from Constants import lmbd # [bps]
from Constants import gmm # dimensionless but related to D/2        
# constants: weights for the enthalpy terms 
from Constants import febase # [kcal/mol]
from Constants import feshift # (dimensionless, usually = 1)

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

PROGRAM      = "chreval.py"  # name of the program

# tests
TEST0        = False  # for testing of objects
TEST1        = False  # for testing of objects
TEST2        = False  # for testing of objects
TEST5        = False  # test the Kahan algorithm

# debugging options 
SHOWMAIN     = False # Debugging in main()
CHECK_glink  = False # Checks the different Link objects at at (i,j)
DEBUG_Trace  = False # Generally for checking Trace()
DEBUG_minFE  = False # Generally for checking Calculate()
DEBUG_M      = False # additional checks on M-loop routines
DEBUG_I      = False # additional checks on I-loop routines
MONITOR      = False # for checking outputs, particularly M-loop results
CHECK_ILOOP  = False  # for checking N2 vs N4 results
STOP_AT_FIND = False  # stops program when encouters condition

# 2. special local global function

# This is not used so much now that I have introduced GetOpts.py,
# but it may still be useful at the very beginning of the program and
# possibly in other parts.

def usage():
    print "USAGE: %s -f file.heat " % PROGRAM
#



# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ################################################################
# ################################################################



    






####################################################################
######################  MINIMUM FREE ENERGY  #######################
####################################################################

# Calculate() is the main machine of this program
class Calculate:
    def __init__(self, cl):
        flag_debug =  False
        if not cl.set_GetOpts:
            print "ERROR: options are not properly configured"
            sys.exit(1)
        #
        
        self.flnm    = cl.f_heatmap[0] # elements are a list now
        self.fe      = FreeEnergy(cl)
        self.T       = cl.T # temperature is undefined!
        
        # output option
        self.p_all_1D      = cl.p_all_1D
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


        self.dGrange = cl.dGrange
        # The search for suboptimal structures will range between
        # dGmin and dGmin + dGrange.
        
        # settings related to the PET clusters
        
        self.PET_wt0 = cl.PETwt                  # PET wt scale
        self.ctcf_tthresh = self.fe.ctcf_tthresh # tandem threshold
        self.ctcf_cthresh = self.fe.ctcf_cthresh # convergent threshold
        self.pssbl_ctcf   = self.fe.pssbl_ctcf   # found ctcfs in input file
        self.edge_ctcf    = self.fe.edge_ctcf    # edge ctcf from input file
        
        # set up the heat map "enthalpy"
        self.hv = []
        self.hv, self.ctcf_setv, self.N = self.fe.read_heat(self.flnm)
        
        self.pssbl_ctcf   = self.fe.pssbl_ctcf   # found ctcfs in input file
        self.edge_ctcf    = self.fe.edge_ctcf    # edge ctcf from input file
        
        # make a complete list of all possibilities without redundancies
        self.all_ctcf = {}
        self.all_ctcf.update(self.ctcf_setv)
        for pc in self.pssbl_ctcf.keys():
            self.all_ctcf.update({pc: self.pssbl_ctcf[pc]})
        self.all_ctcf.update(self.edge_ctcf)
        
        if flag_debug:
            print "ctcf_setv:        ", self.ctcf_setv
            print "pssbl_ctcf:       ", self.pssbl_ctcf
            print "edge_ctcf:        ", self.edge_ctcf
            print "_________________________________________"
            print "all_ctcf:         ", self.all_ctcf
            sys.exit(0)
        #
        
        # initialize the free energy matrix
        self.dG = []
        self.dG = self.fe.mtools.initialize_matrix(self.dG, self.N, 0.0)
        
        
        # build map layout
        self.smap       = Map(self.N)
        # this is the core.
        
        self.iloop      = LGroup()
        self.V_iloop_mx = INFINITY
        self.n_iloop_mx = 21
        self.wt_iloop   = 0.20
        # I am not really so sure that this method of search is
        # particularly useful anymore. I tried using this when I was
        # developing some speed up techniques with RNA structure
        # prediction in vsfold6, and it appeared to be helpful in
        # searching dsRNA structures. However, using it here, I found
        # that it was not very effective. In general, the best speed
        # up is coming from look up at each point (i,j) as in
        # new_find_best_I() rather than this ratty contraption.
        # nevertheless, I still don't have the courage to kill it.
        
        
        # set a lower/upper bound for the I-loop search with
        # new_find_best_I().
        self.span       = 8
        # Although new_find_best_I() should work just fine, I still
        # find that it sometimes misses degenerate structures. These
        # are typically insignificant, and so I do not think this will
        # particularly matter, but I cannot figure out exactly why an
        # item that is clearly on the list is missed presently.  This
        # needs at least some looking into.
        
        
        # Trace back variables
        self.opt_ss_seq = []
        self.maxlayers = self.N/2 + 10
        # The maxlayers is used in traceback routines in both this
        # package traceback_mFE() and in class Trace().
        
        # In general, I found that maxlayers = 10 was
        # sufficient. However, when using this program in very large
        # sets of data, I discovered some cases where get_traces
        # failed because maxlayer = N/6 was too small.
        
        # Ultimately, it is not clear where the limit should be. In
        # principle, there could be a solution where there are N/2
        # bonds. Then, if you have 2 beads, N/2 = 1, maxlayers > 1
        # really means you have a problem! Therefore, just to be safe,
        # I set this to N/2 -- the absolute theoretical maximum.
        
        # It seems kind of absurd that N/2 be required for a N > 100
        # bead search, but because of the recursion, it is possible
        # (in principle). The only issue I foresee is that I do not
        # know the upper limit to the number of recursion steps before
        # the program complains. In general, with multibranch looping
        # (which is far more entropically favorable) this limit should
        # never be necessary, but it is truly unclear how many levels
        # are required, so this is what I am forced to
        # conclude. Obviously, if the N/2 limit is exceeded, something
        # is DEFINITELY wrong.
        
        
        # pseudoknot parameters
        self.scan_ahead   = 10
        self.pk_thresh    =  4
        self.ctcf_tthresh = self.fe.mtools.ctcf_tthresh
    #
    
    
    
    
        
    # The top of the list need not be an 'M', 'K', 'I' or 'B'!!!  So this
    # makes sure that the structure that is selected is one of
    # these. It is used in I-loop searches; i.e., find_best_I() and
    # related structures. So the purpose is to find a structure with a
    # single closing point at (k,l), not a pMBL or a J-loop. The main
    # question I have now is whether this should also include a
    # PK. Currently (160911) by scanNN_MIB(), scan_narrow(),
    # find_best_I(), localscan_for_I()
    def filter_BIKMSW_from_glink(self, k, l):
        #
        V   = INFINITY
        cp  = []
        ctp = ''
        pks = []
        wps = []
        if DEBUG_I:
            print "-------------"
        for m in range(0, len(self.smap.glink[k][l].lg)):
            ctpx = self.smap.glink[k][l].lg[m].motif[0].ctp
            if DEBUG_I:
                print "(%2d,%2d)[%d]: %8.3f [%s]" % (k,l, m, self.smap.glink[k][l].lg[m].Vij, ctpx)
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
                #
                #
        if DEBUG_I:
            print ".", cp
        #
        return V, cp, ctp, pks, wps 
    #
    
    
    # The top of the list need not be an 'M', 'I' or 'B'!!!  So this
    # makes sure that the structure that is selected is one of
    # these. It is used in I-loop searches; i.e., find_best_I() and
    # related structures. So the purpose is to find a structure with a
    # single closing point at (k,l), not a pMBL or a J-loop. The main
    # question I have now is whether this should also include a
    # PK. Currently (160911) by scanNN_MIB(), scan_narrow(),
    # find_best_I(), localscan_for_I()
    def filter_BIM_from_glink(self, k, l):
        #
        V   = INFINITY
        cp  = []
        ctp = ''
        if DEBUG_I:
            print "-------------"
        for m in range(0, len(self.smap.glink[k][l].lg)):
            ctpx = self.smap.glink[k][l].lg[m].motif[0].ctp
            if DEBUG_I:
                print "(%2d,%2d)[%d]: %8.3f [%s]" % (k,l, m, self.smap.glink[k][l].lg[m].Vij, ctpx)
            if ctpx == 'B' or ctpx == 'I' or ctpx == 'M':
                # 160601wkd: M is also closed at the ends!!!!
                if self.smap.glink[k][l].lg[m].Vij < V:
                    for vv in self.smap.glink[k][l].lg[m].motif:
                        cp +=  vv.get_base() # original: cp = (k,l)
                    V  = self.smap.glink[k][l].lg[m].Vij
                    ctp = ctpx
                #
            #
        if DEBUG_I:
            print ".", cp
        #
        return V, cp, ctp
    #


    # The top of the list is _rarely_ likely to be a 'J'.  The purpose
    # of this function is to find potential motif at (i',j') from (i,j)
    # such that i < i', j' < j. This routine is only used by
    # scanNN_MIB().
    def filter_J_from_glink(self, k, l):
        # 
        V   = INFINITY
        cp  = []
        ctp = ''
        for m in range(0, len(self.smap.glink[k][l].lg)):
            ctpx = self.smap.glink[k][l].lg[m].motif[0].ctp
            if ctpx == 'J':
                if self.smap.glink[k][l].lg[m].Vij < V:
                    for n in range(0, len(self.smap.glink[k][l].lg[m].motif)):
                        cp += self.smap.glink[k][l].lg[m].motif[n].get_base()
                        # print cp
                    #print len(self.smap.glink[k][l].lg[m].motif)
                    # sys.exit(0)
                    V  = self.smap.glink[k][l].lg[m].Vij
                    ctp = ctpx
                #
            #
        #
        return V, cp, ctp
    #
    
    def scanNN_MIB(self, i, j):
        Vbest = INFINITY
        cp    = []
        ctp    = ''
        
        V1   = V2   = V3   = Vbest
        cp1  = cp2  = cp3  = cp
        ctp1 = ctp2 = ctp3 = ctp
        
        
        if (i+1) >= (j-1):
            Vp1m1 = INFINITY
        else:
            V1, cp1, ctp1 = self.filter_J_from_glink(i+1, j-1)        
            Vp1m1 = V1
        #    
        
        
        if (i+1) >= (j):
            Vp1m0 = INFINITY
        else:
            V2, cp2, ctp2 = self.filter_J_from_glink(i+1, j  )        
            Vp1m0 = V2
        #
        
        if i >= (j-1):
            Vp0m1 = INFINITY
        else:
            V3, cp3, ctp3 = self.filter_J_from_glink(i  , j-1)        
            Vp0m1 = V3
        #
        
        
        
        #
        if DEBUG_I:
            print "(%2d,%2d): V[i+1,j  ](%8.3f), V[i+1,j-1](%8.3f), V[i  ,j-1](%8.3f)" \
                % (i, j, Vp1m0, Vp1m1, Vp0m1)
        #
        
        tg = ''
        # find and save the minimum energy out of the set
        # of FE values (if it exists)
        if Vp1m1 < Vp0m1 and Vp1m1 < Vp1m0:
            # in this case, both Vp0m1 and Vp1m0 are eliminated 
            tg = '1'
            
            if Vp1m1 < INFINITY:
                Vbest = V1
                cp    = cp1
                ctp   = ctp1
            else:
                Vbest = INFINITY
            #
            
            if DEBUG_I:
                if len(self.smap.glink[i+1][j-1].lg) > 0:
                    print "best V[i+1,j-1] = dG(%2d,%2d)[%s][%8.3f] <- (%2d,%2d)" \
                        % (i+1, j-1, ctp, Vp1m1, i, j)
                #
                else:
                    print "best V[i+1,j-1], (%2d,%2d) <- (%2d,%2d)   empty" % (i+1, j-1, i, j)
                #
            #
        elif Vp1m0 < Vp0m1:
            # can also satisfy case Vp1m1 == Vp1m0, but Vp0m1 is eliminated
            tg = '2'
            
            if Vp1m0 < INFINITY:
                Vbest = V2
                cp = cp2
                ctp = ctp2
                
                # test for additional caess
                if Vp1m1 == Vbest: # this case is still possible
                    cp += cp1
                #
            else:
                Vbest = INFINITY
            #
            
            
            if DEBUG_I:
                if len(self.smap.glink[i+1][j  ].lg) > 0:
                    print "best V[i+1,j  ] = dG(%2d,%2d)[%s][%8.3f] <- (%2d,%2d)" \
                        % (i+1, j  , ctp, Vp1m0, i, j)
                    #
                else:
                    print "best V[i+1,j  ], (%2d,%2d) <- (%2d,%2d)   empty" % (i+1, j  , i, j)
                #
            #
        else:
            # can also satisfy case Vp1m1 == Vp0m1 or Vp1m0 == Vp0m1
            tg = '3'
            
            if Vp0m1 < INFINITY:
                Vbest = V3
                cp    = cp3
                ctp   = ctp3
                
                # test for additional possible acases
                if Vp1m1 == Vbest:
                    cp += cp1
                if Vp1m0 == Vbest:
                    cp += cp2
                #
            else:
                Vbest = INFINITY
            #
            
            
            if DEBUG_I:
                if len(self.smap.glink[i  ][j-1].lg) > 0:
                    print "best V[i  ,j-1] = dG(%2d,%2d)[%s][%8.3f] <- (%2d,%2d)" \
                        % (i  , j-1, ctp, Vp0m1, i, j)
                #
                else:
                    print "best V[i  ,j-1], (%2d,%2d) <- (%2d,%2d)   empty" % (i  , j-1, i, j)
                # 
                
                if Vp0m1 == INFINITY:
                    # this is not so common, but it can happen when
                    # searching for I loops at (i+1,j-1) in the case
                    # where (i +1, j-1) satisfies j = i+3 and
                    # (i+1,j-1) has a bound interaction 'B'.
                    print "found p0m1 -> 1000"
                    # sys.exit(0)
                #
            #
        #
        
        # ###########   should probably remove  ################
        # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        flag_p0m1 = False
        if flag_p0m1:
            if i == (j-1):
                # this handles some rare cases when searching for I loops;
                # i.e., at (i+1,j-1)
                V4, cp4, ctp4 = self.filter_BIM_from_glink(i  , j  )        
                if V4 < Vbest:
                    Vbest = V4
                    # print "V4: ", cp # trace problems with cp
                    
                    if len(cp) > 0:
                        cp += cp4
                    else:
                        cp = cp4
                    ctp = ctp4
                
                    if DEBUG_I:
                        print "used 4th option: ", Vbest, cp, ctp
                        if STOP_AT_FIND:
                            sys.exit(0)
                    #
                #
            #
        #
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # ###########   should probably remove  ################
        
        if len(cp) > 0:
            if DEBUG_I:
                print "scanNN_MIB: ij(%2d,%2d) ctp(%s), V=%12.4f, cp=%s" % (i, j, ctp, Vbest, cp)
            cp = list(set(cp))
        #
        return Vbest, cp, ctp
    #
    
    def scan_narrow(self, i, j,
                    best_dGI,
                    best_cpI,
                    best_ctpI,
                    flag_hitI,
                    best_dGJ,
                    best_cpJ,
                    best_ctpJ,
                    flag_hitJ):
        if DEBUG_I:
            print "first part: going through a small part of the I loop"
        #
        # chose a small local zone to scan ([i,i+5],[j-5,j+1])
        j_min = max(i + 1, j - self.span)
        for l in range(j_min, j + 1):
            i_max = min(i + self.span, l)
            for k in range(i, i_max):
                # print "ij = (%2d,%2d), kl = (%2d,%2d)" % (i, j, k, l)
                if k == i and l == j:
                    if DEBUG_I:
                        print "debug_new_find_best_I: skipping kl(%d,%d) = ij(%d,%d)" % (k, l, i, j)
                    continue
                #
                if len(self.smap.glink[k][l].lg) == 0:
                    continue
                #
                if k >= j:
                  continue
                #
                
                ctp = self.smap.glink[k][l].lg[0].motif[0].ctp
                if self.hv[k][l] > 0.0 or ctp == 'K':
                    # search the list in glink for the best I,J
                    V, cp, ctp, pks, wps = self.filter_BIKMSW_from_glink(k, l)
                    
                    
                    if i <= k and l <= j:
                        if V <= best_dGJ and V < INFINITY:
                            # found an optimal J-loop
                            if DEBUG_I:
                                print "p1 J:(%2d,%2d) -> (%2d,%2d)[%s] %8.3f >= %8.3f, " \
                                    % (i, j, k, l, ctp, V, best_dGJ), cp
                            #
                            flag_hitJ = True
                            if V == best_dGJ:
                                best_cpJ += cp
                            else:
                                best_cpJ  = cp
                            best_dGJ = V
                            best_ctpJ = ctp
                        #
                    #
                    if i < k and l < j:
                        if V <= best_dGI and V < INFINITY:
                            # found an optimal I-loop
                            if DEBUG_I:
                                print "p1 I:(%2d,%2d) -> (%2d,%2d)[%s] %8.3f >= %8.3f, " \
                                    % (i, j, k, l, ctp, V, best_dGI), cp
                            #
                            flag_hitI = True
                            if V == best_dGI:
                                best_cpI += cp
                            else:
                                best_cpI  = cp
                            best_dGI = V
                            best_ctpI = ctp
                        #
                    #
                #
            #
        #
        return best_dGI, best_cpI, best_ctpI, flag_hitI, best_dGJ, best_cpJ, best_ctpJ, flag_hitJ
    #
                  
    
    # I-loop calculation: (Note, I would still prefer to make this
    # into a separate class, but the best placement presently is
    # here.)
    def new_find_best_I(self, i, j):
        if DEBUG_I:
            print "new_find_best_I(%d,%d)" % (i, j)
        best_dGI = INFINITY
        best_cpI = []
        best_ctpI = ''
        flag_hitI = False
        best_dGJ = best_dGI
        best_cpJ = best_cpI
        best_ctpJ = best_ctpI
        flag_hitJ = False
        
        # 160720wkd: So far, I have not been able to completely rid
        # myself of part 1, though I have been able to reduce
        # self.span to only 4 nt. Presently, all I can say is that,
        # somehow, a small search region at the borders is still
        # required, though I still don't understand why. It should be
        # able to solve this from part 3.
        flag_scan_narrow = True
        if flag_scan_narrow:
            best_dGI, best_cpI, best_ctpI, flag_hitI, \
                best_dGJ, best_cpJ, best_ctpJ, flag_hitJ \
                = self.scan_narrow(i, j,
                                   best_dGI,
                                   best_cpI,
                                   best_ctpI,
                                   flag_hitI,
                                   best_dGJ,
                                   best_cpJ,
                                   best_ctpJ,
                                   flag_hitJ)
        #
        
        if DEBUG_I:
            print "1st part I: ", flag_hitI, best_ctpI, best_dGI, best_cpI
            print "1st part J: ", flag_hitJ, best_ctpJ, best_dGJ, best_cpJ
            
            print "third part: going through nearest neighbors"
        #
        V, cp, ctp = self.scanNN_MIB(i, j)
        if V <= best_dGJ and V < INFINITY:
            # found an optimal J-loop
            if DEBUG_I:
                print "p3 J:(%2d,%2d) -> {%s}[%s] %8.3f >= %8.3f " \
                    % (i, j, cp, ctp, V, best_dGJ)
            #
            flag_hitJ = True
            best_cpJ  = cp
            best_dGJ  = V
            best_ctpJ  = ctp
        #

        V, cp, ctp = self.scanNN_MIB(i+1, j-1)
        if V <= best_dGI and V <  INFINITY:
            # found an optimal I-loop
            if DEBUG_I:
                print "p3 I:(%2d,%2d) -> {%s}[%s] %8.3f >= %8.3f, " \
                    % (i, j, cp, ctp, V, best_dGI)
            #
            flag_hitI = True
            best_cpI  = cp
            best_dGI  = V
            best_ctpI  = ctp
        #
            
        if flag_hitI:
            if DEBUG_I:
                cp = best_cpI
                print "new_find_best_I-[I]: (%2d,%2d) -> {%s}[%8.3f][%s]" \
                    % (i, j, cp, best_dGI, best_ctpI) 
            #
        else:
            best_dGI = INFINITY
        #
        if flag_hitJ:
            if DEBUG_I:
                cp = best_cpJ
                print "new_find_best_I-[J]: (%2d,%2d) -> {%s}[%8.3f][%s]" \
                    % (i, j, cp, best_dGJ, best_ctpJ)
            #
        else:
            best_dGJ = INFINITY
        #
        if len(best_cpI) > 1:
            best_cpI = list(set(best_cpI))
        if len(best_cpJ) > 1:
            best_cpJ = list(set(best_cpJ))
        
        return best_dGI, best_cpI, best_dGJ, best_cpJ
    #                
    
    # I-loop calculation: (Note, I would still prefer to make this
    # into a separate class, but the best placement presently is
    # here.)
    def find_best_I(self, i, j):
        if DEBUG_I:
            print "find_best_I(%d,%d)" % (i, j)
        best_dGI = 0.0
        best_cpI = []
        best_ctpI = ''
        flag_hitI = False
        best_dGJ = best_dGI
        best_cpJ = best_cpI
        best_ctpJ = best_ctpI
        flag_hitJ = False
        for l in range(i+1,j+1):
            for k in range(i, l):
                if k == i and l == j:
                    if DEBUG_I:
                        print "debug_find_best_I: skipping kl(%d,%d) = ij(%d,%d)" % (k, l, i, j)
                    continue
                #
                if len(self.smap.glink[k][l].lg) == 0:
                    continue
                #
                
                ctp = self.smap.glink[k][l].lg[0].motif[0].ctp
                if self.hv[k][l] > 0.0 or ctp == 'K' or ctp == 'S':
                    # search the list in glink for the best I,J
                    V, cp, ctp, pks, wps = self.filter_BIKMSW_from_glink(k, l)
                    
                    if i <= k and l <= j:
                        if V <= best_dGJ:
                            # found an optimal J-loop
                            if DEBUG_I:
                                print "(%2d,%2d) -> (%2d,%2d)[%s] %8.3f >= %8.3f, " \
                                    % (i, j, k, l, ctp, V, best_dGJ), cp
                            #
                            flag_hitJ = True
                            if V == best_dGJ:
                                best_cpJ += cp
                            else:
                                best_cpJ  = cp
                            best_dGJ = V
                            best_ctpJ = ctp
                        #
                    #
                    if i < k and l < j:
                        if V <= best_dGI:
                            # found an optimal I-loop
                            if DEBUG_I:
                                print "(%2d,%2d) -> (%2d,%2d)[%s] %8.3f >= %8.3f, " \
                                    % (i, j, k, l, ctp, V, best_dGI), cp
                            #
                            flag_hitI = True
                            if V == best_dGI:
                                best_cpI += cp
                            else:
                                best_cpI  = cp
                            best_dGI = V
                            best_ctpI = ctp
                        #
                    #
                #
            #
        #
        
        if flag_hitI:
            if DEBUG_I:
                cp = best_cpI
                print "find_best_I-[I]: (%2d,%2d) -> {%s}[%8.3f][%s]" \
                    % (i, j, cp, best_dGI, best_ctpI) 
            #
        else:
            best_dGI = INFINITY
        #
        if flag_hitJ:
            if DEBUG_I:
                cp = best_cpJ
                print "find_best_I-[J]: (%2d,%2d) -> {%s}[%8.3f][%s]" \
                    % (i, j, cp, best_dGJ, best_ctpJ)
            #
        else:
            best_dGJ = INFINITY
        #
        if len(best_cpI) > 1:
            best_cpI = list(set(best_cpI))
        if len(best_cpJ) > 1:
            best_cpJ = list(set(best_cpJ))
        return best_dGI, best_cpI, best_dGJ, best_cpJ
    #                
    
    
    # M-loop calculation: (Note, I would still prefer to make this a
    # separate class, but the best placement presently is here.)
    def find_best_M(self, i, j, lptyp):
        # best_dGMk = 0.0
        best_dGMk = INFINITY
        best_Mk = 0
        if DEBUG_M:
            print "enter find_best_M(%2d,%2d)[%s]" % (i, j, lptyp)
        for k in range(1, j-i-1):
            
            if DEBUG_M:
                print "find_best_M: ij(%2d,%2d)(%2d), (%2d,%2d)[%8.3f]|(%2d,%2d)[%8.3f] " \
                    % (i, j, k, i, i+k, self.dG[i][i+k], i+k+1, j, self.dG[i+k+1][j])
            #
            #if self.dG[i][i+k] < 0.0 and self.dG[i+k+1][j] < 0.0:
            if self.dG[i][i+k] < INFINITY and self.dG[i+k+1][j] < INFINITY:
                V = self.dG[i][i+k] + self.dG[i+k+1][j]
            else:
                continue
            #
            n1 = len(self.smap.glink[i][i+k].lg) 
            n2 = len(self.smap.glink[i+k+1][j].lg)
            # !!!!!! M !!!!!!
            if V < best_dGMk and n1 > 0 and n2 > 0:
                # found some real M-loop
                best_Mk = k
                best_dGMk = V
                #
                if DEBUG_M:
                    if n1 > 0:
                        n1v = [n1, self.smap.glink[i    ][i+k].lg[0].motif[0].get_base()]
                    if n2 > 0:
                        n2v = [n2, self.smap.glink[i+k+1][j  ].lg[0].motif[0].get_base()]
                    print "ij = ", i, j,  ", j-i = ", (j-i), "k = ", k, ", n12v:", n1v, n2v, V
                #
            #
        #
        join = []
        if best_Mk > 0:
            k = best_Mk
            # !!!!!! M !!!!!!
            
            sm1 = self.smap.glink[i][i+k].lg[0].motif[0]
            sm2 = self.smap.glink[i+k+1][j].lg[0].motif[0]
            if DEBUG_M:
                print "find_best_M"
                print "base:   ", sm1.get_base(),     " + ", sm2.get_base()
                print "branch: ", sm1.get_branches(), " + ", sm2.get_branches()
            #
            n = len(sm1.get_base())
            j_1 = sm1.get_base()[n-1][1]
            i_2 = sm2.get_base()[0][0]
            # print j_1, i_2
            if j_1 > i_2:
                print "error!! find_best_M"
                print "(%d,%d)|(%d,%d), k = %d" % (i, i+k, i+k+1, j, best_Mk)
                sys.exit(1)
            #
            
            join = []
            for jk in self.smap.glink[i][i+k].lg[0].motif[0].get_base():
                join += [jk]
            for jk in self.smap.glink[i+k+1][j].lg[0].motif[0].get_base():
                join += [jk]
            
            # 160719wkd: I know this may seem confusing, why didn't I
            # include the branches here. Well, because
            # Motif.get_bases() already knows the type and therefore
            # can automatically decide if the output should correspond
            # to an I-loop or a MBL, or if it should be a pointer J or
            # P.
            if DEBUG_M:
                print "join2: ", lptyp, join
                sm1 = self.smap.glink[i][i+k].lg[0].motif[0]
                sm2 = self.smap.glink[i+k+1][j].lg[0].motif[0]
                tp1 = sm1.ctp
                tp2 = sm2.ctp
                print "---------->: ij(%2d,%2d)(%2d)[%s], (%2d,%2d){%s}|(%2d,%2d){%s}[%8.3f] " \
                    % (i, j, k, lptyp, i, i+k, tp1, i+k+1, j, tp2, best_dGMk), join
                print "base:   ", sm1.get_base(),     " + ", sm2.get_base()
                print "branch: ", sm1.get_branches(), " + ", sm2.get_branches()
            #
            
        else:
            best_dGMk = INFINITY
        #
        return best_dGMk, [join]
    #                

    # if the BIKMSW motif has the best FE of the list, then it is added to
    # the list as a potential lookup item.
    def save_best_iloops(self, link):
        flag_debug = False
        n = len(self.iloop.lg)
        if flag_debug:
            loop = link.motif[0]
            print "save_best_iloops: dG(%2d,%2d)[%s] = %8.3f >= %8.3f, %d" \
                % (loop.i, loop.j, loop.ctp, loop.Vij, self.V_iloop_mx, n)
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

    def check_cpX_list(self, i, j, cp_list_n, cp_list_o, tp):
        flag_error = False
        flag_len_error  = False
        flag_loop_error = False 
        len_n = len(cp_list_n)
        len_o = len(cp_list_o)
        if not len_n == len_o:
            flag_len_error = True
        #
        
        for lpo in cp_list_o:
            cp_oI = lpo
            flag_loop_error = True
            
            for lpn in cp_list_n:
                cp_nI = lpn
                
                if DEBUG_I:
                    print tp, 'old: ', cp_oI[0], 'new: ', cp_nI[0]
                #
                
                if cp_nI == cp_oI:
                    flag_loop_error = False
                    break
                #
            #
        #
        if flag_len_error or flag_loop_error:
            if flag_len_error:
                print "ERROR: link lists for I-loops do not match"
            #
            if flag_loop_error:
                print "ERROR: the lists of loops do not match"
            #
            print "new N^2 approach, links(%3d,%3d): " % (i, j), cp_list_n
            print "old N^4 approach, links(%3d,%3d): " % (i, j), cp_list_o
            #
            flag_error = True
            sys.exit(1)
        #
        return flag_error
    #
    
    def test_find_best_I(self, i, j, where):
        check_find_best_I = CHECK_ILOOP
        Vij_I = INFINITY
        cp_I  = None
        Vij_J = INFINITY
        cp_J  = None
        
        if check_find_best_I:
            flag_errorI = True
            flag_errorJ = True
            Vij_nI, cp_nI, Vij_nJ, cp_nJ = self.new_find_best_I(i, j)                        
            Vij_oI, cp_oI, Vij_oJ, cp_oJ = self.find_best_I(i, j)
            # print i, j, Vij_nI, Vij_oI, Vij_nJ, Vij_oJ
            if Vij_nI == Vij_oI:
                if DEBUG_I:
                    if Vij_nI < 0.0:
                        print (i, j), Vij_nI, Vij_oI
                    #
                #
                flag_errorI = self.check_cpX_list(i, j, cp_nI, cp_oI, 'I')
            if Vij_nJ == Vij_oJ:
                if DEBUG_I:
                    if Vij_nJ < 0.0:
                        print (i, j), Vij_nJ, Vij_oJ
                    #
                #
                flag_errorJ = self.check_cpX_list(i, j, cp_nJ, cp_oJ, 'J')
            #
            
            if not (flag_errorI or flag_errorJ):
                Vij_I = Vij_nI
                cp_I  = cp_nI
                Vij_J = Vij_nJ
                cp_J  = cp_nJ
            elif Vij_nI == Vij_oI and Vij_nJ == Vij_oJ:
                # on rare occasions, there can be problems with degeneracy.
                print "WARNING(minFE): find_best_I() incomplete match! ... in %s" % (where)
                
                print "new (N^2): dG(%2d,%2d) = %8.3f------------" % (i, j, Vij_nI)
                print cp_nI
                #
                print "           dG(%2d,%2d) = %8.3f" % (i, j, Vij_nJ)
                print cp_nJ
                #
                print "old (N^4): dG(%2d,%2d) = %8.3f------------" % (i, j, Vij_oI)
                print cp_oI
                #
                print "           dG(%2d,%2d) = %8.3f" % (i, j, Vij_oJ)
                print cp_oJ
                #
                
                
                # well, assign anyway????
                Vij_I = Vij_nI
                cp_I  = cp_nI
                Vij_J = Vij_nJ
                cp_J  = cp_nJ
                sys.exit(1) # presently also exit
                
                
            else:
                # if even the energies do not match, this is definitely trouble
                print "PROBLEMS(minFE): find_best_I() not matching! ... in %s" % (where)
                
                print "new (N^2): dG(%2d,%2d)[I] = %8.3f------------" % (i, j, Vij_nI)
                print cp_nI
                #
                print "           dG(%2d,%2d)[J] = %8.3f" % (i, j, Vij_nJ)
                print cp_nJ
                #
                print "old (N^4): dG(%2d,%2d)[I] = %8.3f------------" % (i, j, Vij_oI)
                print cp_oI
                #
                print "old        dG(%2d,%2d)[J] = %8.3f" % (i, j, Vij_oJ)
                print cp_oJ
                #
                print "terminating run, please check"
                sys.exit(1)
            #
        #
        else:
            Vij_I, cp_I, Vij_J, cp_J = self.new_find_best_I(i, j)
            #                                   Vij_J is the I-loop without the bond at ij
        #
        return Vij_I, cp_I, Vij_J, cp_J
    #

    def get_bondtype(self, hvij):
        bt = 's'
        if hvij > self.ctcf_cthresh:
            bt = 'c' # default "tandem" setting
        elif hvij > self.ctcf_tthresh:
            bt = 't' # default "convergent" setting
        #
        else:
            bt = 's'
        return bt
    #
    
    # Main calculation script for finding the minimum free energy and
    # (in the future) suboptimal structures.
    def minFE(self, T):
        self.T = T
        dGij  = 0.0
        Vij_I = 0.0; cp_I = None; Vij_J = 0.0; cp_J = None
        Vij_M = 0.0; cp_M = None; Vij_P = 0.0; cp_P = None
        for j in range(1, self.N):
            for i in range(j-1, -1, -1):
                
                # single base
                hvij = self.hv[i][j]
                if hvij > 0.0:
                    bondtype = self.get_bondtype(hvij)
                    dGij = self.fe.calc_dG(i, j, self.hv[i][j], T)
                    
                    if DEBUG_minFE:
                        print "assigned B: (%2d,%2d)[%8.3f]" % (i, j, dGij), [(i,j)]
                    #
                    
                    # ######  since 'B' _exists_, must be recorded  ######
                    link_H = Link(i, j, dGij, 'B', bondtype, [(i,j)])
                    self.smap.glink[i][j].add_link(link_H)
                    self.save_best_iloops(link_H)
                    
                    stem_aa, stem_pp = self.find_ifStem(i, j, dGij, bondtype)
                    # print stem_aa, len(stem_aa), stem_pp, len(stem_pp)
                    if len(stem_aa) > 0:
                        self.make_StemMotif(stem_aa)
                    if len(stem_pp) > 0:
                        self.make_StemMotif(stem_pp)
                        
                    
                    # ######  calculate interactions for I, J, M and P  ######
                    
                    # search for the best I,J
                    Vij_I, cp_I, Vij_J, cp_J = self.test_find_best_I(i, j, "ij bound search")
                    # Vij_I, cp_I, Vij_J, cp_J = self.new_find_best_I(i, j)
                    Vij_I += dGij
                    
                    
                    # search for the best M,P
                    Vij_M, cp_M = self.find_best_M(i+1, j-1, 'M') 
                    Vij_M += dGij
                    Vij_P, cp_P = self.find_best_M(i  , j  , 'P')
                    #
                    
                    
                    # ######################################3
                    if MONITOR or DEBUG_minFE:
                        self.monitor_p1m1("M: Vij_M", i, j, Vij_M-dGij)
                    # ######################################3
                    
                    if DEBUG_minFE:
                        print "(%2d,%2d): B(%8.3f), I(%8.3f), M(%8.3f), J(%8.3f), P(%8.3f)" \
                            % (i, j, dGij, Vij_I, Vij_M, Vij_J, Vij_P)
                    #
                    
                    
                    # ######  collect up results for I, J, M and P  ######
                    
                    # it is possible the J exists but I does not.
                    if Vij_J < 0.0:
                        link_J = Link(i, j, Vij_J, 'J', '-', [cp_J[0]])
                        if len(cp_J) > 1:
                            # case for degeneracy
                            for nn in range(1, len(cp_J)):
                                link_J.add_Motif(i, j, Vij_J, 'J', '-', [cp_J[nn]])
                        self.smap.glink[i][j].add_link(link_J)
                        
                        if DEBUG_minFE:
                            print "assigned J: (%2d,%2d)[%8.3f]" % (i, j, Vij_J), cp_J
                        #
                    #
                    
                    # it is possible the J exists but I does not.
                    if Vij_I < 0.0:
                        link_I = Link(i, j, Vij_I, 'I', bondtype, [cp_I[0]])
                        if len(cp_I) > 1:
                            # case for degeneracy
                            for nn in range(1, len(cp_I)):
                                link_I.add_Motif(i, j, Vij_I, 'I', bondtype, [cp_I[nn]])
                        self.smap.glink[i][j].add_link(link_I)
                        
                        if DEBUG_minFE:
                            print "assigned I: (%2d,%2d)[%8.3f]" % (i, j, Vij_I), cp_I
                        #
                        self.save_best_iloops(link_I)
                    #
                    
                    if Vij_M < 0.0:
                        # a multibranch loop exists
                        link_M = Link(i, j, Vij_M, 'M', bondtype, cp_M[0])
                        if len(cp_M) > 1:
                            # case for degeneracy (though currently not used)
                            for nn in range(1, len(cp_M)):
                                link_M.add_Motif(i, j, Vij_M, 'M', bondtype, cp_M[nn])
                        self.smap.glink[i][j].add_link(link_M)
                        
                        if DEBUG_minFE:
                            print "assigned M: (%2d,%2d)[%8.3f]" % (i, j, Vij_M), cp_M
                        #
                        self.save_best_iloops(link_M)
                    #
                    
                    # it is possible the P exists but M does not.
                    if Vij_P < 0.0:
                        link_P = Link(i, j, Vij_P, 'P', '-', cp_P[0])
                        if len(cp_P) > 1:
                            # case for degeneracy (though currently not used)
                            for nn in range(1, len(cp_P)):
                                link_P.add_Motif(i, j, Vij_P, 'P', '-', cp_P[nn])
                        self.smap.glink[i][j].add_link(link_P)
                        if DEBUG_minFE:
                            print "assigned P: (%2d,%2d)[%8.3f]" % (i, j, Vij_P), cp_P
                        #
                    #
                    
                    
                    # ######  sort the data in terms of free energy  #####
                    self.smap.mergeSortLinks(self.smap.glink[i][j].lg)
                    dGbest = self.smap.glink[i][j].lg[0].Vij # !!!!!
                    self.dG[i][j] = dGbest
                    
                    # up to here, we are looking at regions within
                    # (i,j). Now we look at PKs, and we have to expand
                    # the search outside this region. So dGbest is the
                    # best results for a closing point (i,j).
                    
                    
                    # pseudoknot and CTCF-island solutions
                    pkaa, pkpp = self.find_best_PK(i, j)
                    
                    flag_debug_K = False
                    if not len(pkaa) < 1:
                        i_pk = pkaa[0][0]; j_pk = pkaa[0][1]; dGpk = pkaa[1]
                        link_K = Link(i_pk, j_pk, dGpk, 'K', 'sa', [(i,j)], pkaa[2])
                        self.smap.glink[i_pk][j_pk].add_link(link_K)
                        
                        if flag_debug_K:
                            L = self.scan_ahead
                            if j + L >= self.N:
                                L = self.N - j - 1
                            #
                            print "pkaa: ", pkaa
                            print "ij_pk: ", i_pk, j_pk
                            print "ij:    ", i, j, " j + L: ", (j + L)
                            print "dGpk(%8.2f)   dGbest(%8.2f)   ddGpk(%8.2f)" \
                                % (dGpk, dGbest, (dGpk - dGbest))
                        #
                    #
                    if not len(pkpp) < 1:
                        i_pk = pkpp[0][0]; j_pk = pkpp[0][1]; dGpk = pkpp[1]
                        link_K = Link(i_pk, j_pk, dGpk, 'K', 'sp', [(i,j)], pkpp[2])
                        self.smap.glink[i_pk][j_pk].add_link(link_K)
                        
                        if flag_debug_K:
                            L = self.scan_ahead
                            if j + L >= self.N:
                                L = self.N - j - 1
                            #
                            print "pkpp: ", pkpp
                            print "ij_pk: ", i_pk, j_pk
                            print "ij:    ", i, j, " j + L: ", (j + L)
                            print "dGpk(%8.2f)   dGbest(%8.2f)   ddGpk(%8.2f)" \
                                % (dGpk, dGbest, (dGpk - dGbest))
                        #
                    #
                    nlen = 2
                    if flag_debug_K and (len(pkaa) > nlen or len(pkpp) > nlen): # still checking this
                        print ">>> after find_best_PK(), when len(pkaa) and/or len(pkpp) > %d" % nlen
                        if STOP_AT_FIND:
                            print "planned stop"
                            sys.exit(0)
                        #
                    #
                    
                    if self.all_ctcf.has_key((i,j)):
                        # if it doesn't have the key, then forget it!
                        island = self.find_ctcf_islands(i, j, dGbest)                    
                    
                        if len(island) > 0:
                            # if there's nothing there, then don't bother with it!
                            for ff in island:
                                wyspa = ff[0]; join = ff[1]; dGW = ff[2] 
                                link_W = Link(i, j, dGW, 'W', 'wyspa', join, [], wyspa)
                                self.smap.glink[i][j].add_link(link_W)
                            #
                            
                            # self.traceback_mFE(i,j, 0, True)
                        #
                    #
                    self.smap.mergeSortLinks(self.smap.glink[i][j].lg)
                    self.dG[i][j] = self.smap.glink[i][j].lg[0].Vij # !!!!!
                    
                else:
                    Vij_I, cp_I, Vij_J, cp_J = self.test_find_best_I(i, j, "ij unbound search")
                    # Vij_I, cp_I, Vij_J, cp_J = self.new_find_best_I(i, j)
                    Vij_P, cp_P = self.find_best_M(i, j, 'P')
                    #
                    
                    # push down, first J
                    if Vij_J < INFINITY:
                        link_J = Link(i, j, Vij_J, 'J', '-', [cp_J[0]])
                        if len(cp_J) > 1:
                            # case where there is degeneracy
                            for nn in range(1, len(cp_J)):
                                link_J.add_Motif(i, j, Vij_J, 'J', '-', [cp_J[nn]])
                        self.smap.glink[i][j].add_link(link_J)
                    else:
                        self.dG[i][j] = INFINITY
                    #
                    
                    # push down, best value is 'P'
                    if Vij_P < 0.0:
                        self.dG[i][j] = Vij_P
                        link_P = Link(i, j, Vij_P, 'P', '-', cp_P[0])
                        if len(cp_P) > 1:
                            # case where there is degeneracy
                            for nn in range(1, len(cp_P)):
                                link_P.add_Motif(i, j, Vij_P, 'P', '-', cp_P[nn])
                        self.smap.glink[i][j].add_link(link_P)
                    else:
                        self.dG[i][j] = INFINITY
                    #
                    
                    if len(self.smap.glink[i][j].lg) > 0:
                        self.smap.mergeSortLinks(self.smap.glink[i][j].lg)
                        self.dG[i][j] = self.smap.glink[i][j].lg[0].Vij # !!!!!
                    #
                    
                    if DEBUG_minFE:
                        if len(self.smap.glink[i][j].lg) > 0:
                            
                            if   Vij_J < Vij_P:
                                tg = 'J'
                                print "J rules!"
                                self.show_Vs([Vij_J, Vij_P])
                            else:
                                tg = 'P'
                                print "P rules!"
                                self.show_Vs([Vij_J, Vij_P])
                            #
                            for m in range(0, len(self.smap.glink[i][j].lg)):
                                for vk in self.smap.glink[i][j].lg[m].motif:
                                    tp = vk.ctp
                                    cp = vk.get_branches()
                                    Vp = self.smap.glink[i][j].lg[m].Vij
                                    print "assigned %d  dG(%2d,%2d)[%s][%8.3f] -> %s" % (m, i, j, tp, Vp, cp)
                            #
                        else:
                            # didn't write anything
                            if DEBUG_minFE:
                                print "assigned dG(%2d,%2d)  empty" % (i, j)
                            #
                        #
                    #
                    
                    
                    
                    # ######################################3
                    def trans_labels(arg):
                        # a trick to achieve switch/case in C/C++
                        switcher = {
                            'J': ("J: Vij_J", Vij_J),
                            'P': ("P: Vij_P", Vij_P),
                        }
                        return switcher.get(arg, "nothing")
                    if MONITOR:
                        cc = trans_labels(tg)
                        if not (cc == "nothing"):
                            self.monitor_p0m0(cc[0], i, j, cc[1])
                        else:
                            print "ERROR!"
                            sys.exit(1)
                    #
                    # ######################################3
                    
                    #
                #
                if CHECK_glink:
                    # self.check_smap(i, j, 38, 49, True)
                    self.check_smap(i, j, 27, 70, True)
                #
            #
        #
        if DEBUG_minFE:
            print "finished minFE"
            # sys.exit(0)
        return self.dG, self.smap
    #
    
    def show_Vs(self, Vlist):
        s = ''
        label = ['Vij_J', 'Vij_P']
        for k in range(0, len(Vlist)):
            s += "%s: %8.3f\n" % (label[k], Vlist[k])
        print s
    #
    
    
    
    
    
    
    # #####################################################
    # tools for checking the variables for self.dG, etc.
    # #####################################################
    # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    
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
                    tp = vvkn.ctp # conn type
                    bt = vvkn.btp # bond type
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
                        tp   = ij_jnk.ctp # conn type
                        bt   = ij_jnk.btp # bond type
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
    
    def monitor_p1m1(self, label, i, j, test_V):
        # test_V = .... test value
        
        xx = len(self.smap.glink[i+1][j-1].lg)
        if xx > 0:
            for ll in range(0, len(self.smap.glink[i+1][j-1].lg)):
                v = self.smap.glink[i+1][j-1].lg[ll]
                ctp = ''; base = []
                for vk in v.motif:
                    ctp += vk.ctp
                    base += vk.get_base()
                if ll == 0:
                    print "%s(%2d,%2d) = %8.3f, Vp1m1 = %8.3f: V[%s] = %8.3f " \
                        % (label, i, j, (test_V), self.dG[i+1][j-1], ctp, v.Vij), base
                else:
                    print "                                              V[%s] = %8.3f " \
                        % (ctp, v.Vij), base
                #
            #
        else:
            print "%s(%2d,%2d) = %8.3f, Vp1m1 = %8.3f [-]%d" \
                % (label, i, j, (test_V), self.dG[i+1][j-1], xx)
        # 
    #
    
    def monitor_p0m0(self, label, i, j, test_V):
        
        # test_V = .... test value
        xx = len(self.smap.glink[i][j].lg)
        if xx > 0:
            for ll in range(0, len(self.smap.glink[i][j].lg)):
                v = self.smap.glink[i][j].lg[ll]
                if ll == 0:
                    print "%s(%2d,%2d) = %8.3f, Vrslt = %8.3f: V[%s] = %8.3f " \
                        % (label, i, j, (test_V), self.dG[i][j], v.ctp, v.Vij), v.motif[0].get_base()
                else:
                    print "                                              V[%s] = %8.3f " \
                        % (v.ctp, v.Vij), v.motif[0].get_base()
                #
            #
        else:
            print "%s(%2d,%2d) = %8.3f, Vp1m1 = %8.3f [-]%d" \
                % (label, i, j, (test_V), self.dG[i][j], xx)
        # 
    #
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # #####################################################
    # tools for checking the variables for self.dG, etc.
    # #####################################################
    
    
    
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
        flag_debug = False
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
            
        ctp  = self.smap.glink[i][j].lg[0].motif[0].ctp
        btp  = self.smap.glink[i][j].lg[0].motif[0].btp
        V    = self.smap.glink[i][j].lg[0].Vij
        join = self.smap.glink[i][j].lg[0].motif[0].get_branches()
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
            flag_debug_S = False
            if flag_debug_S:
                print "inside traceback_mFE():"
                print "stem Motif: ", self.smap.glink[i][j].lg[0].motif[0].show_Motif()
            #
            
            btp = self.smap.glink[i][j].lg[0].motif[0].btp
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
                    if not self.smap.glink[i_h][j_h].lg[0].motif[0].ctp == 'b':
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
            flag_debug_K = False
            if flag_debug_K:
                print "tracelink_mFE(), K:"
            #
            # the root stem
            vR    = self.smap.glink[i][j].lg[0].motif[0].get_base()
            if flag_debug_K:
                print vR
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
            if flag_debug_K:
                print "pks: ", pks
            for x in pks:
                i_K = x[0][0]; j_K = x[0][1]
                if flag_debug_K:
                    print "ijK: ", i_K, j_K
                self.opt_ss_seq[i_K] = '['
                self.opt_ss_seq[j_K] = ']'
                if show_structure:
                    print "%s[%s](%3d, %3d)[%8.3f]: " % (s, "l", i_K,j_K, x[1])
            #
            if show_structure:
                print "%s---" % s
            if flag_debug_K:
                print "pk results: from traceback_mFE()" 
                print ''.join(self.opt_ss_seq)
                sys.exit(0)
        #
        elif ctp == 'W':
            flag_debug_W = False
            # island
            wyspa = self.smap.glink[i][j].lg[0].motif[0].get_wyspa()
            if flag_debug_W:
                print "wyspa: ", wyspa
            for wk in wyspa:
                # some repetition, but who cares!
                
                self.opt_ss_seq[wk[0]] = '|'
                self.opt_ss_seq[wk[1]] = '|'
            #
            
            # island connections
            jn_w = self.smap.glink[i][j].lg[0].motif[0].get_branches()
            if flag_debug_W:
                print "jn_w: ", jn_w
            for lk in jn_w:
                # whether 1 or 100 elements, they can be processed this way.
                i_W = lk[0]; j_W = lk[1]     # new ij of this subdomain
                if len(self.smap.glink[i_W][j_W].lg) > 0:
                    self.traceback_mFE(i_W, j_W, layer+1, show_structure)
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
    
    def find_ctcf_islands(self, i, j, best_dG):
        flag_debug = False

        if flag_debug:
            print "find_ctcf_islands(%d,%d), best_dG = %8.2f" % (i,j, best_dG)
        
        # 160915wkd: CURRENT STRATEGY
        
        # If I compute every case, I have i2,i3,... iN possible
        # binding points. This translates into N individual, and
        # N!/(r!(N-r)!) subsets of tuples: the "nCr" function on the
        # Casio. This is a lot of cases as N grows!
        
        # We have claimed already that we form these CTCF islands, we
        # show them as being a domain and their frequency is much
        # greater than the singletons.
        
        # Therefore, for the moment, I will work with the heuristic
        # that any region that predicts a group of CTCF interactions
        # at point i or point j and within (i,j) has suffient free
        # energy to bind to the respective i or j all together.

        # Since these singletons can extend from either the i-side or
        # the j-side of the _largest_ domain, and both can contribute
        # significantly to the stability of the domain, I think the
        # only thing we can do is just find all the CTCF connections
        # (k,l) that satisfy i == k OR j == l and split up the
        # structure into these sub-regions. In general, the CTCF
        # connections are significantly stronger than any of the
        # singleton connections.
        
        # If we really want to know all the possible free energy
        # results from all possible patterns (in all the various
        # tuples and alone), then we can add some additional
        # functionality later.
        
        # At least for the relatively simple cases, this will work
        # with no problems. It is also conceivable that we can test
        # the solutions directly at some point. At any rate, since the
        # free energy of binding the CTCFs is typically so much
        # greater than the singletons, I think it is sufficient to
        # simply ASSume that they all bind, presently.
        
        
        keylist = self.all_ctcf.keys()
        vrlps = [] # oVeRLaPS
        zone  = []
        joinW = []
        wyspa = []
        islands = []
        flag_found = False
        
        # first go through the keylist and find common contacts
        for h in keylist:
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
            return islands
        #
        
        wyspa += [(i,j)]
        
        
        # "list(set(vrlps)).sort()" sorts and eliminates redundant
        # items and it doeesn't matter if it is empty
        vrlps = list(set(vrlps))
        vrlps.sort()
        #
        
        if flag_debug:
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
        
        if flag_debug:
            print "zone: ", zone
        #
        
        dGijh  = self.fe.calc_dG(i,  j,  self.hv[i ][j ], self.T, "find_ctcf_islands 1")
        dGh = INFINITY
        
        # now go through the list of zones and compute the singleton
        # contributions
        
        #                  common point (can happen)
        #                       V
        #                         ..................        
        #                       |           ........|
        #                       v          v        v
        #     |..................          |        |       free energy
        #     |........         |          |        |
        #     v       v         v          |        |
        #   ..|.......|.........|..........|........|...
        #     i     ijh1      ijh2       ijh3       j
        #      ------    ------    -------    -----
        #
        #         <-------  loop regions ------>
        
        # naturally, if there are interactions between (jh1,jh2),
        # (jh2,jh3), etc., these would also show up in the growth of
        # the structure
        
        if len(zone) > 0:
            # you would have to calculate L and R separately anyway
            dGh = 0.0 
            for ijh in zone:
                ih = ijh[0]; jh = ijh[1]
                if flag_debug:
                    print "ijh              i(%3d) <- ih(%3d) <- jh(%3d) <- j(%3d): " % (i, ih, jh, j)
                #
                # case where we have an internal CTCF between the 
                if self.all_ctcf.has_key((ih,jh)) and i < ih and jh < j:
                    # need to save this result for the next step
                    wyspa += [(ih,jh)]
                    if flag_debug:
                        print "internal island: i(%3d) <  ih(%3d) <  jh(%3d) <  j(%3d)" % (i, ih, jh, j)
                    #
                #
                if len(self.smap.glink[ih +1][jh-1].lg) > 0:
                    dGh += self.smap.glink[ih+1][jh-1].lg[0].motif[0].Vij
                #
                joinW += [(ih+1,jh-1)]
                
                # 160923wkd: Even if there is nothing inside, this has
                # to be saved so that other parts of the program know
                # what to do with the domains. The problem was that I
                # had two islands
                
                # given: wyspa  = [(38, 43), (40, 43)]
                #    =>   zone  = [(38, 40), (40, 43)]
                #    =>  joinW -> [(39, 39), (41, 42)]
                
                # Although the insides are obviously empty, if you
                # don't indicate this relationship for joinW, then the
                # program just keeps cycling through (38,43) as the
                # reference until the recursion limits kill it.
                
            #
        #
        wyspa.sort()
        if len(wyspa) > 1:
            for ijW in wyspa:
                i_W = ijW[0]; j_W = ijW[1]
                dGh += self.fe.calc_dG(i_W, j_W, self.hv[i_W][j_W], self.T, "find_ctcf_islands 2")
            #
        #
        
        if flag_debug:
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
            if flag_debug:
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
                    print "number of CTCFs > %d, stopping the program for reference" % nzones
                    sys.exit(0)
                #
            #
            #
        #
        #
        
        return islands
    #
    
    def make_StemMotif(self, stem):
        flag_debug = False
        if flag_debug:
            print "make_StemMotif()"
            print "stem: ", stem
        i = stem[0][0]; j = stem[0][1]; dGij = stem[1]
        
        join = []
        btp_special = ''
        btp_general = ''
        for ss in stem[2]:
            join += [ss[0]]
            if ss[2] == 'c' or ss[2] == 'l' or ss[2] == 'r' or ss[2] == 't':
                btp_special = ss[2]
            else:
                btp_general = ss[2]
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
        stmlen = len(join)
        i_t = i; j_t = j
        i_h = i + stmlen - 1; j_h = j - stmlen + 1
        ihh = i_h + 1;        jhh = j_h - 1
        
        
        # I know, it is kind of strange, but whether the stem is
        # parallel or antiparallel, the computation of (i_t,j_t) --
        # the tail -- and (i_h,j_h) -- the head -- is the same.

        # From the diagram below, perhaps one can understood that the
        # entry points and exit points are in the same position on the
        # chain, it is just that the order of the pairing arrangement
        # is changed.
        
        #
        #   antiparallel                parallel
        #     _____                
        #    /     \                            ______
        #    |     |                    j_h    /
        #    |     |                   __._ _./ j_t
        #    |     |                  / _|_|_|___
        #     \   /               ____|/i_t i_h  \
        # i_h  \_/ j_h                |           |
        #      |_|                    |___________|
        # i_t  |_| j_t   
        #                

        # the main difference is that the head of the antiparallel
        # structure can encorporate a motif of type 'B', 'I' or 'M',
        # but the parallel stem motif cannot. This is because the
        # closing points in the parallel stem are at disparate
        # positions (separated by a distance (i_h - i_t) or (j_t -
        # j_h). On the other hand, the closing point (i_h, j_h) on the
        # antiparallel stem is exactly a point where an 'M' motif can
        # attach.
        
        
        if flag_debug:
            print "make_StemMotif(): parallel stem"
            print "pp ij_t (tail): ", i_t, j_t
            print "pp ij_h (head): ", i_h, j_h
            print "pp ijhh (jnct): ", ihh, jhh, ", dijhh: ", (jhh - ihh)
            print "region does not have a defined terminus"
            if i_t == 0 and j_t == 14:
                print "btp(%d,%d) = %s" % (i, j, btp)
                for m in range(0, len(self.smap.glink[i_t][j_t].lg)):
                    print self.smap.glink[i_t][j_t].lg[m].motif[0].show_Motif()
                print "get_bondtype: ", self.get_bondtype(self.hv[i][j])
                print "all_ctcf:     ", self.all_ctcf
                sys.exit(0)
                #
            # sys.exit(0)
        #

        link_S = Link(i_t, j_t, dGij, 'S', btp, join)
        if flag_debug:
            print link_S.motif[0].show_Motif()
            #if btp == 'sp':
            #    sys.exit(0)
        self.smap.glink[i_t][j_t].add_link(link_S)
        
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
    
    
    # see if there is a Motif of type Stem at (i,j)
    def find_ifStem(self, i, j, dGij, btp):
        flag_debug = False 
        pairs_aa = []
        stem_aa = []
        pairs_pp = []
        stem_pp = []

        if flag_debug:
            print "ij:  ", i, j
        
        
        
        # search for an antiparallel connection
        flag_pairs = True
        k = 1
        while flag_pairs:
            iaa = i + k; jaa = j - k
            hvaa = self.hv[iaa][jaa]
            if hvaa > 0.0 and jaa - iaa > 0:
                btpx = self.get_bondtype(hvaa)
                if btpx == 's':
                    btpx = 'sa'
                else:
                    # have to enforce this, even if the prediction is
                    # 't' from get_bondtype()!
                    btpx = 'c'
                #
                dGijaa  = self.fe.calc_dG(iaa,  jaa, hvaa, self.T, "find_ifStem 1")
                pairs_aa += [((iaa, jaa), dGijaa, btpx)]
                if flag_debug:
                    print "ijaa: ", iaa, jaa
                k += 1
            else:
                flag_pairs = False
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
            
            
            # Compute the total free energy of the motif.  In the case
            # of antiparallel stems, they are joined by an
            # antiparallel connector like 'B', 'I', 'K', 'M' or
            # 'W'. Therefore, the _last item_ on the list should NOT
            # be computed from the free energy of a closing loop, as
            # this will result in double counting the closing region.
            dGaa = 0.0
            for k in range(0, len(pairs_aa)-1):
                dGaa += pairs_aa[k][1]
                if flag_debug:
                    print pairs_aa[k][0], dGaa
                #
            kstm = len(pairs_aa) - 1
            iaa = i + kstm; jaa = j - kstm
            if jaa-iaa > 0:
                V, cp, ctp, pks, wps = self.filter_BIKMSW_from_glink(iaa, jaa)
                dGaa += V
                if flag_debug:
                    print "V: ", V
                #
            #
            stem_aa = [(i, j), dGaa, pairs_aa]
            if flag_debug:
                i_t = i;   j_t = j
                i_h = iaa; j_h = jaa
                print "aa ij_t (tail): ", i_t, j_t, dGaa
                print "aa ij_h (tail): ", i_h, j_h
                print "aa ijaa:        ", iaa, jaa
                print stem_aa
                #print "planned exit"
                #sys.exit(0)
        #
        
        # search for an antiparallel connection
        flag_pairs = True
        k = 1
        i_t = i;     j_t = j
        i_h = i + k; j_h = j_t - k
        while flag_pairs:
            # we only search backwards with this because we go over
            # old solutions.
            i_h = i + k; j_h = j_t - k
            ipp = i - k; jpp = j - k
            hvpp = self.hv[ipp][jpp]
            if hvpp > 0.0 and (j_h - i_h) > 0 and ipp >= 0:
                btpx = self.get_bondtype(hvpp)
                if btpx == 's':
                    btpx = 'sp'
                else:
                    # have to enforce this, even if the prediction is
                    # 'c' from get_bondtype()!
                    btpx = 't'
                dGijpp  = self.fe.calc_dG(ipp,  jpp, hvpp, self.T, "find_ifStem 2")
                pairs_pp += [((ipp, jpp), dGijpp, btpx)]
                if flag_debug:
                    print "ijpp: ", ipp, jpp
                k += 1
            else:
                flag_pairs = False
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
            
            i_t = i - stmlenpp + 1; j_t = j
            i_h = i;                j_h = j - stmlenpp + 1
            ihh = i_h + 1;          jhh = j_h - 1
            
            if flag_debug:
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
                if i_t == 0 and j_t == 14:
                    print "btp(%d,%d) = %s" % (i, j, btp)
                    for m in range(0, len(self.smap.glink[i_t][j_t].lg)):
                        print self.smap.glink[i_t][j_t].lg[m].motif[0].show_Motif()
                    print "get_bondtype: ", self.get_bondtype(self.hv[i][j])
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
            
            
            #if len(pairs_pp) > 1:
            #    print "planned exit"
            #    sys.exit(0)
            # appears to sort in ascending order when confronted with this.
            dGpp = 0.0
            for pp in pairs_pp:
                if flag_debug:
                    print pp[1]
                #
                dGpp += pp[1]
            #
            
            # unlike the antiparallel stem, this one has to be ground
            # all the way through.
            
            if  len(self.smap.glink[ihh][jhh].lg) > 0:
                dGpp += self.smap.glink[ihh][jhh].lg[0].Vij
            else:
                if flag_debug:
                    lbound = 6
                    if jhh - ihh > lbound:
                        print "d(jhh - ihh) = %d zone > %d" % ((jhh - ihh), lbound)
                        sys.exit(0)
                    #
                #
            #
            #
            stem_pp = [(i_t, j_t), dGpp, pairs_pp]
            if flag_debug:
                print "pp ij_t (tail): ", i_t, j_t, dGpp
                print "pp ij_h (head): ", i_h, j_h              
                print stem_pp
                #if len(pairs_pp) > 2:
                #    print "pp planned exit"
                #    sys.exit(0)
        #
        #
        #print pairs_aa, pairs_pp
        return stem_aa, stem_pp
    #
    
    def find_best_PK(self, i, j):
        flag_debug_PK = False
        L = self.scan_ahead
        if j + L >= self.N:
            L = self.N - j - 1
        
        dGmin = self.traceback_mFE(i, j, 0)
        # it would probably be faster to use a modified version of
        # traceback to look for interaction points, but it is also a
        # bit more complicated.
        ssv = self.opt_ss_seq
        # provide a list of open locations for pk binding
        
        if flag_debug_PK:
            print "find_best_PK(%d,%d), leading edge(L=%d):  %d" % (i, j, L, j+L)
            print ''.join(self.opt_ss_seq)
        #
        
        best_dG = INFINITY
        best_ijp = ()
        
        best_Spklistaa = []
        last_iaa = i
        best_Spklistpp = []
        last_ipp = j
        flag_find = False
        best_dGaa = 1000
        best_dGpp = 1000

        
        # first, scan to find the best attachment point
        for jp in range(j+L, j, -1):
            # this finds the best hit on each scan
            
            # every time, we scan across j we look for the best result on that scan across i
            
            best_ddGaa = 1000
            best_ddGpp = 1000
            best_daa = ()
            best_dpp = ()
            
            # find best PK along iaa for give jp
            for ip in range(i+1, j):
                iaa = ip
                ipp = j - iaa + i
                #!! print "ijpp = (%3d,%3d), ijaa = (%3d,%3d)" % (ipp, jp, iaa, jp)
                # choose something significant!!!
                hvf = self.hv[iaa][jp]
                if hvf > self.pk_thresh:
                    # originally was "if hvb > self.pk_thresh and hvb < self.ctcf_tthresh:"
                    dGijp = self.fe.calc_dG(iaa, jp, hvf, self.T, "find_best_PK 1")
                    if ssv[iaa] == '.' and ssv[jp] == '.':
                        dGaa_cur = 0.0
                        for aa in best_Spklistaa:
                            dGaa_cur += aa[1]
                        #
                        
                        # find _a_ best on for a given j
                        if last_iaa < iaa and dGijp < best_ddGaa:
                            #!! print "1. ijp(aa): ", iaa, jp
                            best_ddGaa = dGijp
                            best_daa   = (iaa, jp)
                            #
                        elif dGijp < best_dGaa and dGijp < dGaa_cur:
                            #!! print "2. ijp(aa): ", iaa, jp
                            best_ddGaa = dGijp
                            best_daa   = (iaa, jp)
                            best_Spklistaa = []
                            
                        #
                    #
                #
                hvb = self.hv[ipp][jp]
                if hvb > self.pk_thresh:
                    # originally was "if hvb > self.pk_thresh and hvb < self.ctcf_tthresh:"
                    dGijp = self.fe.calc_dG(ipp, jp, hvb, self.T, "find_best_PK 2")
                    
                    if ssv[ipp] == '.' and ssv[jp] == '.':
                        dGpp_cur = 0.0
                        for pp in best_Spklistpp:
                            dGpp_cur += pp[1]
                        #
                        
                        # find _a_ best on for a given j
                        if last_ipp > ipp and dGijp < best_ddGpp:
                            #!! print "1. ijp(pp): ", ipp, jp
                            best_ddGpp = dGijp
                            best_dpp   = (ipp, jp)
                            #
                        elif dGijp < best_dGpp and dGijp < dGpp_cur:
                            #!! print "2. ijp(pp): ", ipp, jp
                            best_ddGpp = dGijp
                            best_dpp   = (ipp, jp)
                            best_Spklistpp = []
                            
                        #
                    #
                #
            #
            if best_ddGaa < 0.0:
                if last_iaa < iaa:
                    i_aa = best_daa[0]; j_aa = best_daa[1]
                    btp = self.get_bondtype(self.hv[i_aa][j_aa])
                    if btp == 's':
                        btp += 'a'
                    else:
                        btp = 'c' 
                    # When hvb < self.ctcf_tthresh and regardless of
                    # whether the bondtype is 't' or 'c', we have to
                    # force the situation (btp = 'c') because this is
                    # a case wherein the stem is an antiparallel stem
                    
                    best_Spklistaa += [(best_daa, best_ddGaa, btp)] 
                    last_iaa      = i_aa
                    flag_find = True
                    if flag_debug_PK:
                        print "last_iaa; ", last_iaa
                        print "best_daa:  ", best_daa, best_ddGaa, dGijp
                        print "best_Spklistaa: ", best_Spklistaa
                        print "ij:  ", i, j
                    #
                #
            #
            if best_ddGpp < 0.0:
                if last_ipp > ipp:
                    i_pp = best_dpp[0]; j_pp = best_dpp[1]
                    btp = self.get_bondtype(self.hv[i_pp][j_pp])
                    if btp == 's':
                        btp += 'p'
                    else:
                        btp = 't' 
                    # When hvb < self.ctcf_tthresh and regardless of
                    # whether the bondtype is 't' or 'c', we have to
                    # force the situation (btp = 't') because this is
                    # a case wherein the stem is an antiparallel stem
                    best_Spklistpp += [(best_dpp, best_ddGpp, btp)] 
                    last_ipp = i_pp
                    flag_find = True
                    if flag_debug_PK:
                        print "last_ipp; ", last_ipp
                        print "best_dpp:  ", best_dpp, best_ddGpp, dGijp
                        print "best_Spklistpp: ", best_Spklistpp
                        print "ij:  ", i, j
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
                ip = i; jp = 0
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
                pklist_aa = [(ip, jp), dGaa, best_Spklistaa]
                
            #
            if len(best_Spklistpp) > 0:
                ip = i; jp = 0
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
        #print "pklist_aa: ", pklist_aa
        #print "pklist_pp: ", pklist_pp
        # sys.exit(0)
        return pklist_aa, pklist_pp
    #
        
#

####################################################################
########################  BACK TRACK TRACE  ########################
####################################################################



# Trace() does the trace back steps for a particular structure.

# At this point, Calculate() has been run (there should be a check for
# that here!!!). The program runFE() only searches for the structure
# with the minimum free energy. However, if we take advantage of the
# additional components in the push down stack LGroup(), we can scan a
# distribution of structures. Ultimately, there should be a sort
# function that organizes the suboptimal structures in some
# progressive logical order.
class Trace:
    def __init__(self, calc):
        # calc: the object Calculate()
        self.N        = calc.N
        self.T        = calc.T
        self.hv       = calc.hv
        self.all_ctcf = calc.all_ctcf
        self.smap     = calc.smap
        self.calc     = calc
        self.dGmin    = INFINITY
        
        
        # the linkage thread variable
        self.lt = [] # [LThread()] # 
        self.debug_Trace = DEBUG_Trace
        self.debug_traces_top = False # find hot spot regions for further search
        self.debug_get_traces = False # search for suboptimal structures
        self.warning_boundary = 10000
        self.maxlayers = self.calc.maxlayers
        # This is used in traceback routines. See comments in the
        # constructor of Calculate()
        
        # Hot Spots
        self.wt_HS = 0.1
        self.V_HS = INFINITY
        self.hotspot = []
        
        
    #
    
    # This should mainly be used for debugging. It displays all the
    # glink data for all the data points in the upper triangle, so it
    # can print out a lot of data. However, it can be useful for
    # detecting problems in the calculated results.
    def disp_allFE(self, smap):
        for j in range(0, self.N):
            for i in range(0,j):
                if len(smap.glink[i][j].lg) > 0:
                    for k in range(0, len(smap.glink[i][j].lg)):
                        ctp  = smap.glink[i][j].lg[k].motif[0].ctp # conn type
                        bdt  = smap.glink[i][j].lg[k].motif[0].btp # bond type
                        Vij  = smap.glink[i][j].lg[k].Vij
                        join = smap.glink[i][j].lg[k].motif[0].get_base()
                        print "(%2d, %2d)[%2d][%s_%s][%8.2f]: " % (i,j, k, ctp, bdt, Vij), join
                else:
                    print "(%2d, %2d): " % (i,j), "empty"
                #
            #
        #
    #
    
    
    
    
    # search for the top contributing regions for suboptimal structure
    def find_HotSpot(self, dGmin, wt = 0.7):
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
        
        j_max = self.N
        j = j_max - 1
        i = 0
        for j in range(2,self.N):
            for i in range(0, j-1):
                if (len(self.smap.glink[i][j].lg) > 0) and (self.hv[i][j] > 0.0):
                    if self.smap.glink[i][j].lg[0].Vij < self.V_HS:
                        V  = self.smap.glink[i][j].lg[0].Vij
                        tp = self.smap.glink[i][j].lg[0].motif[0].ctp
                        self.hotspot += [(i,j, tp, V)]
                        # To help avoid scooping up redundant
                        # neighboring results, this just looks in the
                        # neighborhood around free energy regions that
                        # have some kind of potential bond.
                    #
                #
            #
        #
        self.hotspot = self.ins_sort_hotspot(self.hotspot)
    #
    
    # sort the results from HotSpots
    def ins_sort_hotspot(self, s):
        for i in range(1,len(s)):    
            j = i                    
            while j > 0 and s[j][3] < s[j-1][3]: 
                s[j], s[j-1] = s[j-1], s[j] # syntactic sugar: swap the items
                j=j-1 
            #
        #
        return s
    #
    
    # The top of the list need not be an 'M', 'K', 'I' or 'B'!!!  So this
    # makes sure that the structure that is selected is one of
    # these. It is used in I-loop searches; i.e., find_best_I() and
    # related structures. So the purpose is to find a structure with a
    # single closing point at (k,l), not a pMBL or a J-loop. The main
    # question I have now is whether this should also include a
    # PK. Currently (160911) by scanNN_MIB(), scan_narrow(),
    # find_best_I(), localscan_for_I()
    def get_kref(self, k, l, V, ctp):
        #
        flag_debug = True
        
        kref = 0
        if flag_debug:
            print "get_kref", k, l, V, ctp
        for m in range(0, len(self.smap.glink[k][l].lg)):
            ctpx = self.smap.glink[k][l].lg[m].motif[0].ctp
            btpx = self.smap.glink[k][l].lg[m].motif[0].btp
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
    
    
    # This is in a similar style as the M-loop calculation in
    # Calculate(). However, here, it saves more of the information and
    # it does not care if the search turns up type 'P' or 'J' or any
    # of the multitude of other motif types.  Everything that
    # satisfies the cutoff is fair game.

    # An important thing to remember about "list_M" is that "list_M"
    # is a place holder (a container) for the particular motif
    # object. Therefore, the container motif must point _to_ the motif
    # object.
    def localscan_for_M(self, i, j, list_M):
        debug_localscan_for_M1 = False
        debug_localscan_for_M2 = DEBUG_Trace
        best_Mk = ()
        best_dGMk = 0.0
        flag_M  = False
        best_Ik = ()
        best_dGIk = 0.0
        flag_I  = False
        motif = []
        
        if debug_localscan_for_M1:
            print "enter localscan_for_M(%2d,%2d)" % (i, j)
        for k in range(1, j-i-1):
            
            if debug_localscan_for_M2:
                print "localscan_for_M: ij(%2d,%2d)(%2d), (%2d,%2d)[%8.3f]|(%2d,%2d)[%8.3f] " \
                    % (i, j, k, i, i+k, self.calc.dG[i][i+k], i+k+1, j, self.calc.dG[i+k+1][j])
            #
            V = self.calc.dG[i][i+k] + self.calc.dG[i+k+1][j]
            if V < self.V_HS:
                V = self.calc.dG[i][i+k] + self.calc.dG[i+k+1][j]
                
                n1 = len(self.smap.glink[i][i+k].lg)   # n1 > 0 if Lgroup has data
                n2 = len(self.smap.glink[i+k+1][j].lg) # n2 > 0 if Lgroup has data
                if debug_localscan_for_M2:
                    print "n12: ", n1, n2, (i, j), k
                #
                if n1 > 0 and n2 > 0:
                    flag_M = True
                    # found some real M-loop
                    join = []
                    for jk in self.smap.glink[i][i+k].lg[0].motif[0].get_base():
                        join += [jk]
                    #
                    for jk in self.smap.glink[i+k+1][j].lg[0].motif[0].get_base():
                        join += [jk]
                    #
                    motif_P = Motif(i, j, V, 'P', '-', join)
                    list_M += [(i,j, [motif_P], V)]
                    if V < best_dGMk: 
                        best_Mk = (i,j, [motif_P], V, k)
                        best_dGMk = V
                        #
                        if debug_localscan_for_M2:
                            n1 = [n1, self.smap.glink[i    ][i+k].lg[0].motif[0].get_base()[0]]
                            n2 = [n2, self.smap.glink[i+k+1][j  ].lg[0].motif[0].get_base()[0]]
                            print "ij = ", i, j,  ", j-i = ", (j-i), "k = ", k, ", n12:", n1, n2, V
                        #
                    #
                elif n1 > 0:
                    ctp  = self.smap.glink[i][i+k].lg[0].motif[0].ctp
                    base = self.smap.glink[i][i+k].lg[0].motif[0].get_base()
                    V    = self.smap.glink[i][i+k].lg[0].Vij
                    if debug_localscan_for_M2:
                        print "base: ", base
                    #
                    if not ctp == 'P' and len(base) > 1:
                        base = [base[0]]
                    #
                    # self.get_kref(i, i+k, V, ctp)
                    motif_J = Motif(i, j, V, 'P', '-', base)
                    if debug_localscan_for_M2:
                        if ctp == 'S':
                            print "n1: "
                            print "boundaries: ", i, (i + k)
                            print motif_J.show_Motif()
                            # sys.exit(0)
                        #
                    #
                    list_M += [(i,j, [motif_J], V)]
                    flag_I = True
                    if V < best_dGIk:
                        best_Ik = (i,j, [motif_J], V, k)
                        best_dGIk = V
                        #
                        if debug_localscan_for_M2:
                            n1x = [n1, ctp, motif_J.get_base()]
                            print "ij = ", i, j,  ", j-i = ", (j-i), "k = ", k, ", n12:", n1x, 0, V
                        #
                    #
                elif n2 > 0:
                    ctp  = self.smap.glink[i+k+1][j].lg[0].motif[0].ctp
                    base = self.smap.glink[i+k+1][j].lg[0].motif[0].get_base()
                    V    = self.smap.glink[i+k+1][j].lg[0].Vij
                    if debug_localscan_for_M2:
                        print "base: ", base
                    #
                    if not ctp == 'P' and len(base) > 1:
                        base = [base[0]]
                    #
                    # self.get_kref(i+k+1, j, V, ctp)
                    motif_J = Motif(i, j, V, 'P', '-', base)
                    if debug_localscan_for_M2:
                        if ctp == 'S':
                            print "n2: "
                            print "boundaries: ", (i+k+1), j
                            print motif_J.show_Motif()
                            # sys.exit(0)
                        #
                    #
                    list_M += [(i,j, [motif_J], V)]
                    flag_I = True
                    if V < best_dGIk:
                        best_Ik = (i,j, [motif_J], V, k)
                        best_dGIk = V
                        #
                        if debug_localscan_for_M2:
                            n2 = [n2, ctp, motif_J.get_base()]
                            print "ij = ", i, j,  ", j-i = ", (j-i), "k = ", k, ", n12:", 0, n2, V
                        #
                    #
                else:
                    if debug_localscan_for_M2:
                        print "ij = ", i, j,  ", j-i = ", (j-i), "k = ", k, ", n12:", 0, 0, V, " empty"
                    #
                #
            else:
                continue
            #
        #
        if debug_localscan_for_M1:
            if flag_M:
                print best_dGMk, best_Mk[3]
                k = best_Mk[4]
                join = best_Mk[2][0].get_branches()
                best_dGMk = best_Mk[3]
                tp1 = self.smap.glink[i    ][i+k].lg[0].motif[0].ctp
                tp2 = self.smap.glink[i+k+1][j  ].lg[0].motif[0].ctp
                print "best_dGMk->: ij(%2d,%2d)(%2d), (%2d,%2d){%s}|(%2d,%2d){%s}[%8.3f] " \
                    % (i, j, k, i, i+k, tp1, i+k+1, j, tp2, best_dGMk), join
                #
                #if i == 0 and j == 14:
                #    sys.exit(0)
                #
            else:
                print "no hits for M in localscan_for_M"
            if flag_I:
                print best_dGIk, best_Ik[3]
                k = best_Ik[4]
                join = best_Ik[2][0].get_branches()
                best_dGIk = best_Ik[3]
                tp1 = '-'; tp2 = '-'
                if len(self.smap.glink[i    ][i+k].lg) > 0:
                    tp1 = self.smap.glink[i    ][i+k].lg[0].motif[0].ctp
                if len(self.smap.glink[i+k+1][j  ].lg) > 0:
                    tp2 = self.smap.glink[i+k+1][j  ].lg[0].motif[0].ctp
                print "best_dGIk->: ij(%2d,%2d)(%2d), (%2d,%2d){%s}|(%2d,%2d){%s}[%8.3f] " \
                    % (i, j, k, i, i+k, tp1, i+k+1, j, tp2, best_dGMk), join
                #
                #if i == 0 and j == 14:
                #    sys.exit(0)
                #
            else:
                print "no hits for M in localscan_for_M"
        #
        return list_M
    #                
    
    
    
    
    
    # This is in a similar style as the M-loop calculation in
    # Calculate(). However, here, it saves any specific information
    # about I-loops. Basically, localscan_for_M() is the major player,
    # but to ensure that there are no issues, we add this into the
    # picture.
    
    # An important thing to remember about "list_M" is that "list_M"
    # is a place holder (a container) for the particular motif
    # object. Therefore, the container motif must point _to_ the motif
    # object.
    def localscan_for_I(self, i, j, list_M):
        debug_localscan_for_I = False
        best_dGI = 0.0
        best_cpI = []
        best_ctpI = ''
        flag_hitI = False
        for l in range(i+1,j+1):
            for k in range(i, l):
                ctpkl = self.smap.glink[k][l].lg[0].motif[0].ctp
                if self.hv[k][l] > 0.0 or ctpkl == 'K' or ctpkl == 'W':
                    # 160915wkd: added a case for 'W'
                    
                    # search the list in glink for the best I,J
                    V, cp, ctp, pks, wyspa = self.calc.filter_BIKMSW_from_glink(k, l)
                    # print "%s %8.3f >= %8.3f, " % (ctp, V, self.V_HS), cp
                    if V < self.V_HS:
                        if debug_localscan_for_I:
                            print "%s %8.3f >= %8.3f, " % (ctp, V, self.V_HS), cp
                        #
                        bondtype = self.calc.get_bondtype(self.hv[k][l])
                        if ctp == 'S':
                            # just want to see if anything is found
                            # through this routine. However,
                            # presently, I seem to have this one
                            # (localscan_for_I) turned off.
                            print "found boundtype %s in localscan_for_I"
                            print V, cp, ctp
                            for m in range(0, len(self.smap.glink[k][l].lg)):
                                print self.smap.glink[k][l].lg[m].motif[0].show_Motif()
                            sys.exit(0)
                        #
                        
                        motif = Motif(i, j, V, ctp, bondtype, cp, pks, wyspa)
                        list_M += [(i,j, [motif], V)]
                        
                        if V < best_dGI:
                            # found at least I-loop
                            flag_hitI = True
                            best_cpI  = [cp]
                            best_dGI = V
                            best_ctpI = ctp
                        #
                    #
                #
            #
        #
        if flag_hitI:
            if debug_localscan_for_I:
                cp = best_cpI[0]
                print "best_cpI: ", best_cpI
                print "localscan_for_I: (%2d,%2d) -> {%s}[%8.3f][%s]" \
                    % (i, j, cp, best_dGI, best_ctpI)
                print "localscan_for_I: flag_hitI = True"
                # sys.exit(0)
            #
        else:
            best_dGI = INFINITY
            if debug_localscan_for_I:
                print "localscan_for_I: no hits"
        #
        return list_M
    #                
    
    
    # This is in a similar style as the M-loop calculation in
    # Calculate(). However, here, it saves any specific information
    # about I-loops. Basically, localscan_for_M() is the major player,
    # but to ensure that there are no issues, I also incorporated this
    # into the picture.
    
    # An important thing to remember about "list_M" is that "list_M"
    # is a place holder (a container) for the particular motif
    # object. Therefore, the container motif must point _to_ the motif
    # object.
    def new_localscan_for_I(self, i, j, list_M):
        debug_new_localscan_for_I = False
        best_dGI = 0.0
        best_cpI = [(0,0)]
        best_ctpI = ''
        flag_hitI = False
        for lgk in self.calc.iloop.lg:
            ctp  = lgk.motif[0].ctp
            base = lgk.motif[0].get_base()
            V    = lgk.Vij
            if not ctp == 'P' and len(base) > 1:
                base = [base[0]]
            #
            if not (i <= base[0][0] and base[0][1] <= j):
                continue
            #
            
            # print "ij(%2d,%2d)[%s][%8.3f], base: " % (i, j, ctp, V), base
            if V < self.V_HS:
                if debug_new_localscan_for_I:
                    print "ij(%2d,%2d)[%s][%8.3f], base: " % (i, j, ctp, V), base
                #
                motif = Motif(i, j, V, 'P', '-', base)
                list_M += [(i,j, [motif], V)]
                if V < best_dGI:
                    # found at least I-loop
                    flag_hitI = True
                    best_cpI  = base
                    best_dGI  = V
                    best_ctpI = ctp
                #
            #
        #
        
        
        #
        if debug_new_localscan_for_I:
            if flag_hitI:
                cp = best_cpI[0]
                print "new_localscan_for_I: (%2d,%2d) -> (%2d,%2d)[%8.3f][%s]" \
                    % (i, j, cp[0], cp[1], best_dGI, best_ctpI) 
            #
            else:
                best_dGI = INFINITY
                if debug_new_localscan_for_I:
                    print "new_localscan_for_I: no hits"
                #
            #
            
        #
        return list_M
    #                
    
    def test_localscan_for_I(self, i, j, list_M, where):
        flag_test_localscan_for_I = False
        if flag_test_localscan_for_I:
            print "test_localscan_for_I():"
        #
        
        if flag_test_localscan_for_I:
            flag_error = False
            error_msg  = ''
            list_nM = self.new_localscan_for_I(i, j, list_M)                        
            list_oM = self.localscan_for_I(i, j, list_M)
            
            if flag_test_localscan_for_I:
                print "len(nM, oM): ", len(list_nM), len(list_oM)
            #

            ##### #####
            if len(list_nM) == len(list_oM):
                for k in range(0, len(list_nM)):
                    # print "list_noM[%4d]: " % k, (list_nM[k][0], list_nM[k][1]), (list_oM[k][0], list_oM[k][1])
                    if not (list_nM[k][0] == list_oM[k][0] and list_nM[k][1] == list_oM[k][1]):
                        error_msg  = "mismatch at k=%d: \n" % k
                        error_msg += "new %s\n" % list_nM[k]
                        error_msg += "old %s\n" % list_oM[k]
                        flag_error = True
                        break
                    #
                #
            else:
                error_msg = "new list size (%d) and old list size (%d) are different" \
                    % (len(list_nM), len(list_oM))
                flag_error = True
            #
            
            if not flag_error:
                list_M = list_nM
            else:
                # presently, I don't know what to do if there are
                # problems, so I will just have to ask the program to
                # stop.  Hopefully, it doesn't.
                print "PROBLEMS(minFE): localscan_for_I() from %s" % (where)
                print "                 %s" % error_msg
                n = min(len(list_nM), len(list_oM))
                print "   new(N^2)    old(N^4)"
                for k in range(0,n):
                    i_n = list_nM[k][0]; j_n = list_nM[k][1]
                    i_o = list_oM[k][0]; j_o = list_oM[k][1]
                    print "  n(%3d,%3d)   o(%3d,%3d)" % (i_n, j_n, i_o, j_o)
                sys.exit(1)
            #
        #
        else:
            list_M = self.new_localscan_for_I(i, j, list_M)
        #
        return list_M
    #
    
    
    
    # This operation removes redundant solutions from the list
    # obtained by localscan_for_M(i, j). The algorithm is derived from
    # the program reduce_list(t) that is listed at the end of this
    # package.
    def prune_list_M(self, list_M):
        debug_prune_list_M = False
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
        # sys.exit(0)
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
    def get_traces_top(self, i_top, j_top, ndx, layer, flag_filter):
        #
        
        if self.debug_traces_top:
            print "get_traces_top (%d,%d)[%s]" \
                % (i_top, j_top, self.smap.glink[i_top][j_top].lg[0].motif[0].ctp)
        #
        
        if not (self.hv[i_top][j_top] > 0.0):
            print "WARNING(Trace.get_traces_top()): enthalpy function hv is zero"
            print "                                 at end point (%d,%d)!" % (i_top,j_top)
        #
        
        if self.debug_traces_top:
            print "reference threads"
        #
        
        for k in range(0, len(self.smap.glink[i_top][j_top].lg)):
            
            if len(self.lt) > self.warning_boundary:                
                if self.wt_HS < 0.95:
                    self.adjust_FEboundaries(0.05)
                elif self.wt_HS < 0.99:
                    self.adjust_FEboundaries(0.01)
                elif self.wt_HS < 0.999:
                    self.adjust_FEboundaries(0.001)
                #
                #
            #
            
            V    = self.smap.glink[i_top][j_top].lg[k].Vij
            if flag_filter and V > self.V_HS:
                continue

                
            if self.debug_traces_top:
                print "set new thread"
            self.lt += [LThread(self.N)] # create a new thread
            ndx = len(self.lt) - 1 # assign the particular LThread()
            # 
            
            
            s = self.calc.space(3*layer)
            tp   = self.smap.glink[i_top][j_top].lg[k].motif[0].ctp
            if self.debug_traces_top:
                base = self.smap.glink[i_top][j_top].lg[k].motif[0].get_base()
                join = self.smap.glink[i_top][j_top].lg[k].motif[0].get_branches()
                pk   = self.smap.glink[i_top][j_top].lg[k].motif[0].get_pks()
                print "%s[%s](%3d, %3d)[%8.3f][ndx = %d, layer = %d]: " \
                    % (s, tp, i_top, j_top, V, ndx, layer), base, join, pk
            #
            
            self.get_traces(i_top, j_top, tp, V, ndx, layer)
        #
        if self.debug_traces_top:
            print "pass first part of get_traces_top()"
        #
        # these operations must run in consort to work well.
        list_M = []
        list_M = self.localscan_for_M(i_top, j_top, list_M)
        
        list_M  = self.test_localscan_for_I(i_top, j_top, list_M, "get_traces_top")
        # list_M = self.new_localscan_for_I(i_top, j_top, list_M)
        if self.debug_traces_top:
            print "finished test_localscan_for_I"
        #
        
        list_M = self.ins_sort_hotspot(list_M)
        list_M = self.prune_list_M(list_M)
        if self.debug_traces_top:
            k = 0
            for lM in list_M:
                ctp = lM[2][0].ctp
                if ctp == 'P' or ctp == 'J':
                    print "%3d (%3d, %3d)[%8.2f]" % (k, lM[0], lM[1], lM[3]), \
                        lM[2][0].ctp, "         ", lM[2][0].get_base()
                else:
                    print "%3d (%3d, %3d)[%8.2f]" % (k, lM[0], lM[1], lM[3]), \
                    lM[2][0].ctp, lM[2][0].get_base(), lM[2][0].get_branches()
                #
                k += 1
            #
            # sys.exit(0)
        #
        if self.debug_traces_top:
            print "add_localscan_for_M: enter"
        #
        self.add_localscan_for_M(i_top, j_top, list_M)
        self.lt = self.prune_lt(self.lt)
        # print "add_localscan_for_M: exit"
        # sys.exit(0)
    #
    
    # removes redundant solutions from the list obtained by
    # localscan_for_M(i, j). This is a ratehr simple approach that
    # works on similarity based on lengths, free energy, and other
    # obvious kinds of tests.
    def prune_lt(self, lt):
        debug_prune_lt = False
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
            sys.exit(0)
        #
        return uniq_lt
    #
    
    def add_localscan_for_M(self, i_top, j_top, list_M):
        debug_add_localscan_for_M = self.debug_get_traces
        layer = 0
        k_ref = 0
        if debug_add_localscan_for_M:
            print "add_localscan_for_M(ij_top = (%d,%d))" % (i_top, j_top)
        #
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
                ctpx = self.smap.glink[i][j].lg[m].motif[0].ctp
                btpx = self.smap.glink[i][j].lg[m].motif[0].btp
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
                    sys.exit(0)
                #
                
                if ctp == 'P':
                    print "add_localscan_for_M(): planned exit when ctp = 'P'"
                    # sys.exit(0)
                #
            #
            
            dGij = V
            bondtype = '-'
            # closing at (i_top,j_top), build M-loop/quasi stem
            if i_top < ib and jb < j_top:
                
                # Build M loops from position (i_top,j_top)
                # (exclusively) based on search result obtained at
                # bondaries [ib,jb]
                
                # so the condition (i_top < ib < jb < j_top) is
                # required because we want to _close_ the structure at
                # (i_top,j_top).
                
                if debug_add_localscan_for_M:
                    print "searching M: "
                self.lt += [LThread(self.N)] # create a new thread
                ndx = len(self.lt) - 1 # assign the particular LThread()
                
                
                if self.hv[i_top][j_top] > 0.0:
                    # create the top structure
                    dGij = self.calc.fe.calc_dG(i_top, j_top, self.hv[i_top][j_top], self.T, "add_localscan_for_M")
                    if debug_add_localscan_for_M:
                        print "add_localscan_for_M(): (%2d,%2d)[%8.3f]" \
                            % (i_top, j_top, dGij), [(i_top,j_top)]
                        #
                    #
                
                    # add closing point for the structure
                    bondtype = self.calc.get_bondtype(self.hv[i_top][j_top])
                    # only associated with the closing point 'M' or 'I'!!!
                    if ctp == 'P':
                        self.lt[ndx].add_lsegment((i_top, j_top), ndx, dGij, 'M', bondtype)
                    else:
                        self.lt[ndx].add_lsegment((i_top, j_top), ndx, dGij, 'I', bondtype)
                    #
                else:
                    dGij = V
                    self.lt[ndx].add_lsegment((i_top, j_top), ndx, dGij, ctp, btp)
                #
                
                
                if ctp == 'S':
                    # case for a stem (antiparallel or parallel)
                    self.make_Stemtrace(i, j, ndx, k_ref, layer + 1, ctp, btp, "add_localscan_for_M")
                elif ctp == 'K':
                    # case for a pseudoknot
                    self.make_PKtrace(i, j, ndx, k_ref, layer + 1, ctp, btp, "add_localscan_for_M")
                elif ctp == 'W':
                    # case for 'W'
                    self.make_Islandtrace(i, j, ndx, k_ref, layer + 1, ctp, btp, "add_localscan_for_M")
                else:
                    # all other cases
                    for lk in branches:
                        # new ij
                        i_P = lk[0]; j_P = lk[1]
                        tpP = self.smap.glink[i_P][j_P].lg[0].motif[0].ctp
                        VpP = self.smap.glink[i_P][j_P].lg[0].Vij
                        self.get_traces(i_P, j_P, tpP, VpP, ndx, layer+1)
                    #
                #
            #
            # now build a P-loop between [i_top, j_top] (inclusive)
            # based on search result obtained at bondaries [ib,jb]
            if debug_add_localscan_for_M:
                print "searching P: "
            self.lt += [LThread(self.N)] # create a new thread
            ndx = len(self.lt) - 1 # assign the particular LThread()
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
                    if debug_add_localscan_for_M:
                        print "lk: ", lk, "ij_P: ", i_P, j_P
                    #
                    
                    tpP = self.smap.glink[i_P][j_P].lg[0].motif[0].ctp
                    VpP = self.smap.glink[i_P][j_P].lg[0].Vij
                    self.get_traces(i_P, j_P, tpP, VpP, ndx, layer+1)
                #
            #
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
        flag_debug = False
        if flag_debug:
            print "make_Islandtrace"
        #
        
        self.lt[ndx].add_lsegment((i,j), ndx, 0.0, 'W', 'bgn')
        
        
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
            self.lt[ndx].add_lsegment((i_W,j_W), ndx, dGijW, 'W', 'wyspa')
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
                ctpW = self.smap.glink[i_W][j_W].lg[0].motif[0].ctp
                VpW  = self.smap.glink[i_W][j_W].lg[0].Vij
                if flag_debug:
                    print "ij_W = (%d,%d)[%s], ndx(%d), layer(%d)" % (i_W, j_W, ctpW, ndx, layer)
                self.get_traces(i_W, j_W, ctpW, VpW, ndx, layer+1)
            # if smap is empty, then this is all there is.
        #
        self.lt[ndx].add_lsegment((i,j), ndx, 0.0, 'W', 'end')

        if flag_debug:
            n_wyspa = 5
            if len(wyspa) > n_wyspa:
                print "make_Islandtrace: found islands > %d" % n_wyspa
                print "Island results: from %s" % iprog
                for thr in self.lt[ndx].thread:
                    print thr.disp_lseg()
                if STOP_AT_FIND:
                    sys.exit(0)
                #
            #
        #
        # #####12345
    #
    
    
    # need to add function for 'W'
    def make_PKtrace(self, i, j, ndx, k_ref, layer, ctp, btp, iprog = "get_trace"):
        flag_debug = False
        if flag_debug:
            print "make_PKtrace"
        #
        
        # root stem
        self.lt[ndx].add_lsegment((i,j), ndx, 0.0, 'K', 'bgn')
        v = self.smap.glink[i][j].lg[k_ref].motif[0].get_base()
        i_R  = v[0][0]; j_R = v[0][1]
        ctpR = self.smap.glink[i_R][j_R].lg[k_ref].motif[0].ctp
        btpR = self.smap.glink[i_R][j_R].lg[k_ref].motif[0].btp
        VpR  = self.smap.glink[i_R][j_R].lg[k_ref].Vij
        
        if ctpR == 'I' or ctpR == 'M':
            dGij = self.calc.fe.calc_dG(i_R, j_R, self.hv[i_R][j_R], self.T, "make_PKtrace")
            self.lt[ndx].add_lsegment((i_R,j_R), ndx, dGij, ctpR, btpR)
        #
        
        join = self.smap.glink[i_R][j_R].lg[k_ref].motif[0].get_branches()
        if flag_debug:
            print "join = ", join
        for lk in join:
            # whether 1 or 100 elements, they can be processed this way.
            i_M = lk[0]; j_M = lk[1]     # new ij of this subdomain
            if len(self.smap.glink[i_M][j_M].lg) > 0:
                ctpM = self.smap.glink[i_M][j_M].lg[0].motif[0].ctp
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
            self.lt[ndx].add_lsegment((x[0][0],x[0][1]),ndx,x[1], 'l', x[2])
        #
        self.lt[ndx].add_lsegment((i,j), ndx, 0.0, 'K', 'end')
        
        if flag_debug:
            print "pk results: from %s" % iprog
            for thr in self.lt[ndx].thread:
                print thr.disp_lseg()
            if btpx == 'sp':
                print "parallel"
                # sys.exit(0)
        # #####12345pk
    #
    

    # need to add function for 'W'
    def make_Stemtrace(self, i, j, ndx, kref_t, layer, ctp, btp, iprog = "get_trace"):
        
        # ndx is the index for the particular thread of
        # self.lt[ndx].thread
        
        # kref_t the particular reference to the tail section of the
        # stem
        
        flag_debug = False
        if flag_debug:
            print "make_Stemtrace(called from %s)" % iprog
        #
        
        # root stem
        flag_pp = False
        self.lt[ndx].add_lsegment((i,j), ndx, 0.0, 'S', 'bgn')
        ctpS = self.smap.glink[i][j].lg[kref_t].motif[0].ctp
        btpS = self.smap.glink[i][j].lg[kref_t].motif[0].btp
        VpS  = self.smap.glink[i][j].lg[kref_t].Vij
        join = self.smap.glink[i][j].lg[kref_t].motif[0].get_branches()
        
        # establish the head and tail of the stem
        stmlen = len(join) - 1
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
            for k in range(0, len(join)-1):
                # whether 1 or 100 elements, they can be processed this way.
                i_S = join[k][0]; j_S = join[k][1]     # new ij of this subdomain
                # print "stem ij_S: ", i_S, j_S
                ctpS = 'S'
                bondtype = self.calc.get_bondtype(self.hv[i_S][j_S])
                if bondtype == 's':
                    bondtype += 'a'
                else:
                    if not bondtype == 'sa':
                        # this decision must be enforced in the case of a stem!
                        bondtype = 'c'
                #
                dGij = self.calc.fe.calc_dG(i_S, j_S, self.hv[i_S][j_S], self.T, "make_Stemtrace 1")
                self.lt[ndx].add_lsegment((i_S,j_S), ndx, dGij, 'S', bondtype)
            #
            
            # now add the head
            for m in range(0, len(self.smap.glink[i_h][j_h].lg)):
                # print "motif(ij_h): ", self.smap.glink[i_h][j_h].lg[m].motif[0].show_Motif()
                ctp = self.smap.glink[i_h][j_h].lg[m].motif[0].ctp
                if ctp == 'W':
                    # chr11_111099793_111320552_res5kb.heat satisfied this condition
                    print "aa ij_t (tail): ", i_t, j_t
                    print "aa ij_h (head): ", i_h, j_h
                    print "found an island at the head of antiparallel stem."
                    if STOP_AT_FIND:
                        sys.exit(0)
                    #
                #
                if ctp == 'B' or ctp == 'I' or ctp == 'M' or ctp == 'S' or ctp == 'W':
                    kref_h = m
                    break
                #
            #
            # now add the terminus (the head of the stem)
            ctpS = self.smap.glink[i_h][j_h].lg[kref_h].motif[0].ctp
            btpS = self.smap.glink[i_h][j_h].lg[kref_h].motif[0].btp
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
            # sys.exit(0)
        else:
            if flag_debug:
                print "parallel"
            #
            
            # scan up the antiparallel stem from the tail up the head
            # (the whole stem without a cut off at the top)
            for k in range(0, len(join)):
                # whether 1 or 100 elements, they can be processed this way.
                i_S = join[k][0]; j_S = join[k][1]     # new ij of this subdomain
                if flag_debug:
                    print "stem ij_S: ", i_S, j_S
                #
                ctpS = 'S'
                bondtype = self.calc.get_bondtype(self.hv[i_S][j_S])
                if bondtype == 's':
                    bondtype += 'p'
                else:
                    if not bondtype == 'sp':
                        # this decision must be enforced in the case of a stem!
                        bondtype = 't'
                #
                dGij = self.calc.fe.calc_dG(i_S, j_S, self.hv[i_S][j_S], self.T, "make_Stemtrace 2")
                self.lt[ndx].add_lsegment((i_S,j_S), ndx, dGij, 'S', bondtype)
            #
            
            # now add the head at (i_h + 1, j_h - 1)
            if len(self.smap.glink[ihh][jhh].lg) > 0:
                for m in range(0, len(self.smap.glink[ihh][jhh].lg)):
                    if flag_debug:
                        print "motif(ij_h): ", self.smap.glink[ihh][jhh].lg[m].motif[0].show_Motif()
                    #
                    ctp = self.smap.glink[ihh][jhh].lg[m].motif[0].ctp
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
                    ctph = self.smap.glink[ihh][jhh].lg[kref_h].motif[0].ctp
                    btph = self.smap.glink[ihh][jhh].lg[kref_h].motif[0].btp
                    Vph  = self.smap.glink[ihh][jhh].lg[kref_h].Vij
                    if flag_debug:
                        print "pp: ij_t = (%d,%d)" % (i_t, j_t)
                        print "pp: ij_h = (%d,%d)" % (i_h, j_h)
                        print "pp: ijhh = (%d,%d)[%s], ndx(%d), layer(%d)" % (ihh, jhh, ctph, ndx, layer)
                        if  i_t == 0 and j_t == 14:
                            print "planned exit when evaluating a parallel stem"
                            sys.exit(0)
                        
                    self.get_traces(ihh, jhh, ctph, Vph, ndx, layer+1)
                else:
                    if flag_debug:
                        print "motif(ij_h):  empty" 
        #    
        # sys.exit(0)
        
        self.lt[ndx].add_lsegment((i,j), ndx, 0.0, 'S', 'end')
        
        if flag_debug:
            print "stem results: from %s" % iprog
            for thr in self.lt[ndx].thread:
                print thr.disp_lseg()
            #
            
            #if flag_pp:
            #    print "stop when pp"
            #    sys.exit(0)
        # #####12345stem
    #
    
    
    # find the minimum free energy
    def get_traces(self, i, j, ctp, V, ndx, layer):
        
        # ndx is the index for the particular thread of
        # self.lt[ndx].thread
        
        debug_get_traces = self.debug_get_traces
        if debug_get_traces and layer == 0:
            print "get_traces"
        # 
        if layer > self.maxlayers:
            print "layer(%d) > maxlayer(%d)" % (layer, self.maxlayers)
            print "ERROR: something wrong in the recursion of get_traces!"
            sys.exit(1)
        #
        
        
        #####
        flag_found = False
        kref_t = 0
        for m in range(0, len(self.smap.glink[i][j].lg)):
            ctpx   = self.smap.glink[i][j].lg[m].motif[0].ctp
            Vx    = self.smap.glink[i][j].lg[m].Vij
            if debug_get_traces:
                basex = self.smap.glink[i][j].lg[m].motif[0].get_base()
                joinx = self.smap.glink[i][j].lg[m].motif[0].get_branches()
                pkx   = self.smap.glink[i][j].lg[m].motif[0].get_pks()
                print "layer = %3d: %2d '%s' == ref(%s)?,  %8.3f <= ref(%8.3f)?" \
                    % (layer, m, ctpx, ctp, Vx, V),  basex, joinx, pkx
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
            sys.exit(1)
        #
        
        if ctp == 'B' or ctp == 'I' or ctp == 'M':
            # these are the only structures that form links
            dGij = self.calc.fe.calc_dG(i, j, self.hv[i][j], self.T, "get_traces")
            # ctp is already defined, it should correspond to kref_t
            btp  = self.smap.glink[i][j].lg[kref_t].motif[0].btp
            if debug_get_traces:
                print "get_traces(): (%2d,%2d)[%8.3f]" \
                    % (i, j, dGij), [(i,j)]
            #
            self.lt[ndx].add_lsegment((i,j), ndx, dGij, ctp, btp)
            
            # pseudoknot connections have to be handled in the
            # processing.
        #
        s = self.calc.space(3*layer)
        if debug_get_traces:
            join = self.smap.glink[i][j].lg[kref_t].motif[0].get_branches()
            print "%s[%s](%3d, %3d)[%8.3f][ndx = %d, layer = %d]: " \
                % (s, ctp, i,j, V, ndx, layer), join
        if ctp == 'M':
            if debug_get_traces:
                print "found a case of %s" % ctp
                print "%s------------------" % s
            #
            join = self.smap.glink[i][j].lg[kref_t].motif[0].get_branches()
            for lk in join:
                # new ij
                i_M = lk[0]; j_M = lk[1]
                ctpM = self.smap.glink[i_M][j_M].lg[0].motif[0].ctp
                VpM = self.smap.glink[i_M][j_M].lg[0].Vij
                self.get_traces(i_M, j_M, ctpM, VpM, ndx, layer+1)
            #
            #
        elif ctp == 'P':
            if debug_get_traces:
                print "found a case of %s" % ctp
                print "%s------------------" % s
            #
            join = self.smap.glink[i][j].lg[kref_t].motif[0].get_branches()
            for lk in join:
                i_P = lk[0]; j_P = lk[1]
                ctpP = self.smap.glink[i_P][j_P].lg[0].motif[0].ctp
                VpP = self.smap.glink[i_P][j_P].lg[0].Vij
                self.get_traces(i_P, j_P, ctpP, VpP, ndx, layer+1)
            #
            #
        elif ctp == 'I':
            if debug_get_traces:
                print "found a case of I"
            #
            jn  = self.smap.glink[i][j].lg[kref_t].motif[0].get_branches()
            # new ij!!!
            i_I = jn[0][0]; j_I = jn[0][1]
            VpI, cpI, ctpI, pks, wyspa = self.calc.filter_BIKMSW_from_glink(i_I,j_I)
            self.get_traces(i_I, j_I, ctpI, VpI, ndx, layer+1)
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
            
            if VpJ == 1000:
                print "ERROR(get_traces): could not find a linkage at (%d,%d)" % (i_J, j_J)
                sys.exit(1)
            #
            self.get_traces(i_J, j_J, ctpJ, VpJ, ndx, layer+1)
            #
            
        elif ctp == 'S':
            
            if debug_get_traces:
                print "get_traces: S:"
            #
            
            btp  = self.smap.glink[i][j].lg[kref_t].motif[0].btp
            self.make_Stemtrace(i, j, ndx, kref_t, layer+1, ctp, btp, "get_trace")
            
            #
        elif ctp == 'K':
            
            if debug_get_traces:
                print "get_traces: K:"
            #
            
            btp  = self.smap.glink[i][j].lg[kref_t].motif[0].btp
            #print kref_t
            #print self.smap.glink[i][j].lg[kref_t].motif[0].show_Motif()
            #sys.exit(0)
            self.make_PKtrace(i, j, ndx, kref_t, layer+1, ctp, btp, "get_trace")
            
            #
        elif ctp == 'W':
            
            if debug_get_traces:
                print "get_traces: W:"
            #
            btp  = self.smap.glink[i][j].lg[kref_t].motif[0].btp
            self.make_Islandtrace(i, j, ndx, kref_t, layer+1, ctp, btp, "get_trace")
            
            #
        else:
            if debug_get_traces:
                print "%s---" % s
            #
        #
    #
    
    
    
#





#####################################################################
################  CALCULATE BOLTZMANN DISTRIBUTION  #################
#####################################################################

class Boltzmann:
    def __init__(self, calc):
        self.debug_Boltzmann = False
        self.Z    = 1.0
        self.explist = [] # to ensure that the exp doesn't over extend
        self.shift   = 0  #
        self.rewt    = 1.0
        self.kB   = 0.00198
        self.calc = calc
        self.fe   = calc.fe
        self.N    = calc.N
        self.T    = calc.T
        self.flag_set_Z = False
        
        self.debug = False
        
        # to ensure that Z is calculated before doing other operations.
    #
    
    def calc_Z(self, lt, T):
        kBT = self.kB * T
        
        # Result: $Z = 1 + \sum_k {\exp(-dG_k / kBT) }$
        
        # The method for calculating the partition function may seem a
        # bit odd at first brush. After all, the mathematical formula
        # is really quite straight forward.
        
        # The problem is that floats have a limit in size
        # that can be easily exceeded. So we have first construct the
        # exponents and then analyze them through a subsidiary
        # function KahanSumExp() to add the exponentials together.
        
        # math.exp
        self.explist = []
        for ltk in lt:
            self.explist += [ -ltk.dG / kBT ]
        #
        
        self.Z, self.shift = KahanSumExp(self.explist)
        re_wt = float(2**self.shift)
        self.flag_set_Z = True # calculation of Z is now done!
        
        # 
        if self.shift > -1: # self.debug:
            print "Z = ", self.Z
            print "shift: %4d,  new exp wt: %.1f" % (self.shift, re_wt)

            # print out the first ten arguments of on the partition
            # function exponential term
            #                             #        #
            print "ndx    arg       dG/kBT          dG        p(ndx)"
            for k in range(0, len(self.explist)):
                if k > 10: # could be a huge list
                    break
                print "%2d  %8.2f   %8.2f  %12.2f    %8.3g" \
                    % (k, self.explist[k], (-lt[k].dG/kBT), lt[k].dG, self.calc_p(lt[k], T))
            #
            # sys.exit(0)
        #
    #
    
    
    def set_LThread_TdS(self, lt, T):
        if not self.flag_set_Z:
            print "ERROR(Boltzmann.calc_p()): Z is undefined"
            sys.exit(1)
        #
        
        for m in range(0, len(lt)):
            lt[m].TdS = 0.0
            for tr in lt[m].thread:
                v = tr.ij_ndx
                i = v[0]
                j = v[1]
                # print "ij= ", i, j
                lt[m].TdS += self.fe.TdS(i,j,T)
            #
        #
        if self.debug_Boltzmann:
            dt = DispThreads(self.calc.N)
            for m in range(0, len(lt)):
                ss_m = dt.makeLThreadDotBracket_1b(lt[m], True)
                print m, lt[m].TdS, ss_m, lt[m].p
            sys.exit(0)
        #
        return lt
    #
    
    def set_LThread_p(self, lt, T):
        if not self.flag_set_Z:
            print "ERROR(Boltzmann.calc_p()): Z is undefined"
            sys.exit(1)
        # Note:
        # $p = exp( -dG_k / kBT ) / Z$ 
        # $p = exp( -dG_k / kBT ) / \[ 1 + \sum_k {\exp(-dG_k / kBT)} \]$ 
        
        kBT = self.kB * T
        # math.exp
        if self.shift > 0:
            for k in range(0, len(lt)):
                exponent = - (lt[k].dG / kBT)
                if exponent - self.shift * log(2) < -708.396:
                    lt[k].p = 0.0
                else:
                    lt[k].p = exp( exponent - float(self.shift)*log(2) ) / self.Z
            #
        else:
            for k in range(0, len(lt)):
                exponent = - (lt[k].dG / kBT)
                if exponent - self.shift * log(2) < -708.396:
                    lt[k].p = 0.0 # surely, this is neglectable
                else:
                    lt[k].p = exp( exponent ) / self.Z
                #
            #
        #
        return lt
    #
    
    def calc_p(self, ltk, T):
        if not self.flag_set_Z:
            print "ERROR(Boltzmann.calc_p()): Z is undefined"
            sys.exit(1)
        #
        kBT = self.kB * T
        # math.exp
        if self.shift > 0:
            # print (float(self.shift)*log(2))
            exponent = - (ltk.dG / kBT)
            if exponent - self.shift * log(2) < -708.396:
                p = 0.0 # surely, this is neglectable
            else:
                p = exp( exponent - float(self.shift)*log(2) ) / self.Z
        else:
            exponent = - (ltk.dG / kBT)
            if exponent - self.shift * log(2) < -708.396:
                p = 0.0 # surely, this is neglectable
            else:
                p = exp( exponent ) / self.Z
            #
        #
        return p
    #
#

class Cluster:
    def __init__(self, calc):
        self.calc = calc
        self.N    = calc.N
        # initialize the free energy matrix
        self.clusters = []
        self.clusters = self.calc.fe.mtools.initialize_matrix(self.clusters, self.calc.N, 0.0)
        self.debug = False
    #
    
    def clusterlist(self, ltlist):
        
        wt = self.calc.fe.mtools.normalize_matrix(self.calc.dG, "neg")
        if self.debug:
            print len(ltlist)
        #
        ij_min = 100000.0
        for ltk in ltlist:
            pBoltz = ltk.p
            for tr in ltk.thread:
                if tr.ctp == 'P' or tr.ctp == 'J':
                    print "found a P for ", tr.ij_ndx
                    continue
                v = tr.ij_ndx
                i = v[0]; j = v[1]
                # print v
                self.clusters[i][j] += 1.0*pBoltz*wt[i][j] #  ij
                self.clusters[j][i] += 1.0*pBoltz*wt[i][j] #  ji
                if self.clusters[i][j] < ij_min:
                    ij_min = self.clusters[i][j]
                #
            #
        #
        if ij_min < 0.0:
            shift = - ij_min
            print "encountered positive entropy values"
            print "up-shifting free energy map data by ", shift
            for j in range(0, self.N):
                for i in range(0, self.N):
                    if not self.clusters[i][j] == 0.0:
                        self.clusters[i][j] += shift
                        self.clusters[j][i] += shift
                    #
                #
            #
        #
    #
#





class Manager:
    def __init__(self):
        # important data variables
        self.flnm = ''
        self.T    = -1.0
        self.smap = []
        self.dG   = []
        self.N    = -1
        # objects
        self.calc  = '' # main calculation routine
        self.trace = '' # trace back routines
        self.clust = '' # clustering routines
        self.btz   = '' # boltzmann calculations
        self.dt    = '' # data handling routines
        # other
        self.debug_Manager = SHOWMAIN
        
    #
    
    def runCalculations(self, CL):
        self.calc = Calculate(CL)
        self.T    = self.calc.T
        self.flnm = self.calc.flnm
        self.N    = self.calc.N
        self.dG, self.smap = self.calc.minFE(self.T)
        if self.debug_Manager:
            print "number of iloop entries: ", len(self.calc.iloop.lg)
            for lgk in self.calc.iloop.lg:
                print lgk.motif[0].base, lgk.motif[0].ctp, lgk.Vij, lgk.motif[0].join
            # sys.exit(0)
        #
        
        if self.debug_Manager:
            print self.calc.fe.mtools.disp_fmatrix(self.calc.hv, "enthalpy")
        #
        
        self.trace = Trace(self.calc)
        
        if DEBUG_Trace:
            self.trace.disp_allFE(self.smap)
        #
        
        if self.debug_Manager:
            print "\n\n\n"
            print "Results:"
            print "total free energy"
            print "dG(0,%d)[%s] = %10.3f [kcal/mol]" \
                % (self.calc.N-1,
                   self.smap.glink[0][self.calc.N-1].lg[0].motif[0].ctp,
                   self.smap.glink[0][self.calc.N-1].lg[0].Vij)
            print "\n"
        #
        
        dGmin = INFINITY
        
        if self.debug_Manager:
            print "traceback mFE:"
            dGmin = self.calc.traceback_mFE(0, self.calc.N-1, 0, self.debug_Manager)
            print "B: bound; I: I-loop, M: M-loop"
            print "min FE: %8.2f [kcal/mol]" % dGmin
            print string.join(self.calc.opt_ss_seq, '')
        else:
            dGmin = self.calc.traceback_mFE(0, self.calc.N-1, 0, False)
        #
        
        print "dGmin = %8.3f" % dGmin
        #print string.join(self.calc.opt_ss_seq, '')
        #print "break here"
        #sys.exit(0)
        
        print "looking for hot spots in the free energy"
        self.trace.find_HotSpot(dGmin, 0.65)
        
        if self.debug_Manager:
            print "HotSpots: "
            if len(self.trace.hotspot) > 0:
                for hsk in self.trace.hotspot:
                    print "(%3d,%3d)[%s] %8.3f" % (hsk[0], hsk[1], hsk[2], hsk[3])
                #
            #
        #
        
        print "doing traceback on hot spots"
        flag_filter = True
        if len(self.trace.hotspot) == 0:
            # covers the case where there is only one structure
            print "get_traces_top(hotspots = 0)"
            self.trace.get_traces_top(0, self.calc.N-1, 0, 0, flag_filter)
        #
        
        print "number of hotspots found: ", len(self.trace.hotspot)
        # sys.exit(0)
        
        for hs in self.trace.hotspot:
            i_hs = hs[0]; j_hs = hs[1]
            if self.debug_Manager:
                print "get_traces_top(ijhs = (%d,%d))" % (i_hs, j_hs)
            #
            self.trace.get_traces_top(i_hs, j_hs, 0, 0, flag_filter)
        #
        
        
        self.dt = DispThreads(self.calc.N)
        if self.debug_Manager:
            print "DispThreads"
            for tltk in self.trace.lt:
                print "----"
                self.dt.disp_LThread(tltk)
            #
        #
        print "size of matrix: %d" % self.N
        print "found %d structures:" % len(self.trace.lt)
        
        self.btz = Boltzmann(self.calc)
        self.btz.calc_Z(self.trace.lt, self.T)
        self.trace.lt = self.btz.set_LThread_p(  self.trace.lt, self.T)
        self.trace.lt = self.btz.set_LThread_TdS(self.trace.lt, self.T)
        
        self.clust = Cluster(self.calc)
        self.clust.clusterlist(self.trace.lt)
        
    #

    
    def printResults(self):
        
        # make directory and display and store files
        
        ff   = self.flnm.split('.')
        flhd = ff[0]
        
        # check if the directory already exists
        try:
            if not os.path.isdir(flhd):
                os.mkdir(flhd) # build a subdirectory to store figures in
            else:
                print "WARNING: the directory '%s'" % flhd
                print "         already exists, overwriting files\n" 
        except OSError:
            print "ERROR: problems making %s" % flhd
            sys.exit(1)
        #
        os.chdir(flhd) # move to that directory
        

        # Save information on the secondary structure in long
        # Janusz-Bonieski format.
        ssflnm = flhd + "_summary.txt"
        try:
            fp = open(ssflnm, 'w')
            header = "# summary of structures from %s\n" % flhd
            # thermodynamic parameters
            # entropy parameters from the chreval
            header += "# thermodynamic parameters:\n"
            header += "#   entropy:\n"
            header += "#     local Kuhn length       = %.3g [bps]\n"      % self.calc.fe.xi      
            header += "#     segment length          = %.3g [bps]\n"      % self.calc.fe.seg_len 
            header += "#     lambda (binding dist)   = %.3g [nts^(-1)]\n" % self.calc.fe.lmbd
            header += "#     gmm (SAW parameter)     = %.3g (no units)\n" % self.calc.fe.gmm
            header += "#     Temperature             = %.3g [K]\n"        % self.T
            # constants: weights for the enthalpy terms 
            header += "#   enthalpy:\n"
            header += "#     input data rescaling wt = %.3g\n" % self.calc.fe.rescale_wt
            header += "#     febase                  = %.3g [kcal/mol]\n" % self.calc.fe.base
            header += "#     feshift                 = %.3g\n" % self.calc.fe.shift
            header += "# statistics:\n"
            header += "#   total number of structures:          %8d\n"  % len(self.trace.lt)
            header += "#   final fraction of structures extracted:    %6.3f\n" % self.trace.wt_HS
            header += "#   upper limit of the free energy:         %8.2f\n" % self.trace.V_HS
            header += "#   minimum free energy:                    %8.2f\n" % self.trace.dGmin
            header += "#   tandem CTCF loop threshold:             %8.2f [kcal/mol]\n" % self.calc.fe.ctcf_tthresh
            header += "#   convergent CTCF loop threshold:         %8.2f [kcal/mol]\n" % self.calc.fe.ctcf_cthresh
            header += "# ----\n"
            header += "#\n"
            fp.write(header)
            fp.close()
        except OSError:
            print "ERROR: problems opening %s" % (ssflnm)
            sys.exit(1)
        #
        
        
        k = 1
        file_results = ''
        prnt_results = ''
        print "number of threads obtained: %d" % len(self.trace.lt)
        for thrds in self.trace.lt:
            file_results = "> %s    dG = %8.3f   p = %12.8f\n" % (string.zfill(k, 5), thrds.dG, thrds.p)
            file_results += self.dt.makeLThreadDotBracket_1b(thrds, 2)
            if k < 5:
                prnt_results += "> dG = %8.3f   p = %12.8f\n" % (thrds.dG, thrds.p)
                prnt_results += file_results
                prnt_results += self.dt.makeLThreadDotBracket_1b(thrds, 2)
            #
            if self.calc.p_all_1D:
                outfile = flhd + '_%s.DBN' % string.zfill(k, 5)
                self.dt.printLThreadDotBracket_VARNA(outfile, thrds)
                hmapflnm = flhd + '_%s.heat' % string.zfill(k, 5)
                self.dt.printLThreadHeatMap(hmapflnm, thrds, True)
                #
                ijpair_flnm     = flhd + '_%s.chpair' % string.zfill(k, 5)
                self.dt.print_ijPairList(ijpair_flnm, thrds, True)
                #
                ijsimres_flnm   = flhd + '_%s.simres' % string.zfill(k, 5)
                self.dt.print_SimRNArestraints(ijsimres_flnm, thrds, "slope", True)
                #
            else:
                if k < 50:
                    outfile = flhd + '_%s.DBN' % string.zfill(k, 5)
                    self.dt.printLThreadDotBracket_VARNA(outfile, thrds)
                    hmapflnm = flhd + '_%s.heat' % string.zfill(k, 5)
                    self.dt.printLThreadHeatMap(hmapflnm, thrds, True)
                    #
                    ijpair_flnm     = flhd + '_%s.chpair' % string.zfill(k, 5)
                    self.dt.print_ijPairList(ijpair_flnm, thrds, True)
                    #
                    ijsimres_flnm   = flhd + '_%s.simres' % string.zfill(k, 5)
                    self.dt.print_SimRNArestraints(ijsimres_flnm, thrds, "slope", True)
                    #
                #
            #
            
            # 161025wkd: try moving this step to here. I know this
            # means that the program must open this file each time and
            # deposit the conttents. The reason I do it this way is
            # because a single string of 100k structures could be many
            # million bytes. This may have been the cause of the
            # program crashing with a "killed" statement somewhere
            # between printing out the structures here, and printing
            # out the matrix and the summary file.

            # It was a strange bug because I had just successfully
            # calculated a heatmap of almost twice the size of this
            # one. Therefore, it must be the particular quantity of
            # output perhaps.
            
            # If this work around is successful but having to wait is
            # just too annoying, perhaps we can make this transaction
            # once every 10 times around (or something like that). We
            # must first wait and see if this fixes the issue
            try:
                fp = open(ssflnm, 'a')
                fp.write(file_results)
                fp.close()
            except OSError:
                print "ERROR: problems opening %s" % (ssflnm)
                sys.exit(1)
            #
            k += 1
        #

        if len(self.calc.fe.mtools.pssbl_ctcf) > 0:
            
            file_results = "strong binding sites: %d\n" % (len(self.calc.all_ctcf))
            for rh in self.calc.all_ctcf.keys():
                file_results += "(%4d, %4d), %10.0f\n" % (rh[0], rh[1], self.calc.all_ctcf[rh])
            prnt_results += file_results
        #
        
        # print the results to the terminal
        print prnt_results
        
        try:
            fp = open(ssflnm, 'a')
            file_results = "\n//\n" + file_results
            fp.write(file_results)
            fp.close()
        except OSError:
            print "ERROR: problems opening %s" % (ssflnm)
            sys.exit(1)
        #
        
        try:
            pclusters = self.calc.fe.mtools.disp_fmatrix(self.clust.clusters, "clusters", False)
            fp = open(flhd + "_BDwt.clust", 'w') # Boltzmann distribution weighted
            fp.write(pclusters)
            fp.close()
        except OSError:
            print "ERROR: problems opening %s" % (flhd + "_clustering.txt")
            sys.exit(1)
        #
    #    
#



####################################################################
####################################################################
####################################################################

# ###########################
# ######  MAIN DRIVER  ######
# ###########################

def reduce_list(t):
    # from
    # http://stackoverflow.com/questions/7961363/removing-duplicates-in-lists
    s = []
    for i in t:
        if i not in s:
            s.append(i)
    return s
#

def test0():
    # Test an algorithm for finding unique elements in a simple list
    
    t = [1, 2, 3, 1, 2, 5, 6, 7, 8]
    print "orig: ", t
    print "rev:  ", reduce_list(t)
#

def test1():

    # Test for verifying and demonstrating the syntax for the
    # sorting algorithm for handling groups of links
    smap = Map(30)
    
    i = 3
    j = 29
    
    link = Link(i, j, float(-5), 'B', 's', [(i,j)])
    smap.glink[i][j].add_link(link)
    
    link = Link(i, j, float(-15), 'I', 's', [(4,26)])
    link.add_Motif(i, j, float(-15), 'I', 's', [(5,27)])
    smap.glink[i][j].add_link(link)
    
    link = Link(i, j, float(-25), 'M', 's', [(4, 6), (8, 12), (15, 25)])
    link.add_Motif(i, j, 'P', 's', float(-25), [(4, 6), (8, 12), (15, 25)])
    smap.glink[i][j].add_link(link)
    
    
    print len(smap.glink[i][j].lg)
    for lgk in smap.glink[i][j].lg:
        print "next structure"
        for lpk in range(0, len(lgk.motif)):
            print lgk.motif[lpk].ctp, lgk.motif[lpk].get_base(), lgk.motif[lpk].get_branches()
        #
    #
    
    print "number of elements: ", len(smap.glink[i][j].lg)
    print "before sorting:"
    for lgk in smap.glink[i][j].lg:
        join = lgk.motif[0].get_branches()
        ctp  = lgk.motif[0].ctp
        Vij  = lgk.Vij
        print "%s    %8.2f: " % (ctp, Vij), join
    #
    smap.mergeSortLinks(smap.glink[i][j].lg)
    print "after sorting:"
    for lgk in smap.glink[i][j].lg:
        join = lgk.motif[0].get_branches()
        ctp  = lgk.motif[0].ctp
        Vij  = lgk.Vij
        print "%s    %8.2f: " % (ctp, Vij), join
    #
    
    
#    

def test2():
    debug_test2 = True
    sys.argv = ["chreval_upgrade160714.py", "-f", "test_degeneracy2.heat"]
    print sys.argv
    CL = GetOpts(PROGRAM)
    
    calc = Calculate(CL)
    T    = calc.T
    flnm = calc.flnm
    N    = calc.N
    dG, smap = calc.minFE(T)
    if debug_test2:
        print "number of iloop entries: ", len(calc.iloop.lg)
        for lgk in calc.iloop.lg:
            base = lgk.motif[0].get_base()
            ctp  = lgk.motif[0].ctp
            Vij  = lgk.Vij
            join = lgk.motif[0].get_branches()
            print base, ctp, Vij, join
        # sys.exit(0)
    #
    # Test for verifying and demonstrating the syntax for the
    # sorting algorithm for handling groups of links
    trace = Trace(calc)
    print "traceback mFE:"
    dGmin = trace.traceback_mFE(0, calc.N-1, 0, True)
    print "B: bound; I: I-loop, M: M-loop"
#    

    
def test5():
    global flag_KahanSumExp_USE_SORT
    flag_KahanSumExp_USE_SORT = True
    dG_vr_kBT = [10, 37, 34, 0.1, 0.0004, 34, 37.1, 37.2, 36.9, 709, 710, 711]
    Z, shift = KahanSumExp(dG_vr_kBT)
    print "{0} x 2^{1}".format(Z, shift)
    wt = float(2**shift)
    p = []
    for k in range(0, len(dG_vr_kBT)):
        p += [ exp( dG_vr_kBT[k] - float(shift)*log(2) ) / Z ]
    summation = 0.0
    for pk in p:
        summation += pk
        
    if shift > 1.0: # self.debug:
        print Z
        print "shift: ", shift, "new exp wt: ", wt
        for k in range(0, len(dG_vr_kBT)):
            print "%2d  %8.2f   %12.3g" % (k, dG_vr_kBT[k], p[k])
        print "summation of p: %12.5g" % summation
        # sys.exit(0)
    #
#    


def main(cl):
    CL = GetOpts(PROGRAM)
    
    #
    if SHOWMAIN:
        print cl
        print "number of args: %d" % (len(cl))
    #
    
    if len(cl) < 2:
        emsg = "ERROR: too few arguments"
        CL.error(emsg)
    #
    
    manager = Manager()
    manager.runCalculations(CL)
    manager.printResults()
    
#




if __name__ == '__main__':
    if TEST0:
        test0()
    elif TEST1:
        # for running tests of tools
        test1()
    elif TEST2:
        test2()
    elif TEST5:
        test5()
    else:
        # running the program
        main(sys.argv)

