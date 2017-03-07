#!/usr/bin/env python

# FreeEnergy.py 
#



# Presently, I think the name is rather unoriginal (i.e., "CHRomatin
# EVALuate"). Anyway, what is important is that it works, and maybe we
# can come up with a better name later.


from math import exp, log, ceil, floor
import sys
import os
import string
from MatrixTools import MatrixTools
from GetOpts import GetOpts
from Threads import DispThreads
from Threads import LThread
from Functions import KahanSumExp

# ################################################################
# ######################  Global constants  ######################
# ################################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# coarse-grained resolution factor
from Constants import seg_len
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
# ############  Local global functions and constants  ############
# ################################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# 1. control parameters for the program

PROGRAM      = "FreeEnergy.py"  # name of the program

# tests
TEST0        = False  # for testing of objects

# debugging options 
DEBUG = False  # stops program when encouters condition

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
#####################  FREE ENERGY PARAMETERS ######################
####################################################################

# this provides the free energy parameters
class FreeEnergy:
    def __init__(self, cl):
        self.debug_FreeEnergy = False
        global febase
        
        # input variables
        self.T        = cl.T
        self.N        = -1
        
        # constants: in entropy evaluation
        self.kB       = kB
        self.xi       = xi
        self.lmbd     = lmbd
        self.gmm      = gmm
        self.seg_len  = cl.seg_len
        self.w        = self.seg_len*self.xi/(self.lmbd)**2
        
        # constants: weights for the enthalpy of binding
        # Enthalpy base
        self.base  = cl.dHbase    # dH(ij) = dHbase + ln(dHshift + counts(ij))
        self.shift = cl.dHshift
        if self.debug_FreeEnergy:
            print self.base, self.shift
            # sys.exit()
        #
        
        # PET cluster weight
        self.add_PET_wt          = cl.add_PET_wt
        self.add_PET_wt_to_edges = cl.add_PET_wt_to_edges
        self.PETwt               = cl.PETwt
        
        self.CTCF_scale          = cl.CTCF_scale
        self.ctcf_tthresh        = 0.20*self.CTCF_scale
        self.ctcf_cthresh        = 0.40*self.CTCF_scale
        
        # This started when I was confronted with the somewhat strange
        # data I got from Nenski that had counts in almost every bin
        # and very large numbers.
        
        self.pssbl_ctcf          = {}
        self.edge_ctcf           = {}
        self.from_Nenski         = cl.from_Nenski
        
        self.mtools              = MatrixTools()
        self.mtools.ctcf_tthresh = self.ctcf_tthresh
        self.mtools.ctcf_cthresh = self.ctcf_cthresh
        self.mtools.from_Nenski  = cl.from_Nenski
        
        
        # this re-scales the data by some fraction
        self.rescale_wt          = cl.rescale_wt

        self.allowed_extns       = cl.allowed_extns
    #
    
    # entropy calculation
    def TdS(self, i, j, T):
        if j <= i:
            print "ERROR: i(%d) >= j(%d)!!!" % (i,j)
            sys.exit(1)
        #
        # calculate the CLE
        #
        n = float(j - i + 1)
        ps_n = self.w*n # xi*N/lmbd**2, where N = n * self.seg_len
        # math.log
        tds = self.kB*T*(self.gmm*log(ps_n) - (self.gmm+0.5)*(1.0 - 1.0/ps_n))/self.xi
        return tds
    #
    
    # enthalpy of binding calculation
    def dH(self, v_ij):
        # print v_ij
        # math.log
        dh = self.base - log(self.shift + float(v_ij))
        return dh
    #
    
    def calc_dG(self, i, j, hp, T, local = "undefined"):
        flag_debug = self.debug_FreeEnergy
        if hp <= 0.0:
            print "ERROR(FreeEnergy.calc_dG()): chromatin interaction point at (%d,%d)" % (i, j)
            print "                             is empty. Must be a non-zero value!"
            print "                             call point: %s" % local
            sys.exit(1)
        #
        dGij = self.TdS(i,j, T) + self.dH(hp)
        if flag_debug:
            print "calc_dG: dG(%2d,%2d)[hv(%4.0f)] = %8.3f, [dH(%8.3f),TdS(%8.3f)]" % \
                (i, j, hp, dGij, self.dH(hp), self.TdS(i,j, T))
            #
        #
        return dGij
    #
    
    
    # ####################################################################
    # it might work to put this next function and its partner in
    # MatrixTools. It seems like most outputs will be matrix data, so
    # it makes sense to generate similar style outputs and to have one
    # program that reads them and processes them accordingly.
    # ####################################################################
    # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    
    # This reads in the heatmap data and converts it into something
    # that can be used for calculating the enthalpy with this data.
    def read_heat(self, flnm):
        
        # read in the matrix data
        gmtrx    = self.mtools.read_MatrixFile(flnm, self.allowed_extns, self.rescale_wt)
        self.N   = gmtrx.length
        hv       = gmtrx.heatmap
        clusters = gmtrx.clusters
        self.pssbl_ctcf = self.mtools.pssbl_ctcf # potential CTCFs
        self.edge_ctcf  = self.mtools.edge_ctcf  # the edge CTCF
        
        # ####  Distinguishing between singletons and PET clusters  ####
        
        # From here, we have to decide about PET clusters. Presumably,
        # when the matrix elements are something like 1 to 10, we can
        # assume that the interactions are probably singletons. The
        # CTCF clusters are mainly distinguished in having a much
        # greater interaction frequency, basically numbers greater
        # than 10.
        
        # With respect to CTCF clusters, these are then divided into
        # convergent structures ('c': the most common and therefore
        # strongly bound structures) and tandem right or tandem left
        # structures ('t': encountered less frequently and not as
        # strongly bound).  The tandem CTCF PET clusters are assumed
        # to range from about 10 to 40. The convergent CTCF PET
        # clusters should have numbers greater than 40.
        
        # Presently, the program only distinguishes the PET cluster by
        # the magnitude of the information. Since the tandem right and
        # left structures have the same (or similar) frequency, these
        # are assigned the label 't' in this section for "tandem". If
        # information can be obtained elsewhere about the PET cluster,
        # this can be assigned later.  In general, this must be
        # deduced from the sequence alone.
        
        # Therefore, presently, all I can do is add a second set of
        # labels that specify the character of the
        # information. Presently, this is 's' for "singleton", 'c' for
        # "convergent" and 't' for "tandem" CTCF PET clusters.
        
        # location of the CTCF clusters 
        cv = [] # nothing present
        
        # here, we can search the file using the rules above.
        
        if self.debug_FreeEnergy:
            print self.mtools.disp_fmatrix(hv, "enthalpy: ")
        #
        
        # #####  READJUSTMENT FOR PET CLUSTER WEIGHTS  #####
        # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        
        # Here, I presume that the end point it tied with CTCF (PET
        # clusters). Since the general intensity of these interactions
        # is 0 to 10, this is just a guess on the weight from the CTCF
        # binding sites, but it comes out only to 4 kcal/mol, so
        # setting it to 100 does not seem that horrible an idea in my
        # opinion.

        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        # 160914wkd: In my opinion, this should not be used anymore,
        # but I leave it here anyway for the moment. At any rate, the
        # reason that was initially put here was because there seemed
        # to be no PET data, but now there generally is, so this
        # section is unnecessary.
        
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
        if self.add_PET_wt:
            # 
            if self.add_PET_wt_to_edges:
                # 160602wkd: Initially, I thought that the location of
                # the edge for the PET cluster was always at (0,N-1).
                # However, it seems my understanding that this was not
                # the correct. There can be some overhang in assigning
                # the 1 kbp range where the PET cluster could be in
                # bin (0,N-2) or bin (0,N-1). Possible other issues
                # exist. Nevertheless, because it adds some additional
                # options for debugging and other things, I decided to
                # leave the option here.
                
                hv[0  ][self.N-1] = self.PETwt
                hv[self.N-1][0  ] = self.PETwt
                cv = {(0, self.N-1) : self.PETwt }
                print "Note: Additional weight added to PET cluster"
                print "      borders at (%d,%d)" % (0, self.N-1)
            else:
                # 160602wkd: This is probably a more accurate
                # assumption, one which allows for the possibility of
                # the edge being located at either (0,N-2) or (0,N-1),
                # and allows further free play than that.
                
                hv, cv = self.mtools.add_boundaryWt(self.N, hv, self.PETwt)
            # 
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # #####  READJUSTMENT FOR PET CLUSTER WEIGHTS  ##### 
        
        if self.debug_FreeEnergy:
            print disp_fmatrix(hv, "enthalpy: ")
        #
        return hv, cv, self.N
    #
        
# end of class FreeEnergy


class FESetUp:
    # !!!!!!!!! NOTE !!!!!!!!!

    # This is very similar to the class FreeEnergy (above); however, I
    # constructed this here to provide the necessary tools for
    # test3. In test3, I am simply testing the free energy function
    # under rather specific conditions which I am likely to change,
    # and I really don't care about anything else in the constructor
    # of FreeEnergy. In fact, all the other aspects of the FreeEnergy
    # constructor are actually a nuisance. I don't need to build all
    # the usual stuff and therefore, I also don't need to first go
    # through GetOpt. This is only for here for running test3 and
    # under no circumstance should this be used to actually replace
    # FreeEnergy in general.

    # !!!!!!!!!------!!!!!!!!!
    def __init__(self):
        # input variables
        self.T        = 300.0
        self.N        = -1
        
        # constants: in entropy evaluation
        self.kB       = kB
        self.xi       = xi
        self.lmbd     = lmbd
        self.gmm      = gmm
        self.seg_len  = seg_len
        self.w        = self.seg_len*self.xi/(self.lmbd)**2
        
        # constants: weights for the enthalpy of binding
        self.dHbase  = febase
        self.dHshift = feshift
        
        # PET cluster weight
        self.add_PET_wt          = False
        self.add_PET_wt_to_edges = False
        self.PETwt               = 100.0
        
        self.CTCF_scale          = 100.0
        self.ctcf_tthresh        = 0.20*self.CTCF_scale
        self.ctcf_cthresh        = 0.40*self.CTCF_scale
        
        # This started when I was confronted with the somewhat strange
        # data I got from Nenski that had counts in almost every bin
        # and very large numbers.
        
        self.pssbl_ctcf          = {}
        self.edge_ctcf           = {}
        self.from_Nenski         = False
        
        self.mtools              = MatrixTools()
        self.mtools.ctcf_tthresh = self.ctcf_tthresh
        self.mtools.ctcf_cthresh = self.ctcf_cthresh
        self.mtools.from_Nenski  = self.from_Nenski
        
        
        # this re-scales the data by some fraction
        self.rescale_wt          = 1.0

        # other handling matters
        self.allowed_extns       = ["heat", "data"]
#

def test0():
    cl  = FESetUp()
    fe  = FreeEnergy(cl)
    Tds = fe.TdS(0, 1, cl.T)
    
    print "weight    enthalpy       entropy     free energy"
    print "         [kcal/mol]     [kcal/mol]    [kcal/mol]"
    
    for dw in range(0,10,1):
        dh = fe.dH(float(dw))
        print "%4d     %8.3f      %8.3f       %8.3f" % (dw, dh, Tds, (dh + Tds))
    #
#

def main(cl):
    test0()
#

if __name__ == '__main__':
    # running the program
    main(sys.argv)

