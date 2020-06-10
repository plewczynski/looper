#!/usr/bin/env python3

"""@@@

Main Module:   FreeEnergy.py 

Classes:       FreeEnergy
               SetUpFreeEnergy

Author:        Wayne Dawson
creation date: parts 2016 (in chreval), made into a separate object 170426
last update:   200210 (upgraded to python3) 190718
version:       0
FreeEnergy.py 


Purpose:

This defines the basic free energy functions for the program. It is
intended to be adjustable like the ChromatinModule in that it can be
set up for use with whatever system is of interest be that chromatin,
RNA or proteins.


Comments:

Well, mainly that as of 190403, it is still directed to chromatin.

"""

from math import log
import sys
import os
import string

# ################################################################
# ######################  Global constants  ######################
# ################################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# A hardwired infinity. Somewhat hate this, but it seems there is no
# way around this darn thing.
from Constants import INFINITY
from Constants import kB   # [kcal/molK] (Boltzmann constant)

# coarse-grained resolution factor
from ChrConstants import seg_len
# for entropy calculations
from ChrConstants import xi   # [bps]
from ChrConstants import lmbd # [bps]
from ChrConstants import gmm  # dimensionless but related to D/2        
# constants: weights for the enthalpy terms 
from ChrConstants import febase # [kcal/mol]
from ChrConstants import feshift # (dimensionless, usually = 1)


# for chromatin, these are currently not set in GetOpts 
from ChrConstants import minStemLen     # always 1 bp for chromatin
from ChrConstants import max_aa_gap     # max gap between beads antiparallel
from ChrConstants import max_pp_gap     # max gap between beads parallel
from ChrConstants import minLoopLen     # [nt]
from ChrConstants import dGpk_threshold # [kcal/mol]
from ChrConstants import pk_scan_ahead  # [nt]
from ChrConstants import dGMI_threshold # [kcal/mol]
from ChrConstants import set_dangles    # not used 
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# This is set in GetOpts
from ChrConstants import dG_range 


from BasicTools import roundoff

# program processing routines
from SettingsPacket import InputSettings

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
TEST0        = True  # for testing of objects

# debugging options 
DEBUG = False  # stops program when encouters condition

# 2. special local global function
def usage():
    s  =  "USAGE: %s\n" % PROGRAM
    s +=  "       %s test0\n" % PROGRAM
    return s
#
# At present, there is not much reason to use this at all. The only
# thing that can be done with this object independent of the programs
# that use it is to run a test, and that test is not so interesting
# presently. Nevertheless, maybe it will be useful at some point.


# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ################################################################
# ################################################################



####################################################################
#####################  FREE ENERGY PARAMETERS ######################
####################################################################



class PairDef:
    """@
    
    190425: In vsfold, this had a similar role to ptype; it
    establishes the base pair index for AU, GC, and GU pairs (indices
    between 1 and 6) or zero if it is a nWC pair. In the RNA program,
    in particular, I needed to include more detail such as the stem
    length, the stacking free energy, the di-nucleotide/tri-nucleotide
    base pairs, etc. Here I also include 

    PairDef was moved here because we need the FE defined
    early on with RNA programs. The same should be done with chromatin
    to improve the relative relational functionality of the parts.
    
    """
    
    def __init__(self, p, b, c = 'X', dGp = 0.0):
        self.debug_PairDef =  False
        self.pair =  p
        self.btp  =  b
        self.ctp  =  c
        self.dGp  =  dGp
        """@
        
        Presently (190128, updated 190718), for chromatin, I have
        defined "pair" as either 1 or 0 (paired or unpaired). In
        principle, "pair" should specify some number indicating the
        type of contact. For example, perhaps I would indicate a pair
        involved in a weak singleton interaction as "1" and would
        indicate a strong CTCF interaction as "2". This would be
        useful if the CTCF positions are defined. Since that is not
        always the case, the only parameter currently used is paired
        (1) or unpaired (0).
        
        Since btp can either refer to an antiparallel stem or a
        parallel stem, it is essential. With RNA, there are also
        parameters for stem length and (because there is a real
        sequence) the actual bases in the structure.
        
        If the data is stored, and it seems reasonable that it could
        be for chromatin, then dGp would be the precalculated free
        energy for the respective pair associated with position (i,
        j). This, therefore, could easily take the place of self.hv
        and the incessent calculations that are currently done in
        Calculate() and ChromatinModule().
        
        I must emphasize that this class is still very much under
        development. Presently, I think more information probably
        doesn't hurt, particularly if there are stems detected.

        """

    #
    def __str__(self):
        s = "pair(%d), bond(%s), dGp(%8.2f)" % (self.pair, self.btp, self.dGp)
        return s
    #        
    
    def __repr__(self):
        return self.__str__()
    #
#


class old_SetUpFreeEnergy(object):
    # 200313: this should soon be considered obsolete

    """
    # This is used to set up FreeEnergy() when we don't have the
    # specific parameters.  It is used by test programs like test0()
    # in this modeule. .

    # It is not generally for calculation, but it is sometimes useful
    # for tests.
    """
    
    def __init__(self):
        self.source = "SetUpFreeEnergy"
        self.program = "None"
        self.system  = "Chromatin"
        
        # input variables
        self.T        = 310.0 # not used!!!
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
        
        self.dGrange = 10.0
        
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
        
    #
#

class SetUpFreeEnergy(InputSettings):
    """@ 

    This is used to set up FreeEnergy() when we don't have the
    specific parameters.  It is used by test programs like test0() in
    this module.
    
    """
    def __init__(self):
        
        InputSettings.__init__(self, "Chromatin")   # inherit InputSettings
        
        self.source = "SetUpFreeEnergy"

        if self.show_info:
            print ("dinged here SetUpFreeEnergy")
        #
        
    #
#


# this provides the free energy parameters
class FreeEnergy(object):
    """
    Contains the free energy methods for Chromatin structure prediction
    """
    def __init__(self, cl = SetUpFreeEnergy()):
        self.debug_FreeEnergy = False # True # 
        global febase
        
        
        self.source  = cl.source
        self.program = cl.program
        self.system  = cl.system
        
        ############################
        self.molsys  = cl.molsys
        ############################
        
        # chromatin
        self.f_heatmap = "noname.heat"
        self.chrseq    = ""   # e.g., "cccccccccccc"
        self.chrstr    = ""   # e.g., "((((....))))"
        self.cSeq      = None # Seq(self.chrseq)
        
        
        
        # input variables
        self.T        = cl.T # not used!!!
        
        # constants: in entropy evaluation
        self.kB       = kB         # [kcal/molK],  Boltzmann constant
        self.xi       = xi         # [bp] stem Kuhn length
        self.lmbd     = lmbd       # [bp], distance between bases in units of nt
        
        self.gmm      = gmm        # dimensionless self avoiding walk constant
        self.delta    = 2.0        # exponential weight
        self.nu       = 0.5        # excluded volume term
        
        self.zeta     = self.gmm + 0.5  # presently only Gaussian type functions
        self.seg_len  = cl.seg_len      # [bp], segment length (for RNA = 1 nt) 
        self.w        = self.seg_len*self.xi/(self.lmbd)**2
        #                            weight used in TdS
        
        self.kBT      = self.kB * self.T # [kcal/mol] kB*T
        
        # constants: weights for the enthalpy of binding (used with
        # free energy models like chromatin).
        
        self.base  = cl.dHbase    # dH(ij) = dHbase + ln(dHshift + counts(ij))
        self.shift = cl.dHshift
        if self.debug_FreeEnergy:
            print (self.base, self.shift)
            # sys.exit()
        #
        
        # secondary structure stem parameters
        self.minStemLen = minStemLen
        # minimum loop length (for chromatin it is 1, or RNA 3)
        self.max_aa_gap = max_aa_gap 
        self.max_pp_gap = max_pp_gap 
        """@
        
           max_aa_gap (maximum anti parallel bead gap): in joining two
           stems to form a connected anti-parallel stem, max_aa_gap is
           what simply cannot reasonably qualify anymore as a
           "connected stem" because the distance between the stems is
           simply too large for any reasonable coupling to be expected
           to occur.
        
           max_pp_gap (maximum parallel bead gap): in joining two
           stems to form a connected parallel stem.
        
        """
        
        # secondary structure loop parameters
        self.minLoopLen = minLoopLen
        # minimum loop length (for chromatin it is 1, for RNA = 3)
        self.dGMI_threshold = cl.dGMI_threshold # [kcal/mol]
        """threshold for MBL stability"""
        
        # pseudoknot parameters
        self.scan_ahead     = cl.scan_ahead   # default 10
        self.dGpk_threshold = cl.dGpk_threshold # [kcal/mol]
        
        # MBL parameters
        self.dG_NNbranch    = 1.0*(self.T/310.15) # [kcal/mol] currently fixed
        """dG_NNbranch is the cost of putting branches right next to each other"""
        
        """190130: 
        
        Something less than this is basically the same as the thermal
        energy, so it is probably rather weak. Anyway, it seems like
        -0.2 kcal/mol for 'B' is really weak and it seems unlikely
        that the stability of such a structure is really
        significant. At any rate, the program will work even if you
        set the threshold to zero.
        
        """
        
        self.dangles = set_dangles # not used!!!
        
        
        # free energy range (not used)
        self.dGrange = cl.dGrange
        
        
        # PET cluster weight
        self.add_PET_wt          = cl.add_PET_wt
        self.add_PET_wt_to_edges = cl.add_PET_wt_to_edges
        self.PETwt               = cl.PETwt
        
        # weights for selecting out PET clusters 
        self.CTCF_scale          = cl.CTCF_scale
        self.ctcf_tthresh        = 0.20*self.CTCF_scale
        self.ctcf_cthresh        = 0.40*self.CTCF_scale
        
        """@
        
        Note that these parameters are repeated here and in
        HeatMapTools. It would also be possible that FreeEnergy
        inherits HeatMapTools and extracts their parameters from it. I
        am still thinking about it.
        
        """
        
        
        self.N = -1
        self.hv = []    # heat map weight 
        self.btype = [] # PairDef() <=  [i][j] (main unit)
        """@
        
        Both hv and btype define the statistical weight of the pairing
        interaction, where hv is the actual weight and btype[].dGp
        carries the corresponding free energy. Since chromatin at 5
        kbp resolution is basically a pure polymer, I can use
        ptype[].dGp to store the free energy generated from hv. So
        this is now what I do in the program. Nevertheless, the
        heatmap data has to be read in in the form of hv.
        
        self.N, self.hv and self.btype are assigned after FreeEnergy
        is inherited by ChromatinModule. At that point, the input file
        is read and N, hv and btype assigned. To do this, one invokes
        a call to the method assign_btypes(), assuming that the
        program was called through GetOpts. Otherwise, it has to be
        set in some other way.
        
        """
        
        
        
        self.ctcf_setv = {}
        """@
        
        created boundary weight for cases where CTCF is assumed at the
        edges and this data is missing in the input data. I don't
        think that this should be used anymore, or, at best,
        sparingly. However, I leave it here for the moment.
        
        """
        
        self.exist_heatmap = False
    #
    
    # Boltzmann thermal energy
    def g_kBT(self):
        return self.kBT
    #
    
    # entropy calculation
    def TdS(self, i, j, T):
        if j <= i:
            print ("ERROR: i(%d) >= j(%d)!!!" % (i,j))
            sys.exit(1)
        #
        
        # calculate the CLE
        n = float(j - i + 1)
        ps_n = self.w*n # xi*N/lmbd**2, where N = n * self.seg_len
        # math.log
        tds = self.kB*T*(self.gmm*log(ps_n) - (self.gmm+0.5)*(1.0 - 1.0/ps_n))/self.xi
        return tds
    #
    
    def initialize_btype(self, N):
        """
        This has the form and function of ptype in vsfold5, but it allows
        for more options than the original marker.
        """ 
        for j in range(0, N):
            btype_j = []
            for i in range(0, N):
                btype_j += [PairDef(0, '-', 'X', 0.0)]
            #|endfor
            
            self.btype += [btype_j]
        #|endfor
        
    #
    
    
    def add_btype(self, vs):
        """@
        
        Presently, I believe this is not used. The idea here is that
        you can read in some 1D structural sequence (e.g.,
        "..((((....))))...((...A.)).a.") and generate this particular
        mapping as the desired interaction (perhaps even along with
        the pairing weight measured in free energy).  
        
        """
        if len(self.btype) == 0:
            print ("reset btype:")
            self.initialize_btype(vs.N)
        #
        
        for x in vs.BPlist:
            print (x.i, x.j, x.v, x.name, x.contacts)
            if x.v == 'a':
                self.btype[x.i][x.j] = PairDef(1, 'ap', 'B')
            elif x.v == 'p':
                self.btype[x.i][x.j] = PairDef(1, 'pp', 'B')
            #
            
        #
        
    #
    
    
    # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
    # VVVVVVVVVV   FE/Enthalpy based on Weight functions   VVVVVVVVVV 
    # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
    
    # enthalpy of binding calculation
    def dH(self, v_ij):
        # print (v_ij)
        # math.log
        dh = self.base - log(self.shift + float(v_ij))
        return dh
    #
    
    def calc_dG(self, i, j, hp, T, local = "undefined"):
        flag_debug = self.debug_FreeEnergy
        if hp < 0.0:
            print ("ERROR(FreeEnergy.calc_dG()): chromatin interaction point at (%d,%d)" % (i, j))
            print ("                             is empty. Must be greater or equal to zero!")
            print ("                             call point: %s" % local)
            sys.exit(1)
            
        else:
            dGij = self.TdS(i,j, T)
            if hp > 0.0:
                dGij += self.dH(hp)
            #
            
        #
        
        if flag_debug:
            print ("calc_dG: dG(%2d,%2d)[hv(%4.0f)] = %8.3f, [dH(%8.3f),TdS(%8.3f)]" % \
                (i, j, hp, dGij, self.dH(hp), self.TdS(i,j, T)))
        #
        
        return dGij
    #
    
    def initialize_hv(self, N):
        """
        Presently, this is mostly needed when attempting debugging of a
        subsection of a program. It is used to build the weighted
        matrix self.hv that is usually obtained by reading heatmap
        files.
        """ 
        self.hv = []
        for j in range(0, N):
            hv_j = []
            for i in range(0, N):
                hv_j += [0.0]
            #|endfor
            
            self.hv += [hv_j]
            
        #|endfor
        
    #
    
    # identify the properties of the chromatin by the magnitude of the weight
    def get_bondtype(self, hvij):
        bt = 's'
        if hvij > self.ctcf_cthresh:
            bt = 'c' # default "tandem" setting
        elif hvij > self.ctcf_tthresh:
            bt = 't' # default "convergent" setting
            
        else:
            bt = 's'
        #
        
        return bt
    #
    
    
    def read_heatmap_from_file(self, flnm):
        self.hv, self.ctcf_setv, self.N = self.read_heat(flnm)
        # Note: read_heat is inherited from HeatMapTools
        
        # make a bogus sequence for chromatin
        self.sseq = ''
        for k in range(0, self.N):
            self.sseq += 'c'
        #
        
        self.exist_heatmap = True
    #
    
    def build_heatmap_from_vseq(self, vs):
        """@
        
        Note: This was originally called "make_heatmap". However,
        HeatMapTools also has a method make_heatmap that also reads in
        structures (and is inherited by FreeEnergy and, therefore,
        shadowed!). 
        
        The aim of HeatMapTools' make_heatmap is to develop actual
        heatmaps of any type given a collection of input sequences. It
        is called by my_generation.py. Here, in these operations, we
        are only building a map for _one_ specific sequence in the
        form of a heatmap representation. I'm not even sure that most
        of this is necessary except for the fact that assign_btypes()
        needs this self.hv assigned so that it can map out btypes and
        free energy related information.
        
        """
        
        self.N = vs.N
        self.initialize_hv(self.N)
        #print ("size of heatmap %d**2" % self.N)
        
        # secondary structure contacts
        print ("secondary structure contacts:")
        for v in vs.BPlist:
            #print (v)
            self.hv[v.i][v.j] = 3
            print ("(%3d,%3d)  %5d" % (v.i, v.j, self.hv[v.i][v.j]))
        #
        
        # related pseudoknot structure contacts
        print ("pseudoknot contacts:")
        for v in vs.PKlist:
            self.hv[v.i][v.j] = 3
            print ("(%3d,%3d)  %5d" % (v.i, v.j, self.hv[v.i][v.j]))
        #
        
        # multiloop contacts,
        print ("CTCF island contacts")
        for w in vs.MPlist:
            wb = w.i; we = w.j
            self.hv[wb][we] = 50
            print ("(%3d,%3d)  %5d" % (wb, we, self.hv[wb][we]))
            for k in range(0, len(w.contacts)):
                
                ck = w.contacts[k]
                self.hv[wb][ck] = 10
                self.hv[ck][we] = 10
                print ("(%3d,%3d)  %5d" % (wb, ck, self.hv[wb][ck]))
                print ("(%3d,%3d)  %5d" % (ck, we, self.hv[ck][we]))
                
                for l in range(k+1, len(w.contacts)):
                    cl = w.contacts[l]
                    self.hv[ck][cl] = 5
                    print ("(%3d,%3d)  %5d" % (ck, cl, self.hv[ck][cl]))
                #|endfor
                
            #|endfor
            
        #endfor
        
        self.exist_heatmap = True
    #
    
    
    
    def assign_btypes(self):
        """@
        
        This is definitely used by chreval in ChromatinModule to read
        in a heatmap (either heat or eheat). It requires HeatMapTools
        to function properly.
        
        Presently called when the command line arguments call "GetOpts".
        
        """
        debug_assign_btypes = False # True # 
        if debug_assign_btypes:
            print ("Entered: assign_btypes")
        #
        
        if not self.exist_heatmap:
            print ("ERROR: missing heat map")
            print ("       must call read_heatmap_from_file(flnm) ")
            print ("       or build_heatmap_from_vseq(class Vienna) first")
            sys.exit(0)
        #
        
        
        # make a bogus typedefinition for position ij
        self.btype = []
        self.initialize_btype(self.N)
        for j in range(0, self.N):
            for i in range(0, self.N):
                if self.hv[i][j] > 0.0:
                    # 190524 was self.btype[i][j] = PairDef(1, btp, 'B', self.hv[i][j])
                    if   i < j:
                        btp = self.get_bondtype(self.hv[i][j])
                        dGij = self.calc_dG(i, j, self.hv[i][j], self.T)
                        self.btype[i][j] = PairDef(1, btp, 'B', dGij)
                        
                    elif j < i:
                        btp = self.get_bondtype(self.hv[i][j])
                        dGij = self.calc_dG(j, i, self.hv[i][j], self.T)
                        self.btype[i][j] = PairDef(1, btp, 'B', dGij)
                        
                    else:
                        print ("ERROR: for bonds i(%d) = j(%d) is not allowed!" % (i, j))
                        sys.exit(1)
                    #
                    
                else:
                    if   i < j:
                        dGij = self.calc_dG(i, j, self.hv[i][j], self.T)
                        self.btype[i][j] = PairDef(0, '-', 'X', dGij)
                        
                    elif j < i:
                        dGij = self.calc_dG(j, i, self.hv[i][j], self.T)
                        self.btype[i][j] = PairDef(0, '-', 'X', dGij)
                        
                    else:
                        self.btype[i][j] = PairDef(0, '-', 'X', 1000.0)
                    #
                    
                #
                
            #|endfor i
            
        #|endfor j
        
        
        
        
        # make a complete list of all possibilities without redundancies
        self.all_ctcf = {}
        self.all_ctcf.update(self.ctcf_setv)
        for pc in self.pssbl_ctcf.keys():
            self.all_ctcf.update({pc: self.pssbl_ctcf[pc]})
        #
        
        self.all_ctcf.update(self.edge_ctcf)
        
        if debug_assign_btypes:
            print ("ctcf_setv:        ", self.ctcf_setv)
            print ("pssbl_ctcf:       ", self.pssbl_ctcf)
            print ("edge_ctcf:        ", self.edge_ctcf)
            print ("_________________________________________")
            print ("all_ctcf:         ", self.all_ctcf)
            print ("Exiting: assign_btypes")
            
            
            sys.exit(0)
        #
        
        return self.N
        
    #
    
    
    def add_hv(self, vs):
        """
        This is used with the map building function to generate heatmaps
        from constructed sequences
        """
        if len(self.hv) == 0:
            print ("reset hv:")
            self.initialize_hv(vs.N)
        #

        print ("add_hv: len(hv)", len(self.hv))
        
        for x in vs.BPlist:
            print (x.i, x.j, x.v, x.name, x.contacts)
            self.hv[x.i][x.j] = 1.0
        #
        sys.exit(0)
        
    #
    
    # new_cle_make_ptypes
    def build_btype_from_hm(self, hm):
        """@

        This is under construction!!!!! It is not finished!!!!
        
        Builds the enthalpic contribution for stems so these
        contributions don't need to be calculated every time. If the
        the Kuhn length is constant, then the entropic contribution
        from the cle could be calculated here, like it is in vsfold;
        however, we are aiming at a variable Kuhn length, which is
        context dependent ad cannot be decided here. Therefore, this
        section will only evaluate the base pairing weights.
        
        """
        
        debug_build_btype_from_hm = True # False # 
        
        stem_len = 0
        
        btp_test  = 0
        has_stems = False
        has_bps   = False 
                
        if debug_build_btype_from_hm:
            print ("Enter build_btype_from_hm(N = %d)" % self.N)
        #endif
        
        
        # initialize ptype, phstem and ptstem 
        if debug_build_btype_from_hm:  
            print ("initialize btype")
        #endif
        self.initialize_btype(self.N)
        
        if debug_build_btype_from_hm:  
            print ("PairType: scanning for stems.... ")
        #endif
        
        
        # need to think whether it is minLoopLen + 1 or + 2
        for j in range(self.minLoopLen + 1, self.N):
            for i in range(0, j - self.minLoopLen - 1):
                
                hv_test = self.hv[i][j] 
                
                stem_len = 0;
                
                if hv_test > 0:
                   
                    if self.btype[i][j].pair == 0: # not assigned yet
                        if debug_build_btype_from_hm:  
                            print ("ij(%2d,%2d), btype.pair(%d)" \
                                % (i, j, hv_test))
                        #endif
                        
                        """@
                        
                        190517): Here we have to handle both parallel
                        and anti-parallel stems, where with RNA, we
                        only had to consider anti-parallel stems. For
                        the given pair (i,j), we search here for both
                        the anti-parallel type and parallel type
                        stems. 
                        
                        
                        anti-prallel stems    .              parallel stems
                                              .         
                                ___           .         
                               /   \          .          
.                              \___/          .                  _______.
                            p. |___| .q       .               /          \
                               |___|          .               |           |
                               |___|          .               |  i     p  /
                               |___|          .          _____|__.______ .
                               |___|          .               *__\_\_\_\_________  
                            i. |___| .j       .                  q      j
                                              .
                        
                        
                        The coordinates i and j should always be
                        understood as the _input_ and _output_ of the
                        stem. The pair (i,j) represent the motif's
                        coordinates, and then _inside_ of that object
                        defining the motif, the stem is described as
                        pairs. Therefore, we need to implement the
                        object Stem, complete with its definitions to
                        make this fully functional.
                        
                        I think another question I presently have is
                        whether I can carry this approach beyond just
                        a single set of contiguous contacts in the
                        parallel case. In other words, 
                        
                        ..ABCD..EFGH...IJK...L......abcd...efgh...ijk...l...

                        Overall, in the above example, the coordinate
                        i will correspond to 'A', p will correspond to
                        'L', q will correpond to 'a' and j will
                        correspond to 'l'. However, the program will
                        try to divide these up in general.

                        
                        For the anti-parallel case, I am quite clear,
                        but the parallel case, I am not.

                        """
                        
                        # anti parallel stems
                        stem_len = self.findStemLen(S, i, j)
                        
                        if stem_len >= self.minStemLen: # enough to be a stem
                            p = i + stem_len - 1;
                            q = j - stem_len + 1;
                            self.btype[p][q].h2tlen = stem_len;
                            
                            
                            for k in range(0, stem_len - self.minStemLen + 1):
                                # assign the head region
                                #self.btype[p-k][q+k].h2tlen = (stem_len - k)
                                # assign the tail region
                                #self.btype[i+k][j-k].t2hlen = (stem_len - k)
                                
                                # formation of a pair: here, we simply
                                # assume any pair to be type 1.
                                hvijk = self.hv[i+k][j-k]
                                self.btype[i+k][j-k].pair = 1
                                self.btype[i+k][j-k].bpt  = 'sa'
                                self.btype[i+k][j-k].ctp  = 'S';
                                self.btype[i+k][j-k].dGp \
                                    = self.calc_dG(i+k, j-k, hvijk, T,
                                                   "build_btype_from_hm")
                                
                            #|endfor
                            
                            
                            if debug_build_btype_from_hm:
                                s = "S(%2d,%2d)[len(%d)]: " % (i, j, stem_len)
                                for k in range(0, stem_len):
                                    s += " (%2d,%2d)[%d] " \
                                         % (i+k, j-k, self.btype[i+k][j-k].pair) 
                                #|endfor
                                
                                s += "\n"
                                print (s)
                            #endif
                            
                            # 180428: deleted all that connect stem search stuff.
                        #    //  end: if stem_len >= self.minStemLen: #
                    #        //  end: if self.btype[i][j].pair == 0:
                #            //  end: if hv_test > 0:
            #                //  end: for i in range(0, i < j - self.minLoopLen - 2):
        #                    //  end: for j in range(self.minLoopLen+2, self.N):
        
        
        
        
        if debug_build_btype_from_hm:  
            print ("final check" )
        #endif
        
        """@
        
        final check to see if there is any pairing at all
        
        Verifies that at least one part of the RNA sequence has a stem
        that satisfies the rules of Vsfold5 and the input parameters
        from the user.
        
        """
        
        
        for j in range(self.minLoopLen, self.N):
            for i in range(0, j):
                if self.btype[i][j].t2hlen > 1:
                    has_stems = True
                    has_bps   = True # follows logically
                    break
                #
                
            #|endfor
            
        #|endfor
        
        if not has_stems:
            for j in range(self.minLoopLen, self.N):
                for i in range(0, j):
                    if self.btype[i][j].pair != 0:
                        has_bps = True
                        break
                    #
                    
                #|endfor
                
            #|endfor
            
        #
        
        status = 0
        if not (has_stems and has_bps):
            status = 9999;
            print ("ERROR!: no significant pairing found for a\n")
            print ("        minimum stem length %d bps.\n" % (self.minStemLen))
            print ("        Output will be questionable\n")
            sys.exit(1)
            
        else: 
            if debug_build_btype_from_hm:
                print ("minStemLen = %d" % (self.minStemLen))
            #endif
            
        #
        
        for j in range(self.minLoopLen, self.N):
            for i in range(0, j - self.minLoopLen):
                self.btype[j][i] = self.btype[i][j]
            #|endfor
            
        #|endfor
        
        if debug_build_btype_from_hm:
            """@
            
            This is rarely used, in part because it is very cumbersome
            to display. If the sequence length is longer than 40 nt,
            this should probably be avoided because it becomes very
            difficult to read, understand and interpret. 
            
            I leave this stuff around because these programs are still
            being developed and it is not impossible that it will be
            needed for some problems in the future.
            
            In my experience, when I really wanted to see this data,
            my mind was not concerned about calculation results!!!
            Therefore, the buck stops here.
            
            """
            taglist = ["t2hlen", "h2tlen", "ctp", "pair"]
            for tl in taglist:
                print (self.display_btype(tl))
            #
            
            
            print ("Planned exit after evaluating build_btype_from_hm()\n")
            sys.exit(0);  
        #endif
        
        # set values in VsfoldData
        return 0
    #
    
    def findStemLen(self, S, i, j):
        debug_findStemLen = False # True # 
        
        if debug_findStemLen:
            print ("Enter debug_findStemLen(ij(%d,%d)" % (i, j))
        #
        
        p = i; q = j; stack_pq = 0;
        ntype = get_bspr(S, p, q, debug_findStemLen)
        while not ntype == 0:
            p += 1; q -= 1
            if p < self.N and q > 0:
                if q - p > self.minLoopLen:
                   ntype = get_bspr(S, p, q, debug_findStemLen)
                   if debug_findStemLen:
                       print ("pq(%2d, %2d), ntype = %d" % (p, q, ntype))
                    #
                
                
                else:
                    ntype = 0;
                    self.btype[p][q].pair = 0 
                    self.btype[p][q].btp  = '-' 
                    if debug_findStemLen:
                        print ("pq(%2d, %2d), set ntype = %d" % (p, q, ntype))
                    #
                    
                #
                
            else:
                ntype = 0;
            #
            
            stack_pq += 1
            
        #|endwhile
        
        if debug_findStemLen:
            print ("findStemLen(%d,%d) = %d" % (i, j, stack_pq))
        #endif
        
        return stack_pq
    #
    
    
    # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    # AAAAAAAAAA   FE/Enthalpy based on Weight functions   AAAAAAAAAA
    # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA


# end of class FreeEnergy



def test0():
    cl  = SetUpFreeEnergy()
    fe  = FreeEnergy(cl)
    Tds = fe.TdS(0, 1, cl.T)
    
    print ("weight    enthalpy       entropy     free energy")
    print ("         [kcal/mol]     [kcal/mol]    [kcal/mol]")
    
    for dw in range(0,10,1):
        dh = fe.dH(float(dw))
        print ("%4d     %8.3f      %8.3f       %8.3f" % (dw, dh, Tds, (dh + Tds)))
    #
#



def main(cl):
    test0()
#

if __name__ == '__main__':
    """
    run various testing programs
    """
    
    # print (len(sys.argv))
    if len(sys.argv) > 1:
        opt = sys.argv[1]
        if opt == "test0":
            test0()
        else:
            print ("ERROR: unrecognized option (%s)" % opt)
            print (usage())
            sys.exit(1)
        #
    else:
        main(sys.argv)
    #
#
