#!/usr/bin/env python3

"""@@@

Main Module:   SettingsPacket.py 

Classes:       InputSettings

Author:        Wayne Dawson
creation date: 190829
last update:   200603 (various updates, minor bug fixes)
version:       0.0


Purpose:

Works as a packet of default settings for class FreeEnergy and
BranchEntropy. 

Puts the default settings in one place and allows the
user to change them before initializing the computational programs.

"""

import sys
from os import path

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
from FileTools   import FileTools
from FileTools   import getHeadExt

from Seq import Seq

# set up Vienna and gMatrix and other FE parameterS
from MolSystem import genRNASeq    # make fake RNA sequences from structures
from MolSystem import genChrSeq    # make fake Chromatin seq from structures
from MolSystem import sysDefLabels # system dictionary: RNA, Chromatin
from MolSystem import dEXTS        # allowed parameter file extension dictionary

from MolSystem import MolSystem    # system information packet

#from ViennaDataFmt import ViennaDataObj # Vienna parameter representation for RNA
#from ViennaParFiles import ParamData    # Vienna parameter file for RNA
#from Proc_gMatrix import Build_gMatrix  # gMatrix parameter representation for RNA

# system independent constants
from Constants import kB   # [kcal/molK] (Boltzmann constant)
from Constants import T37C # [K] absolute temp at 37 C.
from Constants import T_0C # [K] absolute temp at  0 C.

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PROGRAM = "SettingsPacket.py"


#

class InputSettings(object):
    def __init__(self, system = "Chromatin"):
        
        self.show_info = False # True # 
        
        # initially, we enter in some bogus names and paramters so
        # that they are defined.
        
        self.system  = "unassigned" # left undefined because needs to be checked
        
        # presently, these two are mainly used for referencing if
        # problems occur.
        self.source  = "InputSettings"  
        self.program = "fake"
        self.jobtype = "undefined"
        
        # chromatin
        self.f_heatmap = "noname.heat"
        self.chrseq    = ""   # e.g., "cccccccccccc"
        self.chrstr    = ""   # e.g., "((((....))))"
        self.cSeq      = None # Seq(self.chrseq)
        
        # RNA
        self.f_vseq = "noname.vseq"
        self.rnaseq = ""   # e.g., "GGGGuuuuCCCC"
        self.rnastr = ""   # e.g., "((((....))))"
        self.rSeq   = None # Seq(self.rnaseq)
        
        self.N      = 0    # e.g., len(self.rnaseq)
        
        # ####################################################################
        # #####  Command line control parameters that affect FreeEnergy  #####
        # ####################################################################
        # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
        
        # set the basic free energy parameter set: Turner Energy Rules
        # or Vis Energy Rules
        
        self.molsys       = MolSystem()
        
        self.parFlnm      = "none"
        self.vdf          = None # ViennaDataObj
        self.bm           = None # Build_gMatrix
        
        
        # thermodynamic parameters in entropy and free energy
        # evaluation
        
        self.kB       = kB    # [kcal/molK] Boltzmann constant
        
        # input variables
        self.T        =  T37C # [K] temperature
        
        
        #############################################################
        ####  parameters requiring settings according to system  ####
        #############################################################
        # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
        
        # constants: in entropy evaluation
        # coarse-grained resolution factor (segment length)
        self.seg_len = 0
        
        self.xi       = 0.0  # [bp] stem Kuhn length
        self.lmbd     = 0.0  # bp separation distance in units of mer-to-mer distance
        self.gmm      = 0.0  # dimensionless self avoiding walk constant
        self.delta    = 0.0  # exponental weight (related to excluded volume)
        
        
        
        """@
        
        # next: entropic weight function
        
        With Chromatin, the Kuhn length (xi) is basically ignored
        except for its entropy weight in the global entropy
        function. As a result, I have set the entropic weight to be
        
        self.w = self.seg_len*self.xi/(self.lmbd)**2
        
        and this is used throughout in ChromatinModule.
        
        However, with RNA, we have a little problem because the value
        of xi depends on the proposed length of the stems in this new
        approach (this new approach is also different from vsfold5
        where we could pretend that the Kuhn length was fixed, set it,
        and be done with it). For RNA, I think the entropic weight
        should be
        
        self.w = self.seg_len/(self.lmbd)**2
        
        It is in fact true that I am including the actual xi of the
        stem. Therefore, xi must be specified as an argument in the
        calculation of the global entropy for RNA.
        
        """
        
        self.w        = 0.0
        
        
        
        # secondary structure stem parameters
        self.minStemLen     = 0 # [bp]
        self.max_bp_gap     = 0 # [nt]
        self.max_aa_gap     = 0 # [beads]
        self.max_pp_gap     = 0 # [beads]
        
        
        # secondary structure loop parameters
        self.minLoopLen     = 0    # [nt]
        # minimum loop length (for chromatin it is 1, for RNA = 3)
        self.dGMI_threshold = 0.0 # [kcal/mol]
        """threshold for MBL stability"""
        
        # pseudoknot parameters
        self.minPKloop      = 2      # (always) minimum PK loop (2 nt)
        self.add_PK         = True   # (always) include PK search
        self.scan_ahead     = 10     # default 10 nt
        self.dGpk_threshold = 0.0    #  [kcal/mol]
        """threshold for PK stability"""
        
        
        # other settings
        self.Mg_binding     = False  # Mg binding interactions
        
        self.dangles        = 2      # always 2
        
        
        # For programs usind the DPA, the search for suboptimal
        # structures will range between dGmin and dGmin + dGrange.
        self.dGrange        = 0.0   # kcal/mol
        
        # ###########################################
        # #####   Chromatin spesific paramters  #####
        # ###########################################
        # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
        # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
        
        """@
        
        In general, this are not used so much even with chromatin, but
        they are here for historical reasons, backward complatibility
        and so forth.  """
        
        # constants: weights for the enthalpy of binding
        self.dHbase  = -6.0
        self.dHshift =  1.0
        
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
        
        # these are probably not necessary
        self.pssbl_ctcf          = {}
        self.edge_ctcf           = {}
        self.from_Nenski         = False
        
        # this re-scales the data by some fraction
        self.rescale_wt          = 1.0
        
        
        # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        # ###########################################
        # #####   Chromatin spesific paramters  #####
        # ###########################################
        
        
        self.allowed_extns       = []
        
        if system == "RNA":
            
            # general constants for RNA
            from RConstants import xi               # [bps]
            from RConstants import lmbd             # [bps]
            from RConstants import gmm
            # dimensionless: self avoiding random walk (app D/2)
            from RConstants import delta
            # dimensionless: related to excluded volume
            from RConstants import seg_len         # for RNA = 1
            from RConstants import minStemLen      # minimum stem length
            from RConstants import max_bp_gap      # max gap between bps
            from RConstants import minLoopLen      # minimum loop length
            from RConstants import pk_scan_ahead   # hot lead length for pk
            from RConstants import dGpk_threshold  # threshold FE for pks
            from RConstants import dGMI_threshold  # threshold FE for M-/I-loops
            from RConstants import set_dangles     # dangle parameter (always = 2)
            from RConstants import dG_range        # FE range in suboptimal structures
            
            # default FE potential
            self.set_system(system) # "RNA"
            
            
            self.xi       = xi     # [bp] stem Kuhn length
            self.lmbd     = lmbd   # bp separation distance 
            self.gmm      = gmm    # self avoiding walk constant
            self.delta    = delta  # exponental weight 
            self.seg_len  = seg_len  # [bp], for RNA = 1 nt
            self.w = self.seg_len/(self.lmbd)**2
            
            self.minStemLen     = minStemLen # [bp]
            self.max_bp_gap     = max_bp_gap # [nt]
            
            # secondary structure loop parameters
            self.minLoopLen     = minLoopLen # [nt]
            # minimum loop length (for chromatin it is 1, for RNA = 3)
            self.dGMI_threshold = dGMI_threshold # [kcal/mol]
            """threshold for MBL stability"""
            
            # pseudoknot parameters
            self.scan_ahead     = pk_scan_ahead  # default 10 nt
            self.dGpk_threshold = dGpk_threshold #  [kcal/mol]
            """threshold for PK stability"""
            
            self.dangles        = set_dangles # in general, this should be True
            self.dGrange        = dG_range    # kcal/mol
            
            # file handling matters
            self.allowed_extns  = ["vseq", "seq"]
            
        elif system == "Chromatin":
            
            # general constants for Chromatin
            from ChrConstants import xi              # Kuhn length
            from ChrConstants import lmbd            # binding distance 
            from ChrConstants import gmm             # SAW parameter (D/2)
            from ChrConstants import delta           # exp weight
            from ChrConstants import seg_len         # 5 kb
            from ChrConstants import minStemLen      # minimum stem length
            from ChrConstants import max_aa_gap      # max gap between beads anti parallel
            from ChrConstants import max_pp_gap      # max gap between beads parallel
            from ChrConstants import minLoopLen      # minimum loop length
            from ChrConstants import pk_scan_ahead   # hot lead length for pk
            from ChrConstants import dGpk_threshold  # threshold FE for pks
            from ChrConstants import dGMI_threshold  # threshold FE for M-/I-loops
            from ChrConstants import set_dangles     # dangle parameter (always = 2)
            from ChrConstants import dG_range        # FE range in suboptimal structures
            from ChrConstants import febase # [kcal/mol]
            from ChrConstants import feshift # (dimensionless, usually = 1)
            
            
            # default FE potential
            self.set_system(system) # "Chromatin"
            
            self.xi       = xi       # [bp] stem Kuhn length
            self.lmbd     = lmbd     # ratio "bond distance"/"bead-to-bead distance"  
            self.gmm      = gmm      # self avoiding walk constant
            self.delta    = delta    # exponental weight 
            self.seg_len  = seg_len  # [bp], for Chromatin = 5kbp
            self.w = self.seg_len*self.xi/(self.lmbd)**2
            """@
            
            Note that with Chromatin, xi is not really an important
            parameter for most issues and is there primarily as a
            formatlity. So we can afford to encorporate xi into the
            expression for w. For RNA, xi is a variable quantity, so
            we cannot afford to just bury it into w.
            
            """
            
            self.minStemLen     = minStemLen # [bp]
            self.max_aa_gap     = max_aa_gap # [nt]
            self.max_pp_gap     = max_pp_gap # [nt]
            
            # secondary structure loop parameters
            self.minLoopLen     = minLoopLen # [nt]
            # minimum loop length (for chromatin it is 1, for RNA = 3)
            self.dGMI_threshold = dGMI_threshold # [kcal/mol]
            """threshold for MBL stability"""
            
            # pseudoknot parameters
            self.scan_ahead     = pk_scan_ahead  # default 10 nt
            self.dGpk_threshold = dGpk_threshold #  [kcal/mol]
            """threshold for PK stability"""
            
            self.dangles        = set_dangles # in general, this should be True
            self.dGrange        = dG_range    # kcal/mol
            
            self.dHbase  = febase  # default  -6.0
            self.dHshift = feshift # default   1.0
            
            # file handling matters
            self.allowed_extns       = ["heat", "data"]
        else:
            print ("ERROR(SettingsPacket): unrecognized system (%s)" % system)
            sys.exit(1)
        #
        
        
        
        # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        #############################################################
        ####  parameters requiring settings according to system  ####
        #############################################################
        
        
        # AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        # ####################################################################
        # #####  Command line control parameters that affect FreeEnergy  #####
        # ####################################################################
        
        
        self.set_GetOpts    = False  # essentially fake!
        if self.show_info:
            print ("dinged here InputSettings")
        #
        
    #
    
    def set_system(self, sysName):
        # merely sets the system name, nothing more
        flag_pass = False
        if not sysName in sysDefLabels:
            print ("ERROR: unrecognized system name (%s)" % sysName)
            print ("       allowed names ")
            for snm in sysDefLabels.keys():
                print (snm)
            #|endfor
            
        #
        
        
        self.system = sysName
        self.molsys.set_system(sysName)
        
    #
    
    
    def set_parFlnm(self, parflnm):
        self.parFlnm = parflnm
        if not self.parFlnm == "none":
            
            ft = FileTools()
            print (dEXTS[self.system])
            if not ft.check_ext(self.parFlnm, dEXTS[self.system]):
                # RNA       -> gMtrx or par file
                # Chromatin -> heat or eheat
                
                print ("ERROR: gMatrix free energy parameter files require")
                print ("       the extension \"gMtrx\",")
                print ("       input file: %s" % args.parFlnm)
                sys.exit(1)
            #
            
            if not path.isfile(self.parFlnm):
                print ("ERROR: cannot find gMatrix file (%s)." % self.parFlnm)
                sys.exit(1)
            #
            
            flhd, ext = getHeadExt(self.parFlnm)
            print (flhd, ext)
            print ("xxx molsys.set_useParamFl", self.system, self.parFlnm)

            if ext == "par":
                self.molsys.set_useParamFl(self.system, self.parFlnm)
            elif ext == "gMtrx":
                self.molsys.set_usegMatrix(self.system, self.parFlnm)
            #
            
            if self.show_info:
                print ("input parameter file: ", self.parFlnm)
            #
            
        else:
            print ("ERROR: gMatrix file is not specified")
            sys.exit(1)
        #
        
        #print ("planned exit at InputSettings.set_parFlnm():"); sys.exit(0)
        
    #
    
    
    
    # set the basic free energy parameter set: Turner Energy Rules
    # or Vis Energy Rules or gMatrix
    
    def set_FEparamData(self, ptype = "Turner", paramflnm = "none"):
        
        if self.system == "Chromatin":
            
            self.set_ChrParams(ptype, paramflnm) # heat, or eheat file
            
        elif self.system == "RNA":
            
            self.set_RNAParams(ptype, paramflnm) # gMtrx or par file
            
        else:
            print ("ERROR: system must be defined before assigning FE parameter types")
            sys.exit(1)
        #
        
    #
                                
    def set_ChrParams(self, ptype, paramflnm = "none"): # class str
        """@
        
        With Chromatin, the parameter file is generally obtained from
        a heat or eheat file.
        
        """
        #self.set_parFlnm(paramflnm)
        #self.molsys.set_useChromatin(paramfile)
        self.molsys.set_ParamType(ptype, paramflnm)
        
        if self.molsys.paramType == "setheat":
            self.f_heatmap = self.molsys.parFlnm
        #
            
        
        
    #
    
    def set_RNAParams(self, ptype, paramflnm):
        print ("Sorry, This system is only designed to handle Chromatin")
        print ("       The RNA tools are available by agreement and from")
        print ("       the author and supporters.")
        sys.exit(0)
    #
    
    
    
    
    
    def set_source(self, s):
        self.source = s
    #
    
    def set_program(self, s):
        # job type is a lot more clear about the meaning
        self.program = s
        #print ("jobtype: ", s)
        self.molsys.set_program(s)
        #sys.exit(0)
    #
    
    def set_JobType(self, s):
        # job type is a lot more clear about the meaning
        self.jobtype = s
        #print ("jobtype: ", s)
        self.molsys.set_JobType(s)
        #sys.exit(0)
    #
    
    def set_flnm(self, inflnm):
        if self.system == "Chromatin":
            self.f_heatmap = inflnm
        elif self.system == "RNA":
            self.f_vseq = inflnm
            
        else:
            print ("ERROR: system has not been set properly")
            sys.exit(1)
        #
        
    #
    
    
    def set_sequence(self, mseq):
        if self.show_info:
            print ("set_sequence")
            print ("system ", self.system)
        #
        
        if self.system == "Chromatin":
            
            self.chrseq = mseq
            self.cSeq   = Seq(self.chrseq, self.system)
            self.N      = len(self.chrseq)
            self.molsys.set_mseq(self.chrseq)
        elif self.system == "RNA":
            
            self.rnaseq = mseq
            self.rSeq   = Seq(self.rnaseq, self.system)
            self.N      = len(self.rnaseq)
            self.molsys.set_mseq(self.rnaseq)
        else:
            print ("ERROR: system has not been set properly")
            sys.exit(1)
        #
        
    #
    
    def set_structure(self, struct):
        if self.system == "Chromatin":
            self.chrstr = struct
        elif self.system == "RNA":
            self.rnastr = struct
        else:
            print ("ERROR: system has not been set properly")
            sys.exit(1)
        #
        
        self.molsys.set_mstr(struct)
        
    #
    
    
    def set_T(self, T):
        self.T   =  T
    #
    
    def set_xi(self, xi):
        self.xi       = xi  # [bp]
        
        # as I mentioned above, I think this is WRONG!
        self.w        = self.seg_len*self.xi/(self.lmbd)**2
        
    #
    
    def set_gmm(self, gmm):
        self.gmm      = gmm   # dimensionless self avoiding walk constant
    #
    
    def set_delta(self, delta):
        self.delta    = delta # exponental weight (related to excluded volume)
    #
    
    
    def set_dGMI_threshold(self, dGMI_threshold):
        # minimum loop length (for chromatin it is 1, for RNA = 3)
        self.dGMI_threshold = dGMI_threshold # [kcal/mol]
        """threshold for MBL stability"""
    #
    
    
    def set_add_PK(self, b):
        self.add_PK         = b   # include PK search
    #
    
    def set_scan_ahead(self, i):
        self.scan_ahead     = i  # default 10
    #
    
    
    def set_dGpk_threshold(self, dGpk_threshold):
        self.dGpk_threshold = dGpk_threshold #  [kcal/mol]
        """threshold for PK stability"""
    #
    
    def set_Mg_binding(self, b):
        self.Mg_binding = b
    #
    
    def set_dGrange(self, dGrange):
        # free energy range of the search
        self.dGrange = dGrange
    #
        
    """
    ## VVVVVVVV  Hard Wired!!! VVVVVVVVVVV
    
    ## free energy 
    self.kB       = kB  # [kcal/molK]
    self.lmbd     = lmbd
    self.seg_len  = 1.0 # [nt]
    
    ## secondary structure stem parameters
    self.minStemLen     = minStemLen # [bp]
    self.max_bp_gap     = max_bp_gap # [nt]
    
    ## secondary structure loop parameters
    self.minLoopLen     = minLoopLen # [nt]
    
    ## pseudoknot parameters
    self.minPKloop      = 2 # minimum PK loop (2 nt)
    
    ## other parameters 
    self.dangles = set_dangles
    
    ## AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    """
    
    
    def show_InputSettings(self):
        print ("source:    %s" % self.source)
        print ("program:   %s" % self.program)
        print ("purpose:   %s" % self.molsys.jobtype)
        print ("paramType: %s" % self.molsys.paramType)
        print ("system:  %s" % self.system)
        
        if self.system == "Chromatin":
            print ("f_heatmap: %s" % self.f_heatmap)
            print ("chrseq:    %s" % self.molsys.mseq)
            print ("chrstr:    %s" % self.molsys.mstr)
            print ("cSeq:      %s" % self.molsys.mSeq)
        elif self.system == "RNA":
            
            print ("f_vseq:    %s" % self.f_vseq)
            print ("rnaseq:    %s" % self.molsys.mseq)
            print ("rnastr:    %s" % self.molsys.mstr)
            print ("rSeq:      %s" % self.molsys.mSeq)
        #
        
        print ("seq len:     %d" % self.N)
        
        print ("FE data set:")
        print (self.molsys)
    #
    
#    

def testRNA(ss_seq = ""):
    # some allowed RNA readable sequences
    
    #          0         10        20        30        40        50        60        70        80        90
    #          |         |         |         |         |         |         |         |         |         |
    #ss_seq = ".(((.(((((((.[[...)))))..]].((((...))))..))))).....([)]([)].......(..)..((..)).(..)..........................."
    #ss_seq  = ".((............)).."
    #ss_seq  = ".((((........)))).."
    #ss_seq  = ".((((..(...).)))).."
    #ss_seq   = ".(((.(((.(((....))).))).)))."
    if ss_seq == "":
        # default RNA test sequence
        ss_seq = ".(((.(((......(((....))).......))).)))."
    # 
    
    rnaseq   = genRNASeq(ss_seq)
    
    print ("RNA")
    iSetUpRNA = InputSettings("RNA")
    iSetUpRNA.set_source("main SettingsPacket)")
    iSetUpRNA.set_JobType("evaluation")
    iSetUpRNA.set_program("testRNA")
    # options: cantata, sarabande, sonata. 
    
    iSetUpRNA.set_FEparamData("Turner") # default is "Turner", "RNA"
    
    # set up the actual FE Data according to the settings
    iSetUpRNA.set_sequence(rnaseq)
    iSetUpRNA.set_structure(ss_seq)
    
    iSetUpRNA.show_InputSettings()
    #print ("stop here at 1"); sys.exit(0)
    gmfl = "test3s_1mm+1mmp1+2mm+3mm_3nt-w5x5.gMtrx"
    iSetUpRNA.set_FEparamData("gMatrix", gmfl) 
    iSetUpRNA.show_InputSettings()
    #print ("stop here at 2"); sys.exit(0)
    
    parfl = "ViSparams.par"
    iSetUpRNA.set_FEparamData("ParFile", parfl) # default is "Turner", "RNA"
    iSetUpRNA.show_InputSettings()
    #print ("stop here at 3"); sys.exit(0)
    
    
    iSetUpRNA.set_JobType("sarabande")
    bb_seq = '.'*len(rnaseq)
    iSetUpRNA.set_structure(bb_seq)
    iSetUpRNA.show_InputSettings()
    
    
#    


def testChr(ss_seq = ""):
    # some chromatin test sequences
    
    #          0         10        20        30        40        50        60        70        80        90
    #          |         |         |         |         |         |         |         |         |         |
    #ss_seq = ".(((.(((((((.[[...)))))..]].((((...))))..))))).....([)]([)].......(..)..((..)).(..)..........................."
    #ss_seq  = ".((............)).."
    #ss_seq  = ".((((........)))).."
    #ss_seq  = ".((((..(...).)))).."
    #ss_seq   = ".(((.(((.(((....))).))).)))."
    
    if ss_seq == "":
        # default chromatin test sequence
        ss_seq = ".(((.(((......(((....))).......))).)))."
    #
    
    chrseq   = genChrSeq(ss_seq)
    
    print ("Chromatin")
    iSetUpChr = InputSettings("Chromatin")
    iSetUpChr.set_source("main SettingsPacket)")
    
    # set up the actual FE Data according to the settings
    
    iSetUpChr.set_JobType("evaluation")
    iSetUpRNA.set_program("testChr")
    iSetUpChr.set_FEparamData("genheat") # default is "Turner", "RNA"
    iSetUpChr.set_sequence(chrseq)
    iSetUpChr.set_structure(ss_seq)
    
    
    iSetUpChr.show_InputSettings()
    #print ("stop here at 1b"); sys.exit(0)
    
    
    
    iSetUpChr.set_JobType("chreval")
    # options: chreval.
    heatflnm = "/home/dawson/python/chromatin/test.heat"
    iSetUpChr.set_FEparamData("setheat", heatflnm) # default is "Turner", "RNA"
    
    print ("vv")
    # set up the actual FE Data according to the settings
    iSetUpChr.set_sequence(chrseq)
    bb_seq = '.'*len(chrseq)
    iSetUpChr.set_structure(bb_seq)
    
    iSetUpChr.show_InputSettings()
    #print ("stop here at 2b"); sys.exit(0)
#    

def usage_main():
    print ("%s system [sequence]")
    

def main(cl):
    print (cl)
    if len(cl) < 2:
        print ("ERROR: must specify the system and optional sequence")
        sys.exit(1)
    #
    
    if not cl[1] in sysDefLabels:
        print ("ERROR: first argument must be either RNA or Chromatin")
        sys.exit(0)
    #
    
    struct = ""
    if len(cl) == 3:
        struct = cl[2]
    #
    
    if cl[1] == "RNA":
        
        testRNA(struct)
        
    elif cl[1] == "Chromatin":
        
        testChr(struct)
    #
    
#


if __name__ == '__main__':
    # running the program
    main(sys.argv)
#    

