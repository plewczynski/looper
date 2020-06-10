#!/usr/bin/env python3

"""@@@

Main Module:   GetOpts.py 

Classes:       GetOpts

Author:        Wayne Dawson
creation date: mostly 2016, some developments in 2017.
last update:   200311 minor adjustments for analyze_loops, etc.
version:       0

Purpose:

Reads the command line information for chreval and sets up the
calculations.


Comments:

I try to do things a little bit differently, and introduce a kind of
command line reader for related programs in this suite. I don't know
that this will improve anything particularly, but it is worth a try.

##################################################################
170829(Note): Surprisingly, although this program is called as an
object in main of the initiating program, I don't have to specify
sys.argv within this module for argparse to function properly.
##################################################################

"""


import sys
from FileTools import FileTools
import argparse

# ################################################################
# ######################  Global constants  ######################
# ################################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

from Constants import kB # [kcal/molK] (Boltzmann constant)

from ChrConstants import T37C    # [K] 37C in Kelven
# coarse-grained resolution factor
from ChrConstants import seg_len # bead to bead distance
# for entropy calculations
from ChrConstants import xi      # [bps]
from ChrConstants import lmbd    # [bps]
from ChrConstants import gmm     # dimensionless but related to D/2        
from ChrConstants import delta   # dimensionless exponetial scaling weight
# constants: weights for the enthalpy terms 
from ChrConstants import febase  # [kcal/mol]
from ChrConstants import feshift # (dimensionless, usually = 1)

# locked vvvvvvvvvvvvvvvvvvvvvvvvvvvvv
from ChrConstants import minStemLen      # minimum stem length
from ChrConstants import minLoopLen      # minimum loop length
# locked ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
from ChrConstants import pk_scan_ahead   # hot lead length for pk
from ChrConstants import dGpk_threshold  # threshold FE for pks
from ChrConstants import dGMI_threshold  # threshold FE for M-/I-loops
# locked vvvvvvvvvvvvvvvvvvvvvvvvvvvvv
from ChrConstants import set_dangles     # dangle parameter (always = 2)
# locked ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
from ChrConstants import dG_range        # FE range in suboptimal structures

from CVersion import get_now
from CVersion import get_Version
from CVersion import get_lastUpdate

# program processing routines
from SettingsPacket import InputSettings

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ################################################################
# ################################################################


# parameters for measuring similarity between two structures in terms
# of free energy ddG or T*ddS.
set_TdS_range  = 2.0 # [kcal/mol] similarity in delta(TdS)
set_ddG_range  = 2.0 # [kcal/mol] similarity in delta(dG)

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# #################################################################



class GetOpts(InputSettings):
    """@
    
    Reads the command line arguments for the suit of programs
    associated with chromatin structure prediction and analysis
    
    """
    
    # NOTE: __everything__ is read or assumed to be a STRING here
    # whether it is a number or a string must be decided elsewhere.
    
    def __init__(self, program):

        InputSettings.__init__(self, "Chromatin")   # inherit InputSettings
        
        self.source = "GetOpts"
        self.program = program
        self.system  = "Chromatin"
        
        # we're able to update the information to this extent at least
        self.molsys.set_JobType("prediction")
        
        
        # reading command line using argparse.parser
        self.debug_GetOpts = False
        if self.debug_GetOpts:
            print ("program: ", program)
        #
        
        
        self.f_heatmap   = ''
        self.f_activity  = ''
        self.f_output    = ''
        
        """@
        
        basic comparision of stability
        
        Measure stability of a structure based on the magnitude of the
        first Boltzmann probability term. Presumably, the terms that
        have a very high Boltzmann probability also have the highest
        stability. However, it should be remembered that as the size
        of the loop becomes large, it is more likely that you will
        have multiple structures that are very similar in structure
        and free energy. As a result, I started searching for ways to
        combine similar structures together using various rules such
        as hamming distance, word similarity, variation in delta(TdS)
        and variation in delta(dG).
        
        """
        
        self.basic       = False # the first Boltzmann prob
        
        # hamming distance 
        
        # measure similarity between two structures in terms of
        # hamming distance 
        self.hamming     = False # prob wtd by hamming distance
        
        # similarity distance 
        
        # measure similarity between two structures in terms of
        # word similarity
        self.similarity  = False # prob wtd by similarity
        
        # measure similarity between two structures in terms of
        # difference in entropy
        self.TdS         = False # prob wtd by CLE 
        self.TdS_range   = set_TdS_range # [kcal/mol] similarity in delta(TdS)
        
        # measure similarity between two structures in terms of
        # difference in free energy
        self.ddG         = False # prob wtd by d(dG)
        self.ddG_range   = set_ddG_range # [kcal/mol] similarity in delta(dG)
        
        self.allwts      = True # show all weights
        self.from_Nenski = False
        self.p_all_1D    = False # only print max 50 1D files
        
        
        self.use_eheat   = False # read in *.eheat files instead of *.heat
        
        
        self.parser = None
        args = None
        self.parser = self.setup_parser(program)
        if program == "GetOpts.py":
            
            """@@@
            
            This could have been done in setup_parser, and it worked
            when I did it that way. However, it is generally better to
            define stuff in a way where you don't have the program
            just ignore an undefined variable because it is buried in
            an "if clause".  Therefore, I set it up this way to keep
            everything clearly on the table.
            
            """
            
            self.parser.add_argument('-test', action='store', default=program,
                                     dest='testprog', 
                                     help='For testing %s with different programs.' % program)
            #
            
            args = self.parser.parse_args()
            program = args.testprog
        else:
            args = self.parser.parse_args()
        #
        
        # ############################################################
        # now evaluate the input data and assign appropriate variables
        # ############################################################
        
        # variables ok?
        self.test_cmd_args(self.parser)
        
        # extensions ok?
        self.EXTS = {}
        self.EXTS.update({'chreval.py' : ["heat", "eheat"]}) # allowed for chreval.py
        self.EXTS.update({'analyze_loops.py'    : ["bed", "txt"] })
        # both "bed" and "txt" extensions are allowed for analyze_loops.py
        self.EXTS.update({'assemble_heatmaps_and_CCDs.py'    : ["bed"] })
        
        # This is still in planning. Anyway, the "txt" files lack
        # sufficient information on the orientation of the loops, so I
        # only allow the bed files. In fact, I am not convinced that
        # this is used with GetOpts.
        
        self.EXTS.update({'GetOpts.py' : ["heat", "data", "bed", "txt"] }) # allowed
        
        # all programs presently "data" is not used for chreval but
        # both "bed" (Przemek format) and "txt" (Teresa)
        
        self.check_list(program, self.EXTS)
        self.program       = program
        self.allowed_extns = self.EXTS[program]
        
        if self.debug_GetOpts:
            print (self.allowed_extns)
            print (self.EXTS)
        #
        
        # information on input files
        self.flnm          = [] # full file name
        self.flhd          = [] # header for file
        self.ext           = [] # file extension
        
        
        if self.program == 'chreval.py':
            
            self.f_heatmap = args.f_heatmap
            #print (self.f_heatmap)
            self.check_extensions(self.f_heatmap)
            
        elif self.program == 'analyze_loops.py':
            self.f_activity = args.f_activity
            self.check_extensions(self.f_activity)
            self.f_output   = args.f_output
        elif self.program == 'assemble_heatmaps_and_CCDs.py':
            self.f_activity = args.f_activity
            self.f_output   = args.f_output
            self.check_extensions(self.f_activity)
        elif self.program == 'GetOpts.py':
            self.f_activity = args.f_activity
            self.f_heatmap  = args.f_heatmap
            self.f_output   = args.f_output
            print ("file listings")
            print ("anal:     ", self.f_activity)
            print ("chreval:  ", self.f_heatmap)
            print ("out flnm: ", self.f_output)
            self.check_extensions(self.f_activity)
            self.check_extensions(self.f_heatmap)
        else:
            emsg = "ERROR: undefined program"
            self.error(emsg)
        #
        
        # this sets [flnm], [flhd], and [ext]
        
        # output option
        self.p_all_1D      = args.p_all_1D
        """@@@
        
        In general, should only need to see the first 10 structures in
        the created directory. For a large calculation, it is easy to
        produce thousands of structures. For quite large structure,
        the program can easily generate one hundred thousand
        structures. Therefore, to avoid a explosion of files, the
        default is to print a maximum of 100 structures. This is
        primarily so that the program (or the computer) doesn't crash
        due to the creation of so many files.
        
        """
        
        # !*1 = currently not an adjustable parameter
        
        # Entropy parameters
        self.T             = args.T       # in Kelvin
        self.seg_len       = args.seg_len # [bps] number of bps between beads
        self.xi            = args.xi      # [nt] Kuhn length
        self.lmbd          = lmbd         # [nt / seg_len] 
        self.gmm           = args.gmm     # self avoiding walk parameter
        self.delta         = args.delta   # exponential weight factor
        
        
        self.minStemLen     = minStemLen
        """@
        
        Note minStemLen "minimum stem length" is an artifact of
        vsfoldN. I don't think there should be any reason to change
        this anymore. stiffness should ultimately be decided by the
        best stem length and not stupid blanket parameters like this
        on. I keep it here because the minimum stem length for
        chromatin is 1, where as for RNA, the minimum pair is at least
        a dinucleotide base pair.
        
        """
        
        self.minLoopLen     = minLoopLen
        # minimum loop length (for chromatin it is 1, or RNA 3)
        
        
        # pseudoknot parameters
        self.scan_ahead     = args.leadingEdge  # default 10
        self.dGpk_threshold = args.pk_threshold # default 2 (needs better definition!!!)
        
        
        self.dGMI_threshold = args.mi_threshold # [kcal/mol]
        # threshold for MBL stability
        
        """@ 
        
        190130: Something less than this is basically the same as the
        thermal energy, so it is probably rather weak. Anyway, it
        seems like -0.2 kcal/mol for 'B' is really weak and it seems
        unlikely that the stability of such a structure is really
        significant. At any rate, the program will work even if you
        set the threshold to zero.
        
        """
        
        self.dangles        = set_dangles # Free energy range to be listed
        """@
        
        Actually this has always been 2 with vsfoldN and the option
        should be basically ignored, in my opinion.
        
        """
        
        # Enthalpy constants: weights for the enthalpy terms 
        self.dHbase        = args.dHbase
        self.dHshift       = feshift # [dimensionless] == 1
        
        self.dGrange       = args.dGrange
        # options associated with all programs
        self.add_PET_wt    = args.add_PET_wt
        self.add_PET_wt_to_edges = args.add_PET_wt_to_edges
        if args.add_PET_wt_to_edges:
            self.add_PET_wt = True
        #
        
        self.PETwt         = 100.0  # [relative counts]
        if isinstance(args.PET_wt, float):
            # used "isinstance" because PET_wt is defined as "None" in
            # option input
            
            # print ("setting add_PET_wt")
            self.add_PET_wt = True
            self.PETwt      = args.PET_wt
        #
        
        self.CTCF_scale    = 100.0 
        if isinstance(args.CTCF_scale, float):
            self.CTCF_scale = args.CTCF_scale
        #
        
        """@@@
        
        161029wkd: In general, this should not be changed, but it does
        allow that you can change the range of tandem CTCF and
        convergent CTCF by bumping this parameter up or
        down. Conversely, If (for example) you choose 50, then the
        tandem CTCFs will now appear at 10, the convergent CTCFs at
        20. Likewise, if you chose 200, then the tandem CTCFs will
        appear at 40 and the convergent at 80.
        
        The original use of the above (with PET_range and ref_scale
        (now CTCF_scale)) was a very confusing. One month later, I
        could not remember how it worked other than the maximum value
        always came out to 100 afterwards. Therefore, these are not
        used and instead, I have introduce a weight rescale_wt
        instead. Therefore, now I just make it automatic that entries
        in the input heatmap are reweighted by rescale_wt and recorded
        effectively as integers.
        
        This is part of on-going corrections to the program to adjust
        for different data types.
        
        """
        
        self.rescale_wt    = 1.0 
        if isinstance(args.rescale_wt, float):
            self.rescale_wt  = args.rescale_wt
        #
        
        # Nenski data
        self.from_Nenski   = args.from_Nenski
        
        
        if self.program == "analyze_loops.py" or \
           self.program == "assemble_heatmaps_and_CCDs.py" or \
           self.program == "GetOpts.py":
            # options associated with analyze_loops or
            # assemble_heatmaps_and_CCDs
            
            self.use_eheat     = args.use_eheat # use the *.eheat file
            
            self.basic         = args.basic 
            self.hamming       = args.hamming
            self.similarity    = args.similarity
            self.TdS           = args.TdS
            self.TdS_range     = args.TdS_range
            self.ddG           = args.ddG
            self.ddG_range     = args.ddG_range
            self.allwts        = args.allwts
            
            if self.allwts:
                self.basic      = True  
                self.hamming    = True 
                self.similarity = True 
                self.TdS        = True 
                self.ddG        = True
            elif self.similarity:
                self.basic      = False  
                self.hamming    = False 
                self.similarity = True 
                self.TdS        = False 
                self.ddG        = False
            elif self.hamming:
                self.basic      = False  
                self.hamming    = True 
                self.similarity = False 
                self.TdS        = False 
                self.ddG        = False
            elif self.TdS:
                self.basic      = False  
                self.hamming    = False 
                self.similarity = False 
                self.TdS        = True 
                self.ddG        = False
            elif self.ddG:
                self.basic      = False  
                self.hamming    = False 
                self.similarity = False 
                self.TdS        = False 
                self.ddG        = True
            else:
                self.basic      = True  
                self.hamming    = False 
                self.similarity = False 
                self.TdS        = False 
                self.ddG        = False
            #
            
        #
        
        self.set_GetOpts   = True
    #
    
    
    
    
    
    def check_list(self, program, exts):
        """
        verifies whether program is within the list of keys
        """
        
        flag_has_key = False
        keys = exts.keys()
        for key_k in keys:
            if key_k == program:
                flag_has_key = True
                break
            #
            
        #|endfor
        
        if not flag_has_key:
            emsg = "ERROR: no program of name '%s' on record." % program
            self.error(emsg)
        #
        
        return flag_has_key
    #
    
    
    
    def setup_parser(self, program):
        """
        contains all the command line arguments
        """
        parser = argparse.ArgumentParser()
        
        
        parser.add_argument('-T', action='store', default=310.0,
                            dest='T', type=float,
                            help='Set temperature in Kelvin [K].')
        
        s_xi = 'The Kuhn length. For RNA, this is a variable (default value \
        is %6.2f nt).' % xi
        parser.add_argument('-xi',  action='store', default=xi,
                            dest='xi', type=float,
                            help=s_xi)
        
        s_gamma = 'Self avoiding walk parameter (default value for RNA is \
        %5.2f -- sometimes it is increased to about 2.3).' % gmm
        parser.add_argument('-gamma',  action='store', default=gmm,
                            dest='gmm', type=float,
                            help=s_gamma)
        
        s_delta = 'Self avoiding walk parameter (default value for RNA is \
        %5.2f -- it might vary between 1.0 and 2.5).' % delta
        parser.add_argument('-delta',  action='store', default=delta,
                            dest='delta', type=float,
                            help=s_delta)
        
        parser.add_argument('-seg_len', action='store', default=seg_len,
                            dest='seg_len', type=float,
                            help='Set number of base pairs between beads.')
        
        parser.add_argument('-pkLead',  action='store', default=pk_scan_ahead,
                            dest='leadingEdge', type=int,
                            help='For PK search, how much of a \'hot lead\' you \
                            want to use ahead of the secondary structure search \
                            (default is 10, probably not good -- or at least \
                            expensive -- to request more than 20).')
        
        parser.add_argument('-pkThresh',  action='store', default=dGpk_threshold,
                            dest='pk_threshold', type=float,
                            help='Pseudoknot threshold free energy; the smallest \
                            free energy necessary for formation of a stable \
                            pseudoknot).')
        
        parser.add_argument('-miThresh',  action='store', default=dGMI_threshold,
                            dest='mi_threshold', type=float,
                            help='M-,I-loop  threshold free energy; the smallest \
                            free energy necessary for formation of a stable \
                            M-loop or I-loop).')
        
        # ctcf parameters
        parser.add_argument('-add_PET_wt', action='store_true', default=False,
                            dest='add_PET_wt',
                            help='Turn on default PET cluster weight at edges of the heatmap.')
        
        parser.add_argument('-add_PET_wt_to_edges', action='store_true', default=False,
                            dest='add_PET_wt_to_edges',
                            help='Place PET cluster weights at very edges of the heatmap \
                            (strongly not recommended anymore, it was used for early \
                            heat maps that we were analyzing).')
        
        parser.add_argument('-PET_wt',     action='store', default=None,
                            dest='PET_wt', type=float,
                            help='Set PET cluster weight at the edges of the heatmap (used \
                            in conjuction with the option -add_PET_wt_to_edges).')
        
        parser.add_argument('-rescale_wt',  action='store', default=1.0,
                            dest='rescale_wt', type=float,
                            help='Rescales the input data by some factor > 0 (default 1).')
        
        parser.add_argument('-CTCF_assignment_scale',  action='store', default=100.0,
                            dest='CTCF_scale', type=float,
                            help='Set the assignment range for CTCF sites (default 100), \
                            where the range between 20 and 40 would represent tandem \
                            structures and more than 40 convergent CTCF sites.')
        
        s_dHbase = 'Adjusts the baseline of the enthalpy \
        (default %6.3f [kcal/mol]).' % febase
        parser.add_argument('-dHbase',  action='store', default=febase,
                            dest='dHbase', type=float,
                            help=s_dHbase)
        
        
        parser.add_argument('-dGrange',     action='store', default=dG_range,
                            dest='dGrange', type=float,
                            help='Set the free energy range that is searched relative to \
                            the minimum free energy; e.g., dGrange = [10 kcal/mol] and \
                            minimum FE -30 [kcal/mol] means search for structures with FE \
                            between -20 and -30 [kcal/mol.')
        
        
        parser.add_argument('-HiC', action='store_true', default=False,
                            dest='from_Nenski',
                            help='When using Hi-C data, this has a filtering algorithm \
                            that helps to remove self-ligation NN chain noise.')
        
        parser.add_argument('-printAll1D', action='store_true', default=False,
                            dest='p_all_1D',
                            help='WARNING: this option may demand a lot more disk\
                            space!!!! Creates a VARNA readable file for every \
                            1D structure evaluated. (With the default settings, the \
                            program creates a max of 50 files of 1D structures). \
                            Rarely do we look beyond the first few structures, so \
                            why would we want to look at 10 of thousands of them?')
        
        flag_checkfile = False
        
        
        
        if program == "chreval.py" or program == "GetOpts.py":
            grp_files = parser.add_mutually_exclusive_group()        
            # options unique to the chreval.py program
            if flag_checkfile:
                # requires existence of file of given name
                
                # 190515: as I recall, giving it type=file from the
                # start created some problems and, as a result, I
                # found this way of reading in a file rather awkward
                # to use. Perhaps that cryptic remark is because when
                # the file could not be found, the program crashed. I
                # think this part should be deleted, but I seem
                # reluctant for some reason I cannot explain.
                grp_files.add_argument('-f', nargs=1, type=file,
                                       default=[], dest='f_heatmap',
                                       help='Input heatmap data on frequency \
                                       of contact interactions observed in \
                                       chromatin structure (input file \
                                       requires extension \'heat\' or \'data\').')
            else:
                grp_files.add_argument('-f', nargs=1, default=[],
                                       dest='f_heatmap',
                                       help='Input heatmap data on frequency \
                                       of contact interactions observed in \
                                       chromatin structure (input file \
                                       requires extension \'heat\' or \'data\').')
            #
            
        #
        
        
        if program == "analyze_loops.py" or \
           program == "assemble_heatmaps_and_CCDs.py" or \
           program == "GetOpts.py":
            # options unique to the analyze_loops.py or
            # assemble_heatmaps_and_CCDs.py program
            
            # mutually exclusive options used in anal for weighting
            # the Boltzmann probabilities
            grp_filter = parser.add_mutually_exclusive_group()
            
            # basic
            grp_filter.add_argument('-basic', action='store_true', default=False,
                                    dest='basic',
                                    help='(analyze_loops) Only report the probability of the first \
                                    component in the Boltzmann distribution.')
            
            # similarity
            grp_filter.add_argument('-sim', action='store_true', default=False,
                                    dest='similarity',
                                    help='(analyze_loops) Use similarity measure in weighting \
                                    chromatin structures.')
            
            # hamming distance
            grp_filter.add_argument('-ham', action='store_true', default=False,
                                    dest='hamming',
                                    help='(analyze_loops) Use hamming distance in weighting \
                                    chromatin structures.')
            
            # entropy weight
            grp_filter.add_argument('-TdS', action='store_true', default=False,
                                    dest='TdS',
                                    help='Use entropy difference in weighting chromatin structures.')
            
            grp_filter.add_argument('-TdS_range',     action='store', default=set_TdS_range,
                                    dest='TdS_range', type=float,
                                    help='(analyze_loops) Set the range for dTdS (used with \
                                    option -TdS).')
            
            # free energy
            grp_filter.add_argument('-ddG',  action='store_true', default=False,
                                    dest='ddG',
                                    help='(analyze_loops) Use free energy difference in \
                                    weighting chromatin structures.')
            
            # read in *.eheat files instead of *.heat
            grp_filter.add_argument('-eheat',  action='store_true', default=False,
                                    dest='use_eheat',
                                    help='read eheat files instead of heatmap.')
            
            grp_filter.add_argument('-ddG_range',     action='store', default=set_ddG_range,
                                    dest='ddG_range', type=float,
                                    help='(analyze_loops) Set the range for ddG (used in \
                                    conjuction with option -ddG).')
            
            # provide all measures of difference
            grp_filter.add_argument('-all', action='store_true', default=True,
                                    dest='allwts',
                                    help='(analyze_loops, default) make a list with all \
                                    possibilities included (recommended).')
            
            
            grp_files = parser.add_mutually_exclusive_group()        
            # options unique to the analyze_loops.py program
            if flag_checkfile:
                # requires existence of file of given name
                
                # 190515: as I recall, giving it type=file from the
                # start created some problems and, as a result, I
                # found this way of reading in a file rather awkward
                # to use. Perhaps that cryptic remark is because when
                # the file could not be found, the program crashed. I
                # think this part should be deleted, but I seem
                # reluctant for some reason I cannot explain.
                
                grp_files.add_argument('-ff', nargs='+', type=file, default=[],
                                       dest='f_activity',
                                       help='(analysis programs) Input activity data on \
                                       chromatin structure (requires extension \
                                       \'txt\' or \'bed\')')
            else:
                grp_files.add_argument('-ff', nargs='+', default=[],
                                       dest='f_activity',
                                       help='(analysis programs) Input activity data on \
                                       chromatin structure (requires extension \
                                       \'txt\' or \'bed\')')
            #
            
            parser.add_argument('-o', default="loops_results.dat",
                                dest='f_output',
                                help='(analysis programs) Specifies the name of \
                                the output file for the analysis program.')
            
        #
        
        
        """
        ############################################################
        
        other types of read in formats
        
        parser.add_argument('--three', nargs=3) # three args
        
        parser.add_argument('--all', nargs='*', dest='all') # all
        
        parser.add_argument('--optional', nargs='?') # option
        
        ############################################################
        
        """
        
        parser.add_argument('-v', '--version', action='version',
                            version='%(prog)s {version}'.format(version=get_Version()))
        
        
        return parser # the main part of the parser
        
    #
    
    def test_cmd_args(self, parser):
        flag_pass = True
        # check to make sure that the specified command line entries
        # and respective variables all agree in type and conditions
        # and whether there are undefined terms.
        try:
            self.parser.parse_args()
        except IOError:
            msg = 'could not parse the arguments'
            flag_pass = False # not necessary, but anyway....
            self.parser.error(str(msg))
        #
        
        return flag_pass
    #
        
    
    def check_extensions(self, files):
        if self.debug_GetOpts:
            print (files)
        #
        
        # check for incompatible mismatch of file types.  This is not
        # so easy to write corrections to argparse to handle.
        flag_file_mismatch = False
        if self.program == 'analyze_loops.py' and len(self.f_heatmap) > 0:
            flag_file_mismatch = True
            
        elif (self.program == 'assemble_heatmaps_and_CCDs.py' and 
              len(self.f_heatmap) > 0):
            flag_file_mismatch = True
            
        elif self.program == 'chreval.py' and len(self.f_activity) > 0:
            flag_file_mismatch = True
        #
        
        if flag_file_mismatch:
            emsg =  "ERROR: Incompatible file type, '%s' only accepts \n" \
                    % self.program
            emsg += "       files with option -ff and extension(s) %s\n" \
                    % self.allowed_extns
            self.error(emsg)
        #
        
        if len(files) == 0:
            emsg = "ERROR: no input file specified." 
            self.error(emsg)
        #
        
        if self.program == "chreval.py" and len(files) > 1:
            emsg = "ERROR: %s only accepts one file at a time." % self.program
            self.error(emsg)
        #
        
        ft = FileTools()
        for flnm in files:
            
            if not ft.check_ext(flnm, self.allowed_extns):
                emsg = "       terminating....." 
                self.error(emsg)
            #
            
            self.flnm += [flnm]
            self.flhd += [ft.flhd]
            self.ext  += [ft.ext]
        #|endfor
        
        if self.debug_GetOpts:
            print ("input file information")
            print (self.flnm)
            print (self.flhd)
            print (self.ext)
        #
        
    #

                
    def msg(self):
        self.parser.print_help() # prints help
        sys.exit(0)
    #
                
    def error(self, emsg):
        print (emsg)
        self.parser.print_help() # prints help
        sys.exit(1)
    #
        
# ###################
# ####   MAIN    ####
# ###################



def main():
    PROGRAM = "GetOpts.py"
    cl = GetOpts(PROGRAM)
#    

if __name__ == '__main__':
    main()
#
