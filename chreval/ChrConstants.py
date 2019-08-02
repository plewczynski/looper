#!/usr/bin/env python

"""@@@

program:        ChrConstants.py

Creation Date:  cr 2017.03~ (in various levels of development)
Last Update:    190719
Version:        1.0


Purpose:

Contains constants that are often used specifically in various
Chromatin programs. Having one file that contains them all helps to
keep everything regular.

Comments:

190522: These constants have been moved around in various stages from
various modules. Now I seek to have these constants all specifically
associated with features of chromatin.

"""

# #################################################################
# ###############  General configuration CONSTANTS  ###############
# ###############    settings used in FreeEnergy    ###############
# #################################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv


T37C    = 310.15   # [K] 37C in Kelven

# fundamental coarse grained resolution length
seg_len = 5000.0

# constants: in entropy evaluation
xi       = 5.0    # # [bp] default stem Kuhn length
xi_fs    = 5.0    # [nt] free strand Kuhn length (basically useless)
lmbd     = 0.002  # [bps] binding distance := bps / seg_len 
gmm      = 2.3    # dimensionless but related to D/2
delta    = 2.0    # dimensionless exponetial scaling weight

# constants: weights for the enthalpy terms 
febase   = -6.0
# [kcal/mol], weight parameter, also have used 4.7 kcal/mol 
feshift  =  1.0 # (dimensionless, usually = 1)

"""@@@

Presently, the enthalpy is based on a logarithmic assessment with
respect to the number of observed counts of a particular interaction
between parts (i,j) of the chromatin chain. It is clearly a crude
function to assess this binding free energy, but I have no real
information to work with here. For a given set of data, these two
parameters may surely require tuning. The current test set showed
somewhat reasonable results when setting a baseline of 4 kcal/mol for
the binding free energy of the CTCF clusters. Whether this estimate is
remotely correct or not is currently unknown.

"""

# other stuff relatd to fitting that should have defaults recorded in
# one place

minStemLen  = 1 # [bp], minimum stem length (default is 1)
"""@ 

The minimum stem length is an artifact of vsfoldN. I don't think there
should be any reason to change this anymore. stiffness should
ultimately be decided by the best stem length and not stupid blanket
parameters like this on. I keep it here because the minimum stem
length for chromatin is 1, where as for RNA, the minimum pair is at
least a dinucleotide base pair.

"""
max_bp_gap = 1



minLoopLen  = 1 # [nt], minimum loop length (default is 1)
# minimum loop length (for chromatin it is 1, or RNA 3)

# threshold for MBL branch stability
dGMI_threshold = -0.3  # [kcal/mol] threshold FE for M-/I-loops


# pseudoknot parameters
pk_scan_ahead  = 10   # [nt] the "hot lead" length for pk
dGpk_threshold = -0.3 # [kcal/mol] threshold FE for pks
"""@

about the thresholds

Something less than this is basically the same as
the thermal energy, so it is probably rather weak. Anyway,
it seems like -0.3 kcal/mol for 'B' is really weak and it
seems unlikely that the stability of such a structure is
really significant. At any rate, the program will work even
if you set the threshold to zero.
"""        
        
set_dangles = 2 # dangle parameter (always = 2)

# Actually this has always been 2 with vsfoldN and the option
# should be basically ignored, in my opinion.

dG_range = 10.0 # [kcal/mol], FE range in suboptimal structures

# Maximum free energy difference range of interest in suboptimal
# structures

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# #################################################################

