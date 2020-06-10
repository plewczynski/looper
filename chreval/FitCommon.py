#!/usr/bin/env python3

"""@@@

program:       FitCommon.py

Classes:       

Functions:     is_symmetric 
               has_terml_AU

creation date: 181031
last update:   200414 refactoring
version:       0.0

FitCommon contains a variety of odds and ends (constants, dictionary
definitions, function tools, etc) that are used by a variety of
unrelated programs.

Unforutnately, with programming, sometimes the "most logical" place to
put things is not the "most practical". Suffice it to say that there
is enough in common between various programs that it seemed "more
practical" to put these items here.

"""




"""@@@

BPWeights

Contains some constants and functions that are shared by various
programs. This helps avoid getting stuck in a circular definition of
various functions and constants in some of the most commonly used
modules.

"""
dbp_set = {"WC" : True, "GU" : True, "nWC" : True} # this is from command line

"""@

interpretive dictionaries:

d-        == "dictionary"

mono-nucleotide base pairs
WC        == "Watson-Crick"
cWC       == "canonical Watson-Crick"
nWC       == "non-Watson-Crick"
guWC      == "GU or UG base pair"
WC+GU     == "canonical Watson-Crick with WC-like GU"
WCwgu     == "WC+GU"

mm-       == "mismatch", e.g., "1mm" = one mismatch

di-nucleotide base pairs
WCwgu*nWC == "WC+GU at bp1 and nWC at bp2" --  1mm (when GU is considered WC)
1mm       == "(WCwgu*nWC) or (nWC*WCwgu)" --  WC+GU  and nWC
2mm       == "nWC*nWC" --  all bp are nWC


tri-nucleotide base pairs
1mm       == "(WCwgu*WCwgu*nWC) or (WCwgu*nWC*WCwgu) or (nWC*WCwgu*WCwgu)" -- two WC+GU  and one nWC
----
1bp       == "nWC*(WCwgu)*nWC" -- two nWC bps separated by a WCwgu bp.
2mm       == "(WCwgu*nWC*nWC) or (nWC*nWC*WCwgu)" -- two consequative nWC bps and one WCwgu 
              **note that 1bp + 2mm is the total set of 2 mismatches**
----
3mm       == "nWC*nWC*nWC"; i.e., all bp are nWC and not GU

"""
dfitlevel =  { "iloop" : 0,
               "WC"    : 1, 
               "WC+GU" : 2,
               "1mm"   : 3,
               "2mm"   : 4,   # not used with di-nt, only tri-nt
               "1bp"   : 5,   # only tri-nt
               "3mm"   : 6 }  # only tri-nt

"""@

NOTE 1: iloop and WC have the same tables but not the same base line
constants

NOTE 2: 1mm means that one of the base pairs is WC+GU and the other is
nWC.

NOTE 3: if iloop is called, then WC is not used, and visa versa.

"""

dflag_maxlevel_list = { "3mm"   : False, # only tri-nt
                        "1bp"   : False, # only tri-nt
                        "2mm"   : False, # not used with di-nt, only tri-nt
                        "1mm"   : False,
                        "WC+GU" : False,
                        "WC"    : False,
                        "iloop" : False }


"""@

It seemed logical to put dplxckey and dplxgen here too; however, they
dplxckey requires set values that are best kept within the confines of
ViennaDataFmt. Therefore, both these dictionary labels were moved
there. In general, they should not be required for anything outside of
the Vienna package duplex parameters and the Turner Energy Rules that
they were based on.


# keys assigned as numerical values
dplxckey = { "init" : 0,
             "sym"  : 1,
             "tAU"  : 2 }
# values in dictionary format
dplxgen  = { "init" : stemInit,
             "sym"  : stemsym,
             "tAU"  : terminalAU }

"""


# ViennaDataFmts
base_ndx = {'@' : 0,
            '_' : 0, 
            'A' : 1,
            'C' : 2,
            'G' : 3,
            'U' : 4,
            'T' : 4,
            'X' : 5,
            'K' : 6,
            'I' : 7,
            'Y' : 8}

basepair = {'@@' : 0,
            '__' : 0,
            'CG' : 1,
            'GC' : 2,
            'GU' : 3,
            'UG' : 4,
            'AU' : 5,
            'UA' : 6,
            'XY' : 7,
            'YX' : 7 } 

dWCbps    = { "CG" :   "WC",  
              "GC" :   "WC",
              "AU" :   "WC",
              "UA" :   "WC"
} # canonical Watson-Crick base pairs

dGUbps    = { "GU" :  "guWC",
              "UG" :  "guWC"
} # WC-like GU base pairs

dWCwgubps = { "CG" :  "WC+GU",  
              "GC" :  "WC+GU",
              "GU" :  "WC+GU",
              "UG" :  "WC+GU",
              "AU" :  "WC+GU",
              "UA" :  "WC+GU"
} # WC and WC-like GU base pairs

dstacking = { "__" :  "empty",
              "CG" :  "WC+GU",  
              "GC" :  "WC+GU",
              "GU" :  "WC+GU",
              "UG" :  "WC+GU",
              "AU" :  "WC+GU",
              "UA" :  "WC+GU", 
              "XY" :  "WGnst",
              "YX" :  "WGnst"
}  # for stacking interaction


dnWCbps   = { "AA" :   "nWC",
              "AC" :   "nWC",  
              "AG" :   "nWC",
              # <- AU
              "CA" :   "nWC",
              "CC" :   "nWC",
              # <- CG
              "CU" :   "nWC",
              "GA" :   "nWC",
              # <- GC
              "GG" :   "nWC",
              # <- GU
              # <- UA
              "UC" :   "nWC",
              # <- UG
              "UU" :   "nWC"
} # "nWC" = non-Watson-Crick base pairs



# "forward" dangle
dmismatchf = {
    # __
    "__/__" : "mmf", "__/_A" : "mmf", "__/_C" : "mmf", "__/_G" : "mmf", "__/_U" : "mmf",
    "_A/__" : "mmf", "_A/_A" : "mmf", "_A/_C" : "mmf", "_A/_G" : "mmf", "_A/_U" : "mmf",
    "_C/__" : "mmf", "_C/_A" : "mmf", "_C/_C" : "mmf", "_C/_G" : "mmf", "_C/_U" : "mmf",
    "_G/__" : "mmf", "_G/_A" : "mmf", "_G/_C" : "mmf", "_G/_G" : "mmf", "_G/_U" : "mmf",
    "_U/__" : "mmf", "_U/_A" : "mmf", "_U/_C" : "mmf", "_U/_G" : "mmf", "_U/_U" : "mmf",
    
    # CG
    "C_/G_" : "mmf", "C_/GA" : "mmf", "C_/GC" : "mmf", "C_/GG" : "mmf", "C_/GU" : "mmf",
    "CA/G_" : "mmf", "CA/GA" : "mmf", "CA/GC" : "mmf", "CA/GG" : "mmf", "CA/GU" : "mmf",
    "CC/G_" : "mmf", "CC/GA" : "mmf", "CC/GC" : "mmf", "CC/GG" : "mmf", "CC/GU" : "mmf",
    "CG/G_" : "mmf", "CG/GA" : "mmf", "CG/GC" : "mmf", "CG/GG" : "mmf", "CG/GU" : "mmf",
    "CU/G_" : "mmf", "CU/GA" : "mmf", "CU/GC" : "mmf", "CU/GG" : "mmf", "CU/GU" : "mmf",
    
    # GC
    "G_/C_" : "mmf", "G_/CA" : "mmf", "G_/CC" : "mmf", "G_/CG" : "mmf", "G_/CU" : "mmf",
    "GA/C_" : "mmf", "GA/CA" : "mmf", "GA/CC" : "mmf", "GA/CG" : "mmf", "GA/CU" : "mmf",
    "GC/C_" : "mmf", "GC/CA" : "mmf", "GC/CC" : "mmf", "GC/CG" : "mmf", "GC/CU" : "mmf",
    "GG/C_" : "mmf", "GG/CA" : "mmf", "GG/CC" : "mmf", "GG/CG" : "mmf", "GG/CU" : "mmf",
    "GU/C_" : "mmf", "GU/CA" : "mmf", "GU/CC" : "mmf", "GU/CG" : "mmf", "GU/CU" : "mmf",
    
    # GU
    "G_/U_" : "mmf", "G_/UA" : "mmf", "G_/UC" : "mmf", "G_/UG" : "mmf", "G_/UU" : "mmf",
    "GA/U_" : "mmf", "GA/UA" : "mmf", "GA/UC" : "mmf", "GA/UG" : "mmf", "GA/UU" : "mmf",
    "GC/U_" : "mmf", "GC/UA" : "mmf", "GC/UC" : "mmf", "GC/UG" : "mmf", "GC/UU" : "mmf",
    "GG/U_" : "mmf", "GG/UA" : "mmf", "GG/UC" : "mmf", "GG/UG" : "mmf", "GG/UU" : "mmf",
    "GU/U_" : "mmf", "GU/UA" : "mmf", "GU/UC" : "mmf", "GU/UG" : "mmf", "GU/UU" : "mmf",
    
    # UG
    "U_/G_" : "mmf", "U_/GA" : "mmf", "U_/GC" : "mmf", "U_/GG" : "mmf", "U_/GU" : "mmf",
    "UA/G_" : "mmf", "UA/GA" : "mmf", "UA/GC" : "mmf", "UA/GG" : "mmf", "UA/GU" : "mmf",
    "UC/G_" : "mmf", "UC/GA" : "mmf", "UC/GC" : "mmf", "UC/GG" : "mmf", "UC/GU" : "mmf",
    "UG/G_" : "mmf", "UG/GA" : "mmf", "UG/GC" : "mmf", "UG/GG" : "mmf", "UG/GU" : "mmf",
    "UU/G_" : "mmf", "UU/GA" : "mmf", "UU/GC" : "mmf", "UU/GG" : "mmf", "UU/GU" : "mmf",
    
    # AU
    "A_/U_" : "mmf", "A_/UA" : "mmf", "A_/UC" : "mmf", "A_/UG" : "mmf", "A_/UU" : "mmf",
    "AA/U_" : "mmf", "AA/UA" : "mmf", "AA/UC" : "mmf", "AA/UG" : "mmf", "AA/UU" : "mmf",
    "AC/U_" : "mmf", "AC/UA" : "mmf", "AC/UC" : "mmf", "AC/UG" : "mmf", "AC/UU" : "mmf",
    "AG/U_" : "mmf", "AG/UA" : "mmf", "AG/UC" : "mmf", "AG/UG" : "mmf", "AG/UU" : "mmf",
    "AU/U_" : "mmf", "AU/UA" : "mmf", "AU/UC" : "mmf", "AU/UG" : "mmf", "AU/UU" : "mmf",
    
    # UA
    "U_/A_" : "mmf", "U_/AA" : "mmf", "U_/AC" : "mmf", "U_/AG" : "mmf", "U_/AU" : "mmf",
    "UA/A_" : "mmf", "UA/AA" : "mmf", "UA/AC" : "mmf", "UA/AG" : "mmf", "UA/AU" : "mmf",
    "UC/A_" : "mmf", "UC/AA" : "mmf", "UC/AC" : "mmf", "UC/AG" : "mmf", "UC/AU" : "mmf",
    "UG/A_" : "mmf", "UG/AA" : "mmf", "UG/AC" : "mmf", "UG/AG" : "mmf", "UG/AU" : "mmf",
    "UU/A_" : "mmf", "UU/AA" : "mmf", "UU/AC" : "mmf", "UU/AG" : "mmf", "UU/AU" : "mmf",
    
    # XY
    "X_/Y_" : "mmf", "X_/YA" : "mmf", "X_/YC" : "mmf", "X_/YG" : "mmf", "X_/YU" : "mmf",
    "XA/Y_" : "mmf", "XA/YA" : "mmf", "XA/YC" : "mmf", "XA/YG" : "mmf", "XA/YU" : "mmf",
    "XC/Y_" : "mmf", "XC/YA" : "mmf", "XC/YC" : "mmf", "XC/YG" : "mmf", "XC/YU" : "mmf",
    "XG/Y_" : "mmf", "XG/YA" : "mmf", "XG/YC" : "mmf", "XG/YG" : "mmf", "XG/YU" : "mmf",
    "XU/Y_" : "mmf", "XU/YA" : "mmf", "XU/YC" : "mmf", "XU/YG" : "mmf", "XU/YU" : "mmf"
}


# "forward" dangle
dd53fbps =   {
    # d5
    "__/__" : "d5f", "_A/__" : "d5f", "_C/__" : "d5f", "_G/__" : "d5f", "_U/__" : "d5f",
    "C_/G_" : "d5f", "CA/G_" : "d5f", "CC/G_" : "d5f", "CG/G_" : "d5f", "CU/G_" : "d5f",
    "G_/C_" : "d5f", "GA/C_" : "d5f", "GC/C_" : "d5f", "GG/C_" : "d5f", "GU/C_" : "d5f",
    "G_/U_" : "d5f", "GA/U_" : "d5f", "GC/U_" : "d5f", "GG/U_" : "d5f", "GU/U_" : "d5f",
    "U_/G_" : "d5f", "UA/G_" : "d5f", "UC/G_" : "d5f", "UG/G_" : "d5f", "UU/G_" : "d5f",
    "A_/U_" : "d5f", "AA/U_" : "d5f", "AC/U_" : "d5f", "AG/U_" : "d5f", "AU/U_" : "d5f",
    "U_/A_" : "d5f", "UA/A_" : "d5f", "UC/A_" : "d5f", "UG/A_" : "d5f", "UU/A_" : "d5f",
    "X_/Y_" : "d5f", "XA/Y_" : "d5f", "XC/Y_" : "d5f", "XG/Y_" : "d5f", "XU/Y_" : "d5f",
    
    # d3
    "__/__" : "d3f", "__/_A" : "d3f", "__/_C" : "d3f", "__/_G" : "d3f", "__/_U" : "d3f",
    "C_/G_" : "d3f", "C_/GA" : "d3f", "C_/GC" : "d3f", "C_/GG" : "d3f", "C_/GU" : "d3f",
    "G_/C_" : "d3f", "G_/CA" : "d3f", "G_/CC" : "d3f", "G_/CG" : "d3f", "G_/CU" : "d3f",
    "G_/U_" : "d3f", "G_/UA" : "d3f", "G_/UC" : "d3f", "G_/UG" : "d3f", "G_/UU" : "d3f",
    "U_/G_" : "d3f", "U_/GA" : "d3f", "U_/GC" : "d3f", "U_/GG" : "d3f", "U_/GU" : "d3f",
    "A_/U_" : "d3f", "A_/UA" : "d3f", "A_/UC" : "d3f", "A_/UG" : "d3f", "A_/UU" : "d3f",
    "U_/A_" : "d3f", "U_/AA" : "d3f", "U_/AC" : "d3f", "U_/AG" : "d3f", "U_/AU" : "d3f",
    "X_/Y_" : "d3f", "X_/YA" : "d3f", "X_/YC" : "d3f", "X_/YG" : "d3f", "X_/YU" : "d3f"
}


dd53rbps =   {
    # d5 (reverse direction)
    "__/__" : "d5r", "__/A_" : "d5r", "__/C_" : "d5r", "__/G_" : "d5r", "__/U_" : "d5r",
    "_G/_C" : "d5r", "_G/AC" : "d5r", "_G/CC" : "d5r", "_G/GC" : "d5r", "_G/UC" : "d5r",
    "_C/_G" : "d5r", "_C/AG" : "d5r", "_C/CG" : "d5r", "_C/GG" : "d5r", "_C/UG" : "d5r",
    "_U/_G" : "d5r", "_U/AG" : "d5r", "_U/CG" : "d5r", "_U/GG" : "d5r", "_U/UG" : "d5r",
    "_G/_U" : "d5r", "_G/AU" : "d5r", "_G/CU" : "d5r", "_G/GU" : "d5r", "_G/UU" : "d5r",
    "_U/_A" : "d5r", "_U/AA" : "d5r", "_U/CA" : "d5r", "_U/GA" : "d5r", "_U/UA" : "d5r",
    "_A/_U" : "d5r", "_A/AU" : "d5r", "_A/CU" : "d5r", "_A/GU" : "d5r", "_A/UU" : "d5r",
    "_Y/_X" : "d5r", "_Y/AX" : "d5r", "_Y/CX" : "d5r", "_Y/GX" : "d5r", "_Y/UX" : "d5r",
    
    # d3 (reverse direction)
    "__/__" : "d3r", "A_/__" : "d3r", "C_/__" : "d3r", "G_/__" : "d3r", "U_/__" : "d3r",
    "_G/_C" : "d3r", "AG/_C" : "d3r", "CG/_C" : "d3r", "GG/_C" : "d3r", "UG/_C" : "d3r",
    "_C/_G" : "d3r", "AC/_G" : "d3r", "CC/_G" : "d3r", "GC/_G" : "d3r", "UC/_G" : "d3r",
    "_U/_G" : "d3r", "AU/_G" : "d3r", "CU/_G" : "d3r", "GU/_G" : "d3r", "UU/_G" : "d3r",
    "_G/_U" : "d3r", "AG/_U" : "d3r", "CG/_U" : "d3r", "GG/_U" : "d3r", "UG/_U" : "d3r",
    "_U/_A" : "d3r", "AU/_A" : "d3r", "CU/_A" : "d3r", "GU/_A" : "d3r", "UU/_A" : "d3r",
    "_A/_U" : "d3r", "AA/_U" : "d3r", "CA/_U" : "d3r", "GA/_U" : "d3r", "UA/_U" : "d3r",
    "_Y/_X" : "d3r", "AY/_X" : "d3r", "CY/_X" : "d3r", "GY/_X" : "d3r", "UY/_X" : "d3r"
}




def is_symmetric(seq1, seq2):
    """@

    determines if a sequence and its complement are idential
    (symmetric) or not
    
    """
    is_sym = False
    debug = False
    if debug:
        print ("seq1: ", seq1)
        print ("seq2: ", seq2)
    #
    
    # Note: to reverse a sequence, you can use the following
    # a=''
    # a.join(reversed(seq2))
    
    if seq2 == seq1:
        is_sym = True
    #
    
    return is_sym
#

def has_terml_AU(seq1, seq2):
    """@

    determines if the input duplex has a terminal AU base pair or not
    
    """
    tAU_count = 0
    debug = False
    if debug:
        print ("seq1: ", seq1)
        print ("seq2: ", seq2)
    #
    
    end = len(seq1)-1
    if seq1[0] == 'a' and seq2[0] == 'u':
        tAU_count += 1
    elif seq1[0] == 'u' and seq2[0] == 'a':
        tAU_count += 1
    #
    
    if seq1[end] == 'a' and seq2[end] == 'u':
        tAU_count += 1
    elif seq1[end] == 'u' and seq2[end] == 'a':
        tAU_count += 1
    #
    
    if debug:
        print ("tAU_count = %d" % tAU_count)
        # if tAU_count == 1:
        #    sys.exit(0)
    #
    
    return tAU_count
#

