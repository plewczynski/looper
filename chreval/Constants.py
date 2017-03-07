#!/usr/bin/env python

# #################################################################
# ###############  General configuration CONSTANTS  ###############
# ###############    settings used in FreeEnergy    ###############
# #################################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# These are used to assign parameters in FreeEnergy()

INFINITY = 1000.0
# 161019wkd: I know infinity is not 1000, but I need to set some upper
# bound on the program on the program where generally nothing can be
# found by typical free energy. Generally this appears to be a large
# enough number to be treated as "infinity". I made it adjustable for
# the sake of raising the bar if necessary.

# Originally, I had tried to rid the code completely of these sorts of
# artificial features, but it seems that this is not so easy to do and
# they end up cropping up somewhere in the code whatever I
# do. Perhaps, if the aesthetics get a bit too annoying, another way
# is to run the program with this artifical value and, at the end,
# rescale this upper bound to the maximum (if it is positive) or zero,
# if not. At any rate, whatever number I use, it should be at least as
# positive as the largest positive calculated value in the whole lot.

# fundamental coarse grained resolution length
seg_len = 5000.0

# constants: in entropy evaluation
kB       = 0.00198 # [kcal/mol]
xi       = 5.0
lmbd     = 0.002  #   
gmm      = 2.3
# constants: weights for the enthalpy terms 
febase   = -6.0   #febase   = -4.7
feshift  =  1.0

# Presently, the enthalpy is based on a logarithmic assessment with
# respect to the number of observed counts of a particular interaction
# between parts (i,j) of the chromatin chain. It is clearly a crude
# function to assess this binding free energy, but I have no real
# information to work with here. For a given set of data, these two
# parameters may surely require tuning. The current test set showed
# somewhat reasonable results when setting a baseline of 4 kcal/mol
# for the binding free energy of the CTCF clusters. Whether this
# estimate is remotely correct or not is currently unknown.

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# #################################################################



# #################################################################
# ###############  General configuration CONSTANTS  ###############
# ###############   used in 1D structure notation   ###############
# #################################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv


# used for PKs and parallel stems
# this notation at least works with VARNA
# used for PKs and parallel stems
# this notation at least works with VARNA
num2lpr = {  0 : '(',
             1 : '[',
             2 : '{',
             3 : '<',
             4 : 'A',
             5 : 'B',
             6 : 'C',
             7 : 'D',
             8 : 'E',
             9 : 'F',
            10 : 'G',
            11 : 'H',
            12 : 'I',
            13 : 'J',
            14 : 'K',
            15 : 'L',
            16 : 'M',
            17 : 'N',
#           18 : 'O', # doesn't work beyond N
#           19 : 'P',
#           20 : 'Q',
#           21 : 'R',
#           22 : 'S',
#           23 : 'T',
#           24 : 'U',
#           25 : 'V',
#           26 : 'W',
#           27 : 'X',
#           28 : 'Y',
#           29 : 'Z',
}

# used for PKs and parallel stems
# this notation at least works with VARNA
num2rpr = {  0 : ')',
             1 : ']',
             2 : '}',
             3 : '>',
             4 : 'a',
             5 : 'b',
             6 : 'c',
             7 : 'd',
             8 : 'e',
             9 : 'f',
            10 : 'g',
            11 : 'h',
            12 : 'i',
            13 : 'j',
            14 : 'k',
            15 : 'l',
            16 : 'm',
            17 : 'n',
#           18 : 'o', # doesn't work beyond N
#           19 : 'p',
#           20 : 'q',
#           21 : 'r',
#           22 : 's',
#           23 : 't',
#           24 : 'u',
#           25 : 'v',
#           26 : 'w',
#           27 : 'x',
#           28 : 'y',
#           29 : 'z',
}


lpr2num = {  '(' :  0,
             '[' :  1,
             '{' :  2,
             '<' :  3,
             'A' :  4,
             'B' :  5,
             'C' :  6,
             'D' :  7,
             'E' :  8,
             'F' :  9,
             'G' : 10,
             'H' : 11,
             'I' : 12,
             'J' : 13,
             'K' : 14,
             'L' : 15,
             'M' : 16,
             'N' : 17,
#            'O' : 18,  # doesn't work beyond N
#            'P' : 19, 
#            'Q' : 20, 
#            'R' : 21, 
#            'S' : 22, 
#            'T' : 23, 
#            'U' : 24, 
#            'V' : 25, 
#            'W' : 26, 
#            'X' : 27, 
#            'Y' : 28, 
#            'Z' : 29, 
}

rpr2num = {  ')' :  0,
             ']' :  1,
             '}' :  2,
             '>' :  3,
             'a' :  4,
             'b' :  5,
             'c' :  6,
             'd' :  7,
             'e' :  8,
             'f' :  9,
             'g' : 10,
             'h' : 11,
             'i' : 12,
             'j' : 13,
             'k' : 14,
             'l' : 15,
             'm' : 16,
             'n' : 17,
#            'o' : 18, # doesn't work beyond N
#            'p' : 19,
#            'q' : 20,
#            'r' : 21,
#            's' : 22,
#            't' : 23,
#            'u' : 24,
#            'v' : 25,
#            'w' : 26,
#            'x' : 27,
#            'y' : 28,
#            'z' : 29,
}

# the next two dictionaries can be used to increment between index i
# and i+1 and reference a particular character. For example, say you
# are using "Ee". Then lpr2num[E] = 8 / rpr2num[e] = 8,
# PKrndx[rpr2num[e]] = 6 and so if we increment to 7,
# rpr2num[PKfndx[7]], num2lpr[7] = 'F', num2rpr[7] = 'f'.
PKfndx  = {  0 :  1,  # ['[', ']'],
             1 :  3,  # ['<', '>'],
             2 :  4,  # ['A', 'a'],
             3 :  5,  # ['B', 'b'],
             4 :  6,  # ['C', 'c'],
             5 :  7,  # ['D', 'd'],
             6 :  8,  # ['E', 'e'],
             7 :  9,  # ['F', 'f'],
             8 : 10,  # ['G', 'g'],
             9 : 11,  # ['H', 'h'],
            10 : 12,  # ['I', 'i'],
            11 : 13,  # ['J', 'j'],
            12 : 14,  # ['K', 'k'],
            13 : 15,  # ['L', 'l'],
            14 : 16,  # ['M', 'm'],
            15 : 17,  # ['N', 'n'],
#           16 : 18,  # ['O', 'o'], # doesn't work beyond N
#           17 : 19,  # ['P', 'p'],
#           18 : 20,  # ['Q', 'q'],
#           19 : 21,  # ['R', 'r'],
#           20 : 22,  # ['S', 's'],
#           21 : 23,  # ['T', 't'],
#           22 : 24,  # ['U', 'u'],
#           23 : 25,  # ['V', 'v'],
#           24 : 26,  # ['W', 'w'],
#           25 : 27,  # ['X', 'x'],
#           26 : 28,  # ['Y', 'y'],
#           27 : 29,  # ['Z', 'z'],
}


PKrndx  = {  1 :  0,  # ['[', ']'],
             3 :  1,  # ['<', '>'],
             4 :  2,  # ['A', 'a'],
             5 :  3,  # ['B', 'b'],
             6 :  4,  # ['C', 'c'],
             7 :  5,  # ['D', 'd'],
             8 :  6,  # ['E', 'e'],
             9 :  7,  # ['F', 'f'],
            10 :  8,  # ['G', 'g'],
            11 :  9,  # ['H', 'h'],
            12 : 10,  # ['I', 'i'],
            13 : 11,  # ['J', 'j'],
            14 : 12,  # ['K', 'k'],
            15 : 13,  # ['L', 'l'],
            16 : 14,  # ['M', 'm'],
            17 : 15,  # ['N', 'n'],
#           18 : 16,  # ['O', 'o'], # doesn't work beyond N
#           19 : 17,  # ['P', 'p'],
#           20 : 18,  # ['Q', 'q'],
#           21 : 19,  # ['R', 'r'],
#           22 : 20,  # ['S', 's'],
#           23 : 21,  # ['T', 't'],
#           24 : 22,  # ['U', 'u'],
#           25 : 23,  # ['V', 'v'],
#           26 : 24,  # ['W', 'w'],
#           27 : 25,  # ['X', 'x'],
#           28 : 26,  # ['Y', 'y'],
#           29 : 27,  # ['Z', 'z'],
}


# The final dictionary is for counting the number of items on the
# stack; i.e., the stack pointer.

pointer = {  0 :  0,  # ['(', ')'],
             1 :  0,  # ['[', ']'],
             2 :  0,  # ['{', '}'],
             3 :  0,  # ['<', '>'],
             4 :  0,  # ['A', 'a'],
             5 :  0,  # ['B', 'b'],
             6 :  0,  # ['C', 'c'],
             7 :  0,  # ['D', 'd'],
             8 :  0,  # ['E', 'e'],
             9 :  0,  # ['F', 'f'],
            10 :  0,  # ['G', 'g'],
            11 :  0,  # ['H', 'h'],
            12 :  0,  # ['I', 'i'],
            13 :  0,  # ['J', 'j'],
            14 :  0,  # ['K', 'k'],
            15 :  0,  # ['L', 'l'],
            16 :  0,  # ['M', 'm'],
            17 :  0,  # ['N', 'n'],
#           18 :  0,  # ['O', 'o'], # doesn't work beyond N
#           19 :  0,  # ['P', 'p'],
#           20 :  0,  # ['Q', 'q'],
#           21 :  0,  # ['R', 'r'],
#           22 :  0,  # ['S', 's'],
#           23 :  0,  # ['T', 't'],
#           24 :  0,  # ['U', 'u'],
#           25 :  0,  # ['V', 'v'],
#           26 :  0,  # ['W', 'w'],
#           27 :  0,  # ['X', 'x'],
#           28 :  0,  # ['Y', 'y'],
#           29 :  0,  # ['Z', 'z'],
}

# this references the position
stack   = {  0 :  0,  # ['(', ')'],
             1 :  0,  # ['[', ']'],
             2 :  0,  # ['{', '}'],
             3 :  0,  # ['<', '>'],
             4 :  0,  # ['A', 'a'],
             5 :  0,  # ['B', 'b'],
             6 :  0,  # ['C', 'c'],
             7 :  0,  # ['D', 'd'],
             8 :  0,  # ['E', 'e'],
             9 :  0,  # ['F', 'f'],
            10 :  0,  # ['G', 'g'],
            11 :  0,  # ['H', 'h'],
            12 :  0,  # ['I', 'i'],
            13 :  0,  # ['J', 'j'],
            14 :  0,  # ['K', 'k'],
            15 :  0,  # ['L', 'l'],
            16 :  0,  # ['M', 'm'],
            17 :  0,  # ['N', 'n'],
#           18 :  0,  # ['O', 'o'], # doesn't work beyond N
#           19 :  0,  # ['P', 'p'],
#           20 :  0,  # ['Q', 'q'],
#           21 :  0,  # ['R', 'r'],
#           22 :  0,  # ['S', 's'],
#           23 :  0,  # ['T', 't'],
#           24 :  0,  # ['U', 'u'],
#           25 :  0,  # ['V', 'v'],
#           26 :  0,  # ['W', 'w'],
#           27 :  0,  # ['X', 'x'],
#           28 :  0,  # ['Y', 'y'],
#           29 :  0,  # ['Z', 'z'],
}

# this is mostly used for check
counter = {  0 :  0,  # ['(', ')'],
             1 :  0,  # ['[', ']'],
             2 :  0,  # ['{', '}'],
             3 :  0,  # ['<', '>'],
             4 :  0,  # ['A', 'a'],
             5 :  0,  # ['B', 'b'],
             6 :  0,  # ['C', 'c'],
             7 :  0,  # ['D', 'd'],
             8 :  0,  # ['E', 'e'],
             9 :  0,  # ['F', 'f'],
            10 :  0,  # ['G', 'g'],
            11 :  0,  # ['H', 'h'],
            12 :  0,  # ['I', 'i'],
            13 :  0,  # ['J', 'j'],
            14 :  0,  # ['K', 'k'],
            15 :  0,  # ['L', 'l'],
            16 :  0,  # ['M', 'm'],
            17 :  0,  # ['N', 'n'],
#           18 :  0,  # ['O', 'o'], # doesn't work beyond N
#           19 :  0,  # ['P', 'p'],
#           20 :  0,  # ['Q', 'q'],
#           21 :  0,  # ['R', 'r'],
#           22 :  0,  # ['S', 's'],
#           23 :  0,  # ['T', 't'],
#           24 :  0,  # ['U', 'u'],
#           25 :  0,  # ['V', 'v'],
#           26 :  0,  # ['W', 'w'],
#           27 :  0,  # ['X', 'x'],
#           28 :  0,  # ['Y', 'y'],
#           29 :  0,  # ['Z', 'z'],
}


# This contains that actual base pairing information for each of the
# label types; i.e., '(', '[', etc. It is a useful buffer that permits
# the final sorting of information into various categories of such as
# secondary structure (in BPlist), pseudoknots (in PKlist) and CTCF
# islands (in CTCFlist). Using the above conversion formulae, each of
# these labels ('(', '[', 'A', etc.) is identified with an index (0,
# 1, 2, 3...). Hence, all the various letters are organized so that
# they can be sorted to build the final lists.
Xlist   = {  0 :  [],  # ['(', ')'],
             1 :  [],  # ['[', ']'],
             2 :  [],  # ['{', '}'],
             3 :  [],  # ['<', '>'],
             4 :  [],  # ['A', 'a'],
             5 :  [],  # ['B', 'b'],
             6 :  [],  # ['C', 'c'],
             7 :  [],  # ['D', 'd'],
             8 :  [],  # ['E', 'e'],
             9 :  [],  # ['F', 'f'],
            10 :  [],  # ['G', 'g'],
            11 :  [],  # ['H', 'h'],
            12 :  [],  # ['I', 'i'],
            13 :  [],  # ['J', 'j'],
            14 :  [],  # ['K', 'k'],
            15 :  [],  # ['L', 'l'],
            16 :  [],  # ['M', 'm'],
            17 :  [],  # ['N', 'n'],
#           18 :  [],  # ['O', 'o'], # doesn't work beyond N
#           19 :  [],  # ['P', 'p'],
#           20 :  [],  # ['Q', 'q'],
#           21 :  [],  # ['R', 'r'],
#           22 :  [],  # ['S', 's'],
#           23 :  [],  # ['T', 't'],
#           24 :  [],  # ['U', 'u'],
#           25 :  [],  # ['V', 'v'],
#           26 :  [],  # ['W', 'w'],
#           27 :  [],  # ['X', 'x'],
#           28 :  [],  # ['Y', 'y'],
#           29 :  [],  # ['Z', 'z'],
}



# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# #################################################################
