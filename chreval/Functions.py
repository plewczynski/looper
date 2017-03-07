#!/usr/bin/env python

# Important functions used by several program in this package

from math import exp, log, ceil, floor
from difflib import SequenceMatcher # for similar
import sys


def make_file_heading(name, bgn, end, bp_resolution):
    return "%s_%d_%d_%s" % (name, bgn, end, bp_resolution)
#
          

# Measure the similarity of two strings and return a fractional
# percentage of the similarity.
def similar(a, b):
   return SequenceMatcher(None, a, b).ratio()
#

def hamming_bin(s1, s2):
    """Calculate the Hamming distance between two bit strings"""
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))
#

def hamming_str(str1, str2):
    """Count the # of differences between equal length strings str1 and str2"""
    
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
        #
    #
    return diffs
#



# (IMPORTANT!!!) For the partition function (in my particular case), I
# sometimes reach the limits of the IEEE 754 double-precision
# standard, which occurs when I approach "math.exp(709.783)". So this
# function is needed if have to compute arguments in the exponential
# that are larger than this limit 709.783.  According to Wikipedia,

# https://en.wikipedia.org/wiki/Kahan_summation_algorithm

# using the Kahan algorithm is a fairly accurate compromise for
# computing large numbers like this without losing as much precision
# as would normally be the case, or, worse yet, just having to give
# up. In general, the list should be sorted before this operation is
# done so as to limit the times it must be corrected. In the case of
# the computational program, this is already done.

# KahanSumExp is an implementation of the Kahan summation algorithm in
# python. The solution found in stackexchange in the scientific
# computation

# http://scicomp.stackexchange.com/questions/1122/how-to-add-large-exponential-terms-reliably-without-overflow-errors
#

# 161022wkd: The program can work faster by sorting in ascending
# order. Since I want it in desending order, one way to keep this
# practice is to include "expvalues.reverse().  Anyway, when running
# test5(), the data really does need to be sorted first.
flag_KahanSumExp_USE_SORT = False

# have to set what the maximum value is
doubleMAX = (1.0 + (1.0 - (2 ** -52))) * (2 ** (2 ** 10 - 1))

def KahanSumExp(expvalues):
    if flag_KahanSumExp_USE_SORT:
        expvalues.sort() # gives precision improvement in certain cases
    #
    shift = 0 
    esum = 0.0 
    carry = 0.0 
    for exponent in expvalues:
        if exponent - shift * log(2) > 709.783:
            n = ceil((exponent - shift * log(2) - 709.783)/log(2))
            shift += n
            carry /= 2*n
            esum /= 2*n
        elif exponent - shift * log(2) < -708.396:
            n = floor((exponent - shift * log(2) - -708.396)/log(2))
            shift += n
            carry *= 2*n
            esum *= 2*n

        exponent -= shift * log(2)
        value = exp(exponent) - carry 
        if doubleMAX - esum < value:
            shift += 1
            esum /= 2
            value /= 2
        tmp = esum + value 
        carry = (tmp - esum) - value 
        esum = tmp
    return esum, shift
#
