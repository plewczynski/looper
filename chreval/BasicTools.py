#!/usr/bin/env python3

"""@@@

Main Module:   BasicTools.py
 
Object:        None

Functions:     copyPairList
               initialize_matrix
               roundoff         
               sortPairListWRT_n
               tuple2PairList

Author:        Wayne Dawson
creation date: 190221
last update:   200210 (upgrade to python3)
version:       0.0

Purpose:

This contains some very simple functions that can be used for a
multitude of purposes and don't need to be rewritten every time they
are used.

"""

from copy import deepcopy
import sys
from math import sqrt


def sortPairListWRT_n(listX, comp_n):
    # + listX should be a list of tuples (or a list of lists)

    # + comp_n is the particular element in the tuple (or list)
    
    # sorts an input list (of lists) with respect to a specific
    # component (comp_n). For example, suppose we have the
    # following list of tuples
    
    # listX=[(5,13), (21, 25), (1,18)].
    
    # If we chose comp_n = 0, then the output list becomes
    
    # listX=[(1,18), (5,13), (21, 25)].
    
    # But if we chose comp_n = 1, then the output list becomes
    
    # listX=[(5,13), (1,18), (21, 25)].
    
    
    newlistX = deepcopy(listX)
    newlistX = sorted(newlistX, key=lambda x: x[comp_n])
    # from https://stackoverflow.com/questions/10695139/sort-a-list-of-tuples-by-2nd-item-integer-value
    
    # alternatively, for an in-place sort ... 
    
    # foo = [(list of tuples)]
    # foo.sort(key=lambda x:x[0]) # To sort by first element of the tuple in place
    return newlistX
#


# Makes a copy of a list. Useful when you don't want the original list
# to be disturbed. Recall that Python defines a list as a pointer, so
# if you just do myList2 = myList1, then if you delete/add an item in
# myList2, it will also be deleted/added to myList1. However, if you
# do myList2 = copyList(myList1), then if you alter myList2, it does
# not alter myList1.
def copyList(listX):
    newListX = deepcopy(listX)
    return newListX
#

# NOTE: originally, I wrote the program this way:

# def copyList(listX):
#     newListX = []
#     for lX_k in listX:
#         # we must physically write this because we must make sure
#         # not to make a pointer; e.g., `pklist = self.PKlist`
#         newListX += [lX_k]
#     #
# #

# However, I noticed that this does not necessarily remove the issue
# when it is a list of lists; e.g., [[1,20], [25, 40]] -> [[1,40]]. I
# cannot reproduce what I noticed, but somehow, the array inside the
# copy was also altered in the original and in the copy (for at least
# one instance I encountered). The safest way to avoid problems is to
# make a deep copy when it matters how the list might be modified.

# see the following information for more details
# https://stackoverflow.com/questions/2612802/how-to-clone-or-copy-a-list/




# Start with a list of tuples or even a tuple of tuples and convert
# everything to lists. This is usually needed if variables in the
# tuple list will be changed.
def tuple2List(tupleList):
    newList = []
    for nlk in tupleList:
        # expand the tuple (nlk) and rewrite as a List (a):
        a = []
        for ll in range(0, len(nlk)):
            a += [nlk[ll]]
        #
        newList += [a]
    #
    return newList
#


# round off methods
def roundoff(v, dplace):
    if not isinstance(dplace, int):
        print ("ERROR(roundoff) requires a integer for the decimal place")
        print ("                input: (%s)" % dplace)
        sys.exit(1)
    #
    # print ("dplace = ", dplace)
    
    vr = 0.0
    w = float(10**dplace)
    if v < 0.0:
        vr =  float(int(w*v - 0.5))/w
    else:
        vr =  float(int(w*v + 0.5))/w
    #
    return vr
#
        
    
def initialize_matrix(Var, N, w):
    """@

    examples:
    
    initialize_matrix(v, N, 0.0)
    initialize_matrix(v, N, 1)
    
    v = the matrix variable name
    
    N = the size of the matrix (N x N).
    
    w = the value to be used.
    
    Because the python program doesn't care, the type can be either
    float or integer. It can also be used to assign vectors such as
    M-loops: cc = [(0,0)...]. Because this simply acts as a place
    holder, any sort of object can be inserted, as long as it is the
    same for all matrix elements.
    
    However, this also mean that you have to be very careful what you
    enter into the function, because python does not care what you
    insert into this thing (integers, floats, vectors, objects). If
    you don't pay attention to what goes in, garbage is likely to
    follow ... python doesn't care. Of course, this sort of thing
    could also be done with C++.
    
    """
        
    Var = []
    for j in range(0, N):
        vj = []
        for i in range(0, N):
            vj += [ w ]
        #
        Var += [vj]
    return Var
#


def kickOnAt(i, i_set, j, j_set):
    """@
    
    kick on at
    
    Turns on debugging parameters when the specific indices i and j
    meet the condition that (i,j) = (i_set, j_set), and turns these
    debugging parameters off under other conditions.

    """
    
    flag_on = False
    if i_set == i and j == j_set:
        flag_on = True
    #
    
    return flag_on
#

def kickOnBetween(i, i_set, j, j_set):
    """@
    
    kick on between
    
    Turns on debugging parameters when the specific indices i and j
    meet the condition that i_set <= i and j <= j_set, and turns these
    debugging parameters off when i < i_set and/or j > j_set.

    """
    
    flag_on = False
    if i_set == i and j == j_set:
        flag_on = True
    #
    
    return flag_on
#

def kickOnOutside(i, i_set, j, j_set):
    """@
    
    kick on between
    
    Turns on debugging parameters when the specific indices i and j
    meet the condition that i_set > i and j > j_set, and turns these
    debugging parameters off when i_set <= i and j <= j_set.

    """
    
    flag_on = False
    if i_set == i and j == j_set:
        flag_on = True
    #
    
    return flag_on
#


class CountDown(object):
    """@
    
    This is sometimes useful for debugging or checking output if you
    want the program to go through a section n times before stopping,
    not just once.
    
    """
    
    def __init__(self, set_exit_count):
        self.set_exit_count = set_exit_count
        self.count_index = 0
    #
    
    def check_count(self):
        if self.count_index == self.set_exit_count:
            print ("planned stop when this section has been crossed %d times" \
                % self.count_index)
            sys.exit(0)
        #
    #
    
    def increment_CD(self):
        self.count_index += 1
        self.check_count()
    #
#
    

class Index(object):
    
    """@
    
    191229: This class was developed in RNA/Motif.py but has not been 
    used. Originally, I had the idea to make Map a 1D array
    functioning like vsfold ssIndx where it returns the full number
    without the weird style of the Vienna package where they use i +
    shift[j]. However, I am likely to expand this to duplexes where
    this would not work very well, so I am reluctant to take this
    up. Still, I decided to keep this here.
    
    This is how the numbering works
    
    ij = i + j*(j-1)//2  # the // is for integer division
    
    For example:
      i   j   ij  jshift
      0   1   0     0
      0   2   1     1
      1   2   2     2
      0   3   3     3
      1   3   4     3
      2   3   5     3
      0   4   6     6
      1   4   7     6
      2   4   8     6
      3   4   9     6
      0   5  10    10
      1   5  11    10
      2   5  12    10
      3   5  13    10
      4   5  14    10
      0   6  15    15
      etc.
    
    Note that in this scheme, we must require that i < j. The output
    gij(i,j) checks to make sure i < j < N.

    """
    
    def __init__(self, N):
        self.N = N
        
        self.jshift = []
        for j in range(0,N+1):
            self.jshift += [j*(j-1)//2]
        #
    #
    
    def gij(self, i, j):
        if j <= 0 or i < 0 or i >= self.N or j >= self.N or i == j:
            
            if i < 0:
                print ("ERROR(Motif.ndx): cannot accept negative numbers for j(%d)" % (j))
            #
            
            if j <= 0:
                print ("ERROR(Motif.ndx): cannot accept negative numbers for j(%d)" % (j))
            #
            
            if i >= self.N:
                print ("ERROR(Motif.ndx): cannot accept values of i(%d) >= %d" % (i, self.N))
            #
            
            if j >= self.N:
                print ("ERROR(Motif.ndx): cannot accept values of j(%d) >= %d" % (j, self.N))
            #
            if i == j:
                print ("ERROR(Motif.ndx) cannot accept i(%d) = j(%d)" % (i, j))
            #
            sys.exit(1)
        #
        
        v = 0
        if   i < j:
            v = i + self.jshift[j]
        else:
            v = j + self.jshift[i]
        #
        
        return v
    #
    
    
    
    def get_ij(self, v):
        N = self.N - 1
        v_max = N*(N-1)
        
        if v > v_max:
            print ("ERROR(Motif.get_ij()),  1D index out of range v(%d) = v_max(%d)" % (v, v_max))
            sys.exit(0)
        #
        
        if v < 0:
            print ("ERROR(Motif.get_ij()),  1D index out of range v(%d) < 0" % (v))
            sys.exit(0)
        #
        
        if v <= 1:
            if v == 0:
                return 0, 1 # jshift = 0
            else:
                return 0, 2 # jshift = 1
            #
            
        #
        
        k = int(sqrt(v)) # start somewhere short of the closest possibility
        # print (k)
        cont = True
        while cont:
            #print (k, self.jshift[k], v)
            if k < self.N:
                if self.jshift[k-1] <= v and v < self.jshift[k]:
                    k -= 1
                    cont = False
                else:
                    k += 1
                #
            else:
                k = N
                cont = False
            #
            
        #
        
        j = k
        i = v - self.jshift[k]
        
        return i, j
    
#


def invertDict(fwdDict):
    
    """@
    
    This inverts an existing dictionary with the following
    transformation: {key : tag} -> {tag : key}. This assumes that the
    dictionary is one-to-one, so it will be a serious problem with
    loss of information if that is not the case.
    
    
    method for one-to-one reversal of dictionary explained here
    https://stackoverflow.com/questions/2568673/inverse-dictionary-lookup-in-python
    
    reverse { key : index } to { index : key }
      labels = { "base"    : 0,
                 "height"  : 1,
                 "spread"  : 2,
                 "Tm"      : 3,
                 "slope"   : 4 }
    
    transforms to
    
      invlbls   = { 0 : "base",        
                    1 : "height",
                    2 : "spread",
                    3 : "Tm",
                    4 : "slope" }
    
    Remember that this approach is not very smart, so if you have a
    dictionary like
    
      labels = { "base1" : 0,
                 "base2" : 0 }
    
    transforms to
    
      invlbls   = { 0 : "base1" }
    
    So you willlose information. Therefore, to have a proper
    transformation, you need a one-to-one correspondence between both
    the key and the index.
    
    """    
    invDict = dict((key, tag) for tag,key in fwdDict.items())
    return invDict
#

    

def test0():
    print (roundoff(0.453, 2))
    print (roundoff(-0.453, 2))
#

def test1():
    ndx = Index(20)
    print (ndx.jshift)
    print (ndx.get_ij(0))
    print (ndx.get_ij(1))
    print (ndx.get_ij(2))
    print (ndx.get_ij(3))
    print (ndx.get_ij(37))
    print (ndx.get_ij(110))
#


if __name__ == "__main__":
    test1()
#

