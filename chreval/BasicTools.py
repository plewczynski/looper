#!/usr/bin/env python

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
last update:   190522
version:       0.0

Purpose:

This contains some very simple functions that can be used for a
multitude of purposes and don't need to be rewritten every time they
are used.

"""

from copy import deepcopy


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
        print "ERROR(roundoff) requires a integer for the decimal place"
        print "                input: (%s)" % dplace
        sys.exit(1)
    #
    # print "dplace = ", dplace
    
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


def test0():
    print roundoff(0.453, 2)
    print roundoff(-0.453, 2)
#


if __name__ == "__main__":
    test0()
#

