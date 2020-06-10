#!/usr/bin/env python3

"""@@@

Main Module:     Pair.py

Objects:         Pair
                 SortPair

Functions:       BranchList2list 
                 find_this_ij 
                 is_ap_stem 
                 is_pp_stem
                 rlist2vsPair
                 vsPair2list 

 
Author:          Wayne Dawson
creation date:   170126 (refactored as independent unit 190704)
last update:     200211 (upgraded to python3), 190718
version:         0.1

Purpose:

In an effort to unify these various programs with reusable code, this
part was separated from Motif and is referenced in Motif. By doing
things this way, three diverse packages (chromatin, RNA, and SimRNA)
use this same module.

It does mean that, on the face of it, Pair is over defined for some
modules, particularly the SimRNA restraints. On the other hand, it
means that I can use modules like Restraint from SimRNA without having
to modify anything in the chromatin and RNA packages.

"""


# ####################################################
# #####   Various tools that manage class Pair   #####
# ####################################################

# Some additional tools associated with Pair data that were at least
# used in the past. They don't appear to be used now, but I keep them
# here anyway because I am sure they were used in the Thread, Trace or
# Vienna programs at one time.

def find_this_ij(k_target, Xlist):
    """@@@
    
    This one may be basically specific to the single method
    (scan_for_ss) in Vienna2TreeNode. Anyway, it is able to process
    objects of class Pair (from Vienna). So, for the moment, I leave
    it as a Vienna tool for that reason.

    """
    
    debug_find_this_ij = False # True
    
    
    if debug_find_this_ij:
        print ("enter find_this_ij: look for i = ", k_target)
    #
    i = -1; j = -1; v = '-'
    flag_found = False
    for xk in Xlist:
        i, j, v, c = xk.get_ssPair()
        if debug_find_this_ij:
            print (k_target, i,j)
        #
        if i == k_target:
            if debug_find_this_ij:
                print ("found: ", i, j, v)
            #
            flag_found = True
            break
        #
    #
    if not flag_found:
        if debug_find_this_ij:
            print ("position %d not found!" % k_target)
        #
        i = -1; j = -1
    #
    return (i,j)
#

def is_ap_stem(ndx, Xlist):
    flag_is_ap_stem = False
    i, j, v, c = Xlist[ndx].get_ssPair()
    n = len(Xlist)
    ix = -1; jx = -1
    if ndx < n - 1:
        ix, jx, vx, cx = Xlist[ndx+1].get_ssPair()
        if ix == i+1 and jx == j-1:
            flag_is_ap_stem = True
        #
    #
    return flag_is_ap_stem
#

def is_pp_stem(ndx, Xlist):
    flag_is_pp_stem = False
    i, j, v, c = Xlist[ndx].get_ssPair()
    n = len(Xlist)
    ix = -1; jx = -1
    if ndx < n - 1:
        ix, jx, vx, cx = Xlist[ndx+1].get_ssPair()
        if ix == i+1 and jx == j+1:
            flag_is_pp_stem = True
        #
    #
    return flag_is_pp_stem
#


# ##################################################################
# ########  Fundamental pairing list and associated tools  #########
# ##################################################################

#





class Pair(object):
    def __init__(self):
        self.name = "bp"
        self.ch_i = ''
        self.i = -1
        self.ch_j = ''
        self.j = -1
        self.v = 'a' # default, contact is antiparallel
        # 'a' = antiparallel
        # 'p' = parallel
        # '-' = not applicable; e.g., CTCFs
        self.contacts = []
    #
    
    def __repr__(self):
        return '{}: {}(ch{}{} ch{}{})[{}] {}'.format(self.__class__.__name__,
                                                     self.name,
                                                     self.ch_i,
                                                     self.i,
                                                     self.ch_j,
                                                     self.j,
                                                     self.v,
                                                     self.contacts)
    #
    
    def __cmp__(self, other):
        if hasattr(other, 'i'):
            return self.i.__cmp__(other.i)
        #
    #
    
    
    def put_ssPair(self, i, j, nm = 'bp', v = 'a'):
        self.name = nm
        self.ch_i = 'A'
        self.ch_j = 'A'
        self.i = i
        self.j = j
        self.v = v
        return 0
    #
    
    # this is used!
    def put_ssPair_i(self, i, nm = 'bp', v = 'a'):
        self.name = nm
        self.ch_i = 'A'
        self.i = i
        return 0
    #
    
    # can be used to mark triple helices as well as CTCFs
    def put_contacts(self, i):
        self.contacts += [i]
        return 0
    #
    
    # this is used!
    def put_ssPair_j(self, j, nm = 'bp', v = 'a'):
        if not self.name == nm:
            print ("ERROR(Pair): improperly matched type")
            print ("             (i,j) = (%d,%d)" % (self.i, j))
            print ("             name(%s) != new name(%s) " % (self.name, nm))
        #
        self.ch_j = 'A'
        self.j = j
        return 0
    #
    def get_ssPair(self):
        return self.i, self.j, self.v, self.contacts
    #
    
    
    def disp_Pair(self):
        s = ''
        if len(self.contacts) > 0:
            s  = '(%5d,%5d)[%s]: ' % (self.i, self.j, self.v)
            s += "{%5d, " % self.i
            for c in self.contacts:
                s += "%5d, " % c
            #
            s += "%5d}" % self.j
        else:
            s  = '(%5d,%5d)[%s]' % (self.i, self.j, self.v)
        #
        return s
    #
    
    
    def put_dsPair(self, ch_i, i, ch_j, j, nm = 'bp', v = 'a'):
        self.name = nm
        self.ch_i = ch_i
        self.i = i
        self.ch_j = ch_j
        self.j = j
        self.v = v
        return 0
    #
    
    def put_dsPairi(self, ch_i, i, nm = 'bp', v = 'a'):
        self.name = nm
        self.ch_i = ch_i
        self.i = i
        self.v = v
        return 0
    #
    
    def put_dsPairj(self, ch_j, j, nm = 'bp', v = 'a'):
        self.name = nm
        self.ch_j = ch_j
        self.j = j
        self.v = v
        return 0
    #
    
    
    
    def disp_dsPair(self):
        return  '(%s:%5d,%s:%5d)[%s]' % (self.ch_i, self.i, self.ch_j, self.j, self.v)
    #
#        


class SortPair(object):
    def sortvsList(self, pairlist, comp = 'i'):
        if comp == 'i':
            pairlist =  sorted(pairlist, key=self.getKey_i)
        elif comp == 'j':
            pairlist =  sorted(pairlist, key=self.getKey_j)
        #
        return pairlist
    #
    def getKey_i(self, item):
        return item.i
    #    
    def getKey_j(self, item):
        return item.j
    #    
#

def vsPair2list(pairlist):
    slist = []
    for bpk in pairlist:
        slist += [(bpk.i, bpk.j)]
    #
    return slist
#
def BranchList2list(branchlist):
    slist = []
    for bpk in branchlist:
        slist += [(bpk.i, bpk.j)]
    #
    return slist
#

def rlist2vsPair(rlist, stype = '-', name = "bp"):
    # The variable "name" does not seem to be all that necessary in my
    # experience of using Pair, but I have not deleted that parameter
    # yet. Perhaps it might be useful, for example, in chromatin,
    # there are the ctcfs, so maybe checking the name and the stype
    # might prove useful there in verifying content.
    vsPairList = []
    if stype == '-':
        is_ppstem = False
        for k in range(0, len(rlist)-1):
            i1 = rlist[k  ][0]; j1 = rlist[k  ][1]
            i2 = rlist[k+1][0]; j2 = rlist[k+1][1]
            
            if i1 < i2 and i2 < j1 and j1 < j2:
                is_ppstem = True
                
            elif i1 < i2 and j2 < j1:
                if is_ppstem:
                    print ("ERROR(rlist2vsPair): There seems to be a mixture of ")
                    print ("                     parallel and antiparallel interactions")
                    for vv in rlist:
                        print (vv)
                    #|endfor
                    
                #
                
            #
            
        #|endfor
        
        for k in range(0, len(rlist)):
            iss = rlist[k][0]; jss = rlist[k][1]
            if is_ppstem:
                a = Pair()
                a.put_ssPair(iss, jss, name, 'p')
                vsPairList += [a]
            else:
                a = Pair()
                a.put_ssPair(iss, jss, name, 'a')
                vsPairList += [a]
            #
            
        #|endfor
        
    else:
        if stype == 'p':
            for k in range(0, len(rlist)):
                iss = rlist[k][0]; jss = rlist[k][1]
                a = Pair()
                a.put_ssPair(iss, jss, name, 'p')
                vsPairList += [a]
            #|endfor
            
        else:
            for k in range(0, len(rlist)):
                iss = rlist[k][0]; jss = rlist[k][1]
                a = Pair()
                a.put_ssPair(iss, jss, name, stype)
                vsPairList += [a]
            #|endfor
            
        #
        
    #
    
    return vsPairList
#



