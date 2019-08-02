#!/usr/bin/env python

"""@@@

Main Module:     LoopRecords.py 

classes:         Branch
                 MBLptr

# Author:        Wayne Dawson
# creation date: 180709
# last update:   190701
# version:       0.0

Purpose:

Important processing and storage functions in building and analyzing
loops; used by several program in this package

"""

# These classes are used by in the dynamic programming algorithm part
# of Calculate to find optimal and suboptimal structures.

import sys

from Constants import INFINITY
from Motif     import Branch

class MBLptr(object):
    """
    Presently, this is used by the Calculate module (via RNAModules or
    ChromatinModules [or ProteinModules]) to compute I-loops and
    M-loops.
    """ 
    def __init__(self, i, j):
        self.i    = i
        self.j    = j
        self.nm   = 'X'      # name of the link  same as ctp
        self.dr   = "s"      # pair type direction [ap or pp] same as btp
        self.V    = INFINITY # free energy value of fragment
        self.Q    = []       # the branches
        self.n    = 0        # the number of branches
    #
    def resetBranch(self):
        self.nm   = 'X'      # name of the link
        self.dr   = "s"     # pair type direction [ap or pp]
        self.V    = INFINITY # free energy value of fragment
        self.Q    = []       # the branches
        self.n    = 0        # the number of branches
    #
    
    def addBranch(self, i, j):
        self.Q += [Branch(i,j)]
        self.n = len(self.Q)
    #
    
    def pushBranch(self, br):
        self.Q += [br]
        self.n = len(self.Q)
    #
    
    def getlen(self):
        return len(self.Q)
    #
    
    def setBranch(self, k, i, j, nm = "branch"):
        flag_pass = True
        if len(self.Q) >= k:
            print "ERROR: no index %d in this branch (%s)" % (k, nm)
            sys.exit(1)
        #
        self,Q[k].i = i
        self,Q[k].j = j
        # should actully check if the set branch makes sense, but
        # presently, it is not that bright.
    #
    
    def getBranch(self, k):
        return self.Q[k]
    #
    
    def getBranchlist(self):
        branches = []
        for k in range(0, len(self.Q)):
            branches += [(self.Q[k].i, self.Q[k].j)]
        #
        return branches
    #
    
    def getBranchEnds(self):
        ib = self.i; jb = self.j
        if len(self.Q) > 0:
            ib = self.Q[0].i
            jb = self.Q[self.n-1].j
        #
        return ib, jb    
    #
    
    def getBoundaries(self):
        return self.i, self.j
    #
    
    def pushBranchlist(self, bl, flag_reset = False):
        if flag_reset:
            self.Qbranches = []
        #
        for k in range(0, len(bl)):
            self.Q += [Branch(bl[k][0], bl[k][1])]
        #
        
    #
    
    def __str__(self):
        s = "i,j   (%d,%d)\n" % (self.i, self.j)
        s += "nm    %s\n" % self.nm
        s += "n     %d\n" % self.n
        s += "V     %10.2f\n" % self.V
        s += "Q:    "
        for k in range(0, len(self.Q)):
            s += "(%d,%d)  " % (self.Q[k].i, self.Q[k].j)
        #
        s += "\n"
        return s
    #
    
    def __repr__(self):
        print self.__str__()
    #
#


class MBLHandle(object):
    """
    Handling object for finding best motif at (i,j). Must be set up
    and deleted with each entry and exit of searchForMBL().
    """
    def __init__(self,
                 gid,                # motif group identifier
                 i, j,               # position of MBL
                 dG_lb = -INFINITY): # dG lower bound
        
        # 190103: This class should be built up to handle degeneracy,
        # but presently, I am just using what I have from vsfold5.
        
        # input variables
        self.gid    = gid
        self.i      = i
        self.j      = j
        self.dG_lb  = dG_lb #  dG lower bound
        
        # information to be extracted
        self.debug_MBLHandle = True
        self.best_iMBL_dG  = INFINITY
        self.mbls          = []  # class MBLptr()
        self.Mloop_wt      = 0.0
        
        # note self.fe.btype[i][j].pair defines whether the location
        # has a true binding point.
    #
    
    # ~MBLHandle(): Destructor would have to dispose of MBLptr without
    # destroying the _vital_ data in fe and smap
    
    
    
    def update_mbl(self, new_MBL, p, q, debug, lbl):
        new_iMBL_dG = new_MBL.V
        self.debug_MBLHandle = debug
        
        
        if new_iMBL_dG < self.best_iMBL_dG:
            # print "check iMBL FE";
            if self.dG_lb < new_iMBL_dG:  # getFE_LowerBound()
                
                if self.debug_MBLHandle:
                    print "gid(%d): dG_lb(%8.2f) < new_iMBL_dG(%8.2f) < best_iMBL_dG(%8.2f)" \
                        % (self.gid, self.dG_lb, new_iMBL_dG, self.best_iMBL_dG)
                #
                
                # record the positions and the energy in the M-loop branches. 
                
                # print "new_iMBL_dG = %8.2f" % new_iMBL_dG
                
                # print "check iMBL FE"
                
                if new_iMBL_dG < self.best_iMBL_dG:
                    self.mbls = []
                    # reset the whole thing, we've found a new and
                    # better minimum free energy candidate over
                    # (i,j). This means calling destructors to remove
                    # all the contents in C++, because this is a
                    # complex object with different types of data.
                #
                self.best_iMBL_dG = new_iMBL_dG;
                self.mbls += [MBLptr(self.i,self.j)]
                # degenerate or not, we now initialize mbls as a vector
                self.mbls[0].resetBranch()
                
                self.mbls[0].n  = new_MBL.n
                self.mbls[0].V  = new_iMBL_dG
                self.mbls[0].nm = new_MBL.nm
                self.mbls[0].dr = new_MBL.dr
                for l in range(0, new_MBL.n):
                    self.mbls[0].pushBranch(new_MBL.getBranch(l))
                #
                
                # print "passed check point";
                
                if self.debug_MBLHandle:  
                    print "searchForMBL.update_mbl(%s):" % lbl
                    print "      M-loop(%2d,%2d): new_iMBL_dG = %8.2f" \
                        % (self.i, self.j, self.mbls[0].V)
	            print "mbls     name     index        (  p,  q)   "
                    for kv in range(0, len(self.mbls)):
	                for l in range(0, self.mbls[kv].n):
                            pv = self.mbls[kv].Q[l].i; qv = self.mbls[kv].Q[l].j
                            ctp = self.mbls[kv].nm
                            btp = self.mbls[kv].dr
	                    print " %2d    %s[%5s]    %2d          (%3d,%3d) "  \
                                % (kv, ctp, btp, l + 1, pv, qv)
	                #
                    #
	            print "Mloop_wt = %8.2f" % self.Mloop_wt
	            # sys.exit(0)
                #endif
            #
            else:
                if self.debug_MBLHandle:  
                    s  = "searchForMBL.update_mbl(%s): " % lbl
                    s += "getFE_LowerBound(%8.2f) >= new_iMBL_dG(%8.2f) ?? ==> " \
                         % (self.dG_lb, new_iMBL_dG)
                    s += "lower bound exceeded, skipped recording this info"
                    print s
                #endif      
            #
        #
        else:
            if self.debug_MBLHandle:  
                s  = "searchForMBL.update_mbl(%s): \n" % lbl
                s += "V[i(%d),p(%d)],V[q(%d),j(%d)]: " % (self.i, p, q, self.j)
                s += "new_iMBL_dG[%d](%8.2f) < best_iMBL_dG[%d](%8.2f) ?? ==> " \
                     % (self.gid, new_iMBL_dG, self.gid, self.best_iMBL_dG)
                s += "skipped recording this info\n"
                print s
            #endif      
        #
    #
    
#

def show_MBLHandle(mblh, smap):
    # summarize the results of the search 
    
    # 180413(indeed, Friday the 13th): I am changing this function to
    # find both I-loops and M-loops. I think this strategy will make
    # this subprogram an essential tool for doing the speedup. The
    # other part is simplifying the I-loop technology
    
    
    print "  gid  M-id  branch  (  p,  q)    type             FE        link"
    pv = qv = -1
    i = mblh.i; j = mblh.j
    gid = mblh.gid
    
    for km in range(0, len(mblh.mbls)):
        ctpx = mblh.mbls[km].nm
        btpx = mblh.mbls[km].dr
        Vijx = mblh.mbls[km].V
        mjnx = mblh.mbls[km].getBranchlist()
        print "  %2d   %2d   *%2d*     (%3d,%3d)     %s[%5s]   %8.2f   " \
            % (gid,  0, 0, i, j, ctpx, btpx, Vijx), mjnx
        
        
        for kv in range(0, mblh.mbls[km].n):
            pv = mblh.mbls[km].Q[kv].i; qv = mblh.mbls[km].Q[kv].j
            ctpv = smap.glink[pv][qv].lg[0].motif[0].get_ctp()
            btpv = smap.glink[pv][qv].lg[0].motif[0].get_btp()
            Vijv = smap.glink[pv][qv].lg[0].motif[0].get_Vij()
            mjnv = smap.glink[pv][qv].lg[0].motif[0].get_branches()
            # print mjn
            
            print "             %2d      (%3d,%3d)     %s[%5s]   %8.2f   " \
                % (kv+1, pv, qv, ctpv, btpv, Vijv), mjnv
            #
        #
    #
    print "Mloop_wt = %8.2f" %  mblh.Mloop_wt
#




def main(cl):
    # tests the operations of Branch and MBLptr
    print cl
    b1 = Branch(1,2)
    b2 = Branch(3,4)
    print "branch1: ", b1
    print "branch2: ", b2
    a = MBLptr(0, 10)
    a.pushBranch(b1)
    a.pushBranch(b2)
    print "branch3 add-> ", (5,6)
    a.addBranch(5,6)
    print a
    
    
#


if __name__ == '__main__':
    main(sys.argv)
