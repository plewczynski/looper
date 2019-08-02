#!/usr/bin/env python

"""@@@

Main Module:   Vienna.py 

Classes:       Vstruct
               ViennaData

Functions:

Author:        Wayne Dawson
creation date: ~2014/2015
last update:   190402
version:       0.1

Purpose:

This tool reads the Vienna formatted pairing sequence and converts it
into a set of pair-contacts, or in the case of chromatin, both pair
contacts and ctcf contacts.

Originally, this was intended as a tool for handling Vienna package
data with SimRNA. Now, it is used by some other programs, particularly
chreval.py and its derivatives. It also has been significantly
upgraded to handle and represent a very diverse variety of structures
containing pseudoknot like structures and also multipair interactions.

The two functions that carry this action out are

parse_SecondaryStructure       -- only dot bracket structures .((...)).
parse_fullDotBracketStructure  -- all types of complex dot bracket structures


This package can be tested in the following way

vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
Python 2.7.6 (default, Mar 22 2014, 22:59:56) 
[GCC 4.8.2] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> from Vienna import Vstruct
>>> vs = Vstruct()
>>> vs.parse_SecondaryStructure("(((((....)))))")
(((((....)))))
(    0,   13)
(    1,   12)
(    2,   11)
(    3,   10)
(    4,    9)
0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

"""

import sys
# used for PKs and parallel stems
# this notation at least works with VARNA
from Constants import num2lpr
from Constants import num2rpr
from Constants import lpr2num
from Constants import rpr2num
from Constants import PKfndx
from Constants import PKrndx
from Constants import pointer
from Constants import stack
from Constants import counter
from Constants import Xlist

from BasicTools import sortPairListWRT_n
from BasicTools import copyList
from BasicTools import tuple2List


# fundamental representation of pairing
from Pair import Pair
from Pair import SortPair
from Pair import vsPair2list

# other constants and parameters

from Constants import sysDefLabels # system RNA, Chromatin

debug = False # True # 

"""
tests:
  1 = test parse_SecondaryStructure for a simple
      secondary structure

  2 = test parse_fullDotBracketStructure for more
      complex structures that include pseudoknots. 
"""

TEST = 2 



# error handling 
class MyException(Exception):
    def __init__(self, st):
        self.stmnt = st
    #
    def __str__(self):
        return self.stmnt
    #

    """ @
    Example:
      if not self.set_var:
          mssg = "\nERROR, value of var undefined\n"
          ex=MyException(mssg)
          raise ex
      #
    """
#


def LThread2Pair(lt, title = "undefined"):
    """coverts a list of LThread data to a list of Pair data"""
    
    BPlist = []
    for tr in lt.thread:
        btp = tr.btp    #  s, sp, sa, wyspa, bgn, end, etc.
        if not (btp == 'bgn' or btp == 'end' or btp == 'wyspa'):
            v = tr.ij_ndx
            i = v[0]; j = v[1]
            
            bp = Pair()
            # Vienna is specifically designed for secondary structure
            # of single strand (ss) pairs. Therefore, we use put_ssPair
            bp.put_ssPair(i, j, 'bp', btp)
            BPlist += [bp]
        #
    #
    ordering = SortPair()
    # sorting is not necessarily required but it makes it easier to
    # read at some level.
    BPlist = ordering.sortvsList(BPlist)
    
    return BPlist
#



class Vstruct(SortPair):
    def __init__(self, system = "RNA"):
        self.system  = system
        self.vstr    = '' # structure
        self.vseq    = '' # sequence (optional) 
        self.N       = -1  # sequence length
        self.wt      = -1 # used with my_generation to add specified weights to heatmaps
        # standard secondary structure (ss)
        self.BPlist  = []
        self.BProots = []
        # a both ss and pseudoknot (PK) in same line
        self.PKlist  = []
        self.PKroots = []
        # CTCF islands, or triple helices
        self.MPlist  = []
        
        
        """@
        
        170306wkd: I know the following looks a little strange because
        I import Xlist (etc.) from Constants.py; however, for some
        reason that I cannot explain, the module Xlist is changed when
        these objects are changed, even though self.Xlist should be a
        independent object and Xlist should be just a
        template. Nevertheless, if I use
        
          self.Xlist.update({xl : Xlist[xl]}),
          self.stack.update({sl : stack[sl]}),
          etc.
        
        The module values are updated with previous entries generated
        from prior created Vstruct objects.  To make a long story
        short, finally what I had to do to make sure that self.Xlist
        didn't end up with misinformation was to hard wire write it
        independent of the contents of the module Xlist. So, in the
        end, the keys from the module are used, but not the
        contents. Strange, but anyway.....
        
        """
        
        # reset the key reference dictionary of pointers 
        self.Xlist = {}
        for xl in Xlist.keys():
            self.Xlist.update({xl : [] })
        #
        self.stack = {}
        for sl in stack.keys():
            self.stack.update({sl : 0 })
        #
        self.counter = {}
        for cl in counter.keys():
            self.counter.update({cl : 0 })
        #
        self.pointer = {}
        for pl in pointer.keys():
            self.pointer.update({pl : 0 })
        #
    #endInit
    
    def set_system(self, sname):
        if sysDefLabels.has_key(sname):
            self.system = sname
        else:
            print "ERROR: unrecognized system name (%s)" % sname
            print "       allowed names "
            for snm in sysDefLabels.keys():
                print snm
            #
            sys.exit(1)
        #
    #
    
    def init_vs_from_lt(self, vf):
        debug_init_vs_from_lt = False
        if debug_init_vs_from_lt:
            print "Entered init_vs_from_lt"
        #
        object_name = type(vf).__name__
        if not object_name == "LThread2Vienna":
            print "ERROR(init_vs_from_lt) requires type LThread2Vienna"
            sys.exit(1)
        #
        
        self.reset_Vstruct()
        self.BPlist    = vf.vsBPlist
        self.BProots   = vf.vsBProots 
        self.PKlist    = vf.vsPKlist   # PKs
        self.PKroots   = vf.vsPKroots 
        self.MPlist    = vf.vsMPlist   # e.g., ctcf etc
        self.N         = vf.N         
        self.vstr      = vf.vstr      
        self.vseq      = vf.vseq
        
        
        if debug_init_vs_from_lt:
            self.disp_allLists("from init_vs_from_lt (before corrections)")
        #
        # post base pair list processing
        self.PKlist = self.sortvsList(self.PKlist, 'i')
        self.PKlist = self.assign_PKdirection(self.PKlist)
        self.compress_to_BPlist()
        if debug_init_vs_from_lt:
            self.disp_allLists("from init_vs_from_lt (after adjustments)")
            #print "stop at end of init_vs_from_lt"; sys.exit(0)
        #
        
        """@@@
        
        It turns out, with the structure "(.).(.([)])", the resulting
        BPlist and PKlist from chreval.py is somewhat different from
        when the list is read from the sequence.
        
        ERROR(buildTreeDmns near 10): somehow, core PK is not properly defined
        iz( 0)   < pLmin( 7) < jz( 2)   < qLmax( 9)
        <or>   pLmin( 7) < iz( 0)   < qLmax( 9) < jz( 2)
        buildTreeDmns near 10
        structural sequence
        
        secondary structure
        (    0,    2)[a]
        (    4,   10)[a]
        (    6,    8)[a]   <=== 
        (    7,    9)[a]   <=== these should both be [p]
        secondary structure roots
        (    0,    2)[a]
        (    4,   10)[a]
        
        The output from vs is 
        
        structure sequence: 
        (.).(.([)])
        vsBPlist: 
        (    0,    2)[a]
        (    4,   10)[a]
        (    6,    8)[p]  
        (    6,    8)[p]  <=== evidently double counts but at least both are marked [p]
        vsPKlist: 
        (    7,    9)[p]
        vsMPlist: 
        xcycxcxxyyy
        (.).(.[<]>)
        
        .... So it seems that I will have to run an independent check of
        the structure through Threads or here.
        
        """
    #endMethod
    
    
    def set_Vstruct(self, s):
        self.vstr = s
        self.N = len(s)
        return 0
    #endMethod
    
    def set_Vseq(self, s):
        self.vseq = s
        return 0
    #endMethod
    
    def reset_Vstruct(self):
        self.vstr = '' # input dot-bracket secondary structure [+ PK + CTCF]
        self.vseq = '' # not required, but sometimes useful
        self.N = -1
        # standard secondary structure (ss)
        self.BPlist  = []
        self.BProots = []
        # a both ss and pseudoknot (PK) in same line
        self.PKlist  = []
        self.PKroots = []
        
        # CTCF islands
        self.MPlist = []

        """@
        
        170306wkd: I know the following looks a little strange because
        I import Xlist (etc.) from Constants.py; however, for some
        reason that I cannot explain, the module Xlist is changed when
        these objects are changed, even though self.Xlist should be an
        independent object and Xlist should be just a
        template. Nevertheless, if I use
        
          self.Xlist.update({xl : Xlist[xl]}),
          self.stack.update({sl : stack[sl]}),
          etc.
        
        The module values are updated with previous entries generated
        from prior created Vstruct objects.  To make a long story
        short, finally what I had to do to make sure that self.Xlist
        didn't end up with misinformation was to hard wire write it
        independent of the contents of the module Xlist. So, in the
        end, the keys from the module are used, but not the
        contents. Strange, but anyway.....

        """
        
        # reset the key reference dictionary of pointers 
        self.Xlist = {}
        for xl in Xlist.keys():
            self.Xlist.update({xl : []})
        #
        self.stack = {}
        for sl in stack.keys():
            self.stack.update({sl : 0})
        #
        self.counter = {}
        for cl in counter.keys():
            self.counter.update({cl : 0})
        #
        self.pointer = {}
        for pl in pointer.keys():
            self.pointer.update({pl : 0})
        #
        return 0
    #endMethod
    
    
    def scan_vs(self, i):
        if debug:
            print "scan_vs(%d):" % i
        #
        mm = -1
        if i < self.N:
            s = self.vstr[i]
            # separate '(' and ')' from '[', ']', and '.'
            if not (s == '(' or s == ')' or s == '[' or s == ']' or s == '.'):
                print 'ERROR: improperly defined Fontana formated sequence'
                print '       offending character \'%c\' located at position %d' % (s, i+1)
                sys.exit(1)
            #
            if s == '(':
                # add to the stack
                b = Pair()
                b.put_ssPair_i(i, "bp")
                self.BPlist += [b]
                if debug:
                    print b.disp_Pair()
                    # print 'i = %5d, %s' % (i, s)
                #
                self.stack[0] += 1 
                self.pointer[0] = stack.stack[0]
                self.counter[0] +=1
                mm = i + 1
                self.scan_vs(mm)
            elif s == ')':
                self.pointer[0] = self.find_next_Xpoint_n(self.BPlist, self.pointer[0])
                self.BPlist[self.pointer[0]-1].put_ssPair_j(i, "bp")
                if debug:
                    print self.BPlist[self.pointer[0]-1].disp_Pair()
                    # print 'j = %5d, %s' % (i, s)
                #
                
                self.pointer[0] -= 1
                self.counter[0] -= 1
                mm = i + 1
                self.scan_vs(mm)
            else:
                mm = i
                # ignore pseudoknot brackets [[[...]]], if present.
                # Converts these to unpaired strand data.
                if (s == '[') or (s == ']'):
                    s = '.'
                #
                
                # Either look for next closing bracket in sequence or
                # terminate at the end of the sequence if nothing is
                # found.
                while (s == '.') and (mm < self.N):
                    if debug:
                        print 'i = %5d, %s' % (mm, s)
                    #
                    mm += 1
                    if mm == self.N:
                        break
                    #
                    s = self.vstr[mm]
                #   
                self.scan_vs(mm)
            #
        
        else:
            if debug:
                print 'point = %d' % self.pointer[0]
            #
            if not self.counter[0] == 0:
                case = 0
                if self.counter[0] > 0:
                    case = 0
                else:
                    case = 1
                #
                
                print 'jcount = %d' %  self.counter[0]
                print 'ERROR!!! Fontana notation is not correct!'
                if case == 0:
                    print '         %d too many \'(\' brackets!' % self.counter[0]
                else:
                    print '         %d too many \')\' brackets!' % (-self.counter[0])
                #
                
                print 'input structure:'
                print self.vstr
                sys.exit(1)
            #
            if not self.check_Xlist_n(self.BPlist):
                print 'ERROR!!! Fontana notation is not correct!'
                print '         At least one structure of type \')...(\' was found'
                print 'input structure:'
                print self.vstr
                sys.exit(1)
            #
        #
        
        return mm
    #endMethod
    
    
    def scan_allTypes(self, i, layer):
        debug_bp   = False # True # 
        debug_PK   = False # True # 
        debug_ctcf = False # True # 
        if debug_bp or debug_PK or debug_ctcf:
            print "enter scan_allTypes(%d):" % i, self.pointer[0]
        #
        if layer > self.N:
            # if it does this, something is definitely wrong! The
            # variable layer is mainly used to track the recursion
            # level and make sure that things have not gone weird. It
            # seems like this part of the program works fine now. I
            # have not encoutered a stop in ages. However, the
            # debugging code is still here because I don't know when
            # it might be needed.
            print "ERROR(scan_allTypes): layer(%d) exceeds the sequence length(%d)!" \
                % (layer, self.N)
            print "                      Something is wrong. ... Last call i = %d"   % (i)
            sys.exit(1)
        #
        
        mm = -1
        if i < self.N:
            s = self.vstr[i]
            if not (lpr2num.has_key(s) or rpr2num.has_key(s) or s == '|' or s == '.'):
                print 'ERROR(scan_allTypes): improperly defined Fontana formated sequence'
                print '                      offending character \'%c\' located at' % (s)
                print 'position %d' % (i+1)
                sys.exit(1)
            #
            
            if s == '(': # standard bp argument 
                # add to the stack
                BPndx = lpr2num[s]
                if debug_bp:
                    print "BPndx(%2d) => %s" % (BPndx, s)
                #
                b = Pair()
                b.put_ssPair_i(i, "bp")
                self.Xlist[BPndx] += [b]
                if debug_bp:
                    print b.disp_Pair()
                    print 'i = %5d, %s' % (i, s)
                    print "element(%2d) => counter(%2d), pointer(%s), stack(%2d)" \
                        % (BPndx, self.counter[0], self.pointer[0], self.stack[0])
                #
                
                mm = i + 1
                self.stack[BPndx]   += 1
                self.pointer[BPndx]  = self.stack[BPndx]
                self.counter[BPndx] += 1
                self.scan_allTypes(mm, layer + 1)
                
            elif s == ')':
                BPndx = rpr2num[s]
                if debug_bp:
                    print "BPndx(%2d) => %s" % (BPndx, s)
                    flag_details = False
                    if flag_details:
                        print "list of current contacts"
                        for bp in self.Xlist[BPndx]:
                            print bp.disp_Pair()
                        #
                        print "------"
                    #
                #
                
                self.pointer[BPndx] \
                    = self.find_next_Xpoint_n(self.Xlist[BPndx], self.pointer[BPndx])
                self.Xlist[BPndx][self.pointer[BPndx]-1].put_ssPair_j(i, "bp")
                if debug_bp:
                    print self.Xlist[BPndx][self.pointer[BPndx]-1].disp_Pair()
                    print 'j = %5d, %s' % (i, s)
                    print "element(%2d) => counter(%2d), pointer(%s), stack(%2d)" \
                        % (BPndx, self.counter[0], self.pointer[0], self.stack[0])
                #
                
                self.pointer[BPndx] -= 1
                self.counter[BPndx] -= 1
                mm = i + 1
                self.scan_allTypes(mm, layer + 1)
                
            elif s == '[': # PK
                PKndx = lpr2num[s]
                if debug_PK:
                    print "PKndx(%2d) => %s" % (PKndx, s)
                #
                
                mm = i
                b = Pair()
                b.put_ssPair_i(i, "pk")
                self.Xlist[PKndx] += [b]
                if debug_PK:
                    print b.disp_Pair()
                    print 'i = %5d, %s' % (i, s)
                #
                
                self.stack[PKndx]   += 1 
                self.pointer[PKndx]  = self.stack[PKndx]
                self.counter[PKndx] +=1
                mm = i + 1
                self.scan_allTypes(mm, layer + 1)
                
            elif s == ']':
                PKndx = rpr2num[s]
                if debug_PK:
                    print "PKndx(%2d) => %s" % (PKndx, s)
                #
                
                self.pointer[PKndx] \
                    = self.find_next_Xpoint_n(self.Xlist[PKndx], self.pointer[PKndx], "PK")
                self.Xlist[PKndx][self.pointer[PKndx]-1].put_ssPair_j(i, "pk")
                if debug_PK:
                    print self.Xlist[PKndx][self.pointer[PKndx]-1].disp_Pair()
                    print 'j = %5d, %s' % (i, s)
                #
                
                self.pointer[PKndx] -= 1
                self.counter[PKndx] -= 1
                mm = i + 1
                self.scan_allTypes(mm, layer + 1)
                
                
                
                # ########################################
                # CTCF
                # ########################################
                # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
            elif s == '{':  # xxxx  CTCF  xxxx
                CTCFndx = lpr2num[s]
                if debug_ctcf:
                    print "CTCFndx(%2d) => %s" % (CTCFndx, s)
                #
                
                mm = i
                b = Pair()
                b.put_ssPair_i(i, "ctcf")
                self.Xlist[CTCFndx] += [b]
                if debug_ctcf:
                    print b.disp_Pair()
                    print 'i = %5d, %s' % (i, s)
                #
                
                self.stack[CTCFndx]   += 1 
                self.pointer[CTCFndx]  = self.stack[CTCFndx]
                self.counter[CTCFndx] +=1
                mm = i + 1
                self.scan_allTypes(mm, layer + 1)
                
            elif s == '}':
                CTCFndx = rpr2num[s]
                if debug_ctcf:
                    print "CTCFndx(%2d) => %s" % (CTCFndx, s)
                #
                
                self.pointer[CTCFndx] = self.find_next_Xpoint_n(self.Xlist[CTCFndx],
                                                                self.pointer[CTCFndx], "CTCF")
                self.Xlist[CTCFndx][self.pointer[CTCFndx]-1].put_ssPair_j(i, "ctcf")
                if debug_ctcf:
                    print self.Xlist[CTCFndx][self.pointer[CTCFndx]-1].disp_Pair()
                    print 'j = %5d, %s' % (i, s)
                #
                
                
                self.pointer[CTCFndx] -= 1
                self.counter[CTCFndx] -= 1
                mm = i + 1
                self.scan_allTypes(mm, layer + 1)
                
            elif s == '|':
                CTCFndx = lpr2num['{']
                if debug_ctcf:
                    print "CTCFndx(%2d) => |" % (CTCFndx)
                #
                
                self.Xlist[CTCFndx][self.pointer[CTCFndx]-1].put_contacts(i)
                mm = i + 1
                if debug_ctcf:
                    print 'j = %5d, %s' % (i, s)
                #
                
                self.scan_allTypes(mm, layer + 1)
                # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                # ########################################
                # CTCF
                # ########################################
                
                
                
                # ########################################
                # expanded PK notations
                # ########################################
                # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
            elif lpr2num.has_key(s):
                l_mate = s
                r_mate = num2rpr[lpr2num[s]]
                PKndx = lpr2num[s]
                if debug_PK:
                    print "PKndx(%2d) => %s, search for pattern %s and %s" \
                        % (PKndx, s, l_mate, r_mate)
                    
                    print "left(i = %d) ==> %s" % (i, s)
                    print "lpr2num[s]:                          ", lpr2num[s]
                    print "num2lpr[lpr2num[s]]:                 ", num2lpr[lpr2num[s]]
                    print "PKrndx[lpr2num[s]]:                  ", PKrndx[lpr2num[s]]
                    print "num2lpr[PKfndx[PKrndx[lpr2num[s]]]]: ", \
                        num2lpr[PKfndx[PKrndx[lpr2num[s]]]]
                    print "search for connecting element:          '%s'" \
                        % num2rpr[lpr2num[s]]
                    print "element(%2d) => counter(%2d), pointer(%s), stack(%2d)" \
                        % (PKndx, self.counter[PKndx],
                           self.pointer[PKndx], self.stack[PKndx]) 
                #
                
                mm = i
                b = Pair()
                b.put_ssPair_i(i, "pk")
                self.Xlist[PKndx] += [b]
                if debug_PK:
                    print b.disp_Pair()
                    print 'i = %5d, %s' % (i, s)
                #
                
                self.stack[PKndx]   += 1 
                self.pointer[PKndx]  = self.stack[PKndx]
                self.counter[PKndx] += 1
                mm = i + 1
                self.scan_allTypes(mm, layer + 1)
                
            elif rpr2num.has_key(s):
                l_mate = num2lpr[rpr2num[s]]
                r_mate = s
                PKndx = rpr2num[s]
                ndx = PKrndx[rpr2num[s]]
                if debug_PK:
                    print "right(i = %d) ==> %s" % (i, s)
                    print "rpr2num[s]:                          ", rpr2num[s]
                    print "num2rpr[rpr2num[s]]:                 ", num2rpr[rpr2num[s]]
                    print "PKrndx[rpr2num[s]]:                  ", PKrndx[rpr2num[s]]
                    print "num2rpr[PKfndx[PKrndx[rpr2num[s]]]]: ", \
                        num2rpr[PKfndx[PKrndx[rpr2num[s]]]]
                    print "search for connecting element:          '%s'" \
                        % num2rpr[rpr2num[s]]
                    print "element(%2d) => counter(%2d), pointer(%s), stack(%2d)" \
                        % (PKndx, self.counter[PKndx], \
                           self.pointer[PKndx], self.stack[PKndx]) 
                #
                self.pointer[PKndx] \
                    = self.find_next_Xpoint_n(self.Xlist[PKndx],
                                              self.pointer[PKndx], "PK", False)
                self.Xlist[PKndx][self.pointer[PKndx]-1].put_ssPair_j(i, "pk")
                
                if debug_PK:
                    print self.Xlist[PKndx][self.pointer[PKndx]-1].disp_Pair()
                    print 'j = %5d, %s' % (i, s)
                #
                
                self.pointer[PKndx] -= 1
                self.counter[PKndx] -= 1
                mm = i + 1
                self.scan_allTypes(mm, layer + 1)
                #
                # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                # ########################################
                # expanded PK notations
                # ########################################
            else:
                # exit out the non-interacting point
                mm = i
                while (s == '.') and (mm < self.N):
                    if debug_bp or debug_PK or debug_ctcf:
                        print 'i = %5d, %s' % (mm, s)
                    #
                    
                    mm += 1
                    if mm == self.N:
                        break
                    #
                    
                    s = self.vstr[mm]
                #
                
                self.scan_allTypes(mm, layer + 1)
                
            #
        
        else:
            
            # ######################################################
            # ######  When reach this point, means we are finished
            # ######  searching the sequence
            # ######  ##############################################

            #print "came here!!!"
            
            # secondary structure contacts 
            self.BPlist = self.Xlist[0]
            self.BPlist = self.sortvsList(self.BPlist, 'i')
            # self.BPlist = sorted(self.Xlist[0], key=self.getKey_i)
            self.BProots = self.findroots(self.BPlist, False)
            if debug_bp:
                print "BPlist: "
                for nn in self.BPlist:
                    print nn.disp_Pair()
                #
                print "BProots: "
                for nn in self.BProots:
                    print nn.disp_Pair()
                #
                #print "stop at 0 in scan_allTypes"; sys.exit(0)
            #
            
            # post base pair list processing
            
            # pseudoknot contacts 
            self.PKlist = []
            pkndx = PKrndx.keys()
            for x_k in pkndx:
                self.PKlist += self.Xlist[x_k]
            #
            self.PKlist = self.assign_PKdirection(self.PKlist)
            self.PKlist = self.sortvsList(self.PKlist, 'i')
            self.compress_to_BPlist()
            if debug_bp:
                print "PKlist: "
                for nn in self.PKlist:
                    print nn.disp_Pair()
                #
                print "PKroots: "
                for nn in self.PKroots:
                    print nn.disp_Pair()
                #
                # print "stop at 1 in scan_allTypes (after compress_to_BPlist)"; sys.exit(0)
            #
            
            
            # MP contacts; e.g., CTCFs
            self.MPlist = self.Xlist[2]
            for cl in range(0, len(self.MPlist)):
                self.MPlist[cl].v = '-'
            #
            
            # display the layout of results
            if debug_bp or debug_PK or debug_ctcf:
                print "index   pointer    counter       stack"
                for x_k in self.Xlist.keys():
                    print " %2d      %3d        %3d        %3d" \
                        % (x_k, self.pointer[x_k], self.counter[x_k], self.stack[x_k])
                #
            #
            for x_k in self.Xlist.keys():
                if not self.counter[x_k] == 0:
                    case = 0
                    if self.counter[x_k] > 0:
                        case = 0
                    else:
                        case = 1
                    #
                    print 'counter[%d] = %d' %  (x_k, self.counter[x_k])
                    print 'ERROR!!! Fontana notation is not correct!'
                    if case == 0:
                        print '         %d too many \'%s\' brackets!' \
                            % (self.counter[x_k], num2lpr[x_k], num2rpr[x_k])
                    else:
                        print '         %d too many \'%s\' brackets!' \
                            % (-self.counter[x_k], num2lpr[x_k], num2rpr[x_k])
                    #
                    print 'input structure:'
                    print self.vstr
                    sys.exit(1)
                #
            #
            for x_k in self.Xlist.keys():
                if not self.check_Xlist_n(self.Xlist[x_k]):
                    print "ERROR!!! Fontana notation is not correct!"
                    print "         At least one structure of type \'%s...%s\' was found" \
                        % (num2rpr[x_k], num2lpr[x_k])
                    print "         i.e., the order should be %s...%s." \
                        % (num2lpr[x_k], num2rpr[x_k])
                    print 'input structure:'
                    print self.vstr
                    sys.exit(1)
                #
            #
            if debug_bp or debug_PK or debug_ctcf:
                print "mostly finished scan_allTypes(layer=%2d, ndx=%4d, N=%4d" \
                    % (layer, i, self.N)
                # print "stop at 2 in scan_allTypes"; sys.exit(0)
            #
        #
        
        if debug_bp or debug_PK or debug_ctcf:
            print "exit scan_allTypes(%d):" % i, self.pointer[0]
            # print "stop at 3 (exit from) scan_allTypes"; sys.exit(0)
        #
        
        return mm
    #
    
    
    
    def assign_PKdirection(self, plist): # plist <- self.PKlist
        """#########
        !!!!!
        
        Decides if the stem is parallel or antiparallel. Here are a
        couple of examples of tests that I ran from Threads.py.
        
        example 1
         0         10        20        30        40 
         |    .    |    .    |    .    |    .    |  
        ss_seq  = ".A..B..C..D...........a..b..c..d.."
        result:
        (   1,   22)[p]
        (   4,   25)[p]
        (   7,   28)[p]
        (  10,   31)[p]
        
        example 2:
         0         10        20        30        40 
         |    .    |    .    |    .    |    .    |  
        ss_seq  = ".AAA...BBB...aaa...CCC....bbb....ccc...."
        result:
        (   1,   22)[a]
        (   2,   21)[a]
        (   3,   20)[a]
        (   7,   28)[a]
        (   8,   27)[a]
        (   9,   26)[a]
        (  13,   35)[a]
        (  14,   34)[a]
        (  15,   33)[a]
        
        example 3:
         0         10        20        30        40 
         |    .    |    .    |    .    |    .    |  
        ss_seq  = ".AAA...BBB...aaa...CCC....bbb....ccc...." # same
        ss_seq  = ".(((...BBB...)))...(((....bbb....)))...." # same
        result:
        (   1,   15)[a]
        (   2,   14)[a]
        (   3,   13)[a]
        (   7,   28)[a]
        (   8,   27)[a]
        (   9,   26)[a]
        (  19,   35)[a]
        (  20,   34)[a]
        (  21,   33)[a]
        
        It appears to be able to distinguish consequitive pairs that
        would satisfy a parallel pattern and they do not have to be
        contiguous either (as shown in example 1). However, if it is
        just an odd sort of pseudoknot where you have antiparallel
        stems but a kind of ladder structure, it will not call this
        ladder arrangement "parallel" (assigning the "p" to the
        brackets [a] or [p].
        
        190219: appears to work ok with the more odd or difficult
        examples I could think of.
        
        """
        debug_ad = False # True # 
        if debug_ad:
            print "Enter assign_PKdirection: "
            print "plist: (before)"
            for vv in plist:
                print vv.disp_Pair()
            #
        #
        
        
        # since the lists are ordered, we can check through
        # self.BPlist and identify parallel stems by the simple rule
        # of parallel-ness: ibp1 < ibp2 < jbp1 < jbp2.

        for k in range(0, len(self.BPlist)-1):
            ibp1 = self.BPlist[k  ].i; jbp1 = self.BPlist[k  ].j
            ibp2 = self.BPlist[k+1].i; jbp2 = self.BPlist[k+1].j
            if debug_ad:
                print "ibp1(%2d) < ibp2(%2d) < jbp1(%2d) < jbp2(%2d)" \
                    % (ibp1, ibp2, jbp1, jbp2)
            #
            if ibp1 < ibp2 and ibp2 < jbp1 and jbp1 < jbp2:
                # parallel relationship. This is the nearest
                # neighboring Pair and it is in BPlist (the dominant
                # one). Therefore, it should have priority and be
                # assigned properly.
                self.BPlist[k  ].v = 'p'
                self.BPlist[k+1].v = 'p'
            #
        #
        
        if debug_ad:
            print "BPlist: (after checking)"
            for vv in self.BPlist:
                print vv.disp_Pair()
            #
            print "PKlist: (before checking)"
            for vv in self.PKlist:
                print vv.disp_Pair()
            #
            
            # print "stop at 0.1 in assign_PKdirection", sys.exit(0)
        #
        
        
        #print self.PKlist
        
        # pre-search for parallel stems before 
        dd = {}
        n = len(self.BPlist) - 1
        kr = 0
        bptemp = []
        while kr < len(self.BPlist):
            ibp  = self.BPlist[kr  ].i; jbp  = self.BPlist[kr  ].j
            if debug_ad:
                print "ibp(%2d,%2d), kr = %d" % (ibp, jbp, kr)
            #
            
            for kt in range(0, len(self.PKlist)):
                ipk = self.PKlist[kt].i; jpk = self.PKlist[kt].j
                if debug_ad:
                    print "ibp(%2d) < ipk(%2d) < jbp(%2d) < jpk(%2d)" \
                        % (ibp, ipk, jbp, jpk)
                #
                if ibp < ipk and jbp < jpk and ipk < jbp:
                    if dd.has_key((ibp, jbp)):
                        dd[(ibp, jbp)] += [(ipk, jpk)]
                    else:
                        dd.update({(ibp, jbp) : [(ipk, jpk)]})
                        bptemp += [(ibp, jbp)]
                        
                    #
                elif jbp < ipk:
                    break
            #
            kr += 1
        #
        
        if debug_ad:
            print "dd results(before)"
            for ddk in bptemp:
                print ddk, dd[ddk] 
            #
        #
        
        # Now we have to filter the initial data so that only parallel
        # stems are examine.
        
        for ddk in bptemp:
            pkpptest = dd[ddk]
            
            for kt in range(0, len(pkpptest)-1):
                ipp1 = pkpptest[kt  ][0]; jpp1 = pkpptest[kt  ][1]
                ipp2 = pkpptest[kt+1][0]; jpp2 = pkpptest[kt+1][1]
                
                if not (ipp1 < ipp2 and jpp1 < jpp2 and ipp2 < jpp1):
                    # if it is not parallel (ipp1 < ipp2 < jpp1 <
                    # jpp2), then we are not interested at this
                    # point. We are only interested in cases where it
                    # is definitely a parallel stem.
                    dd[ddk] = []
                    break
                #
            #
        #
        
        # delete results that don't satisfy the requirement
        #print self.disp_allLists()
        for k in range(0, len(bptemp)):
            if dd.has_key(bptemp[k]):
                #print "bptemp: ", bptemp, dd[bptemp[k]]
                if len(dd[bptemp[k]]) == 0: 
                    del dd[bptemp[k]]
                #
            #
        #
                       
        #print self.disp_allLists()
        #sys.exit(0)
        for ddk in dd.keys():
            i = ddk[0]; j = ddk[1]

            # find the corresponding Pair class object in BPlist
            for k in range(0, len(self.BPlist)):
                if i == self.BPlist[k].i and j == self.BPlist[k].j:
                    self.BPlist[k].v = 'p'
                #
            #
            for pkk in dd[ddk]:
                # print pkk
                ipk = pkk[0]; jpk = pkk[1]
                # find the corresponding Pair class object in PKlist
                for k in range(0, len(self.PKlist)):
                    if ipk == self.PKlist[k].i and jpk == self.PKlist[k].j:
                        self.PKlist[k].v = 'p'
                    #
                #
            #
        #
                
            
        if debug_ad:
            print "dd results(after)"
            for ddk in dd.keys():
                print ddk, dd[ddk] 
            #
            print "BPlist:"
            for vv in self.BPlist:
                print vv.disp_Pair()
            #
            print "PKlist:"
            for vv in self.PKlist:
                print vv.disp_Pair()
            #
            
            # print "stop at 0.2 in assign_PKdirection", sys.exit(0)
        #
                        
        
        # Note: "plist" is the current PK list, not all bps
        
        plist = self.sortvsList(plist, 'i')
        # plist.sort(key=self.getKey_i)
        for k in range(0, len(plist)-1):
            i_kp0 = plist[k  ].i; j_kp0 = plist[k  ].j
            i_kp1 = plist[k+1].i; j_kp1 = plist[k+1].j
            if debug_ad:
                print "a) ij_kp0[=%2d](%2d,%2d) ... ij_kp1[=%2d](%2d,%2d)" % \
                    (k, i_kp0, j_kp0, k+1, i_kp1, j_kp1);
            #
            if i_kp0 <  i_kp1 and i_kp1 < j_kp0 and j_kp0 < j_kp1:
                """
                First, if it is a parallel stem, then it must be that it looks
                like a PK with the linkage to the right side.
                
                "...A.....B.........a........B..." 
                    ^     ^         ^        ^
                 i_kp0  i_kp1     i_kp1    j_kp0
                
                """
                    
                
                test_km1 = False; test_kp2 = False
                if k == 0 and len(plist) == 2:
                    # we only have these two pairs so there is no k-1
                    # or k+2 case to test this against.
                    test_km1 = True
                    test_kp2 = True
                #
                """
                We are _given_ the following
                
                i_kp0 < i_kp1 < j_kp0 < j_kp1
                
                ============================================
                this case should be accepted
                    k=0
                       i_kp1  i_kp2      j_kp2      j_kp1   
                         v      v          v          v            
                      ..AB......[[........]].........ab...
                        A                            A             
                      i_kp0                        j_kp0
                
                    i_kp0 < i_kp1 < j_kp0 < j_kp1 (given)
                    i_kp1 < i_kp2 < j_kp2 < j_kp1 **
                
                
                this case should be rejected
                    k=1
                       i_km1  i_kp1      j_kp1     j_km1   
                        v       v          v         v            
                      ..AB......[[........]].........ab...
                         A       ^        ^           A             
                       i_kp0   i_kp2    j_kp2       j_kp0
                    
                    i_km1 < i_kp0 < j_km1 < j_kp0
                    i_kp0 < i_kp1 < j_kp1 < j_kp0 fails here
                    i_kp1 < i_kp2 < j_kp2 < j_kp1 **
                
                ============================================
                this case should be rejected
                    k=1
                      i_km1   i_kp1     j_kp1       j_km1   
                        v       v         v           v            
                      ..[[......AB........ab.........]]...
                         A       ^         ^         A             
                       i_kp0   i_kp2     j_kp2     j_kp0
                
                    
                    i_km1 < i_kp0 < j_kp0 < j_kp1 **
                    i_kp0 < i_kp1 < j_kp1 < j_kp0 fails here
                    i_kp1 < i_kp2 < j_kp1 < j_kp2 
                
                this case should be accepted
                    k=2
                       i_km1   i_kp1     j_kp1     j_km1   
                         v       v         v         v            
                      ..[[......AB........ab.........]]...
                                A         A             
                              i_kp0     j_kp0
                
                    i_km1 < i_kp0 < j_kp0 < j_kp1 **
                    i_kp0 < i_kp1 < j_kp0 < j_kp1
                
                
                ============================================
                this case should be accepted
                    k=0
                       i_kp1       j_kp1 i_kp2      j_kp2   
                         v           v     v          v            
                      ..AB..........ab.....[[........]]...
                        A           A             
                      i_kp0       j_kp0
                    
                    i_kp0 < i_kp1 < j_kp0 < j_kp1 (given)
                    i_kp1 < j_kp1 < i_kp2 < j_kp2 **
                            -------------
                
                    
                ============================================
                this case should be accepted
                    k=2
                       i_km1     j_km1   i_kp1        j_kp1   
                         v        v         v           v            
                      ..[[........]].......AB..........ab.
                                           A           A             
                                         i_kp0       j_kp0
                    
                    i_km1 < j_km1 < i_kp0 < j_kp0 **
                    i_kp0 < i_kp1 < j_kp0 < j_kp1 (given)
                            -------------
                
                ============================================
                this case should be rejected
                    k=0
                       i_kp1     i_kp2   j_kp1        j_kp2   
                         v        v        v            v            
                      ..AA........BB.......aa..........bb.
                        A                   A             
                      i_kp0               j_kp0
                    
                    i_kp0 < i_kp1 < j_kp1 < j_kp0 fails here
                            -------------
                    i_kp1 < i_kp2 < j_kp1 < j_kp2 
                
                this case should be rejected
                    k=1
                       i_km1   i_kp1,2    j_km1      j_kp1,2   
                        v         vx        v          xv            
                      ..AA........BB.......aa..........bb.
                         A                 A             
                       i_kp0             j_kp0
                    
                    i_km1 < j_kp0 < i_kp0 < j_km1 **
                            -------------
                    i_kp0 < i_kp1 < j_kp0 < j_kp1 (given)
                    i_kp1 < i_kp2 < j_kp2 < j_kp1 **
                            -------------
                
                this case should be rejected
                    k=2
                       i_km1     i_kp1   j_km1      j_kp1   
                         v         v       v           v            
                      ..AA........BB.......aa..........bb.
                                  A                     A             
                                i_kp0                 j_kp0
                    
                    i_km1 < j_kp0 < i_km1 < j_kp0 
                    i_kp0 < i_kp1 < j_kp1 < j_kp0 fails here
                            -------------
                
                """
                if k > 0:
                    i_km1 = plist[k-1].i; j_km1 = plist[k-1].j
                    """
                    We are _given_ i_kp0 < i_kp1 < j_kp0 < j_kp1
                    
                    k=1,2 
                    
                    Both cases fail the above condition so they don't
                    end up here
                    
                    k=3 (first case satisfied)
                         i_kp0  i_kp1        j_kp0       j_kp1
                           V      v           V            v
                      ..[[[[......ABCD........]]]].........abcd...  (i_km1 < i_kp0 and j_kp0 < j_km1) <===
                          ^                    ^                    (i_kp0 < i_kp1  <  jkp0  < j_kp1)[*]
                        i_km1                j_km1                  should be rejected
                    
                    k=4
                               i_kp0,1,2                j_kp0,1,2         
                                  Vvv                      Vvv           
                      ..[[[[......ABCD........]]]].........abcd...  (i_km1 < i_kp0  < j_km1 < j_kp0)
                           ^                  ^                     (i_kp0 < i_kp1  < j_kp0 < j_kp1)[*]
                         i_km1              j_km1                   (i_kp1 < i_kp2  < j_kp1 < j_kp2)
                                                                    should be accepted
                    
                    k=5
                                i_kp0,1,2                j_kp0,1,2         
                                   Vvv                      Vvv            
                      ..[[[[......ABCD........]]]].........abcd...  (i_km1 < i_kp0  < j_km1 < j_kp0)
                                  ^                        ^        (i_kp0 < i_kp1  < j_kp0 < j_kp1)[*]
                                i_km1                    j_km1      (i_kp1 < i_kp2  < j_kp1 < j_kp2)
                                                                    should be accepted
                    
                    k=6
                                 i_kp0,1                   j_kp0,1         
                                    Vv                       Vv            
                      ..[[[[......ABCD........]]]].........abcd...  (i_km1 < i_kp0  < j_km1 < j_kp0)
                                   ^                        ^       (i_kp0 < i_kp1  < j_kp0 < j_kp1)[*]
                                 i_km1                    j_km1     should be accepted
                    
                    """
                    if debug_ad:
                        print "b) ij_kp0[=%2d](%2d,%2d) ... ij_km1[=%2d](%2d,%2d)" % \
                            (k, i_kp0, j_kp0, k-1, i_km1, j_km1);
                    #
                    """
                    This is the former method, I think it is not so secure
                    
                    k       i_km1 < i_kp0    i_kp0 < j_kp1     j_kp0 < j_km1    j_km1 < j_kp1   ijkp01
                    3             T                T                 T                 T       --> a
                    4             T                T                 F                 T       --> p
                    5             T                T                 F                 T       --> p
                    6             T                T                 F                 T       --> p
                    #if not (i_km1 < i_kp0 and i_kp0 < j_kp1 and j_kp0 < j_km1 and j_km1 < j_kp1):
                    #    test_km1 = True
                    ##
                    """
                    
                    if (not (i_km1 < i_kp0 and j_kp0 < j_km1) or # <=== (i_km1 < i_kp0 and j_kp0 < j_km1) 
                        (j_km1 < i_kp0) or 
                        (i_km1 < i_kp0  and j_km1 < j_kp0 and i_km1 < j_kp0)):
                        test_km1 = True
                    #
                    
                    
                #
                if debug_ad:
                    print "k= %d, len(plist) = %d" % (k, len(plist))
                #
                if k + 2 < len(plist):
                    i_kp2 = plist[k+2].i; j_kp2 = plist[k+2].j
                    if debug_ad:
                        print "c) ij_kp0[=%2d](%2d,%2d) ... ij_kp1[=%2d](%2d,%2d) ... ij_kp2[=%2d](%2d,%2d)" % \
                            (k, i_kp0, j_kp0, k+1, i_kp1, j_kp1, k+2, i_kp2, j_kp2);
                    #
                    """
                    This is the former method, I think it is not so secure
                    
                    k       i_kp0 < i_kp1    i_kp1 < i_kp2     j_kp0 < j_kp2     j_kp2 < j_kp1  ijkp01
                    3             T                T                 T                 F       --> p
                    4             T                T                 T                 F       --> p
                    5             T                T                 T                 F       --> p
                    #if not (i_kp0 < i_kp1 and i_kp1 < i_kp2 and j_kp0 < j_kp2 and j_kp2 < j_kp1):
                    #    test_kp2 = True
                    ##
                    """
                    if (not(i_kp1 < i_kp2 and j_kp2 < j_kp1) or # <=== (i_kp1 < i_kp2 and j_kp2 < j_kp1) 
                        (j_kp1 < i_kp2) or 
                        (i_kp1 < i_kp2  and j_kp1 < j_kp2 and i_kp2 < j_kp1)):
                        test_kp2 = True
                    #
                #
                if debug_ad:
                    print "test_km1(%s), test_kp2(%s)" % (test_km1, test_kp2)
                #
                if test_km1 or test_kp2:
                    plist[k  ].v = 'p'
                    plist[k+1].v = 'p'
                #
            #
        #
        if debug_ad:
            print "PKlist: (after)"
            for vv in plist:
                print vv.disp_Pair()
            #
            print "BPlist: (after)"
            for vv in self.BPlist:
                print vv.disp_Pair()
            #
            print "Exit assign_PKdirection: "
            #print "stop at end of assign_PKdirection"; sys.exit(0)
        #
        
        return plist
    #
    #########
    
    
    def findroots(self, prlist, debug = False):
        debug_findroots = debug
        if debug_findroots:
            print "Enter findroots:"
            print "prlist:"
            for vv in prlist:
                print vv.disp_Pair()
            #
        #
        nprlist = []
        for pr_k in prlist:
            # make sure we don't just make a pointer
            nprlist += [pr_k]
        #
        
        k = 0
        while k < (len(nprlist) - 1):
            i_ref = nprlist[k].i; j_ref = nprlist[k].j; tp_ref = nprlist[k].v
            if debug_findroots:
                print "ij_ref: ", nprlist[k].disp_Pair()
            #
            
            ll = 0
            while ll < len(nprlist):
                i_t = nprlist[ll].i; j_t = nprlist[ll].j; tp_t = nprlist[ll].v
                if i_t == i_ref and j_t == j_ref:
                    ll += 1
                    continue
                #
                if 0: #debug_findroots:
                    print "k = %2d, ll = %2d: i_ref(%2d) < i_t(%2d) and j_t(%2d) < j_ref(%2d)" \
                        % (k, ll, i_ref,  i_t, j_t, j_ref)
                #
                delete_t = False
                if not (tp_t == 'p' and tp_ref == 'p'):
                    # still not sure if this is the only condition or
                    # if this uniquely defines the problem of a
                    # parallel stem, but this is what I can currently
                    # think of to make it so the PK roots contain all
                    # these elements.
                    if i_ref < i_t and j_t < j_ref:
                        if (tp_t == 'a' and tp_ref == 'a'):
                            delete_t = True
                        #
                    #
                #
                if delete_t:    
                    if debug_findroots:
                        print "deleting ", nprlist[ll].disp_Pair()
                    #
                    del nprlist[ll]
                else:
                    ll += 1
                #
            #
            k+= 1
            
        #
        if debug_findroots:
            print "end of findroots: vvvvvvvvvvvvv"
            print "nprlist: "
            for bproot in nprlist:
                print bproot.disp_Pair()
            #
            print "prlist: "
            for bproot in prlist:
                print bproot.disp_Pair()
            #
            print "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
            #print "stop at the exit point from findroots"; sys.exit(0)
        #
        return nprlist
    #
    
    def test_entwined(self, pkroots, ssroots, debug_test_entwined = False):
        """
        Consider this example (a structure in Threads.py)
        
           ".((....AA...))..BB....aa....CC...bb....cc.."
        
        In this routine, ssroots and pkroots leave us with the
        following paired down structure
        
           ".((....AA...))..BB....aa....CC...bb....cc.."
        
             |                  |                   |
             V                  V                   V
        
           ".(.....A.....)..B......a....C.....b.....c.."
        
        and this routine isolates the PK roots
        
           "................B...........C.....b.....c.."
        
        So the output from this routine reference to the PK roots
        corresponding to Bb and Cc.
        
        Likewise, following the same reasoning, this next structure
        is processed in the following way.
        
           ".((....AA...))..BB....CC....aa...bb....cc.." 
        
             |                  |                   |
             V                  V                   V
        
           ".(.....A.....)..B.....C......a....b.....c.."
        
             |                  |                   |
             V                  V                   V
        
           "................B.....C...........b.....c.."
        """
        
        debug_test_entwined = False # True # 
        stop_at_end         = False # True # 
        if debug_test_entwined:
            print "Enter test_entwined: (before)"
            print "ssroots: "
            for vv in ssroots:
                print vv.disp_Pair()
            #
            print "pkroots: "
            for vv in pkroots:
                print vv.disp_Pair()
            #
            print "------"
        #
        ss_pssbl = []  # possible secondary structure roots
        for pklk in pkroots:
            ss_pssbl += [pklk]
        #
        kp = 0
        while kp < len(ss_pssbl):
            ss_entwined = False
            i_p = ss_pssbl[kp].i; j_p = ss_pssbl[kp].j; v_p = ss_pssbl[kp].v
            kr = 0
            for kr in range(0, len(ssroots)):
                # compare the secondary structure roots with the 
                i_r = ssroots[kr].i; j_r = ssroots[kr].j; v_r = ssroots[kr].v
                """
                 ...[.....(.......].......)..     ...(.....[.......).......]..
                    i_p   i_r     j_p     j_r        i_r   i_p     j_r     j_p
                """
                if debug_test_entwined:
                    if (i_p < i_r and j_p < j_r and i_r < j_p):
                        # i_p < i_r < j_p < j_r
                        print "1 i_p(%2d) < i_r(%2d) & j_p(%2d) < j_r(%2d) & i_r(%2d) < j_p(%2d), v_p(%s) v_r(%s)" \
                            % (i_p, i_r, j_p, j_r, i_r, j_p, v_p, v_r)
                    elif (i_r < i_p and j_r < j_p and i_p < j_r):
                        # i_r < i_p < j_r < j_p
                        print "2 i_r(%2d) < i_p(%2d) & j_r(%2d) < j_p(%2d) & i_p(%2d) < j_r(%2d), v_p(%s) v_r(%s)" \
                            % (i_r, i_p, j_r, j_p, i_p, j_r, v_p, v_r)
                    elif (i_r < i_p and j_p < j_r): 
                        print "3 i_r(%2d) < i_p(%2d) < j_p(%2d) < j_r(%2d), v_p(%s) v_r(%s)" \
                            % (i_r, i_p, j_p, j_r, v_p, v_r)
                    elif (i_p < i_r and j_r < j_p): 
                        print "4 i_p(%2d) < i_r(%2d) < j_r(%2d) < j_p(%2d), v_p(%s) v_r(%s)" \
                            % (i_p, i_r, j_r, j_p, v_p, v_r)
                    elif (j_r < i_p): 
                        print "5 i_r(%2d) < j_r(%2d) < i_p(%2d) < j_p(%2d), v_p(%s) v_r(%s)" \
                            % (i_r, j_r, i_p, j_p, v_p, v_r)
                    elif (j_p < i_r):
                        print "6 i_p(%2d) < j_p(%2d) < i_r(%2d) < j_r(%2d), v_p(%s) v_r(%s)" \
                            % (i_p, j_p, i_r, j_r, v_p, v_r)
                    elif (i_r == i_p and j_r == j_p):
                        print "7 ij_p(%2d,%2d) <=> ij_r(%2d,%2d), v_p(%s) v_r(%s)" \
                            % (i_p, j_p, i_r, j_r, v_p, v_r)
                    else:
                        print "don't know what is going on"
                        print "8 ij_p(%2d,%2d) <??> ij_r(%2d,%2d), v_p(%s) v_r(%s)" \
                            % (i_p, j_p, i_r, j_r, v_p, v_r)
                        sys.exit(0)
                    #
                #
                if ((i_p < i_r and j_p < j_r and i_r < j_p) or
                    (i_r < i_p and j_r < j_p and i_p < j_r)):
                    if debug_test_entwined:
                        print "deleting ss_pssbl: ", ss_pssbl[kp].disp_Pair()
                    #
                    del ss_pssbl[kp] # pk meshed with ssroots
                    ss_entwined = True
                    break
                    #
                #
            #
            if not ss_entwined:
                kp += 1
            #
        #
        
        if debug_test_entwined:
            print "check overlaps with BPlist:  vvvvvvvvvvvvvv"
        #
        kv = 0
        while kv < len(ss_pssbl):
            ic = ss_pssbl[kv].i; jc = ss_pssbl[kv].j
            if debug_test_entwined:
                print "ijc(%2d,%2d)" % (ic, jc)
            #
            flag_del = False
            for bpk in self.BPlist:
                iv = bpk.i; jv = bpk.j
                if    iv < ic and jv < jc and ic < jv:
                    if debug_test_entwined:
                        print "overlap with region iv(%2d) < ic(%2d) < jv(%2d) < jc(%2d)" \
                            % (iv, ic, jv, jc)
                    #
                    del ss_pssbl[kv]
                    flag_del = True
                    break
                elif ic < iv and jc < jv and iv < jc:
                    if debug_test_entwined:
                        print "overlap with region ic(%2d) < iv(%2d) < jc(%2d) < jv(%2d)" \
                            % (ic, iv, jc, jv)
                    #
                    del ss_pssbl[kv]
                    flag_del = True
                    break
                #
            #
            if not flag_del:
                kv += 1
            #
        #
        
        if debug_test_entwined:
            print "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
            print "Exit test_entwined: (after)"
            print "ssroots: "
            for vv in ssroots:
                print vv.disp_Pair()
            #
            print "pkroots: "
            for vv in pkroots:
                print vv.disp_Pair()
            #
            print "ss_pssbl: "
            for nn in ss_pssbl:
                print nn.disp_Pair()
            #
            print "------"
            if stop_at_end:
                print "requested to stop at the end of test_entwined()"
                sys.exit(0)
            #
        #
        return ss_pssbl
    #

    
    def compress_to_BPlist(self):
        """
        Compress_to_BPlist: organizes BPlist, and PKlist so that it fits a
        kind of standard format of secondary structure and PK
        notation. It ialso generates two new important pieces of
        information: BProots and PKroots. These latter lists can be
        used to identify the base of each secondary structure domain
        and the arrangement of the pseudoknots that hangs over them.
    
    
    Example 1:
                  0         10        20        30        40 
                  |    .    |    .    |    .    |    .    |  
       ss_seq  = ".A..B..C..D...........a..b..c..d.."
       
       result:
       new BPlist:
       (   1,   22)[p]
       new BProots:
       (   1,   22)[p]
       new PKlist:
       (   4,   25)[p]
       (   7,   28)[p]
       (  10,   31)[p]
       new PKroots:
       (   4,   25)[p]
       (   7,   28)[p]
       (  10,   31)[p]
    
        [Note that for parallel stems, this operation does not really
        do all that much to reduce the size of the set of roots of
        parallel stems. It seems like this could be handled, as it is
        just difference of 1 between (i1,j1) and (i2,j2)
        
        i.e., (i2 - i1) = 1 and (j2 - j1) = 1                           
    
        However, presently, the output will always include except the
        root stem.]
    
    Example 2:
                  0         10        20        30        40 
                  |    .    |    .    |    .    |    .    |  
       ss_seq  = ".AAA...BBB...CCC....aaa...bbb....ccc...."
       
       result:
       new BPlist:
       (   1,   22)[a]
       (   2,   21)[a]
       (   3,   20)[a]
       new BProots:
       (   1,   22)[p]
       new PKlist:
       (   7,   28)[a]
       (   8,   27)[a]
       (   9,   26)[a]
       (  13,   35)[a]
       (  14,   34)[a]
       (  15,   33)[a]
       new PKroots:
       (   7,   28)[a]
       (  13,   35)[a]
    
    Example 3:
                  0         10        20        30        40 
                  |    .    |    .    |    .    |    .    |  
       ss_seq  = ".AAA...BBB...aaa...CCC....bbb....ccc...." # same
       ss_seq  = ".(((...BBB...)))...(((....bbb....)))...." # same
   
       result:
       new BPlist:
       (   1,   15)[a]
       (   2,   14)[a]
       (   3,   13)[a]
       (  19,   35)[a]
       (  20,   34)[a]
       (  21,   33)[a]
       new BProots:
       (   1,   15)[a]
       (  19,   35)[a]
       new PKlist:
       (   7,   28)[a]
       (   8,   27)[a]
       (   9,   26)[a]
       new PKroots:
       (   7,   28)[a]

        """
        debug_cbpl  = False # True # 
        stop_at_end = False # True # 
        if debug_cbpl:
            print "enter compress_to_BPlist:"
        #  
        
        # First, we need to build direct copies of self.PKlist and
        # self.BProots and sort the lists with respect to the first
        # element.
        ssroots = copyList(self.BProots)
        ssroots = self.sortvsList(ssroots, 'i')
        #ssroots = sorted(ssroots, key=self.getKey_i)
        
        pklist = copyList(self.PKlist)
        # pklist is already sorted. 
        pkroots = self.findroots(pklist, debug_cbpl) 
        pkroots = self.sortvsList(pkroots, 'i')
        #self.PKroots = sorted(pkroots, key=self.getKey_i)
        
        # print out the current list of ss-domains and PK regions
        if debug_cbpl:
            print "Before: vvvvvvvvvvvvvvvvvv"
            print "base pair lists: "
            print "BPlist:"
            for sspr_k in self.BPlist:
                print sspr_k.disp_Pair()
            #
            print "root lists: "
            print "ssroots:"
            for pkpr_k in self.BProots:
                print pkpr_k.disp_Pair()
            #
            print "pseudoknot lists"
            print "PKlist:"
            for pkpr_k in self.PKlist:
                print pkpr_k.disp_Pair()
            #
            print "pkroots:"
            for pkpr_k in self.PKroots:
                print pkpr_k.disp_Pair()
            #
            print "^^^^^^^^^^^^^^^^^^ :Before"
            #print "stop at 0 in compress_to_BPlist"; sys.exit(0)
        #
        
        """
        Second, we need to isolate the regions where the various PK root
        entries are free of any secondary structure root
        entanglement. The resulting PK root entries may still be an
        entangled list of PKs, but after this operation, the remaining
        root PKs are independent of any secondary structure roots. In
        the first operation, we prune out all PK roots that are
        entangled with the root secondary structure.
        
        In this first step, we scan the list of secondary structure
        against the PK list and look for regions where there are
        clearly additional possible roots that may also be
        equivalently represented as secondary structure.
        
        """
        
        ss_pssbl = self.test_entwined(pkroots, ssroots, debug_cbpl)
        
        """
        ".((....AA...))..BB....aa....CC...bb....cc.."
        
        In this routine, ssroots and pkroots leave us with the
        following paired down structure
        
           ".((....AA...))..BB....aa....CC...bb....cc.."
           
             |                  |                   |
             V                  V                   V
           
           ".(.....A.....)..B......a....C.....b.....c.."
           
             |                  |                   |
             V                  V                   V
           
           "................B...........C.....b.....c.."
           
        In this next step, Bb and Cc are examined and, in this
        example, this algorithm will dig out Bb.
           
             |                  |                   |
             V                  V                   V
           
           "................B.................b........"
           
           ss_pssbl:
           (  16,  34)[a]
           (  28,  40)[a]
           revised ss_psbbl:
           (  16,  34)[a]
        """
       
        kp = 0
        while kp < (len(ss_pssbl)-1):
            i_p = ss_pssbl[kp].i; j_p = ss_pssbl[kp].j
            # remember, this is an ordered list!!!!
            kn = 0
            while kn < len(ss_pssbl):
                # compare the secondary structure roots with the 
                i_n = ss_pssbl[kn].i; j_n = ss_pssbl[kn].j
                """
                ...[.....(.......].......)..     ...(.....[.......).......]..
                   i_p   i_n     j_p     j_n        i_n   i_p     j_n     j_p
                """
                if ((i_p < i_n and j_p < j_n and i_n < j_p) or
                    (i_n < i_p and j_n < j_p and i_p < j_n)):
                    if debug_cbpl:
                        if   (i_p < i_n and j_p < j_n and i_n < j_p):
                            print "i_p(%2d) < i_n(%2d) < j_p(%2d) < j_n(%2d)" \
                                % (i_p, i_n, j_p, j_n)
                        elif (i_n < i_p and j_n < j_p and i_p < j_n):
                            print "i_n(%2d) < i_p(%2d) < j_n(%2d) < j_p(%2d)" \
                                % (i_n, i_p, j_n, j_p)
                        else:
                            print "problems! cbpl"
                            print "ij_n(%2d,%2d), ij_p(%2d,%2d)" % (i_n, j_n, i_p, j_p)
                    del ss_pssbl[kn] # pk meshed with ssroots
                else:
                    kn += 1
                #
            #
            kp += 1
        #
        if debug_cbpl:
            print "revised ss_pssbl: "
            for nn in ss_pssbl:
                print nn.disp_Pair()
            #
            #print "stop at 1 in compress_to_BPlist"; sys.exit(0)
        #
        
        """
        What have we done? Suppose that the following PK roots were
        extracted in the first step
        
           A.....B.......a...C.....b.......c
        
        This is the first result we would find in ss_pssbl.
        
        Because Bb is entangled between both Aa and Cc, this is 
        restructured to the following
        
           A.....B.......a...C.....b.......c -> (.....B.......)...(.....b.......)
        """
        pk2ssroots = copyList(ss_pssbl)
        pk2ssroots = self.sortvsList(pk2ssroots, 'i')
        #pk2ssroots = sorted(pk2ssroots, key=self.getKey_i)
        
        if debug_cbpl:
            print "PKroots:"
            for vv in pkroots:
                print vv.disp_Pair()
            #
            print "pk2ssroots:"
            for vv in pk2ssroots:
                print vv.disp_Pair()
            #
        #
        
        kv = 0
        for pk2ss in pk2ssroots:
            it = pk2ss.i; jt = pk2ss.j
            while kv < len(pkroots):
                ipkr = pkroots[kv].i; jpkr = pkroots[kv].j
                if debug_cbpl:
                    print "ijt(%2d,%2d) vs ijpkr(%2d,%2d)" % (it, jt, ipkr, jpkr)
                #
                if it == ipkr and jt == jpkr:
                    if debug_cbpl:
                        print "delete: ", pkroots[kv]
                    #
                    del pkroots[kv]
                else:
                    kv += 1
                #
            #
        #
        
        del self.PKroots
        self.PKroots  = copyList(pkroots)
        
        if debug_cbpl:
            print "------"
            print "pkroots: "
            for nn in self.PKroots:
                print nn.disp_Pair()
            #
            print "revised ss_pssbl: "
            for nn in ss_pssbl:
                print nn.disp_Pair()
            #
            print "pk2ssroots: (should be the same as ss_pssbl)"
            for nn in pk2ssroots:
                print nn.disp_Pair()
            #
            print "------"
            
            #print "stop at 2 in compress_to_BPlist"; sys.exit(0)
        #
        
        # now we have to rebuild the pk and ss structures from the
        # root information we have acquired and transfer this information. 
        
        
        BPgroup = copyList(self.BPlist)
        # We don't bother to sort this list because we may end up
        # adding to it from pklist and (moreover) the order of the
        # BPgroup list doesn't matter in the operations being done
        # here.
        for k in range(0,len(pk2ssroots)):
            i_pkr = pk2ssroots[k].i; j_pkr = pk2ssroots[k].j; v_pkr = pk2ssroots[k].v 
            if debug_cbpl:
                print "k  = %2d, ij_pkr(%2d,%2d)[%s]" % (k, i_pkr, j_pkr, v_pkr)
            #
            kl = 0
            while kl < len(pklist):
                i_pkt = pklist[kl].i; j_pkt = pklist[kl].j; v_pkt = pklist[kl].v
                if debug_cbpl:
                    print "kl = %2d, ij_pkt(%2d,%2d)[%s]" % (kl, i_pkt, j_pkt, v_pkt)
                    #sys.exit(0)
                #
                if v_pkr == 'a' and v_pkt == 'a':
                    if debug_cbpl:
                        print "v_pkr = %s, v_pkt = %s" % (v_pkr, v_pkt)
                        if i_pkr <= i_pkt and j_pkt <= j_pkr:
                            print "add to BPgroup: ", pklist[kl]
                        #
                    #
                    if i_pkr <= i_pkt and j_pkt <= j_pkr:
                        BPgroup += [pklist[kl]]
                        del pklist[kl]
                    else:
                        kl += 1
                    #
                elif i_pkr == i_pkt and j_pkt == j_pkr:
                    if debug_cbpl:
                        print "ij_pkr(%2d,%2d) <=> ij_pkt(%2d,%2d)" \
                            % (i_pkr, j_pkr, i_pkt, j_pkt)
                    #
                    BPgroup += [pklist[kl]]
                    del pklist[kl]
                else:
                  kl += 1
                #
            #
        #
        
        # Now we sort this list (after possibly updating it)
        
        BPgroup = self.sortvsList(BPgroup, 'i')
        #BPgroup.sort(key=self.getKey_i)
        
        # update the PKlist and rest the initial self.BPlist and
        # self.PKlist, finally reconstruct these two lists.
        PKgroup = copyList(pklist)
        PKgroup = self.sortvsList(PKgroup, 'i')
        #PKgroup.sort(key=self.getKey_i)
        del self.BPlist
        del self.PKlist
        
        self.BPlist = copyList(BPgroup)
        self.PKlist = copyList(PKgroup)
        self.BProots = self.sortvsList(self.findroots(self.BPlist), 'i')
        #self.BProots = sorted(self.findroots(self.BPlist, False), key=self.getKey_i)
        
        
        """
        This first sequence does exactly what I want it to do....
        
        ss_seq  = ".ABCD....EFG...efg..((..[[.)).]]..HIJ...hij....abcd.."
        > python Threads.py
              ...
        updated base pair lists: 
        new BPlist:
        (    1,   47)[p]
        (    9,   15)[p]
        (   20,   28)[a]
        (   21,   27)[a]
        (   34,   40)[p]
        new BProots:
        (    1,   47)[p]
        (    9,   15)[p]
        (   20,   28)[a]
        (   34,   40)[p]
        updated pseudoknot lists: 
        new PKlist:
        (    2,   48)[p]
        (    3,   49)[p]
        (    4,   50)[p]
        (   10,   16)[p]
        (   11,   17)[p]
        (   24,   31)[a]
        (   25,   30)[a]
        (   35,   41)[p]
        (   36,   42)[p]
        new PKroots:
        (    2,   48)[p]
        (    3,   49)[p]
        (    4,   50)[p]
        (    9,   15)[p]
        (   10,   16)[p]
        (   11,   17)[p]
        (   24,   31)[a]
        (   34,   40)[p]
        (   35,   41)[p]
        (   36,   42)[p]
        ^^^^^^^^^^^^^^^^^^ :After
        
        The following sequence is a bit more nasty and it required
        some fixes to test_entwined to get it to work a little
        better. The current output follows.
        
        ss_seq  = ".((((((.A.AAAAA.<BC..))))).((((.>bc....a..aaaaa....))))..)."
        
        After: vvvvvvvvvvvvvvvvvvv
        updated base pair lists: 
        new BPlist:
        (    1,   57)[a]
        (    2,   25)[a]
        (    3,   24)[a]
        (    4,   23)[a]
        (    5,   22)[a]
        (    6,   21)[a]
        (   27,   54)[a]
        (   28,   53)[a]
        (   29,   52)[a]
        (   30,   51)[a]
        new BProots:
        (    1,   57)[a]
        updated pseudoknot lists: 
        new PKlist:
        (    8,   46)[a]
        (   10,   45)[a]
        (   11,   44)[a]
        (   12,   43)[a]
        (   13,   42)[a]
        (   14,   39)[a]
        (   16,   32)[p]
        (   17,   33)[p]
        (   18,   34)[p]
        new PKroots:
        (    8,   46)[a]
        (   16,   32)[p]
        (   17,   33)[p]
        (   18,   34)[p]
        ^^^^^^^^^^^^^^^^^^ :After
        
        """
        if debug_cbpl:
            print "After: vvvvvvvvvvvvvvvvvvv"
            print "updated base pair lists: "
            print "new BPlist:"
            for bpg in self.BPlist:
                print bpg.disp_Pair()
            #
            print "new BProots:"
            for bpr in self.BProots:
                print bpr.disp_Pair()
            #
            print "updated pseudoknot lists: "
            print "new PKlist:"
            for pkg in self.PKlist:
                print pkg.disp_Pair()
            #
            print "new PKroots:"
            for pkr in self.PKroots:
                print pkr.disp_Pair()
            #
            print "^^^^^^^^^^^^^^^^^^ :After"
            print "Exit compress_to_BPlist:"
            #print "stop at 3 (exit from) compress_to_BPlist"; sys.exit(0)
            if stop_at_end:
                sys.exit(0)
            #
        #  
        return 0
    #
    
    
    
    def expand_island(self, CTCFisland):
        MPlist = [CTCFisland]
        n_k = len(CTCFisland.contacts)
        # print n_k
        if n_k == 1:
            for l in CTCFisland.contacts:
                b1 = Pair()
                b1.put_ssPair(CTCFisland.i, l, "ctcf", '-')
                b2 = Pair()
                b2.put_ssPair(l, CTCFisland.j, "ctcf", '-')
                MPlist += [b1,b2]
            #
        #
        elif n_k > 1:
            b = Pair()
            i = CTCFisland.i
            j = CTCFisland.contacts[0]
            b.put_ssPair(i, j, "ctcf", '-')
            MPlist += [b]
            for l in range(1,len(CTCFisland.contacts)):
                for k in range(0,l):
                    if l > k + 1:
                        continue
                    #
                    b = Pair()
                    i = CTCFisland.contacts[k]
                    j = CTCFisland.contacts[l]
                    # print i, j
                    b.put_ssPair(i, j, "ctcf", '-')
                    MPlist += [b]
                #
            #
            b = Pair()
            n = len(CTCFisland.contacts)
            i = CTCFisland.contacts[n-1]
            j = CTCFisland.j
            b.put_ssPair(i, j, "ctcf", '-')
            MPlist += [b]
        #
        
        return MPlist
    #
    
    
    """
    Since the pair elements are derived from the order in which they
    appear on the stack, after you encounter a "mate", you search back
    through the stack until you find one unassigned element. So
    essentially, you first put a bunch of elements (e.g., 'X') on the
    stack (of course, assuming you have a string of them), now, you
    have just encountered the partner of this example ('x') and you
    look at its position on the stack using p (i.e., the
    pointer[rpr2num['x']] position last incremented). From here of
    course, you know that the position on the stack must decrease by
    exactly the number of 'X' items you put on the stack (more to the
    point, it BETTER match!!!)
    
    For a more concrete example, suppose you have the structure
    
       ".X.XXXXX......x..xxxxx.."
    
    where element X has the index 27
    
    we reach the top of the stack with counter[27] = pointer[27] = stack[27] = 6
    
    self.Xlist[27] should now look like this:
    
      ( 1, -1)
      ( 3, -1)
      ( 4, -1)
      ( 5, -1)
      ( 6, -1)
      ( 7, -1)
    
    Now we go through this list looking for the -1. In this case, it
    is at the top of the list, so it is the first element. So we
    assign pp = 5 (because the list starts from 0 so the last element
    is 5).
    
     ( 1, -1)
     ( 3, -1)
     ( 4, -1)
     ( 5, -1)
     ( 6, -1)
     ( 7, 14)
    
    As we find the remaining elements, pointer[27] decreases. So after
    6 iterations, we obtain the correct assignment of elements on the
    stack.
    
     ( 1, 21)
     ( 3, 20)
     ( 4, 19)
     ( 5, 18)
     ( 6, 17)
     ( 7, 14)
     ".X.XXXXX......x..xxxxx.."
    
    Personally, I think the looping should not occur, because that
    means that the stack is out of order, but maybe this search can be
    employed to detect problems in either building the stack or
    closing it. So I have left this somewhat incongruous feature
    within the tool for the moment (since it is still under
    development). Better safe than sorry as they say.

    """
    
    def find_next_Xpoint_n(self, plist, p, nm = "BP", debug = False):
        if debug:
            print "find_next_Xpoint_n, p = %d, name(%s)" % (p, nm)
            print "plist:"
            for pk in plist:
                print pk.disp_Pair()
            #
        #
        pp = p
        j = plist[p-1].j
        if debug:
            print "pp(%3d), j(%3d)" % (pp, j)
        #
        while j > -1 and pp >= 0:
            pp -= 1
            j = plist[pp-1].j
            if debug:
                print "pp(%3d), j(%3d)" % (pp, j)
            #
        #
        if debug:
            print "result: pp = %d" % pp
        #
        return pp
    #
    
    
    def check_Xlist_n(self, plist):
        pass_X = True
        for b in plist:
            if b.j < 0:
                pass_X = False
            #
        #
        return pass_X
    #
    def print_vstr(self):
        print self.vstr
    #
    def print_vseq(self):
        
        if not self.vseq == "":
            print self.vseq
        #
        #print self.vseq
        
    #
    
    
    def print_Xlist_n(self, plist):
        for b in plist:
            print b.disp_Pair()
        #
    #

    # only used for secondary structure
    def parse_SecondaryStructure(self, r_ss, rseq = ""):
        self.reset_Vstruct()
        self.set_Vstruct(r_ss)
        if not rseq == "":
            self.set_Vseq(rseq)
        #
        
        self.scan_vs(0)
        self.print_vstr()
        self.print_Xlist_n(self.BPlist)
        return 0
    #
    
    def disp_allLists(self, sttmnt = "Results:"):
        print sttmnt
        print "structural sequence"
        self.print_vseq()
        self.print_vstr()
        print "secondary structure"
        self.print_Xlist_n(self.BPlist)
        print "secondary structure roots"
        self.print_Xlist_n(self.BProots)
        #
        if len(self.PKlist) > 0:
            print "pseudoknot linkages: "
            self.print_Xlist_n(self.PKlist)
            print "pseudoknot linkage roots"
            self.print_Xlist_n(self.PKroots)
        #
        if len(self.MPlist) > 0:
            print "CTCF connects"
            self.print_Xlist_n(self.MPlist)
        #
        if len(self.MPlist) > 0:
            islands = []
            for cl in self.MPlist:
                islands += [self.expand_island(cl)]
            #
            print "CTCF breakdown:"
            kk = 1
            for island_k in islands:
                print "island(%d):" % kk
                for ii in island_k:
                    print ii.disp_Pair()
                #
                kk += 1
            #
        #
    #
    
    
    # can do RNA or chromatin with all types of structures
    def parse_fullDotBracketStructure(self, r_ss, rseq = "", print_result=False):
        debug = False
        self.reset_Vstruct()
        self.set_Vstruct(r_ss)
        if not rseq == "":
            self.set_Vseq(rseq)
        #
        
        self.scan_allTypes(0, 0)
        if print_result:
            if debug:
                self.disp_allLists("Results at finishing Vienna")
                print "Exit parse_fullDotBracketStructure"
            else:
                self.disp_allLists()
            #
        #
        return 0
    #
#



# reads in Vienna package data created by the SimRNA_trafl2pdbs
# processing function
class ViennaData(object):
    #
    def __init__(self):
        # secondary structure created from SimRNA_trafl2pdbs processing
        self.ss_flnm = ''
        self.ss_data = []
        self.n_structs = -1  
    #
    
    ################################
    ########   Functions    ########  
    ################################
    
    def reset_ViennaData(self):
        self.ss_flnm = ''
        self.ss_data = []
        self.n_structs = -1  
        return 0
    #
    
    # reading in trajectory data (*.trafl) from SimRNA run
    def get_ViennaData(self, fl):
        self.ss_flnm = fl
        try: 
            input = open(self.ss_flnm, 'r')
        except:
            print 'ERROR: %s does not exist!' % self.ss_flnm
            sys.exit(1)
        #
        while True:
            string = input.readline().strip()        # ((((.....)))) data
            if not string: break
            self.ss_data += [string]
        #
        input.close()
        self.n_structs   = len(self.ss_data)
        if self.n_structs == 0:
            mssg = '\nERROR: no data read\n'
            raise MyException(mssg)
        return 0
    #
    
    # reading in trajectory data (*.trafl) from SimRNA run
    def read_SimRNA_ss(self, data):
        
        # print "read_SimRNA_ss(): \n%s" % data
        
        self.ss_data = data
        self.n_structs   = len(self.ss_data)
        if self.n_structs == 0:
            mssg = '\nERROR: no data read\n'
            raise MyException(mssg)
        # print self.n_structs, self.ss_data
        return 0
    #
#

def test1(cl):
    vs = Vstruct()
    ss_seq = "(((((....)))))"
    if len(cl) > 1:
        ss_seq = cl[1]
    #
    vs.parse_SecondaryStructure(ss_seq)
    # output 
    # (((((....)))))
    # (    0,   13)
    # (    1,   12)
    # (    2,   11)
    # (    3,   10)
    # (    4,    9)
    # 0
#   

def test2(cl):
    vs = Vstruct()
    #ss_seq = "((((.....))))...((....))."
    #ss_seq = "{(((((.[[...)))))..]]...|..([)]([)]...}"
    #ss_seq =  "{(((((.A.AAAAA.....)))))......a..aaaaa...|..([)]([)]...}....{.....}"
    #ss_seq =  "{(((((.A.AAAAA.CB..)))))..bc..a..aaaaa...|..([)]([)]...}....{.....}"
    #ss_seq =  "{((((((.A.AAAAA......))))).BBBB........a..aaaaa....bbbb..).|..([)]([)].ABC.abc.}....{.....}"
    #ss_seq =  "{((((((.A.AAAAA......))))).((((........a..aaaaa....))))..).|..([)]([)].ABC.abc.}....{.....}"
    ss_seq  =  "{((((((.A.AAAAA.<BC..))))).((((.>bc....a..aaaaa....))))..).|..([)]([)]...}....{.A.B.C..a.b.c.....}"
    #ss_seq =  "{((((((.A.AAAAA.<BC..))))).((((.>bc.DE.a..aaaaa..de))))..).|..([)]([)]...}....{.A.B.C..a.b.c.....}"
    ss_seq  =  "((((...AAAA...))))...BBBB....aaaa...((((...bbbb..))))..."
    if len(cl) > 1:
        ss_seq = cl[1]
    #
    vs.parse_fullDotBracketStructure(ss_seq, True)
    # output 
    # |(((((.[[...)))))..]]...|..([)]([)]...|
    # secondary structure
    # (    1,   16)
    # (    2,   15)
    # (    3,   14)
    # (    4,   13)
    # (    5,   12)
    # (   27,   29)
    # (   31,   33)
    # pk connects
    # (    7,   20)
    # (    8,   19)
    # (   28,   30)
    # (   32,   34)
    # CTCF connects
    # (    0,   24)
    # (    0,   38)
    # (   24,   38)
#


def main(cl):
    if TEST == 0 or TEST == 1:
        test1(cl)
    elif TEST == 2:
        test2(cl)
    #
#

# Main
if __name__ == '__main__':
    main(sys.argv)
#
