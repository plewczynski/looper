#!/usr/bin/env python

# package name: Vienna
# creater:      Wayne K Dawson
# creation date: ~2014/2015
# last update:   170126


# Originally, this was intended as a tool for handling Vienna package
# data with SimRNA. Now, it is used by some other programs,
# particularly chreval.py. It also has been significantly upgraded in
# for using with chreval such that it is able to represent a very
# diverse variety of structures containing pseudoknot like structures.

# This tool reads the Vienna formatted pairing sequence and converts
# it into a set of pair-contacts, or in the case of chromatin, both
# pair contacts and ctcf contacts.

# This can be tested in the following way
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# Python 2.7.6 (default, Mar 22 2014, 22:59:56) 
# [GCC 4.8.2] on linux2
# Type "help", "copyright", "credits" or "license" for more information.
# >>> import Vienna
# >>> vs = Vienna.Vstruct()
# >>> vs.convert_Vstruct("(((((....)))))")
# (((((....)))))
# (    0,   13)
# (    1,   12)
# (    2,   11)
# (    3,   10)
# (    4,    9)
# 0
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

import sys
debug = False
# debug = True

TEST = 2

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



# error handling 
class MyException(Exception):
   def __init__(self, st):
      self.stmnt = st
   #
   def __str__(self):
      return self.stmnt
   # Example:
   # if not self.set_var:
   #    mssg = "\nERROR, value of var undefined\n"
   #    ex=MyException(mssg)
   #    raise ex
#


#
class Pair:
   def __init__(self):
      self.name = ''
      self.i = -1
      self.j = -1
      self.v = 'a' # default, contact is antiparallel
      # 'a' = antiparallel
      # 'p' = parallel
      # '-' = not applicable; e.g., CTCFs
      self.contacts = []
   #

   def __repr__(self):
      return '{}: {} {} {} {} {}'.format(self.__class__.__name__,
                                         self.name,
                                         self.i,
                                         self.j,
                                         self.v,
                                         self.contacts)
   #
   
   def __cmp__(self, other):
      if hasattr(other, 'i'):
         return self.i.__cmp__(other.i)
      #
   #
   
   
   def put_Pair(self, i, j, nm = 'bp', v = 'a'):
      self.name = nm
      self.i = i
      self.j = j
      self.v = v
      return 0
   #
   
   # this is used!
   def put_Pair_i(self, i, nm = 'bp'):
      self.name = nm
      self.i = i
      return 0
   #
   
   # can be used to mark triple helices as well as CTCFs
   def put_contacts(self, i):
      self.contacts += [i]
      return 0
   #
   
   # this is used!
   def put_Pair_j(self, j, nm = 'bp'):
      if not self.name == nm:
         print "ERROR(Pair): improperly matched type"
         print "             (i,j) = (%d,%d)" % (self.i, j)
         print "             name(%s) != new name(%s) " % (self.name, nm)
      self.j = j
      return 0
   #
   def get_Pair(self):
      return self.i, self.j
   #
   
   def disp_Pair(self):
      s = ''
      if len(self.contacts) > 0:
         s  = '(%5d,%5d)[%s]: ' % (self.i, self.j, self.v)
         s += "{%5d, " % self.i
         for c in self.contacts:
            s += "%5d, " % c
         s += "%5d}" % self.j
      else:
         s  = '(%5d,%5d)[%s]' % (self.i, self.j, self.v)
      return s
   #       
#        



class Vstruct:
   def __init__(self):
      self.vstr='' # structure
      self.N = -1  # sequence length
      self.wt = -1 # used with my_generation to add specified weights to heatmaps
      # standard secondary structure (ss)
      self.BPlist=[]
      # a both ss and pseudoknot (PK) in same line
      self.PKlist= []
      # CTCF islands
      self.CTCFlist = []
      # reference
      self.Xlist = {}
      for xl in Xlist.keys():
         self.Xlist.update({xl : Xlist[xl]})
      #
      self.stack = {}
      for sl in stack.keys():
         self.stack.update({sl : stack[sl]})
      #
      self.counter = {}
      for cl in counter.keys():
         self.counter.update({cl : counter[cl]})
      #
      self.pointer = {}
      for pl in pointer.keys():
         self.pointer.update({pl : pointer[pl]})
      #
   #
   def set_Vstruct(self, s):
      self.vstr = s
      self.N = len(s)
      return 0
   #
   #
   def reset_Vstruct(self):
      self.vstr='' # input dot-bracket secondary structure [+ PK + CTCF]
      self.N = -1
      # standard secondary structure (ss)
      self.BPlist = []
      
      # pseudoknot (PK) branches in the same line
      self.PKlist  = []
      
      # CTCF islands
      self.CTCFlist = []
      
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
   # 
   def scan_vs(self, i):
      if debug:
         print "scan_vs(%d):" % i
      #
      mm = -1
      if i < self.N:
         s = self.vstr[i]
         if not (s == '(' or s == ')' or s == '[' or s == ']' or s == '.'):
            print 'ERROR: improperly defined Fontana formated sequence'
            print '       offending character \'%c\' located at position %d' % (s, i+1)
            sys.exit(1)
         #
         if s == '(':
            # add to the stack
            b = Pair()
            b.put_Pair_i(i, "bp")
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
            self.BPlist[self.pointer[0]-1].put_Pair_j(i, "bp")
            if debug:
               print self.BPlist[self.pointer[0]-1].disp_Pair()
               # print 'j = %5d, %s' % (i, s)
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
            # Either look for next closing bracket in sequence or
            # terminate at the end of the sequence if nothing is
            # found.
            while (s == '.') and (mm < self.N):
               if debug:
                  print 'i = %5d, %s' % (mm, s)
               mm += 1
               if mm == self.N:
                  break
               s = self.vstr[mm]
               
            self.scan_vs(mm)
            #
         #
      #
      else:
         if debug:
            print 'point = %d' % self.pointer[0]
         if not self.counter[0] == 0:
            case = 0
            if self.counter[0] > 0:
               case = 0
            else:
               case = 1
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
         if not self.check_Xlist_n(self.BPlist):
            print 'ERROR!!! Fontana notation is not correct!'
            print '         At least one structure of type \')...(\' was found'
            print 'input structure:'
            print self.vstr
            sys.exit(1)
                
      return mm
   #
   
   
   def scan_ctcf(self, i, layer):
      debug_bp   = False
      debug_PK   = False
      debug_ctcf = False
      if debug_bp or debug_PK or debug_ctcf:
         print "scan_ctcf(%d):" % i, self.pointer[0]
      #
      if layer > self.N:
         print "layer(%d) exceeds the sequence length(%d)" % (layer, self.N)
         sys.exit(1)
      #
      
      mm = -1
      if i < self.N:
         s = self.vstr[i]
         if not (lpr2num.has_key(s) or rpr2num.has_key(s) or s == '|' or s == '.'):
            print 'ERROR: improperly defined Fontana formated sequence'
            print '       offending character \'%c\' located at position %d' % (s, i+1)
            sys.exit(1)
         #
         if s == '(': # standard bp argument 
            # add to the stack
            BPndx = lpr2num[s]
            if debug_bp:
               print "BPndx(%2d) => %s" % (BPndx, s)
            #
            b = Pair()
            b.put_Pair_i(i, "bp")
            self.Xlist[BPndx] += [b]
            if debug_bp:
               print b.disp_Pair()
               # print 'i = %5d, %s' % (i, s)
            mm = i + 1
            self.stack[BPndx]   += 1
            self.pointer[BPndx]  = self.stack[BPndx]
            self.counter[BPndx] += 1
            self.scan_ctcf(mm, layer + 1)
         elif s == ')':
            BPndx = rpr2num[s]
            if debug_bp:
               print "BPndx(%2d) => %s" % (BPndx, s)
               flag_details = False
               if flag_details:
                  print "list of current contacts"
                  for bp in self.Xlist[BPndx]:
                     print bp.disp_Pair()
                  print "------"
               #
            #
            self.pointer[BPndx] = self.find_next_Xpoint_n(self.Xlist[BPndx], self.pointer[BPndx])
            self.Xlist[BPndx][self.pointer[BPndx]-1].put_Pair_j(i, "bp")
            if debug_bp:
               print self.Xlist[BPndx][self.pointer[BPndx]-1].disp_Pair()
               # print 'j = %5d, %s' % (i, s)
            self.pointer[BPndx] -= 1
            self.counter[BPndx] -= 1
            mm = i + 1
            self.scan_ctcf(mm, layer + 1)
         elif s == '[': # PK
            PKndx = lpr2num[s]
            if debug_PK:
               print "PKndx(%2d) => %s" % (PKndx, s)
            #
            mm = i
            b = Pair()
            b.put_Pair_i(i, "pk")
            self.Xlist[PKndx] += [b]
            if debug_PK:
               print b.disp_Pair()
               print 'i = %5d, %s' % (i, s)
            #
            self.stack[PKndx]   += 1 
            self.pointer[PKndx]  = self.stack[PKndx]
            self.counter[PKndx] +=1
            mm = i + 1
            self.scan_ctcf(mm, layer + 1)
         elif s == ']':
            PKndx = rpr2num[s]
            if debug_PK:
               print "PKndx(%2d) => %s" % (PKndx, s)
            #
            self.pointer[PKndx] = self.find_next_Xpoint_n(self.Xlist[PKndx], self.pointer[PKndx], "PK")
            self.Xlist[PKndx][self.pointer[PKndx]-1].put_Pair_j(i, "pk")
            if debug_PK:
               print self.Xlist[PKndx][self.pointer[PKndx]-1].disp_Pair()
               print 'j = %5d, %s' % (i, s)
            self.pointer[PKndx] -= 1
            self.counter[PKndx] -= 1
            mm = i + 1
            self.scan_ctcf(mm, layer + 1)
            
            
            
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
            b.put_Pair_i(i, "ctcf")
            self.Xlist[CTCFndx] += [b]
            if debug_ctcf:
               print b.disp_Pair()
               print 'i = %5d, %s' % (i, s)
            #
            self.stack[CTCFndx]   += 1 
            self.pointer[CTCFndx]  = self.stack[CTCFndx]
            self.counter[CTCFndx] +=1
            mm = i + 1
            self.scan_ctcf(mm, layer + 1)
         elif s == '}':
            CTCFndx = rpr2num[s]
            if debug_ctcf:
               print "CTCFndx(%2d) => %s" % (CTCFndx, s)
            #
            self.pointer[CTCFndx] = self.find_next_Xpoint_n(self.Xlist[CTCFndx], self.pointer[CTCFndx], "CTCF")
            self.Xlist[CTCFndx][self.pointer[CTCFndx]-1].put_Pair_j(i, "ctcf")
            if debug_ctcf:
               print self.Xlist[CTCFndx][self.pointer[CTCFndx]-1].disp_Pair()
               print 'j = %5d, %s' % (i, s)
            self.pointer[CTCFndx] -= 1
            self.counter[CTCFndx] -= 1
            mm = i + 1
            self.scan_ctcf(mm, layer + 1)
            #
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
            self.scan_ctcf(mm, layer + 1)
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
               print "num2lpr[PKfndx[PKrndx[lpr2num[s]]]]: ", num2lpr[PKfndx[PKrndx[lpr2num[s]]]]
               print "search for connecting element:          '%s'" % num2rpr[lpr2num[s]]
               print "element(%2d) => counter(%2d), pointer(%s), stack(%2d)" \
                  % (PKndx, self.counter[PKndx], self.pointer[PKndx], self.stack[PKndx]) 
            #
            mm = i
            b = Pair()
            b.put_Pair_i(i, "pk")
            self.Xlist[PKndx] += [b]
            if debug_PK:
               print b.disp_Pair()
               print 'i = %5d, %s' % (i, s)
            #
            
            self.stack[PKndx]   += 1 
            self.pointer[PKndx]  = self.stack[PKndx]
            self.counter[PKndx] += 1
            mm = i + 1
            self.scan_ctcf(mm, layer + 1)
            
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
               print "num2rpr[PKfndx[PKrndx[rpr2num[s]]]]: ", num2rpr[PKfndx[PKrndx[rpr2num[s]]]]
               print "search for connecting element:          '%s'" % num2rpr[rpr2num[s]]
               print "element(%2d) => counter(%2d), pointer(%s), stack(%2d)" \
                  % (PKndx, self.counter[PKndx], self.pointer[PKndx], self.stack[PKndx]) 
            #
            self.pointer[PKndx] = self.find_next_Xpoint_n(self.Xlist[PKndx], self.pointer[PKndx], "PK", False)
            self.Xlist[PKndx][self.pointer[PKndx]-1].put_Pair_j(i, "pk")
            if debug_PK:
               print self.Xlist[PKndx][self.pointer[PKndx]-1].disp_Pair()
               print 'j = %5d, %s' % (i, s)
            self.pointer[PKndx] -= 1
            self.counter[PKndx] -= 1
            mm = i + 1
            self.scan_ctcf(mm, layer + 1)
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
               mm += 1
               if mm == self.N:
                  break
               s = self.vstr[mm]
               
            self.scan_ctcf(mm, layer + 1)
            #
         #
      #
      else:
         # standard bps (secondary structure)
         self.BPlist = self.Xlist[0]
         
         # PKs
         self.PKlist = []
         pkndx = PKrndx.keys()
         for x_k in pkndx:
            self.PKlist += self.Xlist[x_k]
         #
         self.PKlist = self.assign_PKdirection(self.PKlist)
         self.compress_to_BPlist()
         
         
         # CTCF contacts
         self.CTCFlist = self.Xlist[2]
         for cl in range(0, len(self.CTCFlist)):
            self.CTCFlist[cl].v = '-'
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
               print "         i.e., the order should be %s...%s." % (num2lpr[x_k], num2rpr[x_k])
               print 'input structure:'
               print self.vstr
               sys.exit(1)
            #
         #
      #
      return mm
   #


   #########
   def getKey_i(self, item):
      return item.i
   #
   # !!!!!
   # decide if the stem is parallel or an
   def assign_PKdirection(self, plist):
      debug_ad = False
      plist.sort(key=self.getKey_i)
      for k in range(0, len(plist)-1):
         i_kp0 = plist[k  ].i; j_kp0 = plist[k  ].j
         i_kp1 = plist[k+1].i; j_kp1 = plist[k+1].j
         if i_kp0 <  i_kp1 and i_kp1 < j_kp0 and j_kp0 < j_kp1:
            # the general condition "...(.....[.....).....]..." is
            # satisfied
            test_km1 = False; test_kp2 = False
            if k > 0:
               i_km1 = plist[k-1].i; j_km1 = plist[k-1].j
               #     ii   i       jj      j
               #     kk   k       kk      k
               #     m    p        m      p
               #     10   1       01      1 
               # ..((((...ABCD....))))....abcd...
               if not (i_km1 < i_kp0 and i_kp0 < j_kp1 and j_kp0 < j_km1 and j_km1 < j_kp1):
                  test_km1 = True
               #
            if k + 2 < len(plist):
               #      i   ii         j      jj
               #      k   kk         k      kk
               #          pp                pp
               #      0   12         0      21 
               # ..ABCD...((((....abcd....))))...
               i_kp2 = plist[k+2].i; j_kp2 = plist[k+2].j
               if not (i_kp0 < i_kp1 and i_kp1 < i_kp2 and j_kp0 < j_kp2 and j_kp2 < j_kp1):
                  test_kp2 = True
            if test_km1 or test_kp2:
               plist[k  ].v = 'p'
               plist[k+1].v = 'p'
            #
         #
      #
      return plist
   #
   #########
   
   def compress_to_BPlist(self):
      debug_cbpl = False
      
      # need to build real copies of self.BPlist and self.PKlist
      
      prlist1 = []
      for pl1 in self.BPlist:
         # make sure we don't just make a pointer
         prlist1 += [pl1]
      #
      
      prlist2 = []
      for pl2 in self.PKlist:
         # make sure we don't just make a pointer
         prlist2 += [pl2]
      #
      prlist1.sort(key=self.getKey_i)
      prlist2.sort(key=self.getKey_i)
      
      BPgroup = []
      PKgroup = []
      
      if debug_cbpl:
         print "Before: ------------------"
         print "BPlist:"
         for pl1 in prlist1:
            print pl1.disp_Pair()
         #
         print "PKlist:"
         for pl2 in prlist2:
            print pl2.disp_Pair()
         #
      #
      
      for p1 in range(0, len(prlist1)):
         i_1 = prlist1[p1].i; j_1 = prlist1[p1].j
         if debug_cbpl:
            print "(1) insert ij_1 ", i_1, j_1
         #
         BPgroup += [prlist1[p1]]
         if p1 < len(prlist1)-1:
            if (prlist1[p1+1].i == i_1 + 1) and (prlist1[p1+1].j == j_1 - 1):
               # just find the shortest point on the stem and start from
               # there.
               continue
         #
         p2 = 0
         while p2 < len(prlist2):
            i_2 = prlist2[p2].i; j_2 = prlist2[p2].j; btp_2 = prlist2[p2].v
            if debug_cbpl:
               print "len(prlist2); ", len(prlist2)
               print "i_1(%d) < i_2(%d) < j_1(%d) < j_2(%d)" % (i_1, i_2, j_1, j_2)
            #
            
            # (1)      i_1  i_2   j_1  j_2        (2)      i_2  i_1   j_2  j_1
            #        ..(....[.....)....]..               ..[....(.....]....)..
            
            #if (i_1 < i_2 and i_2 < j_1 and j_1 < j_2) or (i_2 < i_1 and i_1 < j_2 and j_2 < j_1):
            if (i_1 < i_2 and i_2 < j_1 and j_1 < j_2 or btp_2 == 'p'):
               if debug_cbpl:
                  print "(2) remove ij_2 ", i_2, j_2
               #
               PKgroup += [prlist2[p2]]
               del prlist2[p2]
               #p2 -= 1
            else:
               p2 += 1
            #
         #
      #
      BPgroup = prlist1 + prlist2
      BPgroup.sort(key=self.getKey_i)
      PKgroup.sort(key=self.getKey_i)
      del self.BPlist
      del self.PKlist
      
      self.BPlist = BPgroup
      self.PKlist = PKgroup
      
      if debug_cbpl:
         print "After: ------------------"
         print "new BPlist:"
         for gp in self.BPlist:
            print gp.disp_Pair()
         #
         print "new PKlist:"
         for pk in self.PKlist:
            print pk.disp_Pair()
         #
         sys.exit(0)
      #
      return 0
   #
   
   def expand_island(self, CTCFisland):
      CTCFlist = [CTCFisland]
      n_k = len(CTCFisland.contacts)
      # print n_k
      if n_k == 1:
         for l in CTCFisland.contacts:
            b1 = Pair()
            b1.put_Pair(CTCFisland.i, l, "ctcf", '-')
            b2 = Pair()
            b2.put_Pair(l, CTCFisland.j, "ctcf", '-')
            CTCFlist += [b1,b2]
         #
      #
      elif n_k > 1:
         b = Pair()
         i = CTCFisland.i
         j = CTCFisland.contacts[0]
         b.put_Pair(i, j, "ctcf", '-')
         CTCFlist += [b]
         for l in range(1,len(CTCFisland.contacts)):
            for k in range(0,l):
               if l > k + 1:
                  continue
               b = Pair()
               i = CTCFisland.contacts[k]
               j = CTCFisland.contacts[l]
               # print i, j
               b.put_Pair(i, j, "ctcf", '-')
               CTCFlist += [b]
            #
         #
         b = Pair()
         n = len(CTCFisland.contacts)
         i = CTCFisland.contacts[n-1]
         j = CTCFisland.j
         b.put_Pair(i, j, "ctcf", '-')
         CTCFlist += [b]
      #

      return CTCFlist
   #
   
   
   
   # Since the pair elements are derived from the order in which they
   # appear on the stack, after you encounter a "mate", you search
   # back through the stack until you find one unassigned element. So
   # essentially, you first put a bunch of elements (e.g., 'X') on the
   # stack (of course, assuming you have a string of them), now, you
   # have just encountered the partner of this example ('x') and you
   # look at its position on the stack using p (i.e., the
   # pointer[rpr2num['x']] position last incremented). From here of
   # course, you know that the position on the stack must decrease by
   # exactly the number of 'X' items you put on the stack (more to the
   # point, it BETTER match!!!)
   
   # For a more concrete example, suppose you have the structure 
   
   # ".X.XXXXX......x..xxxxx.."
   
   # where element X has the index 27
   
   # we reach the top of the stack with counter[27] = pointer[27] = stack[27] = 6
   
   # self.Xlist[27] should now look like this:
   
   # ( 1, -1)
   # ( 3, -1)
   # ( 4, -1)
   # ( 5, -1)
   # ( 6, -1)
   # ( 7, -1)
   
   # Now we go through this list looking for the -1. In this case, it
   # is at the top of the list, so it is the first element. So we
   # assign pp = 5 (because the list starts from 0 so the last element
   # is 5).
   
   # ( 1, -1)
   # ( 3, -1)
   # ( 4, -1)
   # ( 5, -1)
   # ( 6, -1)
   # ( 7, 14)
   
   # As we find the remaining elements, pointer[27] decreases. So
   # after 6 iterations, we obtain the correct assignment of elements
   # on the stack.
   
   # ( 1, 21)
   # ( 3, 20)
   # ( 4, 19)
   # ( 5, 18)
   # ( 6, 17)
   # ( 7, 14)
   # ".X.XXXXX......x..xxxxx.."
   
   # Personally, I think the looping should not occur, because that
   # means that the stack is out of order, but maybe this search can
   # be employed to detect problems in either building the stack or
   # closing it. So I have left this somewhat incongruous feature
   # within the tool for the moment (since it is still under
   # development). Better safe than sorry as they say.
   
   def find_next_Xpoint_n(self, plist, p, nm = "BP", debug = False):
      if debug:
         print "find_next_Xpoint_n, p = %d, name(%s)" % (p, nm)
         print "plist:"
         for pk in plist:
            print pk.disp_Pair()
      pp = p
      j = plist[p-1].j
      if debug:
         print "pp(%3d), j(%3d)" % (pp, j)
      while j > -1 and pp >= 0:
         pp -= 1
         j = plist[pp-1].j
         if debug:
            print "pp(%3d), j(%3d)" % (pp, j)
      #
      if debug:
         print "result: pp = %d" % pp
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
   

   def print_Xlist_n(self, plist):
      for b in plist:
         print b.disp_Pair()
      #
   #
    
   def convert_Vstruct(self, s):
      self.reset_Vstruct()
      self.set_Vstruct(s)
      self.scan_vs(0)
      self.print_vstr()
      self.print_Xlist_n(self.BPlist)
      return 0
   #
   
   def convert_CTCFstruct(self, s, debug=False):
      self.reset_Vstruct()
      self.set_Vstruct(s)
      self.scan_ctcf(0, 0)
      if debug:
         print "sequence"
         self.print_vstr()
         print "secondary structure"
         self.print_Xlist_n(self.BPlist)
         #
         if len(self.PKlist) > 0:
            print "pk connects"
            self.print_Xlist_n(self.PKlist)
         #
         if len(self.CTCFlist) > 0:
            print "CTCF connects"
            self.print_Xlist_n(self.CTCFlist)
         #
         if len(self.CTCFlist) > 0:
            islands = []
            for cl in self.CTCFlist:
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
      return 0
   #

#


# reads in Vienna package data created by the SimRNA_trafl2pdbs
# processing function
class ViennaData:
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
   vs.convert_Vstruct(ss_seq)
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
   if len(cl) > 1:
      ss_seq = cl[1]
   vs.convert_CTCFstruct(ss_seq, True)
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

