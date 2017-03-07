#!/usr/bin/python
#   Program:        FileTools
#   Author:         Wayne Dawson
#   Version:        0.0
#   Creation Date:  140221 (derived from SimRNATools c 140515)
#   Last Update:    161118

# the main tool that I use from this is 

import re
import sys

debug = False
# debug = True

# Hopefully some tools that can be useful check_ext(testfile, extset, program=''):

#############################
###   class Definitions   ###
#############################

class FileTools:
   def __init__(self):
      self.flnm = ''
      self.flhd = ''
      self.ext  = ''
   #
   
   
   
   # 161118wkd: This was evidently an early tool I wrote and I don't
   # think this one is so useful now. In fact, if you only need to
   # divide out the data items from an initially indeterminate number
   # of spaces (' ') in a string (st), then it is better to use
   # st.split().  Indeed, it is better to write
   # "st.strip().split()". This operation removes all ' ' characters
   # and the '\n' character as well.  So st.strip().split() is a
   # better strategy and then, if you want to extract some other
   # things out of it, then do so in the resulting list.
   def clean_split(self, st):
      ssv = st.split(' ')
      vlist = []
      # remove blank spaces
      for ssv_k in ssv:
         if not ssv_k == '':
            vlist += [ssv_k]
      return vlist
   #
   
   
   # generates a hollerith
   
   # 161118wkd: This does work, and I have used it a lot in the
   # past. However, I found this can also be done by using the package
   # "string" and, for a hollerith string of up to 5 zeros and a
   # positive integer k, you can generate this with the command
   # "string.zfill(k, 5)". Hence, for k = 25, string.zfile(k=25,5)
   # leads to "00025".
   def hollerith(self, n, slen):
      sn = str(n)
      hn = ''
      for i in range(0,slen-len(sn)):
         hn += '0'
      hn += sn
      return hn
   #
   
   # checks the input file <testfile> extension <extset>.  The
   # variable <extset> can either be a single string (e.g., 'trafl')
   # or a list of strings (e.g., ['trafl', 'dGrnk']).
   def old_check_ext(self, testfile, allowed_exts):
      #
      self.flnm = testfile
      ss = testfile.split('.')
      self.ext = ss[len(ss)-1]
      if len(self.ext) == 1:
         print "ERROR: file (%s) requires an extension: " % (self.flnm), allowed_exts
         sys.exit(1)
      #
      self.flhd = testfile[:len(testfile)-len(self.ext)-1]
      

      flag_pass = False
      if type(allowed_exts) is list:
         # print 'more than one extension found'
         for ext in allowed_exts:
            testext='.' + ext + '$'
            p = re.compile(testext)
            a = p.findall(testfile)
            if (len(a) > 0):
               flag_pass = True
      else:
         # print 'only one extension found'
         testext='.' + allowed_exts + '$'
         p = re.compile(testext)
         a = p.findall(testfile)
         if (len(a) > 0):
            flag_pass = True
      #
      
      if not flag_pass:
         xtnsn = ''
         if type(allowed_exts) is list:
            xtnsn = '\'*.%s\'' % allowed_exts[0]
            for i in range(1,len(allowed_exts)):
               xtnsn += ' or \'*.%s\'' % allowed_exts[i]
            xtnsn += '.'
            if len(allowed_exts) == 1:
               print "ERROR: file must have the extension %s" % xtnsn
            else:
               print "ERROR: file must have one of these extensions %s" % xtnsn
         else:
            xtnsn = '\'*.%s\'' % allowed_exts
            print "ERROR: file must have the extension %s" % xtnsn
      return flag_pass
   #



   # checks the input file <testfile> extension <allowed_exts>.  The
   # variable <allowed_exts> can either be a single string (e.g., 'trafl')
   # or a list of strings (e.g., ['trafl', 'dGrnk']).
   def check_ext(self, testfile, allowed_exts, program=''):
      #
      self.flnm = testfile
      splf = testfile.split('.')
      self.ext = splf[len(splf)-1]
      if len(self.ext) == 1:
         print "ERROR: file (%s) should have one of the following extensions: " % (self.flnm), allowed_exts
         sys.exit(1)
      #
      self.flhd = testfile[:len(testfile)-len(self.ext)-1]
      
      flag_pass = False
      if type(allowed_exts) is list:
         # print 'more than one extension found'
         for ext in allowed_exts:
            if ext == self.ext:
               flag_pass = True
               break
            #
         #
      else:
         if allowed_exts == self.ext:
            flag_pass = True
         #
      #
      
      if not flag_pass:
         xtnsn = ''
         print "[program: %s]:" % program
         if type(allowed_exts) is list: # need to know if it is a list or a string
            xtnsn = '\'*.%s\'' % allowed_exts[0]
            for i in range(1,len(allowed_exts)):
               xtnsn += ' or \'*.%s\'' % allowed_exts[i]
            xtnsn += '.'
            if len(allowed_exts) == 1:
               print "ERROR: file must have the extension %s" % xtnsn
            else:
               print "ERROR: file must have one of these extensions %s" % xtnsn
            #
            print "       input file '%s' --> (extension: '%s')" % (self.flnm, self.ext)
         else:
            xtnsn = '\'*.%s\'' % allowed_exts
            print "ERROR: file must have the extension %s" % xtnsn
            print "       input file '%s' --> (extension: '%s')" % (self.flnm, self.ext)
      return flag_pass
   #

# 
