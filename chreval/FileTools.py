#!/usr/bin/python
"""@@@

  Program:        FileTools.py

  Classes:        FileTools

  Functions:      getHeadExt

  Author:         Wayne Dawson
  Version:        0.0
  Creation Date:  140221 (derived from SimRNATools c 140515)
  Last Update:    190409

I eventually found that Python has a lot of useful file manipulation
tools already available, so this is only a marginally useful set of
tools. The main tools that I still use from this module are

getHeadExt() 

and

check_ext(testfile, [ext1, ext2, ...], program_name):

"""

import re
import sys

debug = False # True # 



def getHeadExt(testfile):
    debug_getHeadExt = False
    """@

    output: file_header + _extension
    
    Obviously, this assumes that you don't think multiple dots in your
    private extension definition are ok. If this is not your way of
    thinking, then this will only read the stuff following the last
    dot.

       >>> flnm = "test.txt"
       >>> flhd, ext = getHeadExt(flnm)
       >>> print flhd, ext
       test txt
       
    
    """
    if debug_getHeadExt:
        print "getHeadExt: ", testfile
    #
    splf = testfile.split('.')
    ext = ''
    flhd = testfile
    if len(splf) > 1:
        ext = splf[len(splf)-1]
        flhd = testfile[:len(testfile)-len(ext)-1]
    #
    return flhd, ext
#

#############################
###   class Definitions   ###
#############################

class FileTools:
    def __init__(self):
        self.flnm = ''
        self.flhd = ''
        self.ext  = ''
    #
    
    def check_ext(self, testfile, allowed_exts, program=''):
        """@
        
        Checks the input file <testfile> extension <allowed_exts>.  The
        variable <allowed_exts> can either be a single string (e.g.,
        'trafl') or a list of strings (e.g., ['trafl', 'dGrnk']). This
        version of check_ext() should be the one that is used from now
        on, though the old_check_ext() works basically the same.
        
        Example code:
        
        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        >>> flnm = "file.trafl"
        >>> ft = FileTools()
        >>> file_passed = ft.check_ext(flnm, ['trafl', 'dGrnk'], "myprogram")
        True
        >>> flnm = "file.txt"
        >>> file_passed = ft.check_ext(flnm, ['trafl', 'dGrnk'], "myprogram")
        [program: myprogram]:
        ERROR: file must have one of these extensions '*.trafl' or '*.dGrnk'.
               input file 'file.txt' --> extension: 'txt')
        >>> if not file_passed:
        >>>     print "ERROR: terminated"; sys.exit(1)
        >>> #
        $
        
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        
        output: True/False # True when extension matches expected input file type
        
        """
        
        #
        self.flnm = testfile
        splf = testfile.split('.')
        self.ext = splf[len(splf)-1]
        if len(self.ext) == 1:
            print "ERROR: file (%s) should have one of the following extensions: " \
                % (self.flnm), allowed_exts
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
                #
                xtnsn += '.'
                if len(allowed_exts) == 1:
                    print "ERROR: file must have the extension %s" % xtnsn
                else:
                    print "ERROR: file must have one of these extensions %s" % xtnsn
                #
                print "       input file '%s' --> (extension: '%s')" \
                    % (self.flnm, self.ext)
            else:
                xtnsn = '\'*.%s\'' % allowed_exts
                print "ERROR: file must have the extension %s" % xtnsn
                print "       input file '%s' --> (extension: '%s')" \
                    % (self.flnm, self.ext)
            #
        #
        return flag_pass
    #
    
    # ##########################################################
    # ######  Obsolete code that I cannot seem to delete  ######
    # ##########################################################
    # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

    
    def clean_split(self, st):
        """@
        
        161118wkd: This was evidently an early tool I wrote and I
        don't think this one is so useful now. In fact, if you only
        need to divide out the data items from an initially
        indeterminate number of spaces (' ') in a string (st), then it
        is better to use
        
           >>> st.split()  
        
        Indeed, it is better to write
        
           >>> st.strip().split() 
        
        This operation removes all ' ' characters and the '\n'
        character as well.  So st.strip().split() is a better strategy
        and then, if you want to extract some other things out of it,
        then do the same procedure -- i.e., strip() -- on each item in
        the resulting list.
        
        """
        ssv = st.split(' ')
        vlist = []
        # remove blank spaces
        for ssv_k in ssv:
            if not ssv_k == '':
                vlist += [ssv_k]
            #
        #
        return vlist
    #
    
    def hollerith(self, n, slen):
        """@
        
        generates a hollerith
        
        161118wkd: This does work, and I have used it a lot in the
        past. However, I found this task can also be done using the
        package "string". For example, to produce a hollerith string
        with up to 5 zeros and a positive integer k, you can generate
        this same output with the following command
        
           >>> k = 25
           >>> string.zfill(k, 5) 
           00025
        
        """
        
        sn = str(n)
        hn = ''
        for i in range(0,slen-len(sn)):
            hn += '0'
        #
        hn += sn
        return hn
    #
    
    def old_check_ext(self, testfile, allowed_exts):
        """@ 
        
        This program is now obsolete!
        
        Checks the input file <testfile> extension <extset>.  The
        variable <extset> can either be a single string (e.g.,
        'trafl') or a list of strings (e.g., ['trafl', 'dGrnk']). The
        difference between the old version (here) and the new version
        (above) is simply that I shifted from analyzing with a regular
        expression to just comparing the arguments in the list. It
        seems that the new way works, but I left this old approach
        here anyway. It is slightly more cumbersome. It is basically
        obsolete, and I no longer use it and so I don't tend to update
        it.
        
        """
        
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
                #
            #
            
        else:
            # print 'only one extension found'
            testext='.' + allowed_exts + '$'
            p = re.compile(testext)
            a = p.findall(testfile)
            if (len(a) > 0):
                flag_pass = True
            #
        #
        
        if not flag_pass:
            xtnsn = ''
            if type(allowed_exts) is list:
                xtnsn = '\'*.%s\'' % allowed_exts[0]
                for i in range(1,len(allowed_exts)):
                    xtnsn += ' or \'*.%s\'' % allowed_exts[i]
                #
                xtnsn += '.'
                
                if len(allowed_exts) == 1:
                    print "ERROR: file must have the extension %s" % xtnsn
                else:
                    print "ERROR: file must have one of these extensions %s" % xtnsn
                #
                
            else:
                xtnsn = '\'*.%s\'' % allowed_exts
                print "ERROR: file must have the extension %s" % xtnsn
            #
        #
        return flag_pass
    #
    
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # ##########################################################
    # ######  Obsolete code that I cannot seem to delete  ######
    # ##########################################################
# 
