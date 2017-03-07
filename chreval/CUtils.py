#!/usr/bin/env python

#
import sys



class String:
    def __init__(self, span = 4):
        if span > 20:
            print "WARNING: hollerith string can have up to %d zeros" % span
        elif span <= 1:
            print "ERROR: hollerith range should be greater than 1"
            sys.exit(1)
        #
        self.hspan = span
    #
    
    # make a hollerith expression
    def hollerith(self, n):
        try:
            nt = int(n)
        except ValueError:
            print "ERROR(StringU.hollerith): input '%s' variable not an integer" % n
            sys.exit(1)
        #
        if n < 0:
            print "ERROR(StringU.hollerith): integer %d must be greater than or equal to zero." % n
            sys.exit(1)
        #
        d = 0
        v = 10**d
        while nt > 1:
            if 1*v < nt and nt < 10*v:
                break
            nt /= 10
            d+= 1
        # print d
        s = ''
        for k in range(0, (self.hspan-d)):
            s += '0'
        s += '%d' % int(n)
        return s
    #
    
    def h2i(self, h):
        k = 0
        lh = len(h)
        while h[k] == '0' and k < lh:
            k += 1
        #
        s = h[k:]
        try:
            i = int(s)
        except ValueError:
            print "ERROR(StringU.h2i): input hollerith '%s' does not yield an integer" % h
            sys.exit(1)
        #
        return i
    #
    
    
    
    # this is _very_ strict about what is a float
    def isFloat(self, value):
        flag_isFloat = False
        flag_isNumber = False
        s = str(value)
        a = s.strip().split('.')
        # print a, len(a), a[0].isdigit()

        if not len(a) == 0:
            a0 = a[0]
            if a0[0] == '-':
                a0 = a0[1:]
            #
            
            if len(a) == 2:
                if a0.isdigit() and a[1].isdigit():
                    flag_isFloat = True
                    flag_isNumber = True
                #
            elif len(a) == 1 and a0.isdigit():
                flag_isNumber = True
            #
        #
        if not flag_isNumber:
            print "Warning(float): '%s' is not a number" % value
        #
        return flag_isFloat
    #
    
    
    # this is _very_ strict about what is an integer
    def isInt(self, value):
        flag_isInt = False
        flag_isNumber = False
        
        s = str(value)
        a = s.strip().split('.')
        # print a
        if not len(a) == 0:
            a0 = a[0]
            if a0[0] == '-':
                a0 = a0[1:]
            #
            
            if len(a) == 1 and a0.isdigit():
                flag_isInt = True
                flag_isNumber = True
            #
            elif len(a) == 2:
                if a0.isdigit() and a[1].isdigit():
                    flag_isNumber = True
                #
            #
        #
        
        if not flag_isNumber:
            print "Warning(int): '%s' is not a number" % value
        #
        return flag_isInt
    #
    
    
    def field(self, num, length):
        s = str(num)
        if not (self.isFloat(num) or self.isInt(num)):
            s = ''
        #
        
        n = len(s)
        if n < length:
            for i in range(n, length):
                s += ' '
            #
        #
        return s
    #
            
#
