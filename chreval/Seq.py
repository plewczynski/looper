#!/usr/bin/env python3

"""@@@

   Program:        Seq
   Author:         Wayne Dawson
   Version:        0.0
   Creation Date:  170804
   Last Update:    200211 (upgraded to python3) 191016


Purpose: an RNA sequence object with a variety of functions to
manipulate it, if desired. The main value in using this interface is
that it ensures that the sequence is actually RNA and not some other
beast with strange stuff on it. Once that information is obtained,
then a few different manipulations can be done such as reversing the
sequence, and getting the complementary sequence (starting from either
the 3' or 5' end,

"""
    
import sys


"""@

191028; it seems that this stuff isn't used here; however, it may be
useful information for some other projects. I think this was
originally used to save some intermediate RNA sequence information and
I used the random number to identify it and the read it. With SimRNA,
there are some hardwired things that I probably had pull some tricks
in order to do what I wanted to do. Eventually, this reduced to this
simple purpose, but these little fragment remained. I keep it here
because I may remember some use for it in the future.

Example: how to make a random number tag using random.

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
import random # random numbers

tag = str(int(random.random()*100000))  # add a random number tag
tmpflnm_head ='temp_' + tag
tmpflnm      = tmpflnm_head + '.trafl'
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

"""

    

debug = False # True # 


class Seq:
    def __init__(self, sq, system = "RNA"):
        self.seq = sq

        if system == "RNA":
            if not self.check_is_RNAseq():
                print ("ERROR: not an RNA sequence")
                print ("       sequence: %s" % self.seq)
                sys.exit(1)
            #
        #
    #
    
    def  __str__(self):
        return self.seq
    #
    
    def __repr__(self):
        return __str__()
    #
    
    
    def setSeq(self, s):
        self.seq = s
    #
    
        
    def reverse(self):
        return ''.join(reversed(self.seq))
    #
    
    def check_is_RNAseq(self):
        """check if sequence satisfies the basic requirement of being
        considered 'RNA'."""
        is_RNAseq = True # innocent until proven guilty
        
        sq = list(self.seq.upper()) # convert string to list
        for k in range(0, len(sq)):
            if (sq[k] == 'A' or
                sq[k] == 'C' or
                sq[k] == 'G' or
                sq[k] == 'U' or
                sq[k] == 'T'):
                # presently, these are the only allowed options, and
                # they are case sensitive!!!
                continue
            else:
                print ("unrecognized symbol (%s) at position %d of seq" % (sq[k], k))
                print ("sequence is not RNA or is formatted incorrectly")
                is_RNAseq = False
                #break; # can complain more than once
            #
        #
        return is_RNAseq
    #
    
    
    
    def getCRNA_3p(self):
        """complementary RNA sequence starting from 3 prime end"""
        sq = list(self.seq.upper()) # convert string to list
        for k in range(0, len(sq)):
            if sq[k] == 'A':
                sq[k] = 'U'
            elif sq[k] == 'C':
                sq[k] = 'G'
            elif sq[k] == 'G':
                sq[k] = 'C'
            elif sq[k] == 'U':
                sq[k] = 'A'
            elif sq[k] == 'T':
                sq[k] = 'A'
            else:
                print ("unrecognized symbol (%s) at position %d of seq" % (sq[k], k))
                print ("sequence is not RNA or is formatted incorrectly")
                sq = ''
                break;
            #
        #
        newseq = ''.join(sq) # compresses sq into a string
        return newseq 
    #
    
    
    def getCRNA_5p(self):
        """complementary RNA sequence starting from 5 prime end"""
        a = self.getCRNA_3p()
        return ''.join(reversed(a))
    #
    
    def getCRNA(self):
        """default choice for complementary RNA"""
        return self.getCRNA_3p()
    #
#

    


def setCRNA(cl):
    #print (len(cl))
    #print (cl)
    if len(cl) > 2:
        print ("only can take one argument (the sequence)")
        sys.exit(1)
    #
    a = Seq(cl[1])
    b = a.getCRNA()
    print (b)
#

# Main
if __name__ == '__main__':
   setCRNA(sys.argv)
#

