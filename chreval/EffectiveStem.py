#!/usr/bin/env python3

"""@@@

program:       EffectiveStem.py

Classes:       EffectiveStem

Functions:     


creation date: 200511
last update:   200511
version:       0.0

EffectiveStem allows sharing the assessment of every stem in a
universal way. There may need to be things like sequence
dependence. Presently, that is not considered.

"""

import sys


class EffectiveStem(object):
    def __init__(self, max_aa_gap = 2, max_pp_gap = 2):
        
        
        
        # secondary structure stem parameters
        self.max_aa_gap   = max_aa_gap # anti-parallel stems
        self.max_pp_gap   = max_pp_gap # parallel stems
        # default parameters are for chromatin

        self.old_approach = False
        # for testing and comparing with the former method of Vienna
        # and Vienna2TreeNode.
        
        self.debug_ef = False
    #
    
    
    
    def is_connected_aaStem(self, slen1, ph1, qh1, slen2, pt2, qt2):
        if self.old_approach:
            # Some older evaluation schemes in Vienna and
            # Vienna2TreeNode were strict and did not permit any
            # effective stem of any kind
            
            return False
        #
        
        """@
        
        *** all input parameters are integers
        
        For the moment, I use my standard L/2 rule but I add the twist
        that both stems should be at least as long as the longest
        interval that separates them
        
            |<-slen1->|      |<-slen2->|
                    ph1 .... pt2     
             _ _ _ _ _ /    \ _ _ _ _ _
            |_|_|_|_|_|      |_|_|_|_|_|
                    qh1\..../qt2
            
        For example
        
        (1)
                        ...      
             _ _ _ _ _ /    \ _ _ 
            |_|_|_|_|_|      |_|_|
                       \..../
        
        pt2-ph1-1 = 3 and qh1-qt2-1 = 4, therefore dn = 4 
        slen1 = 6     and slen2 = 3
        ==> pair_gap = 4
        
        Whereas slen1 >= 4, slen2 = 3 and not (slen2 >=
        pair_gap). Therefore, this is not a connected stem
        
        (2)
                        ...      
             _ _ _ _ _ /    \ _ _ _ _ _
            |_|_|_|_|_|      |_|_|_|_|_|
                       \..../
        
        pt2-ph1-1 = 3 and qh1-qt2-1 = 4, therefore dn = 4 
        slen1 = 6     and slen2 = 6
        ==> pair_gap = 4
        
        Since slen1 >= 4 and slen2 >= 4, this is a connected stem
        
        (3)
        
               slen1
              ->| |<-|<----- slen2 ------->|
                                        
                 _ /\ _ /\ _ /\ _ /\ _ /\ _ 
                |_|  |_|  |_|  |_|  |_|  |_|
                   \/   \/   \/   \/   \/
        
        pt2-ph1-1 = 1 and qh1-qt2-1 = 1, therefore dn = 1 
        slen1 = 2     and slen2 = 10
        ==> pair_gap = 1
        
        Since slen1 >= 1 and slen2 >= 1, this is a connected stem
        
        """
        
        is_connected = False
        dp12 = (pt2 - ph1 - 1)
        dq12 = (qh1 - qt2 - 1)
        pair_gap = 0
        
        # find the largest gap (difference)
        if dp12 < dq12:
            pair_gap = dq12
        else:
            pair_gap = dp12
        #
        
        if self.debug_ef:
            print ("Enter is_connected_aaStem:")
            print ("      slen1(%3d), ph1(%3d), qh1(%3d)" % (slen1, ph1, qh1))
            print ("      slen2(%3d), ph2(%3d), qh2(%3d)" % (slen2, pt2, qt2))
            
            print ("dp12 = %d, dq12 = %d ==>> pair_gap(%d) < max_aa_gap(%d)?" \
                   % (dp12, dq12, pair_gap, self.max_aa_gap))
        #
        
        # test 1: Is the gap itself less that the cutoff?
        if pair_gap > self.max_aa_gap:
            pair_gap = self.max_aa_gap
        else:
            is_connected = True
        #
        
        
        # test 2: if passed 1, then are both respective stem lengths
        # at least as long as the maximum gap itself?
        if slen1 >= pair_gap and slen2 >= pair_gap and is_connected:
            is_connected = True
        else:
            is_connected = False
        #
        
        if self.debug_ef:
            print ("slen1(%d) >= %d && slen2(%d) >= %d" \
                   % (slen1, pair_gap, slen2, pair_gap))
            print ("result: is_connect = ", is_connect)
        #
        
        return is_connected
    #
    
    
    
    
    def is_connected_ppStem(self, slen1, ph1, qh1, slen2, pt2, qt2):
        if self.old_approach:
            # Some older evaluation schemes in Vienna and
            # Vienna2TreeNode were strict and did not permit any
            # effective stem of any kind
            
            return False
        #
        """@
        
        For chromatin, this is basically just an empty function to
        passify the torturous scanning routines in
        Vienna2TreeNode. Presently, this is only used by chromatin, so
        it means that all operations involving connected stems never
        happen because it is immediately resumed that it is not
        connected. Perhaps I need to think about this, but anyway, for
        the moment, I will follow the same plan established for
        is_connected_asStem.
        
        """
        is_connected = False
        dp12 = float(pt2 - ph1 - 1)
        dq12 = float(qt2 - qh1 - 1)
        pair_gap = 0
        
        
        
        # find the largest gap (difference)
        if dp12 < dq12:
            pair_gap = dq12
        else:
            pair_gap = dp12
        #
        
        
        if self.debug_ef:
            print ("Enter is_connected_ppStem:")
            print ("      slen1(%3d), ph1(%3d), qh1(%3d)" % (slen1, ph1, qh1))
            print ("      slen2(%3d), ph2(%3d), qh2(%3d)" % (slen2, pt2, qt2))
            
            print ("dp12 = %d, dq12 = %d ==>> pair_gap(%d) < max_pp_gap(%d)?" \
                   % (dp12, dq12, pair_gap, self.max_pp_gap))
        #
        
        
        # test 1: Is the gap itself less that the cutoff?
        if pair_gap > self.max_pp_gap:
            pair_gap = self.max_pp_gap
        else:
            is_connected = True
        #
        
        # test 2: if passed 1, then are both respective stem lengths
        # at least as long as the maximum gap itself?
        if slen1 >= pair_gap and slen2 >= pair_gap and is_connected:
            is_connected = True
        else:
            is_connected = False
        #
        
        if self.debug_ef:
            print ("slen1(%d) >= %d && slen2(%d) >= %d" \
                   % (slen1, pair_gap, slen2, pair_gap))
            print ("result: is_connect = ", is_connect)
        #
        
        return is_connected
    #
    
#

def test0():
    
    es = EffectiveStem(8, 8)
    # test the connect stem functions
    slen1 = 1; ph1 =  8; qh1 = 46;
    slen2 = 4; pt2 = 10; qt2 = 45
    print (es.is_connected_aaStem(slen1, ph1, qh1, slen2, pt2, qt2))
    
    slen1 = 4; ph1 =  3; qh1 = 46;
    slen2 = 4; pt2 = 10; qt2 = 45
    print (es.is_connected_aaStem(slen1, ph1, qh1, slen2, pt2, qt2))
    
    slen1 = 1; ph1 =  8; qh1 = 43;
    slen2 = 4; pt2 = 10; qt2 = 45
    print (es.is_connected_ppStem(slen1, ph1, qh1, slen2, pt2, qt2))
    
#



def main(cl):
    print (cl)
    test0()
    
#

if __name__ == '__main__':
    # running the program
    main(sys.argv)
#    
