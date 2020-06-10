#!/usr/bin/env python3

"""@
Program:    Xi_L

Classes:    Xi_Stem

Functions:  

author:        Wayne Dawson
creation date: 2014 (approximate: originally from vsfold cpp package)
last update:   200210 (upgrade to python3)


Purpose: 

Purpose: I was originally trying to develop a new approach with
variable Kuhn length and wanted to be able to compute the Kuhn length
of a given RNA stem based upon the length of the stem. This initiated
the idea to look at the solution for the worm-like-chain model and the
freely-jointed-polymer-chain model. 

For ssRNA stems and for dsRNA/dsDNA, this is used to evaluate the
ACTUAL Kuhn length relative to the MAXIMUM Kuhn length of the
material.  It was necessary to add this because the Kuhn length
appears to be largely linear in character for very short sequences (as
typically found in ssRNA), but it grows to some maximum value (here
xi_m) when the sequences become very long.  So basically, I've
observed that when you have typical ssRNA, the Kuhn length is
proportional to the stem length, but when you have very long
dsRNA/dsDNA strands, the Kuhn length tends toward a maximum (xi_m).


Comments:

   1) The main function that is used in this package is f_wlc, because
      it is arguably the most accurate description of the
      interaction. However, the other stuff is still here for
      comparison.


   2) This requires that a maximum be entered Xi_stem.put_xi_m(x), and
      then Kuhn lengths are measured relative to this input using
      various types of functions: the "worm-like chain" (wlc) model
      (probably the most correct), the "freely-jointed polymer chain"
      (fjpc) model (to some extent correct, but it has been shown not
      to satisfied polymer stretching as well as the wlc), and a
      integrated exponential function resembling at filling of a
      capacitor (which is introduced mainly for comparison).

"""

# global variables and imported objects
debug = False
#debug = True
from math import cosh
from math import sinh
from math import exp
from math import sqrt

import sys


# this is the main meat of this module.
class Xi_stem:
    
    def __init__(self):
        self.debug = debug
        
        self.xi_m  = 150.0  # [bps] default upper bound
        
        """@
        
        In this object, xi_m is the LIMIT on the Kuhn length.  It is
        NOT the length of the chain, it is the MAXIMUM possible Kuhn
        length. Hence, when you are setting this value, you should be
        thinking about dsRNA or dsDNA where you would get the maximum
        length, NOT the length of the stem.
        
        """
        self.xi_fs = 4.0    # nt
        
        """@
        
        In this object, this is the lower bound for the Kuhn length.
        This is because you can chose sequences that are shorter than
        4 bp (for example), but the Kuhn length at this real is just
        equal to the length of the strand because it is just way too
        short.
        
        """
        
        self.w     = 3.0    # [dimensionaless scaling factor for fjpc]
        
        """@
        
        The weight "w" is required for the
        freely-jointed-polymer-chain model because the equation is
        based on a one dimensional solution.  I think that it is, to
        some extent, a mistake to have solved the expression in this
        way, as it is not clear that this is actually legitimate.
        Nevertheless, I see no point in questioning Grosberg and
        Klokov at this point.  It may be that a one dimensional
        solution is allowed and I am just being pedantic and wrong.
        
        Therefore, in the meantime, I will just weight the equation as
        though it was solved one dimensionally and can be expanded to
        three dimensions by multiplying by a factor of three, and be
        done with it.
        
        """
    #
    

    # for testing, no other use.
    def hello(self):
        print ("hello from class Xi_stem")
    #
    
    def xi_r(self, x):
        """@
        
        This calculates the ratio of some stem length divided by the
        maximum Kuhn length.
        
        """
        return float(x)/self.xi_m
    #

    def f_fjpc(self, slen):
        """@
        
        The Freely jointed polymer chain (fjpc) solution.  This is the
        Langevin equation
        [http://en.wikipedia.org/wiki/Langevin_equation], or so it was
        claimed by Grosberg and Klokov at least. The basic form here
        is tanh(x)-1/x.

        """
        wx = self.w*self.xi_r(slen)
        """@
        
        The wx is a bit of an adaptation that I had to resort to
        because it was not really making sense when the one
        dimensional equation was being compared with the three
        dimensional worm-like-chain model over a range of Kuhn lengths
        of dsRNA/dsDNA. 
        
        I do not exactly have a derivation for the "why?" part of it,
        it seems that self.w ~ D, where D is the dimensionality. So if
        this is 3 dimensions, then self.w = 3 (the default value
        presently). I think it can be left in default because, for the
        most part, we don't really rely on the fjpc model for
        parameters, we depend on the wlc.

        """
        coth_wx = cosh(wx)/sinh(wx)
        xi = self.xi_m*(coth_wx - 1.0/wx);
        return xi
    #
    
    def f_exp(self, slen):  
        """@
        
        This _integrated_ exponential function is used mainly for
        comparison with the "freely jointed polymer chain" (fjpc) and
        "worm like chain" (wlc) models.  It resembles the capacitor
        equation where the filling of charge is linear for small
        values and gradually tapers off to some maximum.

        """
        xi = self.xi_m*(1.0 - exp(-self.xi_r(slen)))
        return xi
    #
    
    
    def f_wlc(self, slen):
        """calculates the Kuhn length based on the Worm-like-chain function"""
        if (self.debug == True):
            print ('f_wlc: xi_m = %8.3f nt' % self.xi_m)
        #
        
        """@
        
        Since the persistence length is roughly half the Kuhn length,
        we have to double the value to make it generate the proper
        length.
        
        """
        xr = 2.0*self.xi_r(slen)

        """@
        
        Note that the equation is already adjusted to account for the
        fact that we are using Kuhn length instead of persistence
        length.

        """
        xi = self.xi_m*(1 - (1.0/xr)*(1- exp(-xr)))
        return xi
    #
    
    
    def put_xi_m(self, x):
        """set the maximum Kuhn length"""
        self.xi_m = x
    #
    
    
    def get_xi_m(self):
        """read the maximum Kuhn length"""
        return self.xi_m
    #
    
     
    def put_xi_fs(self, x):
        """set the minimum Kuhn length"""
        self.xi_fs = x
    #
    
    def put_w(self, x):
        """@
        
        Set the weight factor on the "freely jointed polymer chain"
        (fjpc) model.  Basically, if the chain is in three dimensions,
        it should be 3, two dimensions 2 and one dimension 1.

        """
        self.w = x
    #
    
    
    def get_stemlen_from_xi(self, xi):
        """@
        
        xi = xi_m { 1 - (xi_m/(2x) [ 1 - exp(-2x/xi_m) ] }
        
        I determine x by expanding the exponential to 3rd order
        
        xi = x - (2/3)(x**2/xi_m)
        
        rearranging, we get an expression in powers of x
        
        x**2 - 3xi_m/2 + 3xi_m*xi/2 = 0
        
        the solution is:
        
        x = 3xi_m/4 + [(3xi_m/4)**2 - 3*xi_m*xi/2]**(1/2)
        
        """
        b  = 3.0*self.xi_m/4.0
        ac = 3.0*self.xi_m*xi/2.0
        slen = b - sqrt(b**2 - ac)
        return slen
    #
#    

def main(cl):
    print (cl)
    a = Xi_stem()
    a.put_xi_m(150.0)
    for k in range(1, 20):
        xi = a.f_wlc(k)
        print ("%2d     %8.2f" % (k, xi))
    #

    for k in range(1, 20):
        slen = a.get_stemlen_from_xi(float(k))
        print ("%2d     %8.2f" % (k, slen))
    #
    
#
    


if __name__ == '__main__':
    # running the program
    main(sys.argv)
#
