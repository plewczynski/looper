#!/usr/bin/env python3

"""@

Program         gammafcn.py

Classes:        GammaFcn

Functions:      main

Author:         Wayne Dawson
Creation date:  ~2014 (approximate period)
Last Update:    200210 upgrade to python3


Purpose: 

Does simple 1D integration using Simpson's rule. 


Comments:

    1) This was originally taken from my vsfold5 package and
       translated into python.

    2) To access and run these scripts

       > from gammafcn import GammaFcn
       > w = GammaFcn()
       > w.hello()   # for example
       > w.rx_functSmpsn(1., 5., 2.0, 1.75, 0.5)

"""
from math import pow
from math import log
from math import exp
import math
import sys
import argparse



class GammaFcn:
    
    def __init__(self):
        self.debug = False # True #
        
        self.set_zeta = False
        self.set_dnw  = False
        
        
        self.zeta  = 1.0   # zeta = [G(g+3/d)/G(g+1/d)]^d/2
        
        self.dnw1  = 1.0   # dnw1 = delta*(1-nu)
        
        self.dnw2  = 0.0   # dnw2 = delta*(2*nu-1)
        
        self.w_bar = 1.0   # w_bar = (1-nu)*(delta*gamma+1)/z
        
        """@
        
          symbol       meaning               default
         g = gamma   self avoiding wt           1.75
         d = delta   exponental wt              2
         n = nu      solvent wt                 0.5
         z = zeta    [G(g+3/d)/G(g+1/d)]^d/2    1
        
        """
        
        self.tol      = 0.0001 # default tolerance
        
        
        """180216(chinese new year, wd): 
        
        About the value used in the tolerance. Integration of these
        functions is rather expensive (even using Simpson's rule), so
        the above tolerance was evidently a compromise between speed
        and accuracy (evidently tested 070911wd based on last
        entry). Certainly, if there is an intrinsic gamma function,
        this function should be used.
        
        """
    #
    
    def hello(self):
        print("hello from GammaFcn")
    #
    
    def set_tolerance(self, f):
        """set the desired tolerance degree of accuracy."""
        try:
            self.tol = float(f)
        except(ValueError):
            print ("ERROR: input \'%s\' is not a float variable")
            sys.exit(1)
        #
    #
            
    # the gamma function
    def gamma(self, x, w): 
        F = exp(-x)*pow(x,w);
        return F;
    #
    
    # the slow local CLE function
    def f_cwt(self, x, delta, gmm, nu):
        """
        delta    d     ==> 2
        gamma    g     ==> 1.75
        nu       n     ==> 1/2
        """
        
        F = 0.0
        """@
        
        NOTE: this function is assumed to be weighted by
        kBT*(g+1/d)/(D xi) in other programs, so this weight is not
        included here!!!!!
        
        """
        if ((0.99999 < x) and (x < 1.00001)):
            F = 0.0
        else:
            if (((1.99999 < delta) and (delta < 2.00001)) 
                and ((0.49999 < nu) and (nu < 0.500001))):
                F = (log(x)/(1.0-x)+1.0);
                """F = log(x)/(1.0-x) + 1.0;"""
            else: 
                F = pow(x,self.dnw2)*(1.0-pow(x,self.dnw1))*self.w_bar
                F = F + log(x)
                F = F/(1.0-x)
                """@
                
                F = dnw1*{log(x)+pow(x,dnw2)*[1-pow(x,dnw)]}/(1.0-x);
                
                The expression does not require omega-bar because the
                sole purpose of that expression is to match the
                boundary conditions so that the second term comes out
                to be 1.
                
                Rhis seems to have some problems for nu != 1/2 and
                delta != 2.
                
                """
        return F;
    #
    
    
    
    # 
    def calc_zeta(self, gmm, delta):
        """zeta = [Gamma(gmm + 3/delta)/Gamma(gmm+ 1/delta)]^delta/2 """
        z = 0.0
        if ((1.99999 < delta) and (delta < 2.00001)):
            z = gmm+0.5;
       
        else:
            z = self.gmm_functSmpsn(gmm+3.0/delta);
            z = z/self.gmm_functSmpsn(gmm+1.0/delta);
            z = pow(z,delta/2.0);
        #
        return z;
    #
    
    
    
    def calc_omega_bar(self, gmm, delta):
        wb = 1;
        if ((1.99999 < delta) and (delta < 2.00001)):
            wb = gmm+0.5;
            
        else:
            if (self.set_zeta == False):
                self.zeta = self.calc_zeta(gmm,delta)
                
            wb = (gmm + 1.0/delta)/self.zeta;
        #
        self.w_bar = wb;
        return wb;
    #
    
    def rx_functSmpsn(self, xi_1, xi_2, delta, gmm, nu):
        
        """@
        
        This is the more recent integration method to calculate the
        local entropy. It uses Simpson's rule, which is maybe not the
        fastest way, but much faster and more reliable than the
        original way that I used before where I used a fixed dx.

        """
        
        
        
        # initial variables
        x = 0.0
        F = 0.0;  Fo = 0.0;
        deltaF = 0.0
        tolF = 0.0
        flagm = False
        k = 0
        
        
        self.zeta = self.calc_zeta(gmm, delta)
        self.set_zeta = True;
        self.w_bar = (gmm + 0.5)/self.zeta;
        if self.debug:
            print ('w_bar %12.6f\n' % self.w_bar)
        #
        if ((delta < 2.0) or (delta > 2.0) or (nu < 0.5) or (nu > 0.5)):
            self.dnw1 = delta*(1.0-nu)
            self.dnw2 = delta*(2.0*nu-1.0)
            if self.debug:
                print ('changed nu to %7.4f\n' % (nu))
            #
        #
        if self.debug:
            print ('nu: %7.4f, ' % nu)
            print ('weights: dnw1 = %10.3f, ' % self.dnw1)
            print ('dnw2 = %10.3f\n' % self.dnw2)
        #
        
        self.set_dnw = True;
        if (xi_1 > xi_2):
            tmp = xi_1
            flagm = True
            xi_1 = xi_2
            xi_2 = tmp
        #
        
        x =  xi_2/xi_1
        count = 2
        dx = (xi_2 - xi_1)/count
        
        S_even = 0.0
        S_edge = self.f_cwt(xi_1, delta, gmm, nu)  
        + self.f_cwt(xi_2, delta, gmm, nu)
        S_odd = self.f_cwt(xi_1+dx, delta, gmm, nu)
        F = (S_edge + 4.0*S_odd)*(dx/3.0)
        
        
        if (F > Fo):
            deltaF = F - Fo
        else:
            deltaF = Fo - F
        #
        
        
        if (F > 0):
            tolF = F*self.tol
        else:
            tolF = -F*self.tol
        #
        
        
        while ((deltaF > tolF) or (dx > 0.1)):
            count = count*2;
            dx = (xi_2 - xi_1)/count;
            Fo = F;
            S_even = S_even + S_odd;
            S_odd = 0.0;
            
            
            
            
            k = 1
            while (k < count/2): 
                x = xi_1 + (2*k-1)*dx
                S_odd  = S_odd + self.f_cwt(x, delta, gmm, nu)
                k = k + 1
            #
            F = (S_edge + 2.0*S_even + 4.0*S_odd)*(dx/3.0);
            
            # ###############
            # exit conditions
            
            if (F > Fo):
                deltaF = F - Fo
            else:
                deltaF = Fo - F
            #
            
            
            if (F > 0):
                tolF = F*self.tol
            else:
                tolF = -F*self.tol
            # ################
        #
        
        if (flagm == True):
            tmp = xi_1
            xi_1 = xi_2;
            xi_2 = tmp;
            F = -F;
        #
        
        if self.debug:
            print ('f_Simpson(%10.3f, %10.3f, %10.3f)\n' % (xi_1, xi_2, F))
        #
        return F;
    #
    
    
    def gmm_functSmpsn(self, z):
        """@
        
        * Compute the gamma function
        * 
        * int_{0}^{\infty} exp( -x ) x^{z-1} dx
        
        """
        x = 0.0;
        
        F = 0.0; Fo = 0.0;
        deltaF = 1000.0;
        tolF = 0.0;
        
        k = 0;
        
        xi_1 = 0.0;
        xi_2 = 10000.;
        w    = z - 1;
        
        x =  xi_2;
        count = 2;
        dx = (xi_2 - xi_1)/count;
        
        S_even = 0.0;
        S_edge = self.gamma(xi_1, w) + self.gamma(xi_2, w);
        S_odd = self.gamma(xi_1+dx, w);
        F = (S_edge + 4.0*S_odd)*(dx/3.0);
        if (F > Fo):
            deltaF = F - Fo
        else:
            deltaF = Fo - F
        #    
        
        if (F > 0):
            tolF = F*self.tol
        else:
            tolF = -F*self.tol
        #
        
        
        while (deltaF > tolF) or (dx > 0.1):
            count = count*2
            dx = (xi_2 - xi_1)/count;
            Fo = F;
            S_even = S_even + S_odd;
            S_odd = 0.0;
            
            
            #for (k = 1; k < count/2; k++) {
            k = 1
            while (k < count/2): 
                x = xi_1 + (2*k-1)*dx
                S_odd = S_odd + self.gamma(x, w)
                k = k + 1
            #
            F = (S_edge + 2.0*S_even + 4.0*S_odd)*(dx/3.0);
            
            # ###############
            # exit conditions
            
            if (F > Fo):
                deltaF = F - Fo
            else:
                deltaF = Fo - F
            #
            
            if (F > 0):
                tolF = F*self.tol
            else:
                tolF = -F*self.tol
            # ################
        #
        
        if self.debug:
            print ('f_Simpson(%10.3f, %10.3f, %10.3f)\n' % (xi_1, xi_2, F))
        return F;
    #
#





def main(cl):
    # print (cl)
    
    parser = argparse.ArgumentParser(description="calculate various local CLE functions")
    
    parser.add_argument("-delta", default=2.0,
                        dest='delta',
                        help="weight on the exponential function")
    parser.add_argument("-gamma", default=1.5,
                        dest='gamma',
                        help="weight function parameter: ----------------------- \
                        $p(r) = C r^{2 gamma} exp( - theta r^delta ) dr$")
    parser.add_argument("-nu", default=0.5,
                        dest='nu',
                        help="exponental weight on the rms end-to-end distance.")
    parser.add_argument("-len", default=10,
                        dest='seqLen',
                        help="the sequence length; for the local CLE function ")
    parser.add_argument("-tolerance", default=1.e-4,
                        dest='tolerance',
                        help="set integration tolerance (very precice values\
                        not recommended, default 1e-4)")
    
    args = parser.parse_args()
    try:
        delta  = float(args.delta)
    except(ValueError):
        print ("entered delta (%s) is not a float" % args.delta)
        sys.exit(1)
    #
    try:
        gmm    = float(args.gamma)
    except(ValueError):
        print ("entered gamma (%s) is not a float" % args.gamma)
        sys.exit(1)
    #
    try:
        length = int(args.seqLen)
    except(ValueError):
        print ("entered length (%s) is not a integer" % args.length)
        sys.exit(1)
    #
    try:
        nu     = float(args.nu)
    except(ValueError):
        print ("entered gamma (%s) is not a float" % args.nu)
        sys.exit(1)
    #
    try:
        tol    = float(args.tolerance)
    except(ValueError):
        print ("entered tolerance (%s) is not a float" % args.tolerance)
        sys.exit(1)
    #
    
    gf = GammaFcn()
    gf.hello()
    gf.set_tolerance(tol)
    
    print ("length of sequence    local CLE")
    w = gf.rx_functSmpsn(1.0, length, delta, gmm, nu)
    print ("   %6d          %12.6f" % (length, w))
    print ("gmmfnct(%9.5f):   %9.6f" % (gmm, gf.gmm_functSmpsn(gmm)))
    print ("zeta:                 %9.6f" % gf.calc_zeta(gmm, delta))
    
#


if __name__=='__main__':
    main(sys.argv)
#

