#!/usr/bin/env python3

"""@@@

program:       LocalEntropy.py

Classes:       

Functions:     calc_dGloc 
               calc_dGlocal_exact

creation date: 200414
last update:   200414
version:       0.0

FitCommon contains a variety of odds and ends (constants, dictionary
definitions, function tools, etc) that are used by a variety of
unrelated programs.

Unforutnately, with programming, sometimes the "most logical" place to
put things is not the "most practical". Suffice it to say that there
is enough in common between various programs that it seemed "more
practical" to put these items here.

"""

import sys
from collections   import OrderedDict

from BasicTools  import roundoff

# #######################
# Thermodynamic constants
# #######################

from Constants  import kB
from Constants import T37C  # [K] absolute temp at 37 C.
from Constants import T_0C  # [K] absolute temp at  0 C.


from gammafcn   import GammaFcn
from Xi_L       import Xi_stem
from RConstants import xi_max




# At present, there is not much reason to use this at all. The only
# thing that can be done with this object independent of the programs
# that use it is to run a test, and that test is not so interesting
# presently. Nevertheless, maybe it will be useful at some point.

xi_wt = GammaFcn()
x_stem = Xi_stem()  
# Maximum Kuhn length (xi_max) and temperature (T).  It's easier to
# just set this here. Presently, xi_max is hard wired, but this could
# be changed if necessary.
x_stem.put_xi_m(xi_max)

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ################################################################
# ################################################################


# compute the local CLE contribution for a given specific strand length.
def calc_dGloc(slen):
    # Maximum Kuhn length (xi_m) and temperature (T).
    # It's easier to just set this here
    xi_m = 150.0
    T = 37.0+273.15
    x_stem.put_xi_m(xi_m)  # sets the maximum Kuhn length
    xi = x_stem.f_wlc(float(slen))   
    xi_v = xi
    if (xi < 10.0):
        xi_v = 10.0
    #
    
    # print ("xi = %8.2f" % xi_v)
    
    dGxi = kB*T*((2.25)/3.0)*xi_wt.rx_functSmpsn(1., xi_v, 2.0, 1.75, 0.5)
    dGxi = float(slen)/xi_v*dGxi
    return dGxi
#

####################################################################
#####################  FREE ENERGY PARAMETERS ######################
####################################################################

class LocalEntropy(object):
    def __init__(self):
        self.xi2stm  = OrderedDict() # stem length corresponding to Kuhn length xi
        self.dGl_stm = {}            # local entropic contribution to free energy
        self.dGl_fs  = 0.0           # local entropic contribution to free energy
    #
#

xi2stm  = OrderedDict() # stem length corresponding to Kuhn length xi
dGl_stm = {}            # local entropic contribution to free energy
dGl_fs  = 0.0           # local entropic contribution to free energy







# compute the local CLE contribution for a given specific strand length.
def calc_dGlocal_exact(T, slen, delta, gmm, nu, xi_mn = 4.0):
    """@
    
    calculates the local entropy for any delta, gamma, nu and T with a
    minimum Kuhn length defind by xi_mn

        T      = set temperature (in Kelvin)
        slen   = upper limit of integration ("strand length" or 
                 "stem length", depending on the context)
        delta  = exponential weight (gaussian is 2)
        gmm    = self avoiding walk parameter (default 1.75)
        nu     = excluded volume weight (default 0.5)
        xi_mn  = minimum stem Kuhn length 

    """
    
    xi = x_stem.f_wlc(float(slen))
    # the actual Kuhn length of stem of length slen. 
    
    xi_v = xi
    if (xi < xi_mn):
        xi_v = xi_mn
    #
    
    # print ("xi = %8.2f" % xi_v)
    xi_wt.set_tolerance(1.0e-4)
    dGxi = kB*T*((2.25)/3.0)*xi_wt.rx_functSmpsn(1., xi_v, 2.0, 1.75, 0.5)
    dGxi = (float(slen)/xi_v)*dGxi
    return dGxi
#



def setup_dGl_stm_xpress(T):
    
    """@
    
    preweighted local entropy values for default xi_fs = 5 [nt]. This
    data can be obtained by using the following command
    
          > python RFreeEnergy.py test2
    
    Currently, the size of this list seems somewhat extreme; most of
    the time, this should not require calculation of stems that exceed
    20 bp (in effective length). Nevertheless, I have included numbers
    up to a little over 100 to ensure that a single hairpin of 100 bps
    can be fitted with this model. There is no limit ultimately, but
    this should handle most of the stems that will ever be requested
    by any reasonable operator.
    
    I guess if one tries to calculate a stem that is longer than 100
    nt, maybe the program should really calculate the actual Kuhn
    length directly and issue a warning that the data set does not
    contain stems of that length. Alternatively, another way would be
    to continue the table entering values in intervals of 10 bps when
    the effective stem length exceeds 100 and interpolate the data
    inbetween. At any rate, values over 100 have to be calculated the
    long way using calc_dGlocal_exact() or set up from setup_dGl_stm()
    with some input option. 
    
    Since this calculation is very expensive to do directly, and the
    calculation is highly parameter dependent (resting on the default
    value for xi_fs, gmm, and delta), it makes a lot of sense to use
    this approach instead of calculating this every time. (See
    comments on setup_dGl_stm().
    
    It is important to remember that if we specify a maximum Kuhn
    length, it is also possible to specify an upper limit on the Kuhn
    length. Presently, I don't have it working, but it would be
    possible to compute a second dictionary with a dGl_xi that
    extrapolates stem length from the entered xi value.
    
    """
    
    # xi_min = 5.0 [nt] (defined in calc_dGlocal_exact)
    
    #            stem length   dGl_stm
    #               [nt]      [kcal/mol]
    dGl_stm.update({  1  :   0.14791167 })
    dGl_stm.update({  2  :   0.29582333 })
    dGl_stm.update({  3  :   0.44373500 })
    dGl_stm.update({  4  :   0.59164666 })
    dGl_stm.update({  5  :   0.73955833 }) # << end of linear scaling, xi_min =  5 bp
    dGl_stm.update({  6  :   1.02070810 })
    dGl_stm.update({  7  :   1.31904379 })
    dGl_stm.update({  8  :   1.63082171 })
    dGl_stm.update({  9  :   1.95344024 })
    dGl_stm.update({ 10  :   2.28500883 }) # << end of linear scaling, xi_min = 10 bp
    dGl_stm.update({ 11  :   2.62410537 })
    dGl_stm.update({ 12  :   2.96963013 })
    dGl_stm.update({ 13  :   3.32071330 })
    dGl_stm.update({ 14  :   3.67665383 })
    dGl_stm.update({ 15  :   4.03666336 })
    dGl_stm.update({ 16  :   4.40067606 })
    dGl_stm.update({ 17  :   4.76809495 })
    dGl_stm.update({ 18  :   5.13858014 })
    dGl_stm.update({ 19  :   5.51184044 })
    dGl_stm.update({ 20  :   5.88762432 })
    dGl_stm.update({ 21  :   6.26571292 })
    dGl_stm.update({ 22  :   6.64591453 })
    dGl_stm.update({ 23  :   7.02806025 })
    dGl_stm.update({ 24  :   7.41200041 })
    dGl_stm.update({ 25  :   7.79760174 })
    dGl_stm.update({ 26  :   8.18474504 })
    dGl_stm.update({ 27  :   8.57332325 })
    dGl_stm.update({ 28  :   8.96323979 })
    dGl_stm.update({ 29  :   9.35440729 })
    dGl_stm.update({ 30  :   9.74674637 })
    dGl_stm.update({ 31  :  10.14018476 })
    dGl_stm.update({ 32  :  10.53465640 })
    dGl_stm.update({ 33  :  10.93010079 })
    dGl_stm.update({ 34  :  11.32646236 })
    dGl_stm.update({ 35  :  11.72368996 })
    dGl_stm.update({ 36  :  12.12173637 })
    dGl_stm.update({ 37  :  12.52055793 })
    dGl_stm.update({ 38  :  12.92011419 })
    dGl_stm.update({ 39  :  13.32036759 })
    dGl_stm.update({ 40  :  13.72128320 })
    dGl_stm.update({ 41  :  14.12282847 })
    dGl_stm.update({ 42  :  14.52497303 })
    dGl_stm.update({ 43  :  14.92768850 })
    dGl_stm.update({ 44  :  15.33094829 })
    dGl_stm.update({ 45  :  15.73472751 })
    dGl_stm.update({ 46  :  16.13900276 })
    dGl_stm.update({ 47  :  16.54375207 })
    dGl_stm.update({ 48  :  16.94895477 })
    dGl_stm.update({ 49  :  17.35459136 })
    dGl_stm.update({ 50  :  17.76064348 })
    dGl_stm.update({ 51  :  18.16709377 })
    dGl_stm.update({ 52  :  18.57392583 })
    dGl_stm.update({ 53  :  18.98112414 })
    dGl_stm.update({ 54  :  19.38867400 })
    dGl_stm.update({ 55  :  19.79656149 })
    dGl_stm.update({ 56  :  20.20477338 })
    dGl_stm.update({ 57  :  20.61329712 })
    dGl_stm.update({ 58  :  21.02212078 })
    dGl_stm.update({ 59  :  21.43123302 })
    dGl_stm.update({ 60  :  21.84062304 })
    dGl_stm.update({ 61  :  22.25028054 })
    dGl_stm.update({ 62  :  22.66019571 })
    dGl_stm.update({ 63  :  23.07035921 })
    dGl_stm.update({ 64  :  23.48076211 })
    dGl_stm.update({ 65  :  23.89139586 })
    dGl_stm.update({ 66  :  24.30225233 })
    dGl_stm.update({ 67  :  24.71332371 })
    dGl_stm.update({ 68  :  25.12460253 })
    dGl_stm.update({ 69  :  25.53608167 })
    dGl_stm.update({ 70  :  25.94775427 })
    dGl_stm.update({ 71  :  26.35961378 })
    dGl_stm.update({ 72  :  26.77165390 })
    dGl_stm.update({ 73  :  27.18386860 })
    dGl_stm.update({ 74  :  27.59625209 })
    dGl_stm.update({ 75  :  28.00879880 })
    dGl_stm.update({ 76  :  28.42150339 })
    dGl_stm.update({ 77  :  28.83436071 })
    dGl_stm.update({ 78  :  29.24736584 })
    dGl_stm.update({ 79  :  29.66051400 })
    dGl_stm.update({ 80  :  30.07380064 })
    dGl_stm.update({ 81  :  30.48722133 })
    dGl_stm.update({ 82  :  30.90077185 })
    dGl_stm.update({ 83  :  31.31444809 })
    dGl_stm.update({ 84  :  31.72824612 })
    dGl_stm.update({ 85  :  32.14216214 })
    dGl_stm.update({ 86  :  32.55619248 })
    dGl_stm.update({ 87  :  32.97033360 })
    dGl_stm.update({ 88  :  33.38458209 })
    dGl_stm.update({ 89  :  33.79893465 })
    dGl_stm.update({ 90  :  34.21338811 })
    dGl_stm.update({ 91  :  34.62793937 })
    dGl_stm.update({ 92  :  35.04258548 })
    dGl_stm.update({ 93  :  35.45732354 })
    dGl_stm.update({ 94  :  35.87215080 })
    dGl_stm.update({ 95  :  36.28706455 })
    dGl_stm.update({ 96  :  36.70206220 })
    dGl_stm.update({ 97  :  37.11714123 })
    dGl_stm.update({ 98  :  37.53229921 })
    dGl_stm.update({ 99  :  37.94753377 })
    dGl_stm.update({100  :  38.36284263 })
    dGl_stm.update({101  :  38.77822358 })
    dGl_stm.update({102  :  39.19367448 })
    dGl_stm.update({103  :  39.60919324 })
    dGl_stm.update({104  :  40.02477785 })
    dGl_stm.update({105  :  40.44042637 })
    dGl_stm.update({106  :  40.85613690 })
    dGl_stm.update({107  :  41.27190759 })
    dGl_stm.update({108  :  41.68773667 })
    
    # xi_min = 10.0 [nt] (defined in calc_dGlocal_exact)
    """
    dGl_stm.update({ 1 :     0.23339035})
    dGl_stm.update({ 2 :     0.46678069})
    dGl_stm.update({ 3 :     0.70017104})
    dGl_stm.update({ 4 :     0.93356138})
    dGl_stm.update({ 5 :     1.16695173})
    dGl_stm.update({ 6 :     1.40034207})
    dGl_stm.update({ 7 :     1.63373242})
    dGl_stm.update({ 8 :     1.86712276})
    dGl_stm.update({ 9 :     2.10051311})
    dGl_stm.update({10 :     2.33390345}) << end of linear scaling, xi_min = 10 bp
    dGl_stm.update({11 :     2.62410537}) 
    dGl_stm.update({12 :     2.96963013})
    dGl_stm.update({13 :     3.32071330})
    dGl_stm.update({14 :     3.67665383})
    dGl_stm.update({15 :     4.03666336})
    dGl_stm.update({16 :     4.40067606})
    dGl_stm.update({17 :     4.76809495})
    dGl_stm.update({18 :     5.13858014})
    dGl_stm.update({19 :     5.51184044})
    dGl_stm.update({20 :     5.88762432})
    dGl_stm.update({21 :     6.26571292})
    dGl_stm.update({22 :     6.64591453})
    dGl_stm.update({23 :     7.02806025})
    """
    
    for kl in dGl_stm.keys():
        # rescaling for different temperature
        dGl_stm[kl] *= (T/T37C)
    #
    
    return dGl_stm
#







def setup_dGl_stm(T, xi_mn, stem_len_max):
    """@
    
    Setup program for pre-calculation of the local entropy. This is
    integration and is very costly, so it is better to calculate a
    list of values and interpolate between them. When non-Gaussian
    parameters are requested (delta != 2 and/or nu != 1/2), this
    program must be called to determine the proper values for the
    local entropy.
    
    In general, this is what should be used, and if a larger maximum
    Kuhn length is required, it should be set up with such values
    using setup_dGl_stm_xpress. This presently requires an independent
    calculation to extend the list and, if xi_mn is changed, it must
    also be included in the calculation at 310.15 K. Perhaps I should
    go to 100 bp for idiot calculations where one tests such numbers
    because it will crash if the program does not detect a value that
    it knows.  
    
    Presently, I have selected 24 bp and a minimum of 5 bp (where
    everything less than 5 bp is scaled linearly) because it makes the
    most sense from an aesthetic point of view. Most RNA that we
    encounter in real life probably is not a lot more than a Kuhn
    length of 12 bp or so because of the breaks in continuity and
    effective continuity, but it is true that a stem that is
    contiguous over 100 bps should have access to Kuhn lengths
    corresponding to that.
    
    """
    global dGl_stm
    
    
    if xi_mn < 4.0:
        print ("ERROR(setup_dGl_stm): set value for free strand ")
        print ("Kuhn length (%8.3f) is too small (should be >= 4.0 [nt].")
        sys.exit(1)
    #
    
    print ("setting up the local variable stem Kuhn length")
    
    delta = 2.0
    gmm = 1.75
    nu = 0.5
    
    # calculate stem weights
    for k in range(1, int(stem_len_max + 0.5)):
        # integer just in case someone includes an idiotic value
        
        # print (k)
        length = float(k)
        if length < xi_mn:
            length = xi_mn
        #
        
        """@
        
        calc_dGlocal_exact(T, length, delta, gmm, nu, xi_mn)
        
        T      = temperature (in Kelvin)
        length = upper limit of integration (length)
        delta  = exponential weight (gaussian is 2)
        gmm    = self avoiding walk parameter (default 1.75)
        nu     = excluded volume weight (default 0.5)
        xi_mn  = minimum stem Kuhn length 
        
        This tool should not be necessary to use in most reasonable
        cases. The only issue presently is for cases where the Kuhn
        length can be longer than 24 nt or if one requires a different
        value for the minimum Kuhn length. For most real RNA cases,
        this would not be likely.
        
        Anyway, one problem that this tool presently has is that it
        does not consider temperature and, because it is expensive to
        do these integral calculations, it also relies on already
        fitted parameters and these values must be added manually.
        
        """
        dG = calc_dGlocal_exact(T, length, delta, gmm, nu, xi_mn)
        if float(k) < xi_mn:
            dG *= (float(k)/xi_mn)
        #
        
        dGl_stm.update({k : dG})
    #
    
    
    print ("finished building local entropy weights")
    return dGl_stm
#


def setup_dGl_fs_xpress(T):
    
    """@
    
    preweighted local entropy values for default xi_fs = 5 [nt]. This
    data can be obtained by using the following command
    
          > python RFreeEnergy.py test2
    
    To obtain other values we would have to use the long way of
    setup_dGl_stm, but presently, this too sits on the default
    value for xi_fs, gmm, and delta, so it makes a lot of sense to use
    this one presently. (See comments on setup_dGl_stm().)

    """
    # xi_min = 4.0 [nt] (defined in calc_dGlocal_exact)
    # 0.59164666
    dGl_fs =  0.49023745    # 0.59164666 # [kcal/mol]
    dGl_fs *= (T/310.15)
    #
    
    return dGl_fs
#

def setup_dGl_fs(T, xi_fs):
    """@
    
    Setup program for pre-calculation of the local entropy. This is
    integration and is very costly, so it is better to calculate a
    list of values and interpolate between them. When non-Gaussian
    parameters are requested (delta != 2 and/or nu != 1/2), this
    program must be called to determine the proper values for the
    local entropy.
    
    """
    global dGl_fs
    
    if xi_fs < 4.0:
        print ("ERROR(setup_dGl_fs): set value for free strand ")
        print ("Kuhn length (%8.3f) is too small (should be >= 4.0 [nt].")
        sys.exit(1)
    #
    print ("setting up the local free strand Kuhn length")
    delta = 2.0
    gmm = 1.75
    nu = 0.5
    
    # calculate free strand weights for a segment of length xi_fs
    dGl_fs = calc_dGlocal_exact(T, xi_fs, delta, gmm, nu, xi_fs)
    
    print ("finished building local entropy weights")
    return dGl_fs
#

def interpolate_dGl_stm(length):
    global dGl_stm
    l_floor = int(length)
    l_plus1 = l_floor + 1
    # print (l_floor, l_plus1)
    dl = length - float(l_floor)
    ddG = dGl_stm[l_plus1] - dGl_stm[l_floor]
    # print (ddG)
    dG = dGl_stm[l_floor] + ddG*dl
    # print (dG)
    return dG
#


def interpolate_dGl_fs(length):
    global dGl_fs
    global xi_fs
    return float(length)*dGl_fs/xi_fs
#



def setup_xi2stm():
    
    for slen_k in range(1, 40):
        xi = x_stem.f_wlc(float(slen_k))
        xi2stm.update({  int(100*roundoff(xi, 2)) :   float(slen_k) })
    #
    return xi2stm
#

def interpolate_xi2stm(xi):
    global xi2stm
    debug_interpolate_xi2stm = False # True # 
    """@
    
    Unfortunately, it is not so easy to find the inverse from the
    dictionary, and I basically concluded after putting this together
    that the estimation function x_stem.xxxx(xi) achieves about the
    same accuracy as this interpolation function.

    Anyway, this is another way to get the stem length from an input
    value of xi with similar accuracy as the theoretical estimation
    carried out to second order.
    
    """
    
    xi_cbp = int(100.0*roundoff(xi, 2))
    xi100_keys = list(xi2stm.keys())
    # It doesn't look like this would be a problem if we wrote
    #
    #     xi100_keys = xi2stm.keys().
    #
    # However, in python3, this might still be better to expressed as
    # an actual list, as that was the intent in the python2 version.
    xi100k_init = 0
    xi100k_next = 0
    dxik_wt = 0.0

    flag_found = False
    
    for xi100k_new in xi100_keys:
        if debug_interpolate_xi2stm:
            print ("xi100k_init(%5d) <= xi_cbp(%5d) < xi100k_new(%5d)" \
                % (xi100k_init, xi_cbp, xi100k_new))
        #

        if xi100k_init <= xi_cbp and xi_cbp < xi100k_new:
            xi100k_next = xi100k_new
            dxik_wt = float(xi_cbp - xi100k_init)/float(xi100k_next - xi100k_init)
            flag_found = True
            break
        else:
            xi100k_init = xi100k_new
        #
    #
    
    # print (xi100k_init, xi100k_next, dxik_wt)
    
    if flag_found:
        stm_floor = float(xi2stm[xi100k_init])
        stm_plus1 = float(xi2stm[xi100k_next])
        dstm = float(stm_plus1 - stm_floor)
        slen = stm_floor + dstm*dxik_wt
        if debug_interpolate_xi2stm:
            print ("stm_floor(%6.1f) stm_plus1(%6.1f) slen(%8.2f)" \
                % (stm_floor, stm_plus1, dstm))
        #
    else:
        n = len(xi100_keys)-1
        slen = float(xi2stm[xi100k_init])
    #
    return slen
#







