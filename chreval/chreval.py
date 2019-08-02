#!/usr/bin/env python

"""@@@

Main Module:   chreval.py (CHRomatin EVALuation program for heatmaps)
Classes:       Manager
Author:        Wayne Dawson
creation date: mostly 2016 a little bit in 2017 up to March
last update:   190111
version:       0

chreval.py (CHRomatin EVALuation program for heatmaps)

Chreval is a dynamic programming method for determining the most
probably chromatin structure arrangement as well as the distribution
of chromatin structure arrangements as a function of free energy
inside of loops

    command line example:
    > chreval.py -f chr17_55674644_55751596_res5kb.heat -add_PET_wt 

Let the file chr17_55674644_55751596_res5kb.heat be called
"chrN_x_y_res5kb.heat" for short. The program requires that
"chrN_x_y_res5kb.heat" consist of a matrix of integer values that
indicate the instances of an interaction. The input file
"chrN_x_y_res5kb.heat" must have the extension "heat". The output
consists of a directory with the same name as the input file minus the
extension ".heat"; i.e., the directory name "chrN_x_y_res5kb". This
directory contains (in separate files) all the calculated loop
structures. These files have the extension "DBN" and can be read by
the 3rd party program varna. Additionally, Chevral deposits two
additional files: chrN_x_y_res5kb_BDwt.clust contains a matrix with
the Boltzmann probabilities for different interactions and
chrN_x_y_res5kb_summary.txt contains a shorthand list of the secondary
structures.

Presently, I think the name is rather unoriginal (i.e., "CHRomatin
EVALuate"). Anyway, what is important is that it works, and maybe we
can come up with a better name later.

"""

from math import exp, log
import sys
import os
import string

# main tool objects
# from Functions import KahanSumExp
from GetOpts import GetOpts
from Cluster import Cluster
from Cluster import ClustData

# Free energy parameters
#from FreeEnergy import FreeEnergy

# Motif object representation
from Motif import Motif # a core object of these building programs
from Motif import Link
from Motif import LGroup
from Motif import Map # main map for all the FE and structures

from Calculate import Calculate
from Trace     import Trace
from Boltzmann import Boltzmann


# LThread object representation
from LThread import DispLThread
from LThread import LThread
from LThread import LThread2Vienna

from Vienna2TreeNode import LThread2DotBracket

# Other objects and tools
from FileTools   import getHeadExt
from ChPair      import ChPairData
from ChPair      import LThread2ChPair
from SimRNATools import SimRNAData

# ################################################################
# ######################  Global constants  ######################
# ################################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# for entropy calculations
# constants: weights for the enthalpy terms 

# A hardwired infinity. Somewhat hate this, but it seems there is no
# way around this darn thing.
from Constants import INFINITY

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ################################################################
# ################################################################


# ################################################################
# ###########  global functions/objects and constants  ###########
# ################################################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

# 1. control parameters for the program

PROGRAM      = "chreval.py"  # name of the program

# tests
TEST0        = False  # for testing of objects
TEST1        = False  # for testing of objects
TEST2        = False  # for testing of objects

# debugging options 
SHOWMAIN     = False # Debugging in main()
DEBUG_Trace  = False # Generally for checking Trace()

# 2. special local global function

def usage():
    print "USAGE: %s -f file.heat " % PROGRAM
#
# The function usage() is not used so much now that I have introduced
# GetOpts.py, but it may still be useful at the very beginning of the
# program and possibly in other parts.



# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ################################################################
# ################################################################

####################################################################
######################  MINIMUM FREE ENERGY  #######################
####################################################################


####################################################################
########################  BACK TRACK TRACE  ########################
####################################################################


#####################################################################
################  CALCULATE BOLTZMANN DISTRIBUTION  #################
#####################################################################

class Manager:
    def __init__(self):
        # important data variables
        self.flnm = ''
        self.T    = -1.0
        self.smap = []
        self.dG   = []
        self.N    = -1
        # objects
        self.calc  = None # main calculation routine
        self.trace = None # trace back routines
        self.clust = None # clustering routines
        self.btz   = None # boltzmann calculations
        self.dt    = None # data handling routines
        self.lt2db = None # LThread2DotBracket tool
        # other
        self.debug_Manager = SHOWMAIN
        
    #
    
    def runCalculations(self, CL):
        self.calc = Calculate(CL)
        self.T    = self.calc.T
        self.flnm = self.calc.flnm
        self.N    = self.calc.N
        self.dG, self.smap = self.calc.minFE(self.T)
        if self.debug_Manager:
            print "number of iloop entries: ", len(self.calc.iloop.lg)
            for lgk in self.calc.iloop.lg:
                print lgk.motif[0].base, lgk.motif[0].ctp, lgk.Vij, lgk.motif[0].join
            #
            #print "stop at 0 in runCalculations"; sys.exit(0)
        #
        
        if self.debug_Manager:
            print "disp_fmatrix:"
            # show the input matrix
            print self.calc.fe.disp_fmatrix(self.calc.fe.hv, "enthalpy")
        #
        
        self.trace = Trace(self.calc)
        
        if DEBUG_Trace:
            print "display of all FE values"
            self.trace.disp_allFE(self.smap)
            #print "stop at 1 in runCalculations"; sys.exit(0)
        #
        
        if self.debug_Manager:
            print "\n\n\n"
            print "Results:"
            print "total free energy"
            print "dG(0,%d)[%s] = %10.3f [kcal/mol]" \
                % (self.calc.N-1,
                   self.smap.glink[0][self.calc.N-1].lg[0].motif[0].ctp,
                   self.smap.glink[0][self.calc.N-1].lg[0].Vij)
            print "\n"
            #print "stop at 2 in runCalculations"; sys.exit(0)
        #
        
        dGmin = INFINITY
        
        if self.debug_Manager:
            print "traceback mFE:"
            dGmin = self.calc.fe.traceback_mFE(0, self.calc.N-1, 0, self.debug_Manager)
            print "B: bound; I: I-loop, M: M-loop"
            print "min FE: %8.2f [kcal/mol]" % dGmin
            print string.join(self.calc.opt_ss_seq, '')
            #print "stop at 3 in runCalculations"; sys.exit(0)
        else:
            dGmin = self.calc.fe.traceback_mFE(0, self.calc.N-1, 0, False)
        #
        
        print "dGmin = %8.2f" % dGmin
        print string.join(self.calc.fe.opt_ss_seq, '')
        
        #print "stop at 4 in runCalculations"; sys.exit(0)
        
        flag_filter = True
        self.trace.setup_HotSpot(dGmin, 0.65)
        hs = [(0, self.calc.N-1,
               self.smap.glink[0][self.calc.N-1].lg[0].motif[0].ctp,
               self.smap.glink[0][self.calc.N-1].lg[0].motif[0].Vij)]
        #
        
        # print hs
        self.trace.get_traces_top(hs, 0, 0, flag_filter)
        
        if self.debug_Manager:
            #self.calc.show_smap_xy(0,self.calc.N-1)
            #print "stop at 4 in runCalculations";  sys.exit(0)
            
            print "chreval: vvvvvvv after get_traces_top"
            a = DispLThread(self.calc.N)
            a.disp_LThread(self.trace.lt)
            self.calc.show_smap_xy(0,32)
            #print "stop at 5 in runCalculations, after get_traces_top()"; sys.exit(0)
        #endif
        
        
        self.dt = DispLThread(self.calc.N)
        if self.debug_Manager:
            print "DispLThread"
            self.dt.disp_LThread(self.trace.lt)
        #
        
        print "size of matrix: %d" % self.N
        print "found %d structures:" % len(self.trace.lt)
        
        # new corrections
        self.lt2db = LThread2DotBracket(self.calc.N, self.calc.fe)

        
        self.btz = Boltzmann(self.calc)
        self.btz.calc_Z(self.trace.lt, self.T)
        self.trace.lt = self.btz.set_LThread_p(  self.trace.lt, self.T)
        self.trace.lt = self.btz.set_LThread_TdS(self.trace.lt, self.T)
        
        self.clust = Cluster(self.calc)
        self.clust.clusterlist(self.trace.lt)
        self.clust.cpiflist(self.trace.lt)
        
    #
    
    
    def printResults(self):
        """main tool to print out the results of the calculation"""

        
        # make directory and display and store files
        
        flhd, ext = getHeadExt(self.flnm)
        
        # check if the directory already exists
        try:
            if not os.path.isdir(flhd):
                os.mkdir(flhd) # build a subdirectory to store figures in
            else:
                print "WARNING: the directory '%s'" % flhd
                print "         already exists, overwriting files\n" 
        except OSError:
            print "ERROR: problems making %s" % flhd
            sys.exit(1)
        #
        os.chdir(flhd) # move to that directory
        
        
        # Save information on the secondary structure in long
        # Janusz-Bonieski format.
        ssflnm = flhd + "_summary.txt"
        try:
            fp = open(ssflnm, 'w')
            header = "# summary of structures from %s\n" % flhd
            # thermodynamic parameters
            # entropy parameters from the chreval
            header += "# thermodynamic parameters:\n"
            header += "#   entropy:\n"
            header += "#     local Kuhn length       = %.3g [bps]\n"      \
                      % self.calc.fe.xi      
            header += "#     segment length          = %.3g [bps]\n"      \
                      % self.calc.fe.seg_len
            header += "#     lambda (binding dist)   = %.3g [nts^(-1)]\n" % self.calc.fe.lmbd
            header += "#     gmm (SAW parameter)     = %.3g (no units)\n" % self.calc.fe.gmm
            header += "#     Temperature             = %.3g [K]\n"        % self.T
            # constants: weights for the enthalpy terms 
            header += "#   enthalpy:\n"
            header += "#     input data rescaling wt = %.3g\n" % self.calc.fe.rescale_wt
            header += "#     febase                  = %.3g [kcal/mol]\n" % self.calc.fe.base
            header += "#     feshift                 = %.3g\n" % self.calc.fe.shift
            header += "# statistics:\n"
            header += "#   total number of structures:          %8d\n"  % len(self.trace.lt)
            header += "#   final fraction of structures extracted:    %6.3f\n" \
                      % self.trace.wt_HS
            header += "#   upper limit of the free energy:         %8.2f\n" \
                      % self.trace.V_HS
            header += "#   minimum free energy:                    %8.2f\n" \
                      % self.trace.dGmin
            header += "#   tandem CTCF loop threshold:             %8.2f [kcal/mol]\n" \
                      % self.calc.fe.ctcf_tthresh
            header += "#   convergent CTCF loop threshold:         %8.2f [kcal/mol]\n" \
                      % self.calc.fe.ctcf_cthresh
            header += "# ----\n"
            header += "#\n"
            fp.write(header)
            fp.close()
        except OSError:
            print "ERROR: problems opening %s" % (ssflnm)
            sys.exit(1)
        #
        
        
        k = 1
        is_chromatin = True
        flag_no_header = False
        save_simRNA   = True
        file_results = ''
        prnt_results = ''
        print "output structures:"
        print "number of threads obtained: %d" % len(self.trace.lt)
        for thrds in self.trace.lt:
            
            file_results = "> %s    dG = %8.3f   p = %12.8f\n" \
                           % (string.zfill(k, 5), thrds.dG, thrds.p)

            if 0:
                file_results += self.dt.makeLThreadDotBracket_1b(thrds, 2, is_chromatin)
            else:
                file_results += self.lt2db.getLThread2DotBracket(thrds, 1, is_chromatin)
                # originally self.dt.getLThread2DotBracket(thrds, 1, True)
            #
            
            if k < 5:
                prnt_results += file_results
            #
            if self.calc.p_all_1D:
                outfile = flhd + '_%s.DBN' % string.zfill(k, 5)
                self.dt.printLThreadDotBracket_VARNA(outfile, thrds, is_chromatin)
                hmapflnm = flhd + '_%s.heat' % string.zfill(k, 5)
                self.dt.printLThreadHeatMap(hmapflnm, thrds, flag_no_header)
                
                chpair_flnm     = flhd + '_%s.chpair' % string.zfill(k, 5)
                chsimres_flnm   = flhd + '_%s.simres' % string.zfill(k, 5)
                chdt = LThread2ChPair(thrds, flhd)
                chdt.print_ChPairData(chpair_flnm)
                srdt = SimRNAData()
                srdt.ChPair2SimRes(chdt, ['N~N'])
                srdt.print_SimRNArestraints(chsimres_flnm, "slope", save_simRNA)
                #
            else:
                if k < 50:
                    outfile = flhd + '_%s.DBN' % string.zfill(k, 5)
                    self.dt.printLThreadDotBracket_VARNA(outfile, thrds)
                    hmapflnm = flhd + '_%s.heat' % string.zfill(k, 5)
                    self.dt.printLThreadHeatMap(hmapflnm, thrds, is_chromatin)
                    #
                    chpair_flnm     = flhd + '_%s.chpair' % string.zfill(k, 5)
                    chsimres_flnm   = flhd + '_%s.simres' % string.zfill(k, 5)
                    chdt = LThread2ChPair(thrds, flhd)
                    chdt.print_ChPairData(chpair_flnm)
                    #print "-------"
                    srdt = SimRNAData()
                    srdt.ChPair2SimRes(chdt, ['N~N'])
                    #print "print_SimRNArestraints"
                    srdt.print_SimRNArestraints(chsimres_flnm, "slope", save_simRNA)
                    #print "next"
                    #
                #
            #
            
            # 161025wkd: try moving this step to here. I know this
            # means that the program must open this file each time and
            # deposit the conttents. The reason I do it this way is
            # because a single string of 100k structures could be many
            # million bytes. This may have been the cause of the
            # program crashing with a "killed" statement somewhere
            # between printing out the structures here, and printing
            # out the matrix and the summary file.
            
            # It was a strange bug because I had just successfully
            # calculated a heatmap of almost twice the size of this
            # one. Therefore, it must be the particular quantity of
            # output perhaps.
            
            # If this work around is successful but having to wait is
            # just too annoying, perhaps we can make this transaction
            # once every 10 times around (or something like that). We
            # must first wait and see if this fixes the issue
            try:
                fp = open(ssflnm, 'a')
                fp.write(file_results)
                fp.close()
            except OSError:
                print "ERROR: problems opening %s" % (ssflnm)
                sys.exit(1)
            #
            k += 1
        #

        if len(self.calc.fe.pssbl_ctcf) > 0:
            
            file_results = "strong binding sites: %d\n" % (len(self.calc.fe.all_ctcf))
            for rh in self.calc.fe.all_ctcf.keys():
                file_results += "(%4d, %4d), %10.0f\n" % (rh[0], rh[1], self.calc.fe.all_ctcf[rh])
            prnt_results += file_results
        #
        
        try:
            fp = open(ssflnm, 'a')
            file_results = "\n//\n"
            fp.write(file_results)
            fp.close()
        except OSError:
            print "ERROR: problems opening %s" % (ssflnm)
            sys.exit(1)
        #
        
        # print the results to the terminal
        print prnt_results
        
        
        try:
            pclusters = self.calc.fe.disp_fmatrix(self.clust.clusters, "clusters", False)
            fp = open(flhd + "_BDwt.clust", 'w') # Boltzmann distribution weighted
            fp.write(pclusters)
            fp.close()
        except OSError:
            print "ERROR: problems opening %s" % (flhd + "_BDwt.clust")
            sys.exit(1)
        #
        try:
            pcpif = self.calc.fe.disp_fmatrix(self.clust.cpif, "cpif", False)
            fp = open(flhd + "_BDwt.cpif", 'w') # Boltzmann distribution weighted
            fp.write(pcpif)
            fp.close()
        except OSError:
            print "ERROR: problems opening %s" % (flhd + "_BDwt.cpif")
            sys.exit(1)
        #
    #    
#



####################################################################
####################################################################
####################################################################

# ###########################
# ######  MAIN DRIVER  ######
# ###########################

def reduce_list(t):
    # from
    # http://stackoverflow.com/questions/7961363/removing-duplicates-in-lists
    s = []
    for i in t:
        if i not in s:
            s.append(i)
    return s
#

def test0():
    # Test an algorithm for finding unique elements in a simple list
    
    t = [1, 2, 3, 1, 2, 5, 6, 7, 8]
    print "orig: ", t
    print "rev:  ", reduce_list(t)
#

def test1():

    # Test for verifying and demonstrating the syntax for the
    # sorting algorithm for handling groups of links
    smap = Map(30)
    
    i = 3
    j = 29
    
    link = Link(i, j, float(-5), 'B', 's', [(i,j)])
    smap.glink[i][j].add_link(link)
    
    link = Link(i, j, float(-15), 'I', 's', [(4,26)])
    link.add_Motif(i, j, float(-15), 'I', 's', [(5,27)])
    smap.glink[i][j].add_link(link)
    
    link = Link(i, j, float(-25), 'M', 's', [(4, 6), (8, 12), (15, 25)])
    link.add_Motif(i, j, 'P', 's', float(-25), [(4, 6), (8, 12), (15, 25)])
    smap.glink[i][j].add_link(link)
    
    
    print len(smap.glink[i][j].lg)
    for lgk in smap.glink[i][j].lg:
        print "next structure"
        for lpk in range(0, len(lgk.motif)):
            print lgk.motif[lpk].ctp, lgk.motif[lpk].get_base(), lgk.motif[lpk].get_branches()
        #
    #
    
    print "number of elements: ", len(smap.glink[i][j].lg)
    print "before sorting:"
    for lgk in smap.glink[i][j].lg:
        join = lgk.motif[0].get_branches()
        ctp  = lgk.motif[0].ctp
        Vij  = lgk.Vij
        print "%s    %8.2f: " % (ctp, Vij), join
    #
    smap.mergeSortLinks(smap.glink[i][j].lg)
    print "after sorting:"
    for lgk in smap.glink[i][j].lg:
        join = lgk.motif[0].get_branches()
        ctp  = lgk.motif[0].ctp
        Vij  = lgk.Vij
        print "%s    %8.2f: " % (ctp, Vij), join
    #
    
    
#    

def test2():
    debug_test2 = True
    sys.argv = ["chreval.py", "-f", "test_degeneracy2.heat"]
    print sys.argv
    CL = GetOpts(PROGRAM)
    
    calc = Calculate(CL)
    T    = calc.T
    flnm = calc.flnm
    N    = calc.N
    dG, smap = calc.minFE(T)
    if debug_test2:
        print "number of iloop entries: ", len(calc.iloop.lg)
        for lgk in calc.iloop.lg:
            base = lgk.motif[0].get_base()
            ctp  = lgk.motif[0].ctp
            Vij  = lgk.Vij
            join = lgk.motif[0].get_branches()
            print base, ctp, Vij, join

        #sys.exit(0)
    #
    # Test for verifying and demonstrating the syntax for the
    # sorting algorithm for handling groups of links
    trace = Trace(calc)
    print "traceback mFE:"
    dGmin = trace.traceback_mFE(0, calc.N-1, 0, True)
    print "B: bound; I: I-loop, M: M-loop"
#    

    


def main(cl):
    CL = GetOpts(PROGRAM)
    
    #
    if SHOWMAIN:
        print cl
        print "number of args: %d" % (len(cl))
    #
    
    if len(cl) < 2:
        emsg = "ERROR: too few arguments"
        CL.error(emsg)
    #
    
    manager = Manager()
    manager.runCalculations(CL)
    manager.printResults()
    
#




if __name__ == '__main__':
    if TEST0:
        test0()
    elif TEST1:
        # for running tests of tools
        test1()
    elif TEST2:
        test2()
    else:
        # running the program
        main(sys.argv)
    #
#
