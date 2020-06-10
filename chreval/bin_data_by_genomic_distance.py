#!/usr/bin/env python3



import sys
from collections import OrderedDict
import argparse

from BasicTools import invertDict
from FileTools  import getHeadExt

lbls = {}

lbls.update({"chr"    :  0 }) # chromosome index 1 2 3, X
lbls.update({"begin"  :  1 }) # beginning of the measured segment  
lbls.update({"end"    :  2 }) #"ending of the measured segment      
lbls.update({"ctcf1"  :  3 }) # direction of the first CTCF (L/R)   
lbls.update({"ctcf2"  :  4 }) # direction of the second CTCF (L/R) 
lbls.update({"nPET"   :  5 }) # number of PET counts 
lbls.update({"cmplx1" :  6 }) # complexity (general) 
lbls.update({"cmplx2" :  7 }) # complexity (definition A) 
lbls.update({"cmplx3" :  8 }) # complexity (definition B) 
lbls.update({"len"    :  9 }) # genomic distance
lbls.update({"dG_0"   : 10 }) # free energy of first structure
lbls.update({"TdS_0"  : 11 }) # entropy of first structure
lbls.update({"dGbar"  : 12 }) # average FE per segment dG_0/len
lbls.update({"TdSbar" : 13 }) # average  S per sgement TdS_0/len
lbls.update({"p_max"  : 14 }) # prob of the min FE structure
# similarity and hamming: see Functions
# p_sim    sim  cmplx  ->  scmplx
# p_ham    ham  cmplx  ->  hcmplx
# p_dTdS   dTdS cmplx  ->  ecmplx
# p_ddG    ddG  cmplx  ->  gcmplx
lbls.update({"p_sim"  : 15 }) # similarity index
lbls.update({"scmplx" : 16 }) # numb iterations
lbls.update({"p_ham"  : 17 }) # hamming distance
lbls.update({"hcmplx" : 18 }) # numb iterations
lbls.update({"p_dTdS" : 19 }) # calculate from d(TdS) prob
lbls.update({"ecmplx" : 20 }) # numb iterations
lbls.update({"p_ddG"  : 21 }) # calculate from d(dG) prob
lbls.update({"gcmplx" : 22 }) # numb iterations
lbls.update({"open"   : 23 }) # from ATAC-seq
lbls.update({"active" : 24 }) # from chromHMM and/or segway
lbls.update({"A"      : 25 }) # region A
lbls.update({"B"      : 26 }) # region B


rlbls = invertDict(lbls) # reverse to { 0 : "chr" }


fldInfo = {
    "chr"    : "chromosome index 1 2 3, X",
    "begin"  : "beginning of the measured segment",
    "end"    : "ending of the measured segment",
    "ctcf1"  : "direction of the first CTCF (L/R)",
    "ctcf2"  : "direction of the second CTCF (L/R)",
    "nPET"   : "number of PET counts ",
    "cmplx1" : "complexity (general) ",
    "cmplx2" : "complexity (definition A) ",
    "cmplx3" : "complexity (definition B) ",
    "len"    : "genomic distance",
    "dG_0"   : "free energy of first structure",
    "TdS_0"  : "entropy of first structure",
    "dGbar"  : "average FE per segment dG_0/len",
    "TdSbar" : "average  S per sgement TdS_0/len",
    "p_max"  : "prob of the min FE structure",
    "p_sim"  : "similarity index",
    "scmplx" : "numb iterations of similarity",
    "p_ham"  : "hamming distance",
    "hcmplx" : "numb iterations of hamming",
    "p_dTdS" : "calculate from d(TdS) prob",
    "ecmplx" : "numb iterations of d(TdS)",
    "p_ddG"  : "calculate from d(dG) prob",
    "gcmplx" : "numb iterations of d(dG)",
    "active" : "from chromHMM and/or segway",
    "open"   : "from ATAC-seq",
    "A"      : "region A",
    "B"      : "region B"
}


fldUnits = {
    "chr"    : "none",
    "begin"  : "bp",
    "end"    : "bp",
    "ctcf1"  : "none",
    "ctcf2"  : "none",
    "nPET"   : "counts",
    "cmplx1" : "order",
    "cmplx2" : "order",
    "cmplx3" : "order",
    "len"    : "[/5 kbp]",
    "dG_0"   : "[kcal/mol]",
    "TdS_0"  : "[kcal/mol]",
    "dGbar"  : "[kcal/mol/5kbp]",
    "TdSbar" : "[kcal/mol/5kbp]",
    "p_max"  : "[units]",
    "p_sim"  : "[units]",
    "scmplx" : "[n/a]",
    "p_ham"  : "[units]",
    "hcmplx" : "[n/a]",
    "p_dTdS" : "[units]",
    "ecmplx" : "[n/a]",
    "p_ddG"  : "[units]",
    "gcmplx" : "[n/a]",
    "active" : "[units]",
    "open"   : "[units]",
    "A"      : "[units]",
    "B"      : "[units]"
}

allowed = {"begin"  :  1, 
           "end"    :  2, 
           "nPET"   :  5, 
           "cmplx1" :  6, 
           "cmplx2" :  7, 
           "cmplx3" :  8, 
           "len"    :  9, 
           "dG_0"   : 10, 
           "TdS_0"  : 11, 
           "dGbar"  : 12, 
           "TdSbar" : 13, 
           "p_max"  : 14, 
           "p_sim"  : 15, 
           "scmplx" : 16, 
           "p_ham"  : 17, 
           "hcmplx" : 18, 
           "p_dTdS" : 19, 
           "ecmplx" : 20, 
           "p_ddG"  : 21, 
           "gcmplx" : 22 
}

rallowed = invertDict(allowed)

def showOptions():
    vv = list(rallowed)
    for vvk in vv:
        print ("%-8s   %-50s" % (rallowed[vvk], fldInfo[rallowed[vvk]]))
    #
#

# sort the results from HotSpots
def insSortListByIndex(alist, ndx = 0):
    for i in range(1,len(alist)):    
        j = i
        """@
        
        sort according to the component ndx in "alist". For example,
        suppose alist is of the following arrangement:
        
           > alist += [(i, j, mbl, V)] # (i, j, class MBL, dG)
        
        if we want to sort the free energy -- the third component in
        the list (starting from 0 to 3) -- then we call
        
           > alist = insSortListByIndex(alist, 3) 
        
        and the third component will be sorted
        
        """
         
        while j > 0 and alist[j][ndx] < alist[j-1][ndx]: 
            alist[j], alist[j-1] = alist[j-1], alist[j] # syntactic sugar: swap the items
            j=j-1 
        #|endwhile
        
    #|endfor
    
    return alist
#


        
class GenomicOrderParams(object):
    def __init__(self):
        self.flhd  = ''
        self.ext   = ''
        self.adata = OrderedDict()
        self.lenbar = 0.0
        self.lenmx  = 0
        self.lenmn  = 1e9
        self.total  = 0
        self.binlist = {}
        self.binavgL = []
        self.bins    = {}
        self.binsize = 0
        
    #
    
    def showStats(self):
        s = ("overall length [/5 kbp]:\n")
        s += ("min len:   %4d\n"   % self.lenmn)
        s += ("max len:   %4d\n"   % self.lenmx)
        s += ("avg len:  %8.2f" % self.lenbar)
        return s
    #
    
    def showGDrange(self, ndx):
        v = self.binlist[ndx]
        s = ("genomic distance from: %d to %d /5kbp" % (v[0], v[1]))
        return s
    #
    
    
    def readAnalLoopData(self, flnm):
        try:
            fp = open(flnm, 'r')
        except IOError:
            print ("Cannot open %s" % flnm)
            sys.exit(1)
        #
        
        lfp = fp.readlines()
        fp.close()

        self.flhd, self.ext = getHeadExt(flnm)
        
        kdt = 0
        for k in range(0, len(lfp)):
            slfp = lfp[k].strip()
            if slfp[0] == '#':
                continue
            #
            
            sl = slfp.split()
            vv = sl[0][:3]
            ndx = sl[0][3:]
            if len(ndx) == 1:
                ndx = '_' + ndx
            # print (vv, sl[0][3:])
            tag = vv + ndx  + '_' + sl[1].zfill(9) + '_' + sl[2].zfill(9)
            sl[1] = int(sl[1])
            sl[2] = int(sl[2])
            sl[lbls["nPET"]] = int(sl[lbls["nPET"]][2:])
            # cmplx1     cmplx2     cmplx3
            sl[lbls["cmplx1"]] = int(sl[lbls["cmplx1"]][3:])
            sl[lbls["cmplx2"]] = int(sl[lbls["cmplx2"]][3:])
            sl[lbls["cmplx3"]] = int(sl[lbls["cmplx3"]][3:])
            
            #print (tag, sl[lbls["nPET"]], sl[lbls["cmplx1"]], sl[lbls["cmplx2"]], sl[lbls["cmplx3"]])
            
            # len    
            sl[lbls["len"]]    = int(sl[lbls["len"]])
            self.lenbar += sl[lbls["len"]]
            if sl[lbls["len"]] > self.lenmx:
                self.lenmx  = sl[lbls["len"]]
            #
            
            if sl[lbls["len"]] < self.lenmn:
                self.lenmn  = sl[lbls["len"]]
            #
            
            # dG_0      TdS_0     dGbar     TdSbar   p_max     
            sl[lbls["dG_0"]]   = float(sl[lbls["dG_0"]])
            sl[lbls["TdS_0"]]  = float(sl[lbls["TdS_0"]])
            sl[lbls["dGbar"]]  = float(sl[lbls["dGbar"]])
            sl[lbls["TdSbar"]] = float(sl[lbls["TdSbar"]])
            sl[lbls["p_max"]]  = float(sl[lbls["p_max"]])
            
            
            # p_sim    scmplx
            # p_ham    hcmplx
            # p_dTdS   ecmplx
            # p_ddG    gcmplx
            
            sl[lbls["p_sim"]]  = float(sl[lbls["p_sim"]]);  sl[lbls["scmplx"]]  = float(sl[lbls["scmplx"]])
            sl[lbls["p_ham"]]  = float(sl[lbls["p_ham"]]);  sl[lbls["hcmplx"]]  = float(sl[lbls["hcmplx"]])
            sl[lbls["p_dTdS"]] = float(sl[lbls["p_dTdS"]]); sl[lbls["ecmplx"]]  = float(sl[lbls["ecmplx"]])
            sl[lbls["p_ddG"]]  = float(sl[lbls["p_ddG"]]);  sl[lbls["gcmplx"]]  = float(sl[lbls["gcmplx"]])
            
            # active    open         A         B
            sl[lbls["active"]] = float(sl[lbls["active"]])
            sl[lbls["open"]]   = float(sl[lbls["open"]])
            #print (len(sl))
            if len(sl) > 25:
                sl[lbls["A"]]      = float(sl[lbls["A"]])
                sl[lbls["B"]]      = float(sl[lbls["B"]])
            #
            
            self.adata.update({tag : sl})
            kdt += 1
            """
            if kdt == 10:
                sys.exit(0)
            #
            """
        #|endfor
        self.total  = kdt
        self.lenbar = float(self.lenbar)/float(kdt)
        print (self.showStats())
        
    #

    def binnedData(self, binsize):
        if self.lenmx == 0:
            print ("ERROR: data not set")
            sys.exit(1)
        #
        
        self.binsize = binsize
        # reset all the binning 
        self.binlist = {}
        self.bins = {}
        self.binavgL = []
        
        span = self.lenmx + 20
        dspan = span // binsize
        
        bgn = 0; end = dspan
        for k in range(0, binsize):
            self.binlist.update({k : (bgn, end)})
            self.bins.update({k : [] })
            lavg = 0
            n = 0
            for vv in self.adata.keys():
                #print (self.adata[vv])
                length = self.adata[vv][lbls["len"]]
                if bgn <= length and length < end:
                    lavg += length
                    n += 1
                    self.bins[k] += [self.adata[vv]]
                #
                
            #|endfor
            
            if n > 0:
                lavg = float(lavg)/float(n)
            else:
                lavg = 0
            #
            
            self.binavgL += [lavg]
            bgn += dspan; end += dspan
        #|endfor
    #
    
    
    
    def showBinSliceOpen(self, kslice, option):
        
        if kslice >= len(self.binlist):
            print ("ERROR: bins are between 0 and %d" % len(self.binlist))
            sys.exit(1)
        #
        
        dataset = self.bins[kslice]
        dataset = insSortListByIndex(dataset, lbls["open"])
        
        index = str(kslice).zfill(2)
        oflnm = (self.flhd + "_open_" + option + "_" + index  + ".dat")
        print ("making ", oflnm)
        fp = open(oflnm, 'w')
        s =  ("# open %s\n" % option)
        s += ("# " + self.showGDrange(kslice) + '\n')
        s += ("# " + self.showStats() + '\n')
        s += ("# open         %s\n" % fldInfo["open"])
        s += ("# %-10s   %s\n" % (option, fldInfo[option]))
        
        s += ("\n#open         %-15s     span \n" % option)
        s += ("#             %-15s   [/5kbp]\n" % fldUnits[option])
        
        fp.write(s)
        for k in range(0, len(dataset)):
            dt = dataset[k]
            s = ("%8.4f    %8.2f            %4d" % (dt[lbls["open"]], dt[lbls[option]], dt[lbls["len"]]))
            fp.write(s + '\n')
            #print (s)
        #|endfor
        
        fp.close()
        
    #
    
    
    def showBinSliceActive(self, kslice, option):
        if kslice >= self.binsize:
            print ("ERROR: bins are between 0 and %d" % len(self.binlist))
            sys.exit(1)
        #
        
        dataset = self.bins[kslice]
        dataset = insSortListByIndex(dataset, lbls["active"])
        
        index = str(kslice).zfill(2)
        oflnm = (self.flhd + "_active_" + option + "_" + index  + ".dat")
        print ("making ", oflnm)
        fp = open(oflnm, 'w')
        s =  ("# active %s\n" % option)
        s += ("# " + self.showGDrange(kslice) + '\n')
        s += ("# " + self.showStats() + '\n')
        s += ("# active       %s\n" % fldInfo["active"])
        s += ("# %-10s   %s\n" % (option, fldInfo[option]))
        
        s += ("\n#active       %-15s     span \n" % option)
        s += ("#             %-15s   [/5kbp]\n" % fldUnits[option])
        #print (s)
        fp.write(s)
        
        for k in range(0, len(dataset)):
            dt = dataset[k]
            s = ("%8.4f    %8.2f            %4d" % (dt[lbls["active"]], dt[lbls[option]], dt[lbls["len"]]))
            fp.write(s + '\n')
            #print (s)
                
        #|endfor
        
        fp.close()

    #

#

def main(cl):
    print (cl)
    
    parser = argparse.ArgumentParser()
    
    
    parser.add_argument('-numbins', action='store', default=20,
                        dest='nbins', type=int,
                        help='how many bins to divide the data.')
    
    parser.add_argument('-bins', nargs="+", default=[10],
                        dest='kslices',
                        help='which bins to read and display.')
    
    parser.add_argument('-f', action='store', default="test_loops_results_all_190528.dat",
                        dest='inflnm', type=str,
                        help='which data file to read')
    
    parser.add_argument('-options', nargs="+", default=["dGbar"],
                        dest='options',
                        help='which fields to show.')
    
    parser.add_argument('-showOpts', action='store_true',
                        default=False,
                        dest='showOpts', 
                        help='provides a list of options')
    
    args = parser.parse_args()
    
    flnm = args.inflnm
    bins = []
    for k in args.kslices:
        bins += [int(k)]
    #
    
    nbins = int(args.nbins)
    opts  = args.options
    sOpts = args.showOpts
    if sOpts:
        showOptions()
        sys.exit(0)
    #
    
    """
    print (flnm)
    print (bins)
    print (nbins)
    print (opts)
    print (args.showOpts)
    showOptions()
    sys.exit(0)
    """
    
    gop = GenomicOrderParams()
    gop.readAnalLoopData(flnm)
    
    gop.binnedData(nbins)
    for k in bins:
        for opt_l in opts:
            print ("Active")
            gop.showBinSliceActive(k, opt_l)
            print ("Open")
            gop.showBinSliceOpen(k, opt_l)
        #|endfor

    #|endfor
    
        
#



if __name__ == '__main__':
    # running the program
    main(sys.argv)
#

