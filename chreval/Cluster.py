#!/usr/bin/env python3

"""@@@

Main Module:   Cluster.py

Objects:       Cluster
               ClustData

Author:        Wayne Dawson
creation date: parts 2016 (in chreval), made into a separate object 170426
last update:   200210 (upgrade to python3), 170719, 190812
version:       0


Purpose:

This is used to prepare data for heat-maps and other activities like
this. Currently, this program reads or helps generate files with the
extensions "clust","heat", "eheat", and "cpif". Classes in this module
are still used by chreval.py and shuffle_heatmap.py.


Comments:

190719:

    Examining the code, it seems that I have removed any requirement
    these settings make use of the object from class Calculate. Class
    Cluster is used by chreval to generate the heatmap and cpif file,
    but it doesn't seem to require any prior information about heatmap
    files. Therefore, I have now entirely decoupled this from anything
    associated with Calculate. Although Calculate/chreval use this
    package, they don't require any special settings from HeatMapTools
    that were required at an earlier step. In fact, even Clust, which
    is used by Calculate, does not require any settings from
    clust. Therefore, I have decoupled any association of HeatMapTools
    with either of the two packages here.


180709:

   There seem to a variety of data formats that have been created in
   this Chreval package. There is Motif (the most complete), then
   there is LThread/LNode, which is hopelessly less complete and only
   marginally different from Pair. Here I have this ChPair.

   I think this data format was originally developed because Pair or
   LNode were too highly entwined with Calculate to use them for
   writing code to generate 3D structure building information.  The
   remaining utility of this class is for building shuffled data of
   the pairing interactions. The shuffling is used to validate the 3D
   cryo-EM like structure fitting programs. Since this file contains
   only the pairing information, it is fairly easy to to shuffle heat
   maps without changing other information.

"""

from FileTools    import FileTools
from FileTools    import getHeadExt
from BasicTools   import initialize_matrix
from HeatMapTools import HeatMapTools
from LThread      import DispLThread
import sys
import random
import os

# #####################################
# #####  configuration variables  #####
# #####################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
EXTS = ["clust","heat", "eheat", "cpif", "csv"] # extension for the input file
PROGRAM = "Cluster.py" # name of the program
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


class Cluster:
    def __init__(self, calc):
        self.calc = calc
        self.N    = calc.N
        # initialize the free energy matrix
        self.clusters = []
        self.clusters = initialize_matrix(self.clusters, self.N, 0.0)
        self.cpif = []
        self.cpif = initialize_matrix(self.cpif, self.N, 0.0)
        self.debug = False
    #
    
    def clusterlist(self, ltlist):
        
        htools = HeatMapTools()
        wt = htools.normalize_matrix(self.calc.fe.dG, "neg")
        if self.debug:
            print (len(ltlist))
        #
        
        ij_min = 100000.0
        for ltk in ltlist:
            pBoltz = ltk.p
            for tr in ltk.thread:
                
                # ignore linking information
                ctp = tr.ctp
                btp = tr.btp
                if tr.ctp == 'P' or tr.ctp == 'J':
                    # print ("found a P at ", tr.ij_ndx)
                    continue
                #
                
                if ctp == 'K' and (btp == 'bgn' or btp == 'end'):
                    continue
                #
                
                if ctp == 'W' and (btp == 'bgn' or btp == 'end'):
                    continue
                #
                
                if ctp == 'S' and (btp == 'bgn' or btp == 'end'):
                    continue
                #
                
                
                
                v = tr.ij_ndx
                i = v[0]; j = v[1]
                # print (v)
                self.clusters[i][j] += 1.0*pBoltz*wt[i][j] #  ij
                self.clusters[j][i] += 1.0*pBoltz*wt[i][j] #  ji
                if self.clusters[i][j] < ij_min:
                    ij_min = self.clusters[i][j]
                #
                
            #|endfor tr in ltk.thread:
            
        #|endfor tr in ltk.thread:
        
        if ij_min < 0.0:
            shift = - ij_min
            print ("encountered positive entropy values")
            print ("up-shifting free energy map data by ", shift)
            for j in range(0, self.N):
                for i in range(0, self.N):
                    if not self.clusters[i][j] == 0.0:
                        self.clusters[i][j] += shift
                        self.clusters[j][i] += shift
                    #
                #
            #
        #
    #
    
    
    def cpiflist(self, ltlist):
        debug_cpiflist = False
        Nlt = len(ltlist)
        if self.debug:
            print (Nlt)
        #
        
        kk = 0 # structure count
        for ltk in ltlist:
            kk += 1
            if debug_cpiflist:
                # this shows that whole contents of ltk =
                # LThread(). One may need to see all of LThread to
                # understand where (or if) something is wrong.
                disthr = DispLThread(self.N)
                disthr.disp_LThread(ltk)
            #
            
            for tr in ltk.thread:
                # ignore linking information
                ctp = tr.ctp
                btp = tr.btp
                if ctp == 'P' or ctp == 'J':
                    # print ("found a P for ", tr.ij_ndx)
                    continue
                #
                
                if ctp == 'K' and (btp == 'bgn' or btp == 'end'):
                    if debug_cpiflist:
                        # verify that data is handled correctly
                        print (kk, tr.disp_lnode())
                    #
                    
                    continue
                
                elif ctp == 'K':
                    if debug_cpiflist:
                        # verify that data is handled correctly
                        print (kk, tr.disp_lnode())
                    #
                    
                #
                
                if ctp == 'W' and (btp == 'bgn' or btp == 'end'):
                    if debug_cpiflist:
                        # verify that data is handled correctly
                        print (kk, tr.disp_lnode())
                    #
                    
                    continue
                
                elif ctp == 'W':
                    if debug_cpiflist:
                        # verify that data is handled correctly
                        print (kk, tr.disp_lnode())
                    #
                #
                
                if ctp == 'S' and (btp == 'bgn' or btp == 'end'):
                    continue
                #
                
                v = tr.ij_ndx
                i = v[0]; j = v[1]
                # print (v)
                self.cpif[i][j] += 1.0 #  ij
                self.cpif[j][i] += 1.0 #  ji
            #|endfor tr in ltk.thread:
            
            if debug_cpiflist:
                print ("planned exit")
                sys.exit(0)
            #
            
        #|endfor ltk in ltlist:
        
        return 0
        
    #

#


class ClustData:
    """@@@
    
    190719: Looking at this now, I'm not really sure why this had to
    be developed separate from the programs in HeatMapTools;
    particularly the method get_data. Anyway, it is used by
    make_heatmap.py to generate heatmaps.
    
    """
    def __init__(self, use_raw_data):
        self.DEBUG         = False # True
        self.use_raw_data  = use_raw_data
        self.from_Nenski   = False
        self.hm_max = -1000.0
        self.hm_min = 1000.0
        self.wt_range = -1.0
        self.weight = {}
        self.dist   = {}
        self.title = "None"
    #
    
    def set_Nenski(self):
        self.from_Nenski = True
    #
    
    def get_data(self, flnm):
        flhd, ext = getHeadExt(flnm)
        self.title = flhd
        print ("getting data from %s" % flnm)
        htools = HeatMapTools() # default setup GenerateHeatMapTools()
        if self.from_Nenski:
            htools.set_Nenski()
        #
        
        htools.fileformat = ext
        
        
        hm = None
        N = -1
        if ext == "csv":
            gmtrx = htools.read_heatmap(flnm)
        else:
            if self.use_raw_data:
                gmtrx = htools.read_MatrixFile(flnm, EXTS) # EXT = "clust","heat"
            else:
                gmtrx = htools.read_MatrixFile_wt(flnm, EXTS) # EXT = "clust"
            #
            
        #
        
        N        = gmtrx.length
        hm       = gmtrx.heatmap
        clusters = gmtrx.clusters
        
        if self.DEBUG:
            print (htools.disp_fmatrix(hm, "heatmap"))
        #
        self.hm_max = htools.hm_max
        self.hm_min = htools.hm_min
        self.wt_range = htools.wt_range
        return hm, N
    #

    # there could be several different types of heatmaps; e.g.,
    # heatmap, clust, cpif. This provides a general way to build and
    # display all of them.
    def get_distribution(self, hm):
        N = len(hm)
        self.weight = {}
        self.dist = {}
        for j in range(0,N):
            for i in range(0,j):
                if not hm[i][j] == 0:
                    d = j - i
                    if d in self.dist:
                        self.dist[d] += 1
                    else:
                        self.dist.update({d : 1 })
                    #
                    
                    wt = self.wt_range*hm[i][j]
                    if wt in self.weight:
                        self.weight[wt] += 1
                    else:
                        self.weight.update({wt : 1})
                    #
                    
                #
            #|endfor
            
        #|endfor
        
        if self.DEBUG:
            print ("weight: ", len(self.weight))
            print ("dist:   ", len(self.dist))
        #
        
    #
    
    
    
    
    def disp_WeightDistrib(self):
        flag_skip = False
        if self.title == "None":
            flag_skip = True
        #
        
        if len(self.weight) == 0:
            flag_skip = True
        #
        
        s = ''
        if not flag_skip:
            keys = list(self.weight.keys())
            keys.sort()
            
            """@
            
            In python3 keys becomes 
            
            dict_keys([ .... ]). 
            
            As far as I can tell, it looks like the order is according
            the way the items were entered. If we don't need to sort
            the data, then we could just write
            
            for v in keys:
               print v
            #
            
            However, we want the keys sorted, so we will have to work
            around this. The procedure is as follows
            
            keys =  list(self.weight.keys())
            keys.sort()
            
            So we must request the list before we do the sorting.
            
            """
            
            s += "# %s\n" % self.title
            s += "# instances of same matrix element\n"
            s += "#  matrix         number\n"
            s += "# element          of\n"
            s += "#  found        instances\n"
            for ww in keys:
                s += " %8.3f          %3d\n" % (ww, self.weight[ww])
        #
        
        return s
    #
    
    
    def disp_DistDistrib(self):
        flag_skip = False
        if self.title == "None":
            flag_skip = True
        #
        
        if len(self.dist) == 0:
            flag_skip = True
        #
        
        s = ''
        if not flag_skip:
            keys = list(self.dist.keys())
            keys.sort()
            s += "# %s\n" % self.title
            s += "# genomic distance vs counts\n"
            s += "# dist     counts\n"
            for dd in keys:
                s += " %4d        %3d\n" % (dd, self.dist[dd])
            #|endfor
            
        #
        
        return s
    #
        
#    
