#!/usr/bin/env python

import sys
from sys import getsizeof
import numpy as np
from matplotlib import pyplot as plt
from math import log, floor, sqrt
from MatrixTools import MatrixTools
import argparse


# #####################################
# #####  configuration variables  #####
# #####################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
SHOWMAIN = False # for debugging main()
EXT = ["clust","heat", "eheat"] # extension for the input file
PROGRAM = "make_heatmap.py" # name of the program
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


# general usage statement
def usage():
    print "USAGE: %s file" % PROGRAM
    print "       ... where the input file *must* have the extension 'heat'!" 
#


# ######################
# #######  PLOT  #######
# ######################

# summary of workable color schemes that I found. 

# ax.imshow(hm, cmap=plt.cm.ocean_r, interpolation='nearest')

# interpolation = 'nearest': this produces squares instead of
# rounded structures.

# cmap = plt.cm.gray: in this configuration, the minimum (here
# 0.0) is black and the maximum is white. The outcome is very
# difficult to see and surely not a good choice

# gray_r: the 'r' evidently tells the program to "reverse" the
# colors scheme so that 0.0 is white and the maximum is black.

# ocean_r: 'reversed', various shades of blue green from white
# (min) to green (max)

# gist_heat_r: 'reversed', various reds from black (max) to white
# (min) with shades of red in between.

# gist_earth_r: 'reversed', earth tones from black (max) to white
# (min) and shades of red in between. Very similar to heat.

# not really sure at this point if this should be an object, but
# anyway....

class Plot:
    def __init__(self):
        self.debug = False
        self.set_interp = 'nearest'
        self.set_cmap   = plt.cm.gist_heat_r
    #
    
    def plot(self, hm):
    
        print "drawing plot..."
        
        # axes
        fig, ax = plt.subplots()
        
        # show grap with colors varying between white = 0 to black =
        # maximum, with varying degrees of red in between.
        ax.imshow(hm, cmap=plt.cm.gist_heat_r, interpolation='nearest')
        
        ax.set_title('heat map')
        
        # Move left and bottom spines outward by 10 points
        ax.spines['left'].set_position(('outward', 10))
        ax.spines['bottom'].set_position(('outward', 10))
        # Hide the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # Only show ticks on the left and bottom spines
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        
        plt.show()
    #
#

class ClustData:
    def __init__(self, raw_data):
        self.DEBUG         = False # True
        self.use_raw_data  = raw_data
        self.from_Nenski   = False
        self.hm_max = -1000.0
        self.hm_min = 1000.0
        self.wt_range = -1.0
    #
    
    def set_Nenski(self):
        self.from_Nenski = True
    #
    
    def disp_distribution(self, hm):
        N = len(hm)
        weight = {}
        dist = {}
        for j in range(0,N):
            for i in range(0,j):
                if hm[i][j] > 0:
                    d = j - i
                    if dist.has_key(d):
                        dist[d] += 1
                    else:
                        dist.update({d : 1 })
                    #
                    if weight.has_key(hm[i][j]):
                        weight[hm[i][j]] += 1
                    else:
                        weight.update({hm[i][j] : 1})
                    #
                #
            #
        #
        keys = weight.keys()
        keys.sort()
        for ww in keys:
            print "%8.3f   %3d" % (ww, weight[ww])
        #
        keys = dist.keys()
        keys.sort()
        for dd in keys:
            print "%4d        %3d" % (dd, dist[dd])
        #
        return weight
            
        
    def get_data(self, flnm):
        print "getting data from %s" % flnm
        mtools = MatrixTools()
        if self.from_Nenski:
            mtools.set_Nenski()
            
        hm = None
        N = -1
        if self.use_raw_data:
            gmtrx = mtools.read_MatrixFile(flnm, EXT) # EXT = "clust","heat"
        else:
            gmtrx = mtools.read_MatrixFile_wt(flnm, EXT) # EXT = "clust"
        #
        N        = gmtrx.length
        hm       = gmtrx.heatmap
        clusters = gmtrx.clusters
        
        if self.DEBUG:
            print mtools.disp_fmatrix(hm, "heatmap")
        #
        self.hm_max = mtools.hm_max
        self.hm_min = mtools.hm_min
        self.wt_range = mtools.wt_range
        return hm, N
    #
#    


def main(cl):
    if SHOWMAIN:
        print cl
        print "number of args: %d" % (len(cl))
    #
    
    # parse the command line
    parser = argparse.ArgumentParser()
    # set up a positional argument (the first one)
    parser.add_argument("heatmap", default=None,
                        help="file with heatmap data generated from chreval.py")
    # in positional arguments, it seems that the specification
    # "type=file" means that the program reads the file.  Here, I
    # don't much like that, so I have to do this some other way.
    parser.add_argument('-raw', action='store_true', default=False,
                        dest='raw_data',
                        help='Display the raw data without any manipulation.')
    
    parser.add_argument('-Nenski', action='store_true', default=False,
                        dest='nenski_data',
                        help='Display data from Nenski.')
    
    
    
    # assign arguments
    args = parser.parse_args()
    # print args
    flnm = args.heatmap
    raw_data = args.raw_data
    nenski_data = args.nenski_data
    print 'flnm ', flnm
    
    cd = ClustData(raw_data)
    if nenski_data:
        print "Data from Nenski"
        cd.set_Nenski()
    print "get_data:"
    hm, N = cd.get_data(flnm)
    print "plot:"
    pt = Plot()
    pt.plot(hm)
    print "minimum: %8.2f" % (cd.hm_min)
    print "maximum: %8.2f" % (cd.hm_max)
    print "range:   %8.2f" % (cd.wt_range)
    cd.disp_distribution(hm)
    
    
#    

if __name__ == '__main__':
    main(sys.argv)
