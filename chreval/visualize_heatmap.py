#!/usr/bin/env python

"""@@@

Main Module:   visualize_heatmap.py 

Classes:       Plot

Functions:     

Author:        Wayne Dawson
creation date: some time in 2016 (in Dariusz' lab)
last update:   190718
version:       0


Purpose: 

Make a _visual_ 2D plot of text data from heatmap files. The format of
the input file should have one of the following extensions: heat,
eheat, csv with the proper corresponding formats.

Make a _visual_ 2D plot of test data of structures calculated using
cluster and cpif files (i.e., the results of chreval calculations).

"""


import sys
from matplotlib import pyplot as plt

# @@ not used!  vvvvvvvv
import numpy as np 

from math import log
from math import floor
from math import sqrt

from  sys import getsizeof
# @@ not used!  ^^^^^^^^

import argparse # for parsing the command line


from FileTools import FileTools
from FileTools import getHeadExt
from Cluster   import ClustData


# #####################################
# #####  configuration variables  #####
# #####################################
# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
SHOWMAIN = False # for debugging main()
EXT = ["clust","heat", "eheat", "cpif", "csv"] # extension for the input file
PROGRAM = "make_heatmap.py" # name of the program
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


# general usage statement
def usage():
    print "USAGE: %s file" % PROGRAM
    print "       ... where the input file *must* have the extension 'heat'!" 
#

"""
# ######################
# #######  PLOT  #######
# ######################

summary of workable color schemes that I found. 

ax.imshow(hm, cmap=plt.cm.ocean_r, interpolation='nearest')

interpolation = 'nearest': this produces squares instead of
rounded structures.

cmap = plt.cm.gray: in this configuration, the minimum (here
0.0) is black and the maximum is white. The outcome is very
difficult to see and surely not a good choice

gray_r: the 'r' evidently tells the program to "reverse" the
colors scheme so that 0.0 is white and the maximum is black.

ocean_r: 'reversed', various shades of blue green from white
(min) to green (max)

gist_heat_r: 'reversed', various reds from black (max) to white
(min) with shades of red in between.

gist_earth_r: 'reversed', earth tones from black (max) to white
(min) and shades of red in between. Very similar to heat.

not really sure at this point if this should be an object, but
anyway....
"""

class Plot(object):
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
    
    parser.add_argument('-HiC', action='store_true', default=False,
                        dest='hic_data',
                        help='This is for HiC data. HiC data typically \
                        looks rather noisy compared to ChIA-PET data, so \
                        this option essentially filters HiC data to look \
                        more like ChIA-PET data.')
    
    parser.add_argument('-show_wts', action='store_true', default=False,
                        dest='show_weights',
                        help='If you want to see summary statistics, \
                        the output file \'[infile].gdist\' contains a \
                        list of genomic distances and the number of \
                        counts found and \'[infile].weight\' lists how \
                        the data is scaled. Sometimes these statistics\
                        are useful, but they start to fill a directory\
                        quickly, so this option is not recommended for \
                        use as the default.')
    
    # assign arguments
    args = parser.parse_args()
    # print args
    flnm = args.heatmap
    raw_data = args.raw_data
    hic_data = args.hic_data
    show_wts    = args.show_weights
    print 'flnm ', flnm
    flhd, ext = getHeadExt(flnm)
    
    cd = ClustData(raw_data)
    if hic_data:
        print "HiC data"
        cd.set_Nenski()
        
        """@@@
        
        180411wkd: It's a long story why this says Nenski rather than
        HiC. My first encountered with HiC data came from data
        generated at Nenski Institute. I didn't understand what was
        going on, because I was accostomed to ChIA-PET data, so, being
        totally uninformed as usual, I named it "Nenski" because I did
        not know what else to do, but it was clear that it was
        different -- very different from the ChIA-PET data I was
        accostomed to. This is why this option with this name appears
        here. Keep in mind that not all data From Nencki is HiC.

        """
    print "get_data:"
    hm, N = cd.get_data(flnm)
    print "plot:"
    pt = Plot()
    pt.plot(hm)
    
    cd.get_distribution(hm)
    sx = '\n'
    sx += "# minimum: %8.2f\n" % (cd.hm_min)
    sx += "# maximum: %8.2f\n" % (cd.hm_max)
    sx += "# range:   %8.2f\n" % (cd.wt_range)

    if show_wts:
        s_weight = cd.disp_WeightDistrib() + sx
        f = flhd + ".weight"
        fp = open(f, 'w')
        fp.write(s_weight)
        fp.close()
        
        s_gdist   = cd.disp_DistDistrib() + sx
        f = flhd + ".gdist"
        fp = open(f, 'w')
        fp.write(s_gdist)
        fp.close()
    #
#    

if __name__ == '__main__':
    main(sys.argv)
#
