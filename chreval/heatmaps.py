#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""@

190124wkd: This is code provided by Przemek for directly analyzing the
genomic data. In general, I don't really use this so it has little
utility for me, but it is available here in case for some reason I
come across the genomic data.

"""

from __future__ import division
import sys
from sys import getsizeof
import numpy as np
from matplotlib import pyplot as plt
from math import log, floor, sqrt
from scipy import ndimage
import argparse

def round_int(x, prec):
    return prec * ((x + prec/2) // prec)

def wzmocnienie(x, k):
    """Heatmaps amplificaion 
    x - numpy array,
    k - amplification factor"""
    x = (-1)**(floor(k)-1)*(x-1)**k+1
    return x

####
def getKaryotype(fname):
    """returns dictionary with chromosome sizeses e.g.: {'chr13': 115169878, ... } """
    # chr - hs1 1 0 249250621 chr1
    # chr - hs2 2 0 243199373 chr2
    # ....
    # band hs1 p36.33 p36.33 0 2300000 gneg
    # ....
    
    # read lines like above
    # here i[:3] is read from "chr - hs1 1 0 249250621 chr1"
    
    # the operation i.strip().split() produces the following:
    
    #        i[0]   i[1]  i[2]    i[3]   i[4]    i[5]       i[6]
    #  i = ['chr', '-',  'hs1',  '1',   '0',  249250621', 'chr1']
    data = [i.strip().split() for i in open(fname) if i[:3] == 'chr']
    hs = {}
    for i in data:
        # this makes the dictionary {i[6] : int(i[5]) }
        hs[i[6]] = int(i[5])
        # so in the above sample {'chr1' : 249250621 }
    #
    return hs
#

###
def heatmap(data, karirotyp, chromosom, resolution=100000):
    """resolution - number of bp in pixel"""
    bins = np.arange(0,karirotyp[chromosom],resolution)
    hm = np.zeros((len(bins), len(bins), 3))
    n = len(data)
    for k,i in enumerate(data): 
        sys.stdout.write( "\rheatmap main loop: {:.2f} %".format(k/n*100)) # this is a cool output showing the progress
        sys.stdout.flush()
        x = round_int(i[1], resolution)/resolution
        y = round_int(i[5], resolution)/resolution
        hm[x][y][1] += i[6]
    print ()

    m = hm[:,:,1].max()
    print (m)
    hm = hm[:,:,1]/m
    print ("Amplification...")
    hm = wzmocnienie(hm, 200)
    print ("Rotation...")
    hm = ndimage.rotate(hm,-135)
    print ("Trimming... ")
    lx, ly = hm.shape    
    hm = hm[ lx/ 2:, :]
    return hm

def main():
    parser = argparse.ArgumentParser(description="Draw heatmap")
    # the order here defines where the items should be placed on the command line.
    parser.add_argument("karyotype", help="file with karyotype (format as in Circos)")
    parser.add_argument("singletons", help="plik z singletonami")
    parser.add_argument("chromosome", help="np chr22" )
    args = parser.parse_args()
    
    #print (args)
    #print (args.karyotype)
    #print (args.chromosome)
    #print (args.singletons)
    #sys.exit(0)
    
    chromosome = args.chromosome
    print ("Openning karyotype...")
    karyotype = getKaryotype(args.karyotype)
    print ("Read singletons...")
    
    # reads the file; e.g., "chrY.PET1.GM12878.CTCF.singletons_cluster_PET_2_3.txt"

    # for example, in the above file we have
    # chrY	13834733	13834754	chrY	13842661	13842787	1
    # chrY	13833273	13833369	chrY	13841257	13841292	1
    # chrY	59003631	59003695	chrY	59031378	59031492	1

    # this first step reads the contents of the file, where the
    # operation i.strip().split() converts the first row above as
    # follows

    # i = ['chrY', '13834733', '13834754', 'chrY', '13842661', '13842787', '1'], i.e., all strings
    
    singletons = [ i.strip().split() for i in open(args.singletons) ]
    
    # this next step formats elements 1, 2 and 4, 5 to integers.
    singletons = [ [i[0], int(i[1]), int(i[2]), i[3], int(i[4]), int(i[5]), int(i[6]) ] for i in singletons ]
    
    # this final step verifies that the contents of the file have i[0]
    # and i[3] containing the argument 'chrY', according to the above
    # example.
    singletons = [ i for i in singletons if i[0] == chromosome and i[3] == chromosome]

    print ("Heatmap...")
    hm = heatmap(singletons, karyotype, chromosome, 100000)
    fig = plt.figure()
    ax1 = fig.add_subplot('111')
    ax1.imshow(hm, interpolation="none", extent=[0, karyotype[chromosome],karyotype[chromosome]/sqrt(2),0], aspect='auto', origin="upper")
    
    print ("drawing plot...")
    plt.imshow(hm)
    plt.show()
    
if __name__ == '__main__':
    main()
