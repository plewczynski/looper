#!/usr/bin/env python3

"""@@@

Main Module:   shuffle_heatmap.py 

Classes:       

Functions:     shuffle_heatmap

Author:        Wayne Dawson
creation date: 2016/2017
last update:   200211 (upgraded to python3), 190718
version:       0


Purpose: 

Reads in a chpair output(input) file, and shuffles the
information. This was part of an effort to compare randomly generated
heatmaps with the actual experimental heatmaps.


Comments:

In some respects, I still don't think that it would make any
difference. The point of an optimal biopolymer structure is that it
will have a balanace between the maximum entropy for a given set of
contacts. Biology will, in effect, intentionally chose arrangements
that satisfy both (apparent) randomness and effectiveness in gene
function. Moreover, since the free energy is this balance between the
maximum entropy for a given set of contacts, there will be a lot more
possible sets of patterns that are possible. So I cannot think that we
would be able to show that the chromatin patterns we observe and the
randomly generated ones (were the same number of contacts are
generated for a given entropy), would be any different. Maybe I am
wrong, but that is what I think.

190719:

   I really dislike this ChPair as yet another specialized class like
   LNode and Pair, only even less informative. Evidently, it was
   introduced in a very big hurry probably because I didn't feel so
   comfortable with using Pair or LNode or neither data representation
   was independent enough to use at the time for this purpose. At any
   rate, I don't need this highly minimalist class anymore and I will
   eventually eliminate it (replacing it with either Pair or LNode),
   though it is only a minor irritant so this plan is not really high
   on my priority list.

"""


from HeatMapTools import HeatMapTools
from FileTools    import FileTools
from FileTools    import getHeadExt
from ChPair       import ChPairData
from ChPair       import shuffle_ChPairData
from Chromatin2SimRNA import SimRNARestraints #from SimRNATools  import SimRNAData
import sys
import random
import argparse
import os
import string

PROGRAM = "shuffle_heatmap.py"



def shuffle_heatmap(cl):
    
    
    # parse the command line
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-f', nargs=1, default=["test.chpair"],
                        dest='f_chpair',
                        help="Input chpair file name [default 'test.chpair']")
    
    parser.add_argument('-o', nargs=1, default=None,
                        dest='f_header',
                        help="Header for output file [default is file header + index]")
    
    parser.add_argument('-sqlen', nargs=1, default=[-1],
                        dest='sqlen', type=int,
                        help='Set the sequence length: useful when value is not in chpair file.')
    
    parser.add_argument('-Nshuf', nargs=1, default=[1],
                        dest='n_shuffled', type=int,
                        help='Number of shuffled restraints.')
    
    
    input_opt = parser.add_mutually_exclusive_group()    
    input_opt.add_argument('-oSimRNA', action='store_true', default=True,
                           dest='oSimRNA',
                           help="output results as SimRNA restraint file.")
    
    input_opt.add_argument('-oHeat', action='store_true', default=False,
                        dest='oHeat',
                           help="output results as Heat Map file.")
    
    #
    # assign arguments
    args = parser.parse_args()
    # print (args)
    iflnm = args.f_chpair[0]
    iflhd, ext = getHeadExt(iflnm)
    oflhd = iflhd
    #print (args.sqlen)
    #print (args.n_shuffled)
    if not args.f_header == None:
        oflhd = args.f_header[0]
    #
    try:
        N = int(args.sqlen[0])
    except ValueError:
        print ("ERROR: input value for the sequence length must be an integer")
        print ("       entered value: '%s'" % args.sqlen)
        sys.exit(1)
    #
    try:
        Nshuf = int(args.n_shuffled[0])
    except ValueError:
        print ("ERROR: number of shuffled structures must be a POSITIVE integer")
        print ("       entered value: '%s'" % args.n_shuffled)
        sys.exit(1)
    #
    if Nshuf < 0:
        print ("ERROR: number of shuffled structures must be a POSITIVE integer")
        print ("       entered value: %d" % Nshuf)
        sys.exit(1)
    #
    use_SimRNA = args.oSimRNA
    use_heatmap = args.oHeat
    print (iflnm)
    ft = FileTools()
    
    if not ft.check_ext(iflnm, ["chpair"], "shuffle_heatmap"):
        print ("ERROR: input file is flawed")
        sys.exit(0)
    #
    
    print ('input file name:   ', iflnm)
    if use_SimRNA:
        print ('output files will contain SimRNA restraint data')
    else:
        print ('output files will contain heatmap data')
    #
    
    chdt = ChPairData()
    if N > 0:
        print ("N = ", N)
        chdt.set_vsqlen(N) # means N was set as an option
    #
    chdt.read_ChPairFile(iflnm)
    #print (chdt.disp_ChPairData(chdt.data))
    for k in range(1, Nshuf+1):
        oflnm = oflhd + "_x%s.simres" % (string.zfill(k, 3))
        print (oflnm)
        nchdt = shuffle_ChPairData(chdt)
        #print (chdt.disp_ChPairData(nchdt.data))
        srdt = SimRNAData()
        srdt.ChPair2SimRes(nchdt, ['N~N'])
        srdt.print_SimRNArestraints(oflnm, "slope", True) 
    
    #
#


# Main
if __name__ == '__main__':
    shuffle_heatmap(sys.argv)
#
