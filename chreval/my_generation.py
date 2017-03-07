#!/usr/bin/env python

# my_generation.py: generates a heatmap for reading by chreval.py.
# This is intended to be used to build test structures for chreval.py
# to find. It is particularly needed for pseudoknot like contraptions
# where all sorts of things hairy things could go wrong.

from Vienna import Vstruct
from MatrixTools import MatrixTools
import sys
import random
import argparse
import os

PROGRAM = "my_generation.py"

def help_ExFile():
    print "Here is an example of an input file:\n"
    print "     > cat example.ss"
    print "     .ABC..........abc.. 20 "
    print "     ....DE...........de  3 "
    print "     .......((..))...... 10 "
    print "     {.................}  5"
    print "\nThe first two sequences express a parallel stem"
    print ".ABCDE........abcde\n"
    print "The stem is divided into two pieces with weights 20 and 3,"
    print "respectively."
    print "The third sequence is an antiparallel stem of weight 10"
    print ".......((..))......\n"
    print "The fourth sequence is a CTCF \"{..}\" with an assigned weight 5."
    print "{.................}\n"
    print "If no weight is ascribed to the sequence, then a default value"
    print "will be used for the elements in that line. Note also that part"
    print "of the parallel stem overlaps with the CTCF. Hence, more than one"
    print "structure is allowed in this scheme."
    print "\ncommand line:"
    print "\n     > %s -f example.ss" % PROGRAM
#    

def help_ExSeq():
    print "command line:\n"
    print "     > %s -seq \".ABCDE.((..)).abcde\" \"{.................}\"" % PROGRAM
    print "\nThe first sequence expresses a parallel stem"
    print ".ABCDE........abcde\n"
    print "and an antiparallel stem"
    print ".......((..))......\n"
    print "The second sequence is a CTCF."
    print "{.................}\n"
    print "To ascribe weights, an input file is required with the flag '-f'."
    print "Use flag '-hExFile' for more information. Note also that part"
    print "of the parallel stem overlaps with the CTCF. Hence, more than one"
    print "structure is allowed in this scheme."
#    


def add_noise(N, hv):
    for j in range(0, N):
        for i in range(0,j):
            u1 = int(random.uniform(0,1) + 0.5)
            u2 = int(random.uniform(0,1) + 0.5)
            v = u1*u2*int(random.uniform(0,2) + 0.5)
            hv[i][j] = v
            hv[j][i] = v
        #
    #
    return hv
#

def read_SeqFile(iflnm):
    print iflnm
    if not os.path.isfile(iflnm):
        print "ERROR: %s not found" % iflnm
        sys.exit(1)
    #
    fp = open(iflnm, 'r')
    lfp = fp.readlines()
    fp.close()
    ss_seq = []
    ss_wt  = []
    for k in range(0, len(lfp)):
        s = lfp[k].strip().split()
        if len(s) > 0:
            if s[0][0] == '#':
                continue
            if len(s) == 1:
                ss_seq += [s[0]]
                ss_wt  += [-1]
            elif len(s) == 2:
                ss_seq += [s[0]]
                try:
                    wt = int(s[1])
                except ValueError:
                    print "ERROR: input file '%s' must contain a sequence and an integer" % iflnm
                    print "       line(%d): " % k, s
                    sys.exit(1)
                #
                if wt < 1:
                    print "ERROR: input file '%s' must contain a sequence and a POSITIVE integer" % iflnm
                    print "       line(%d): " % k, s
                    sys.exit(1)
                #
                ss_wt += [wt]
            else:
                print "ERROR: input file '%s' has too many entries. I don't understand it." % iflnm
                print "       line(%d): " % k, s
                sys.exit(1)
                
        #
    #
    for k in range(0, len(ss_seq)):
        print ss_seq[k], ss_wt[k] 
    return ss_seq, ss_wt
#
        
def generate(ss_seq, ss_wt, w_noise, oflnm):
    debug_generate = False
    vs = []
    k = 0
    N = len(ss_seq[0])
    
    # setup objects and check input sequence (or sequences)
    for k in range(0, len(ss_seq)):
        vs += [Vstruct()]
        n = len(ss_seq[k])
        if not n == N:
            print "ERROR: input sequence is not the same length as the others"
            m = 1
            for m in range(0, len(ss_seq)):
                print "seq(%2d): %s" % (m, ss_seq[m])
            print "new seq: %s" % ss_seq[k]
            sys.exit(1)
        #
        n_vs = len(vs) - 1
        vs[n_vs].wt = ss_wt[n_vs]
    #
    if debug_generate:
        for k in range(0, len(vs)):
            print "%d: %d" % (k, vs[k].wt)
        #
    #
    
    for k in range(0, len(vs)):
        print "line(%2d): " % k
        vs[k].reset_Vstruct()
        if debug_generate:
            print "Xlist:  ", vs[k].Xlist
            print "BPlist: ", vs[k].BPlist
        #
        vs[k].set_Vstruct(ss_seq[k])
        vs[k].scan_ctcf(0, 0)
        vs[k].print_vstr()
        
        if debug_generate:
            print "Xlist:  ", vs[k].Xlist
            print "BPlist: ", vs[k].BPlist
        #
        flag_other = True
        if flag_other:
            # display for information purposes
            if len(vs[k].BPlist) > 0:
                print "secondary structure"
                vs[k].print_Xlist_n(vs[k].BPlist)
            #
            if len(vs[k].PKlist) > 0:
                print "pk connects"
                vs[k].print_Xlist_n(vs[k].PKlist)
            #
            if len(vs[k].CTCFlist) > 0:
                print "CTCF connects"
                vs[k].print_Xlist_n(vs[k].CTCFlist)
            #
        #
    #
    # sys.exit(0)
    mtools = MatrixTools()
    hv = []
    hv = mtools.initialize_matrix(hv, N, 0)
    
    # add noise if requested
    if w_noise:
        hv = add_noise(N, hv)
    #
    for k in range(0, len(vs)):
        for pair in vs[k].BPlist:
            i = pair.i
            j = pair.j
            if debug_generate:
                print "BPlist: ", vs[k].wt
            #
            if vs[k].wt > 0:
                hv[i][j] = vs[k].wt
            else:
                hv[i][j] = 3
            #
            hv[j][i] = hv[i][j]
        #
        for pair in vs[k].PKlist:
            i = pair.i
            j = pair.j
            if debug_generate:
                print "PKlist: ", vs[k].wt
            #
            if vs[k].wt > 0:
                hv[i][j] = vs[k].wt
            else:
                hv[i][j] = 5
            #
            hv[j][i] = hv[i][j]
        #
        for pair in vs[k].CTCFlist:
            i = pair.i
            j = pair.j
            if debug_generate:
                print "CTCFlist: ", vs[k].wt
            #
            if vs[k].wt > 0:
                hv[i][j] = vs[k].wt
            else:
                hv[i][j] = 100
            #
            hv[j][i] = hv[i][j]
        #
    #
    if debug_generate:
        print hv
    #
    heatmap = mtools.make_heatmap(hv)
    print heatmap
    fp = open(oflnm, 'w')
    fp.write(heatmap)
    fp.close()
#
    



def main(cl):
    
    
    # parse the command line
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-o', nargs=1, default=["test.heat"],
                        dest='f_heatmap',
                        help="name of output heatmap file [default 'test.heat']; the current output is in the older format where CTCF and singleton data are jumbled together on a single heatmap.")
    
    parser.add_argument('-add_noise', action='store_true', default=False,
                        dest='w_noise',
                        help='Request to include noise in the output heatmap.')
    
    parser.add_argument('-hExFile', action='store_true', default=False,
                        dest='help_ExFile',
                        help='shows an example of an input file.')
    parser.add_argument('-hExSeq', action='store_true', default=False,
                        dest='help_ExSeq',
                        help='shows an example of a command line call.')
    
    input_opt = parser.add_mutually_exclusive_group()    
    input_opt.add_argument('-f', nargs=1, default=["test.ss"],
                        dest='f_ss',
                           help="name of input 1D structure file [default 'test.ss']; the file can contain multiple lines (with structures of the same length), but each line _must_ contain only one Fontana-type formatted sequence.")
    
    input_opt.add_argument('-seq',     action='store', nargs='+', default=None,
                        dest='seq', 
                           help='The input structure (in Fontana-like format); e.g. \".ABCDE..abcde\" \"{...........}\" generates a CTCF boundary "{...........}" and a parallel stem of length 5, \".ABCDE..abcde\", where \".....E......e\" overlaps with the CTCF.')
    
    
    #
    # assign arguments
    args = parser.parse_args()
    # print args
    if args.help_ExFile:
        help_ExFile()
        sys.exit(0)
    #
    if args.help_ExSeq:
        help_ExSeq()
        sys.exit(0)
    #
    
    ss_seq = args.seq
    ss_wt = []
    oflnm = args.f_heatmap[0]
    iflnm = args.f_ss[0]
    w_noise = args.w_noise
    print 'input file name:   ', iflnm
    print 'output file name:  ', oflnm

    if ss_seq == None:
        # if ss_seq is 
        ss_seq, ss_wt = read_SeqFile(iflnm)
    else:
        for k in range(0, len(ss_seq)):
            ss_wt += [-1]
        #
    #
    print 'input sequences:   '
    for k in range(0, len(ss_seq)):
        if ss_wt[k] > 0:
            print "%s  %4d" % (ss_seq[k], ss_wt[k])
        else:
            print "%s   default weight" % ss_seq[k]
        #
    #
    generate(ss_seq, ss_wt, w_noise, oflnm)
#


# Main
if __name__ == '__main__':
    main(sys.argv)

