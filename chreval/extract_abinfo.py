#!/usr/bin/env python3

"""@@@

Main Program:  extract_abinfo.py 

Classes:       

Author:        Wayne Dawson
creation date: 200319
last update:   200319
version:       0

Purpose:

Just a trick for displaying the entries found in ab_missing

"""

import sys
import os
import chreval
from FileTools     import FileTools
from FileTools     import getHeadExt

from ChromatinData import Data
from ChromatinData import make_file_heading

from BasicTools    import initialize_matrix
from HeatMapTools  import HeatMapTools

from assemble_heatmaps_and_CCDs import Domain
from assemble_heatmaps_and_CCDs import Anchor
from assemble_heatmaps_and_CCDs import BuildHeatMaps
from assemble_heatmaps_and_CCDs import build_tmpDir
from assemble_heatmaps_and_CCDs import assemble_eheatMaps


PROGRAM = "extract_abinfo.py"
EXT1    = ["bed"]
EXT2    = ["txt", "dat"]

# tag for printing
USAGE = ["\n\nUSAGE: %s -ff [ab_ref_file].bed  -rr ab_missing.txt" % PROGRAM,
         "",
         "Purpose:",
         "Used to match the data in the reference *.bed file with the ",
         "output file \"ab_missing.txt\" (or equivalent). Prints out ",
         "similar information as the original.",
         "",
         "Examples:",
         "  > %s -ff compartments.annotated.x_13032020_done.bed -rr ab_missing.txt" \
         % PROGRAM, 
         "  > %s -ff loops.CTCF.annotated.with_A_and_B.bed -rr ctcf_missing.txt" \
         % PROGRAM
]
#


def usage():
    for uk in USAGE:
        print (uk)
    #|endfor
    
#



def main(cl):
    debug_main = False # True # 
    if len(cl) < 2:
        print ("ERROR: too few arguments")
        usage()
        sys.exit(1)
    #
    
    k = 1
    inbedflnm = []
    inrefflnm = ''
    while k < len(cl):
        # print (cl[k])
        if cl[k] == "-h" or cl[k] == "--help" or cl[k] == "-help":
            usage()
            sys.exit(0)
            
        elif cl[k] == "-ff":
            # read multiple files following the commandline argument
            # -ff until some other argument (presently this can only
            # be the help requests "-h", "--help" or "-help").
            ft = FileTools()
            k += 1
            for m in range(k, len(cl)):
                if debug_main: print ("cl[%d]: %s" % (m, cl[m]))
                if cl[m][0] == "-":
                    # if it encounters another command line argument,
                    # it terminates
                    
                    if debug_main:
                        print ("1", m)
                    #
                    
                    k = m - 1 
                    break
                elif ft.check_ext(cl[m], EXT1):
                    # read in another bed file.
                    
                    if debug_main:
                        print ("2", m)
                    #
                    
                    k = m
                    inbedflnm += [cl[m]]
                else:
                    # I don't think it can ever get here in reality
                    
                    if debug_main:
                        print ("3", m)
                    #
                    
                    k = m 
                    break
                #
                
            #|endfor
                
        elif cl[k] == "-rr":
            # read multiple files following the commandline argument
            # -ff until some other argument (presently this can only
            # be the help requests "-h", "--help" or "-help").
            ft = FileTools()
            k += 1
            if ft.check_ext(cl[k], EXT2):
                if debug_main:
                    print ("2", k)
                #
                
                inrefflnm = cl[k]
            else:
                print ("ERROR: reference file must be of type %s" % EXT2)
                if debug_main:
                    print ("3", m)
                #
                
                k = m 
                sys.exit(1)
            #
            
            if debug_main:
                print ("inrefflnm: ", inrefflnm)
            #
            
        else:
            print ("ERROR: unrecognized flag (%s)." % cl[k])
            usage()
            sys.exit(1)
        #
        
        k += 1
    #
    
    # summary of the input files
    print ("\ninput *.bed file(s) to be processed: ")
    for f in inbedflnm:
        print ("%s" % f)
    #|endfor
    
    print ("reference AB file: %s" % inrefflnm)
    
    #print ("stop at main 1"); sys.exit(0)
    
    # now do the processing ...
    
    
    # read in the CTCF bed file(s) 
    abrefdt = Data()
    for f in inbedflnm:
        if os.path.exists(f):
            abrefdt.readfile(f)
        else:
            print ("ERROR: cannot open %s; the file missing" % f)
            sys.exit(1)
        #
        
    #|endfor
    
    abrefdt.ordered_keylist()
    # ordered_keylist: make an ordered list in terms of chromosome
    # number and region. This was before I learned about "from
    # compilations import OrderedDict", but anyway.
    
    abrefdt.display_Data(False)

    try:
        fp = open(inrefflnm, 'r')
        lfp = fp.readlines()
        fp.close()
    except:
        print ("ERROR: cannot open %s" % inrefflnm)
        sys.exit(1)
    #
    
    flhd, ext = getHeadExt(inrefflnm)
    out_flnm = flhd + "_annotated.txt"
    abdata = []
    for k in range(len(lfp)):
        sv = lfp[k].strip()
        if len(sv) == 0:
            continue
        #
        
        if sv[0] == '#':
            print ("%s" % sv)
        #
        
        svv = sv.split('_')
        if len(svv) == 4:  
            sv_key = abrefdt.makekey(svv[:3])
            abdata += [ sv_key ]
        #
    #
    
    # length of original dataset:  1257
    # missing files:

    with_hollerith = False
    text  = "# original file(s):\n"
    for v in inbedflnm:
        text += "#      %s\n" % v
    #
    text += "# reference file:   %s\n" % inrefflnm
    text += "# number in original dataset: %6d\n" % len(abrefdt.cdata)
    text += "# number of missing items:    %6d\n" % len(abdata)
    text += "\n"
    
    text += abrefdt.disp_cdata_head(with_hollerith)
    for sv_key in abdata:
        text += (abrefdt.disp_cdata_key(sv_key, with_hollerith) + '\n')
    #
    
    fp = open(out_flnm, 'w')
    fp.write(text)
    fp.close()
    #print (text)
    print ("output file: %s" % out_flnm)
    print ("Done")
    
    
#

if __name__ == '__main__':
    # running the program
    main(sys.argv)

#
