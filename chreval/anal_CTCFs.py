#!/usr/bin/python

"""@@@

Main Program:  anal_CTCF.py (ANALyze CTCF)

Classes:       Anchor
               Domain
               Anal

Author:        Wayne Dawson
creation date: mostly 2016 and up to March 2017.
last update:   181227
version:       0

Purpose:

Analyze the results of a chreval calculation on a whole group of these
heat maps (say on a chromasome or even the whole set).

anal_CTCF.py (ANALyze CTCF)

The idea behind this code is to take the set of data from a bed file
containing CTCF or RNA polymerase II data (or both) and assemble it
within the boundaries of a heat map file. For example, take a file
from "loops.CTCF.annotated.bed"

                      CTCF    CTCF    PET   cmplx   cmplx   cmplx     
chr    bgn     end     1       2      cnt    ity      a       b       open            active
chr1	838908	912011	R	L	11	0	4	6	0.0686291944243	0.426918183932
chr1	838908	922335	R	R	5	1	6	8	0.0601364066789	0.397329401752
chr1	838908	1000167	R	L	7	3	15	17	0.0910398799447	0.5333717808
chr1	918286	969271	R	R	5	0	4	4	0.119015396685	0.543963910954
chr1	918286	1000167	R	L	75	1	8	8	0.11802493863	0.65697780926
chr1	918286	1059032	R	R	4	2	12	12	0.09421226891	0.648785755901
chr1	967152	1000167	R	L	52	0	3	3	0.152052097531	0.842677570801
chr1	967152	1059032	R	R	6	1	7	7	0.0937744884632	0.711155855464
chr1	967152	1308124	R	L	7	2	34	34	0.0720704339359	0.62287812489
chr1	1306262	1377332	L	L	7	0	12	12	0.0938229914169	0.663528915154

chr1	1890973	1978919	R	L	41	0	-2	0	0.0051054055898	0.0633229481727
chr1	1890973	2106697	R	L	13	2	-1	3	0.0118067530734	0.228560568133
chr1	1890973	3341691	R	L	4	9	44	49	0.0292889451982	0.212245246836
chr1	2021345	2106697	R	L	11	0	-1	1	0.0131689942825	0.359991564345
chr1	2021345	2316695	R	L	9	2	17	18	0.0530184526833	0.47499238192
chr1	2105324	2316695	L	L	7	1	16	16	0.0687653462396	0.52343509753
chr1	2125307	2316695	R	L	259	0	13	12	0.0720316843271	0.512090622192

chr1	2343625	2480764	R	L	80	0	6	6	0.0746177236235	0.392732920613


There are two different starting points 838908 (with 3 entries) and
918286 (also with 3 entries). Hence, we can group the CTCFs to

( 838908, 1000167) c
( 838908,  922335) tR
( 838908,  912011) c

so, graphically, this first set suggests

-----x------------------|---|--------------------------------x------
     |838908      912011|   |922335                          |1000167      
     |_________c________|   |                                |      
     |_________tR___________|                                |
     |_________c_____________________________________________|

and 

( 918286,  969271) tr
( 918286, 1000167) c 
( 918286, 1059032) tr

and

( 967152, 1000167) c	
( 967152, 1059032) tr	
( 967152, 1308124) c	

and

(1306262, 1377332) tL

so these two data sets pictured graphically are arranged as follows:

                           _______tR______________________________________ .. _________
                          |_______c__________________________                          |
                          |_______tR________                 |                         |
                          |                 |                |                         |
                          |918286           |969271          |                         |1059032
                          |                 _____c___________|____________ .. __________________ .. _________
                          |                |_____tR__________|____________ .. _________                       |
                          |                |_____c___________|                         |                      |
                          |                ||                |                         |                      |
                          |          967152||969271          |             ..          |1059032  ..           |1308124
                          |                ||                |                         |                      |
     |838908              |                ||                |1000167      ..          |1059032               |1308124
-----x------------------|-|-|--------------||----------------x------------ .. ---------x-------- .. ----------x----------- .------------
     |838908      912011|   |922335                          |1000167                                 1306262|                     |1377332
     |_________c________|   |                                |                                               |____________ .. _____|
     |_________tR___________|                                |
     |_________c_____________________________________________|




separation between the short fragments:
int( 0.5 + ( 918286 - 912011) / 5000 ) =  1
int( 0.5 + ( 922335 - 918286) / 5000 ) =  1
int( 0.5 + ( 922335 - 912011) / 5000 ) =  2
int( 0.5 + ( 969271 - 967152) / 5000 ) =  0
int( 0.5 + (1308124 -1306262) / 5000 ) =  0

absolute points on the grid from 838908
int( 0.5 + ( 912011 - 838908) / 5000 ) =  15*
int( 0.5 + ( 918286 - 838908) / 5000 ) =  16*
int( 0.5 + ( 922335 - 838908) / 5000 ) =  17*
---
int( 0.5 + ( 967152 - 838908) / 5000 ) =  26*
int( 0.5 + ( 969271 - 838908) / 5000 ) =  26*
---
int( 0.5 + (1000167 - 838908) / 5000 ) =  32
---
int( 0.5 + (1059032 - 838908) / 5000 ) =  44
---
int( 0.5 + (1306262 - 838908) / 5000 ) =  93*
int( 0.5 + (1308124 - 838908) / 5000 ) =  94*
---
int( 0.5 + (1377332 - 838908) / 5000 ) = 108

where the "*" means that there are several CTCFs that are in very
close proximity with respect to the heatmap.


   The input file is a list of loops formatted according to
   Przemek's approach (with extension 'bed').

   command line example:
   > anal_CTCFs.py -ff loops.CTCF.annotated.bed -dG

   where the file "loops.CTCF.annotated.bed" contains a list of data
   formatted according to Przemek. This has the format

   file format1 (loops.CTCF.annotated.bed):

                      CTCF    CTCF    PET   cmplx   cmplx   cmplx     
chr    bgn     end     1       2      cnt    ity      a       b       open            active
chr1	838908	912011	R	L	11	0	4	6	0.0686291944243	0.426918183932
chr1	838908	922335	R	R	5	1	6	8	0.0601364066789	0.397329401752
chr1	838908	1000167	R	L	7	3	15	17	0.0910398799447	0.5333717808
chr1	918286	969271	R	R	5	0	4	4	0.119015396685	0.543963910954
chr1	918286	1000167	R	L	75	1	8	8	0.11802493863	0.65697780926
chr1	918286	1059032	R	R	4	2	12	12	0.09421226891	0.648785755901

   file format2 (structure.loops.GSE6352.annotated.bed):
chr1	1050000	1190000	241	-1	9	9	0.0529785714286	0.418914285714
chr1	1585000	1650000	80	-2	0	0	0.0553384615385	0.737907692308
chr1	1710000	1840000	154	-2	22	22	0.0699230769231	0.874938461538

   file format2 (structure.TAD_Rao.annotated.bed):
chr1	915000	1005000	0.72251	1	10	10	0.116533333333	0.643166666667
chr1	1030000	1235000	1.1954	-1	12	12	0.0460487804878	0.507936585366
chr1	1255000	1450000	0.9312	0	27	27	0.0741435897436	0.653230769231

   In general, if more than one file is used (with extension "bed"),
   all files should be in the same format.

   These files often have different data related to open, active,
   inactive, compartment A or B, etc., so the results typically
   require more than one file.  Therefore, multiple file entries are
   allowed for both file types; however, in general, for "bed"
   files, there should only be one input file.

"""

import sys
import os
import chreval
from FileTools     import FileTools

from ChromatinData import Data
from ChromatinData import make_file_heading

from BasicTools    import initialize_matrix
from HeatMapTools  import HeatMapTools

#
PROGRAM = "anal_CTCFs.py"
EXTS    = ["bed"]
USAGE = ["USAGE: %s -ff [file].bed  [[file2].bed ... ]" % PROGRAM,
         "",
         "Purpose:",
         "To combine all the information from a *.bed file (or a ",
         "short list of *.bed files) and a directory full of ",
         " *.heat files into a single unified set of files with ",
         "the corresponding names *.eheat (for Extended heat).",
         "",
         "The *.eheat file specifies the details of the CTCFs and",
         "the heatmap and is readable by the program chreval.",
         "",
         "How to Run This Program....",
         "(0) It is advised to work with _one_ complete bed file",
         "    than combine many bed files together. Nevertheless,",
         "    the program will follow instructions (even if they",
         "    are wrong!!!).",
         "",
         "(1) You _must_ be located in the directory where ",
         "    the \'*.heat\' files are located. ",
         "",
         "(2) You specify a \'*.bed\' file (or files) onto which ",
         "    you want to graft the relevant CTCF information into",
         "    the \'*.heat\' file. (NOTE: if the \'*.bed\' file is",
         "    not located in the same directory as the \'*.heat\'",
         "    files, then you must specify the path [e.g., ",
         "    \'../*.bed\'] where the file is located.)",
         "",
         "(3) The program automatically reads the contents of the ",
         "    \'bed\' file and deposits the information about the ",
         "    CTCFs into the coresponding \'eheat\' file. These ",
         "    results are deposited in the subdirectory \'tmp\' ",
         "    (in the current directory where the \'*.heat\' files",
         "    are located). Hopefully, this will help avoid cluttering",
         "    the current directory or corrupting the original files ",
         "    in some unfortunate way.",
         "",
         "(4) This program can take a while to finish, especially ",
         "    if you want to convert \'*.heat\' files for a whole ",
         "    genome. Typical times are 1 hr. Therefore, be patient! "]
#


def usage():
    for uk in USAGE:
        print uk
#

def checkDomainInput(b_bgn, b_end, p_bgn, p_end):
    flag_pass = True
    # domain boundaries:  b_bgn, b_end
    # segment boundaries: p_bgn, p_end
    
    # It should _always_ be the case that these input parameters
    # satisfy the following conditions
    
    # b_bgn <= p_bgn < p_end <= b_end
    
    # The program also checks for other really stupid things.
    
    if b_bgn > 10e20: # no genome is this size!
        print "ERROR: domain boundary base is not properly set"
        print "boundary bgn = ", b_bgn
        flag_pass = False
    #
    
    if b_end < 1: # no genome is this small!
        print "ERROR: domain boundary base is not properly set"
        print "boundary end = ", b_end
        flag_pass = False
    #
    
    if b_end <= b_bgn:
        print "ERROR: domain boundaries don't make sense"
        print "boundary bgn = ", b_bgn
        print "boundary end = ", b_end
        flag_pass = False
    #
    
    if b_bgn > p_bgn:
        print "ERROR: segment boundary starts before domain boundary"
        print "boundary bgn = ", b_bgn
        print "point bgn    = ", p_bgn
        flag_pass = False
    #
    
    if b_end < p_end:
        print "ERROR: segment boundary ends after domain boundary"
        print "boundary bgn = ", b_end
        print "point bgn    = ", p_end
        flag_pass = False
    #
    
    if p_end <= p_bgn:
        print "ERROR: segment boundaries don't make sense"
        print "point bgn = ", p_bgn
        print "point end = ", p_end
        flag_pass = False
    #
    
    if not flag_pass:
        sys.exit(1)
    #
#


class Anchor:
    def __init__(self, chrN, bgn, end, nPET, data):
        self.chrN = chrN
        self.bgn  = bgn
        self.end  = end
        self.nPET = nPET
        self.anchordata = data # a line from the bed file
    #
    def showAnchor(self):
        return self.anchordata
    #
#

class Domain:
    def __init__(self):
        self.gmin = 1e308 # ~inf
        self.gmax = -1
        self.length = 0
        self.anchors = [] # [Anchor_1, Anchor_2, .... Anchor_n]
        self.resolution = float(5000) # default is
        self.mtools              =  HeatMapTools() # default setup GenerateHeatMapTools()
        self.mtools.ctcf_tthresh =  10
        self.mtools.ctcf_cthresh =  40
        self.mtools.from_Nenski  =  False
        self.gmap = []
    #
    
    def showDomain(self):
        s = ''
        for a in self.anchors:
            bgn, end = self.locate_on_Map(a.bgn, a.end)
            s += "%s  %6d  %6d\n" % (a.showAnchor(), bgn, end)
        #
        return s
    #

    def set_resolution(self, resolution):
        if resolution <= 0:
            print "ERROR: resolution must be positive and greater than zero"
            sys.exit(1)
        #
        self.resolution = float(resolution)
    #
    
    def addAnchor(self, aa):
        if aa.bgn < self.gmin:
            self.gmin = aa.bgn
        #
        if aa.end > self.gmax:
            self.gmax = aa.end
        #
        self.anchors += [aa]
    #
    
    def locate_on_Map(self, p_bgn, p_end):
        checkDomainInput(self.gmin, self.gmax, p_bgn, p_end)
        # print "gmin,gmax = %5d,%5d  p_bgn,p_end = %5d,%5d" % (self.gmin, self.gmax, p_bgn, p_end)
        v_bgn = int( 0.5 + (float(p_bgn - self.gmin) / self.resolution ))
        v_end = int( 0.5 + (float((p_bgn - self.gmin) + (p_end - p_bgn)) / self.resolution ))
        return v_bgn, v_end
    #
    
    def make_Map(self, flag_display = True):
        debug_make_Map = False
        self.gmap = []
        b_bgn, b_end = self.locate_on_Map(self.gmin, self.gmax)
        if debug_make_Map:
            print "b_bgn,b_end, gmin,gmax", b_bgn, b_end, self.gmin, self.gmax
        #
        self.gmap = initialize_matrix(self.gmap, b_end + 1, 0)
        for a in self.anchors:
            bgn, end = self.locate_on_Map(a.bgn, a.end)
            flnm = "%s_%d_%d_res5kb.heat" % (a.chrN, a.bgn, a.end)
            if debug_make_Map:
                print flnm
                print "bgn,end, a.bgn,a.end", bgn, end, a.bgn, a.end
            #
            if os.path.exists(flnm):
                if debug_make_Map:
                    print "shift: ", bgn
                #
                gmtrx = self.mtools.read_MatrixFile(flnm, ["heat"], 1.0, "make_Map", flag_display, bgn)
                N        = gmtrx.length
                hv       = gmtrx.heatmap
                clusters = gmtrx.clusters
                #                read_MatrixFile(flnm, allowed_extns, rescale_wt)
                i = 0; j = 0
                for jv in range(0, N):
                    for iv in range(0, N):
                        i = iv + bgn
                        j = jv + bgn
                        self.gmap[i][j] = hv[iv][jv]
                    #
                #
            else:
                print "WARNING: cannot open %s; the file missing" % flnm
                i = bgn
                j = end
                self.gmap[i][j] = a.nPET
            #
        #
        
        # this ensures that every anchor is present
        for a in self.anchors:
            bgn, end = self.locate_on_Map(a.bgn, a.end)
            print "bgn,end = (%5d, %5d)[%5d]" % (bgn, end, a.nPET)
            self.gmap[bgn][end] = a.nPET
            self.gmap[end][bgn] = a.nPET
        #
    #

#


class Anal:
    def __init__(self, data):
        self.debug              = False
        self.data               = data
        self.cdata              = data.cdata
        self.keylist            = data.keylist
        self.datatype           = data.datatype
        self.res                = "res5kb"
        self.range_TddS         = 2.0
        self.range_dG           = 2.0
        self.domain             = []
    #
    
    def get_PET_wt(self, ky, wt):
        ctcf1 = self.cdata[ky].ctcf1
        ctcf2 = self.cdata[ky].ctcf2
        wt_ctcf = wt
        if (ctcf1 == 'R' and ctcf2 == 'R') or (ctcf1 == 'L' and ctcf2 == 'L'):
            wt_ctcf *= 0.25
        elif (ctcf1 == 'L' and ctcf2 == 'R'):
            wt_ctcf *= 0.05
        #
        return wt_ctcf
    #
    
    
    def build_CTCF_map(self, cl):
        debug_build_CTCF_map = False
        missing = []
        
        if debug_build_CTCF_map: 
            print "build_CTCF_map()"
        #
        
        keylist = self.keylist
        if len(self.keylist) == 0:
            # not sure that it really works, but it would seem that if
            # there is some data loaded, then it can actually process
            # it at this step too. I think the program should just
            # stop if it doesn't have a set keylist, but we will see.
            print "setting keylist"
            print "length of the dataset: ", len(self.cdata)
            keylist = self.data.ordered_keylist()
        #
        
        print "length of keylist: ", len(self.keylist)
        if debug_build_CTCF_map: 
            print "keylist:"  # this is ordered
            print keylist
            print "cdata.keys():"  # this is just the dictionary data
            print self.cdata.keys()
        #
        
        print "data.chrmsm_grp.keys():"  # This gives a list of the chromosome used
        print self.data.chrmsm_grp.keys()
        print "listing:"
        for chrN in self.data.chrmsm_grp.keys():
            for segm_k in range(0, len(self.data.chrmsm_grp[chrN].chrsegment)-1):
                chrpos = self.data.chrmsm_grp[chrN].chrsegment[segm_k]
                #
                v = self.data.makekey([chrN, chrpos[0], chrpos[1]])
                nPET   = self.cdata[v].nPET
                if debug_build_CTCF_map: 
                    print v, \
                        self.cdata[v].ctcf1, \
                        self.cdata[v].ctcf2, \
                        self.cdata[v].nPET, \
                        self.cdata[v].state["open"], \
                        self.cdata[v].state["active"]
                #
                
                # list of complexity measure
                cc = ''
                for ccx in self.cdata[v].cmplx:
                    cc += '%4d  ' % ccx
                #
                
                # make a string that indicates the contacts for storage 
                ss = "%5s %10d  %10d  %s  %s  %8g  %s %10.6f  %10.6f  %10.6f  %10.6f" \
                    % (chrN, chrpos[0], chrpos[1], 
                       self.cdata[v].ctcf1, self.cdata[v].ctcf2, 
                       self.cdata[v].nPET,
                       cc, 
                       self.cdata[v].state["open"],
                       self.cdata[v].state["active"],
                       self.cdata[v].state["A"],
                       self.cdata[v].state["B"])
                # print ss
                aa = Anchor(chrN, chrpos[0], chrpos[1], nPET, ss)
                flhd = make_file_heading(chrN, chrpos[0], chrpos[1], self.res)
                
                # #######################################################
                # verify that the files exist or not. If they are
                # missing, then at least make a record of it.
                # #######################################################
                
                # print "evaluating: ", os.path.exists(cl.f_heatmap[0]), cl.f_heatmap[0]
                if not os.path.exists(flhd + ".heat"):
                    missing += [ flhd ]
                    continue
                #
                
                if len(self.domain) == 0:
                    newdmn = Domain()
                    newdmn.addAnchor(aa)
                    self.domain += [newdmn]
                else:
                    for dmnlist in self.domain:
                        gmin = dmnlist.gmin
                        gmax = dmnlist.gmax
                    if (chrpos[0] == gmin):
                        dmnlist.addAnchor(aa)
                    elif (chrpos[1] == gmax):
                        dmnlist.addAnchor(aa)
                    elif ((chrpos[0] <= gmin) and (chrpos[1] >= gmax)):
                        dmnlist.addAnchor(aa)
                    elif ((chrpos[0] >= gmin) and (chrpos[1] <= gmax)):
                        dmnlist.addAnchor(aa)
                    elif ((chrpos[0] <= gmin) and (gmin <= chrpos[1])):
                        dmnlist.addAnchor(aa)
                    elif ((chrpos[0] <= gmax) and (gmax <= chrpos[1])):
                        dmnlist.addAnchor(aa)
                    else:
                        newdmn = Domain()
                        newdmn.addAnchor(aa)
                        self.domain += [newdmn]
                    #
                #
            #
            # print self.data.chrmsm_grp[chrN].show_chr_fragments()
            for kdmn in range(0, len(self.domain)):
                print "cluster: %6d   %10d   %10d" % (kdmn, self.domain[kdmn].gmin, self.domain[kdmn].gmax)
                print "chr        bgn         end    ctcfs      #PET ---- cplxty ----      open        active         A          B"
                print self.domain[kdmn].showDomain()
        #
        print "missing:"
        print missing
        
        try:
            tmpdir = "tmp"
            if not os.path.isdir(tmpdir):
                os.mkdir(tmpdir) # build a subdirectory to store figures in
            else:
                print "WARNING: the directory '%s'" % tmpdir
                print "         already exists, overwriting files\n" 
        except OSError:
            print "ERROR: problems making %s" % tmpdir
            sys.exit(1)
        #
        
        for chrN in self.data.chrmsm_grp.keys():
            # print self.data.chrmsm_grp[chrN].show_chr_fragments()
            for kdmn in range(0, len(self.domain)):
                if debug_build_CTCF_map:
                    print ""
                    print "cluster: %6d   %10d   %10d" % (kdmn, self.domain[kdmn].gmin, self.domain[kdmn].gmax)
                self.domain[kdmn].make_Map(False)
                if debug_build_CTCF_map:
                    print len(self.domain[kdmn].gmap[0])
                #
                dmnflnm = "%s/%s_%d_%d.eheat" % (tmpdir, chrN, self.domain[kdmn].gmin, self.domain[kdmn].gmax)
                print dmnflnm
                self.domain[kdmn].mtools.write_ExtHeatMap(dmnflnm, self.domain[kdmn], "0.0")                
            #
        #
        #
    #
#
    


def main(cl):

    if len(cl) < 2:
        print "ERROR: too few arguments"
        usage()
        sys.exit(1)
    #
    k = 1
    inbedflnm = []
    while k < len(cl):
        # print cl[k]
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
                # print cl[m]
                if cl[m][0] == "-":
                    # print "1", m
                    k = m - 1 
                    break
                elif ft.check_ext(cl[m], EXTS):
                    # print "2", m
                    k = m
                    inbedflnm += [cl[m]]
                else:
                    # print "3", m
                    k = m 
                    break
                #
            #
            # print inbedflnm
        #
        else:
            print "ERROR: unrecognized flag (%s)." % cl[k]
            usage()
            sys.exit(1)
        #
        k += 1
    #
    
    # summary of the input files
    print "\ninput *.bed file(s) to be processed: "
    for f in inbedflnm:
        print "%s" % f
    #

    # now do the processing ... 
    dt = Data()
    
    for f in inbedflnm:
        if os.path.exists(f):
            dt.readfile(f)
        else:
            print "ERROR: cannot open %s; the file missing" % f
            sys.exit(1)
        #
    #
    dt.ordered_keylist() # make an ordered list in terms of chromosome
                         # number and region
    dt.display_Data(False)
    
    an = Anal(dt)
    an.build_CTCF_map(inbedflnm)
    
    
#

if __name__ == '__main__':
    # running the program
    main(sys.argv)

