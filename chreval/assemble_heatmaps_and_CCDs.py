#!/usr/bin/env python3

"""@@@

Main Program:  assemble_heatmaps_and_CCDs.py 
               (originally called anal_CTCFs.py (ANALyze CTCFs))

Classes:       Anchor
               Domain
               BuildHeatMaps

Author:        Wayne Dawson
creation date: mostly 2016 and up to March 2017.
last update:   200316 gone through, checked, commented and structured for stitching
version:       0

Purpose:

Build eheat files out of CCDs stored in the *.bed files and the
respective heatmaps by growing clusters of CTCFs into domains around
central points in the bed and heatmap file data. This tool assembles a
final heat map from the groupings for any set of heat maps that
actually exist. It does not build a predefined region of heat map.


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
   > assemble_heatmaps_and_CCDs.py -ff loops.CTCF.annotated.bed 

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


PROGRAM = "assemble_heatmaps_and_CCDs.py"
EXTS    = ["bed"]
USAGE = ["\n\nUSAGE: %s -ff [file].bed  [[file2].bed ... ]" % PROGRAM,
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
        print (uk)
    #|endfor
    
#

def build_tmpDir(dirname):
    
    # now we begin to stitch together the various heatmaps
    
    try:
        tmpdir = dirname
        if not os.path.isdir(tmpdir):
            os.mkdir(tmpdir) # build a subdirectory to store figures in
        else:
            print ("WARNING: the directory '%s'" % tmpdir)
            print ("         already exists, overwriting files\n" )
        #
        
    except OSError:
        print ("ERROR: problems making %s" % tmpdir)
        sys.exit(1)
    #
    
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
        print ("ERROR: domain boundary base is not properly set")
        print ("boundary bgn = ", b_bgn)
        flag_pass = False
    #
    
    if b_end < 1: # no genome is this small!
        print ("ERROR: domain boundary base is not properly set")
        print ("boundary end = ", b_end)
        flag_pass = False
    #
    
    if b_end <= b_bgn:
        print ("ERROR: domain boundaries don't make sense")
        print ("boundary bgn = ", b_bgn)
        print ("boundary end = ", b_end)
        flag_pass = False
    #
    
    if b_bgn > p_bgn:
        print ("ERROR: segment boundary starts before domain boundary")
        print ("boundary bgn = ", b_bgn)
        print ("point bgn    = ", p_bgn)
        flag_pass = False
    #
    
    if b_end < p_end:
        print ("ERROR: segment boundary ends after domain boundary")
        print ("boundary bgn = ", b_end)
        print ("point bgn    = ", p_end)
        flag_pass = False
    #
    
    if p_end <= p_bgn:
        print ("ERROR: segment boundaries don't make sense")
        print ("point bgn = ", p_bgn)
        print ("point end = ", p_end)
        flag_pass = False
    #
    
    if not flag_pass:
        sys.exit(1)
    #
#


class Anchor(object):
    def __init__(self, chrN, bgn, end, nPET, data):
        self.chrN = chrN
        self.bgn  = bgn
        self.end  = end
        self.nPET = nPET
        self.anchordata = data # effectively, a line from the bed file
    #
    
    def showAnchor(self):
        return self.anchordata
    #
    
    def __str__(self):
        s = "%5s   %8d    %8d   %5d" % (self.chrN, self.bgn, self.end, self.nPET)
        return s
    #
    
    def __repr__(self):
        return self.__str__()
    #
    
#

class Domain(object):
    def __init__(self, chrN, fixed = False):
        self.domain_fixed = fixed
        self.chrN = chrN
        self.gmin = 1e308 # ~inf
        self.gmax = -1
        self.length = 0 # seems it is not used presently (200316)
        
        # I think this "atags" is not necessary, but I'm not quite
        # ready to delete it either.
        self.atags = {}
        self.use_atags = True
        
        self.anchors = [] # [Anchor_1, Anchor_2, .... Anchor_n]
        self.resolution = float(5000) # default is 5 kbp
        self.mtools              =  HeatMapTools() # default setup GenerateHeatMapTools()
        self.mtools.ctcf_tthresh =  10
        self.mtools.ctcf_cthresh =  40
        self.mtools.from_Nenski  =  False
        self.gmap = []
    #
    
    def showDomain(self, include_headers = True):
        s = ''
        if include_headers:
            s = "%s  ---------------------------------------------------------------------------------\n" \
                % (self.showDomainRange())
        #
        
        for a in self.anchors:
            bgn, end = self.locate_on_Map(a.bgn, a.end)
            s += "%s  %6d  %6d\n" % (a.showAnchor(), bgn, end)
        #|endfor
        
        if include_headers:
            # s += "length: %9d\n" % self.length  # length is not used!
            s += "^^^^^^^^^^^^^^^^^^^^^^^^^^^^  ---------------------------------------------------------------------------------\n"
        #
        
        return s
    #
    
    def showDomainRange(self):
        s = "%5s  %9d   %9d" % (self.chrN, self.gmin, self.gmax)
        return s
    #
    
    def set_resolution(self, resolution):
        if resolution <= 0:
            print ("ERROR: resolution must be positive and greater than zero")
            sys.exit(1)
        #
        
        self.resolution = float(resolution)
    #
    
    def addAnchor(self, aa):
        if not self.domain_fixed:
            # in working with existing data, we just expand the boundaries
            # of the given domain.
            if aa.bgn < self.gmin:
                # set the lower bound of the domain -- it grows as we scan
                # the region
                self.gmin = aa.bgn
            #
        
            if aa.end > self.gmax:
                # set the upper bound of the domain -- it grows as we
                # build a region
                self.gmax = aa.end
            #
        #
        
        if self.use_atags:
            if not aa in self.atags:
                self.atags.update({ aa : len(list(self.atags)) }) 
                self.anchors += [aa]
            #
            
        else:
            self.anchors += [aa]
        #
        
    #
    
    def locate_on_Map(self, p_bgn, p_end, debug = False):
        checkDomainInput(self.gmin, self.gmax, p_bgn, p_end)
        
        
        v_bgn = int( 0.5 + (float(p_bgn - self.gmin) / self.resolution ))
        # should be that p_bgn >= self.gmin
        v_end = int( 0.5 + (float((p_bgn - self.gmin) + (p_end - p_bgn)) / self.resolution ))
        # (p_end - p_bgn) defines the size of the region in the heatmap
        # (p_bgn - gmin) defines the start position.
        if debug:
            print ("gmin,gmax = %5d,%5d  p_bgn,p_end = %5d,%5d" \
                   % (self.gmin, self.gmax, p_bgn, p_end))
            print ("v_bgn,end: ", v_bgn, v_end)
        #
        return v_bgn, v_end
    #
    
    def make_Map(self, flag_display = True):
        debug_make_Map = False # True # 
        self.gmap = []
        b_bgn, b_end = self.locate_on_Map(self.gmin, self.gmax, debug_make_Map)
        if debug_make_Map:
            print ("b_bgn,b_end, gmin,gmax", b_bgn, b_end, self.gmin, self.gmax)
        #
        
        self.gmap = initialize_matrix(self.gmap, b_end + 1, 0)
        for a in self.anchors:
            bgn, end = self.locate_on_Map(a.bgn, a.end)
            flnm = "%s_%d_%d_res5kb.heat" % (a.chrN, a.bgn, a.end)
            if debug_make_Map:
                print (flnm)
                print ("bgn,end, a.bgn,a.end", bgn, end, a.bgn, a.end)
            #
            
            if os.path.exists(flnm):
                if debug_make_Map:
                    print ("shift: ", bgn)
                #
                
                # this reads in existing heatmaps
                gmtrx = self.mtools.read_MatrixFile(flnm, ["heat"],
                                                    1.0,
                                                    "make_Map",
                                                    flag_display,
                                                    bgn)
                N        = gmtrx.length
                hv       = gmtrx.heatmap
                clusters = gmtrx.clusters
                
                #read_MatrixFile(flnm, allowed_extns, rescale_wt)
                
                i = 0; j = 0
                for jv in range(0, N):
                    for iv in range(0, N):
                        i = iv + bgn
                        j = jv + bgn
                        self.gmap[i][j] = hv[iv][jv]
                    #|endfor
                    
                #|endfor
                
            else:
                print ("WARNING: cannot open %s; the file missing" % flnm)
                i = bgn
                j = end
                self.gmap[i][j] = a.nPET
            #
            
        #|endfor
        
        # this ensures that every anchor is present
        for a in self.anchors:
            bgn, end = self.locate_on_Map(a.bgn, a.end)
            print ("bgn,end = (%5d, %5d)[%5d]" % (bgn, end, a.nPET))
            self.gmap[bgn][end] = a.nPET
            self.gmap[end][bgn] = a.nPET
        #|endfor
        
    #
    
    def __str__(self):
        return self.showDomain()
    #
    
    def __repr__(self):
        return self.__str__()
    #
    
    
#


class BuildHeatMaps(object):
    def __init__(self, data):
        if not type(data).__name__ == "Data":
            print ("ERROR: input class Data does not belong to module ChromatinData")
            sys.exit(1)
        #
        
        self.debug              = False
        self.ctcfData           = data   # class Data (from ChromatinData)
        self.ctcf_cdata         = data.cdata
        self.ctcfkeylist        = data.keylist
        self.ctcfDatatype       = data.datatype
        # data.datatype: -> "type1" or "type2". I think that only the
        # bed files (type2) are allowed for these problems, but maybe
        # I am wrong.
        self.res                = "res5kb"
        self.domain             = []
        self.ctcf_missing       = []
    #
    
    def get_PET_wt(self, ky, wt):
        """@
        
        Presently, it seems that this function is not used. Maybe it
        assigns a weight if you want it relative to some reference
        weight for convergent, divergent and tandem cases.
        
        """
        ctcf1 = self.ctcf_cdata[ky].ctcf1
        ctcf2 = self.ctcf_cdata[ky].ctcf2
        wt_ctcf = wt
        if (ctcf1 == 'R' and ctcf2 == 'R') or (ctcf1 == 'L' and ctcf2 == 'L'):
            wt_ctcf *= 0.25
        elif (ctcf1 == 'L' and ctcf2 == 'R'):
            wt_ctcf *= 0.05
        #
        
        return wt_ctcf
    #
    
    
    def groupCTCFregionsInBedFiles(self, inBedFlnms, debug = False):
        """@
        
        The main job here is group CTCFs over the range that they form
        a single island (the collection of arches in a typical display
        of CTCF data in Tang et al.). 
        
        This procedure goes through the list of CTCF contact points
        listed in the input bed file (which was read and processes
        using an object of class Data in the module ChromatinData) and
        groups all the common CTCFs that fit within a particular
        island of CTCFs. The resulting groups can then be used to
        stitch together heat maps. 
        
        The initial step in the process is to verify that the given
        heat map indicated in the bed file actually exists. It is not
        clear why at this point, but some of the bed entries do not
        have any corresponding heat map files and are rejected as a
        result. The ctcf_missing information is stored in a handle
        "ctcf_missing". The ctcf_missing files cannot be used to
        stitch together a heatmap, but they might be overlapped by
        other CTCFs in the region that overlap the ctcf_missing
        segment.
        
        While the list of heatmap files is being established, we also
        group the CTCF connections together in an object of class
        Domain (here self.domain).
        
        The size of the particular domain grouping of CTCFs is decided
        by the data itself within the bed file. The borders of that
        domain are therefore only fixed after the entire bed datafile
        is scanned. It would, therefore, be possible that only one
        large kluge sized domain is generated. In the current data
        that I have from bed files, it is no so dense and thoroughly
        entwined, but that is something this approach risks.
        
        """
        debug_groupCTCFregionsInBedFiles = debug
        self.ctcf_missing = []
        
        if debug_groupCTCFregionsInBedFiles: 
            print ("groupCTCFregionsInBedFiles()")
            print (".... combining CCD info from following bed files:")
            for v in inBedFlnms:
                print ("         <- %s" % v)
            #|endfor
            
            #sys.exit(0)
        #
        
        keylist = self.ctcfkeylist
        if len(self.ctcfkeylist) == 0:
            # not sure that it really works, but it would seem that if
            # there is some data loaded, then it can actually process
            # it at this step too. I think the program should just
            # stop if it doesn't have a set keylist, but we will see.
            print ("groupCTCFregionsInBedFiles: setting keylist")
            print ("length of the dataset: ", len(self.ctcf_cdata))
            keylist = self.ctcfData.ordered_keylist()
        #
        
        print ("groupCTCFregionsInBedFiles -- length of keylist: ", len(self.ctcfkeylist))
        if debug_groupCTCFregionsInBedFiles: 
            print ("keylist:")  # this is ordered
            print (keylist)
            print ("cdata.keys():")  # this is just the dictionary data
            print (self.ctcf_cdata.keys())
        #
        
        
        if debug_groupCTCFregionsInBedFiles: 
            print ("data.chrmsm_grp.keys():")  # This gives a list of the chromosome used
            print (self.ctcfData.chrmsm_grp.keys())
            print ("groupCTCFregionsInBedFiles:listing -------------vvvvvvv")
        #
        
        for chrN in self.ctcfData.chrmsm_grp.keys():
            for segm_k in range(0, len(self.ctcfData.chrmsm_grp[chrN].chrsegment)-1):
                
                """@
                
                Here, we scan through the various entries in the bed
                file, verify that the file actually exists, and group
                it with any other files that overlap.  """
                
                chrpos = self.ctcfData.chrmsm_grp[chrN].chrsegment[segm_k]
                
                v = self.ctcfData.makekey([chrN, chrpos[0], chrpos[1]])
                nPET   = self.ctcf_cdata[v].nPET
                if debug_groupCTCFregionsInBedFiles: 
                    print (v, \
                        self.ctcf_cdata[v].ctcf1, \
                        self.ctcf_cdata[v].ctcf2, \
                        self.ctcf_cdata[v].nPET, \
                        self.ctcf_cdata[v].state["open"], \
                        self.ctcf_cdata[v].state["active"])
                #
                
                # list of complexity measure
                cc = ''
                for ccx in self.ctcf_cdata[v].cmplx:
                    cc += '%4d  ' % ccx
                #|endfor
                
                # make a string that indicates the contacts for storage 
                ss = "%5s %10d  %10d  %s  %s  %8g  %s %10.6f  %10.6f  %10.6f  %10.6f" \
                     % (chrN, chrpos[0], chrpos[1], 
                        self.ctcf_cdata[v].ctcf1, self.ctcf_cdata[v].ctcf2, 
                        self.ctcf_cdata[v].nPET,
                        cc, 
                        self.ctcf_cdata[v].state["open"],
                        self.ctcf_cdata[v].state["active"],
                        self.ctcf_cdata[v].state["A"],
                        self.ctcf_cdata[v].state["B"])
                
                
                # print (ss)
                aa = Anchor(chrN, chrpos[0], chrpos[1], nPET, ss)
                flhd = make_file_heading(chrN, chrpos[0], chrpos[1], self.res)
                
                # #######################################################
                # verify that the files exist or not. If they are
                # ctcf_missing, then at least make a record of it.
                # #######################################################
                flnm = flhd + ".heat"
                if debug_groupCTCFregionsInBedFiles: 
                    print ("evaluating: ", os.path.exists(flnm), flnm)
                #
                
                if not os.path.exists(flnm):
                    self.ctcf_missing += [ flhd ]
                    continue
                #
                
                if debug_groupCTCFregionsInBedFiles:
                    print ("top  vvvvvvvvvvv==========================")
                    print ("current aa:")
                    print (aa)
                #
                
                if len(self.domain) == 0:
                    # domain is unassigned
                    newdmn = Domain(chrN)
                    newdmn.addAnchor(aa)
                    self.domain += [newdmn]
                    if debug_groupCTCFregionsInBedFiles:
                        print ("0 making new dmn \"%s\"" % aa)
                        print ("newdmn:\n%s" % newdmn)
                    #
                    
                else:
                    """@
                    
                    We scan through the list of domains we have
                    searching for a domain that satisfies the
                    boundaries condistions we impose on chrpos_ij and
                    g_min,max
                    
                    The reason this works is because the bed file is
                    an ordered list. Therefore, aa always increments
                    in the lower index in an orderly manner.  """
                    
                    k_curr = len(self.domain)-1
                    dmn_k = self.domain[k_curr]
                    gmin = dmn_k.gmin
                    gmax = dmn_k.gmax
                    
                    if dmn_k.use_atags:
                        if aa in dmn_k.atags:
                            print (dmn_k.atags)
                            print ("skip: ", aa)
                            sys.exit(0)
                            continue
                        #
                        
                    #
                    
                    if debug_groupCTCFregionsInBedFiles:
                        print ("current dmn_k:\n%s^^^^^^^^^^^^^^^^==========================" \
                               % dmn_k)
                        #sys.exit(0)
                    #
                    
                    
                    if (chrpos[0] == gmin):
                        
                        # chrpos_i <= gmin -- means chrpos_i just
                        # exactly matches the boundary gmin
                        # irrespective of what chrpos_j turns out
                        # to be.
                        
                        dmn_k.addAnchor(aa)
                        if debug_groupCTCFregionsInBedFiles:
                            print ("1 adding \"%s\" to dmn_k" % aa)
                        #
                        
                    elif (chrpos[1] == gmax):
                        
                        # gmax == chrpos_j -- means chrpos_j just
                        # exactly matches the boundary gmax,
                        # irrespective of what chrpos_i turns out
                        # to be.
                        
                        dmn_k.addAnchor(aa)
                        if debug_groupCTCFregionsInBedFiles:
                            print ("2 adding \"%s\" to dmn_k" % aa)
                        #
                        
                    elif ((chrpos[0] <= gmin) and (gmax <= chrpos[1])):
                        
                        # chrpos_i <= gmin < gmax <= chrpos_j --
                        # means that chrpos_ij overflows both the
                        # current domain boundaries
                        
                        dmn_k.addAnchor(aa)
                        if debug_groupCTCFregionsInBedFiles:
                            print ("3 adding \"%s\" to dmn_k" % aa)
                        #
                        
                    elif ((gmin <= chrpos[0]) and (chrpos[1] <= gmax)):
                        
                        # gmin <= chrpos_i < chrpos_j <= gmax --
                        # chrpos_ij is subsummed within both the
                        # current domain boundaries
                        
                        dmn_k.addAnchor(aa)
                        if debug_groupCTCFregionsInBedFiles:
                            print ("4 adding \"%s\" to dmn_k" % aa)
                        #
                        
                    elif ((chrpos[0] <= gmin) and (gmin <= chrpos[1])):
                        
                        # chrpos_i <= gmin && gmin <= chrpos_j --
                        # means that chrpos_i overflows current
                        # domain boundaries
                        
                        dmn_k.addAnchor(aa)
                        if debug_groupCTCFregionsInBedFiles:
                            print ("5 adding \"%s\" to dmn_k" % aa)
                        #
                        
                    elif ((chrpos[0] <= gmax) and (gmax <= chrpos[1])):
                        
                        # chrpos_i <= gmax && gmax <= chrpos_j --
                        # means that chrpos_j overflows current
                        # domain boundaries
                        
                        dmn_k.addAnchor(aa)
                        if debug_groupCTCFregionsInBedFiles:
                            print ("6 adding \"%s\" to dmn_k" % aa)
                        #
                        
                    else:
                        newdmn = Domain(chrN)
                        newdmn.addAnchor(aa)
                        self.domain += [newdmn]
                        if debug_groupCTCFregionsInBedFiles:
                            print ("7 making new dmn \"%s\"" % aa)
                            print ("newdmn:\n%s" % newdmn)
                        #
                        
                            
                    #
                    
                    
                #
                
            #|endfor
            
            print ("\nResults from grouping CTCFs:")
            for kdmn in range(0, len(self.domain)):
                print ("cluster: %6d   %10d   %10d" \
                       % (kdmn, self.domain[kdmn].gmin, self.domain[kdmn].gmax))
                print ("chr        bgn         end    ctcfs      #PET ---- cplxty ----      open        active         A          B")
                print (self.domain[kdmn].showDomain())
                #print ("stop at 3 in groupCTCFregionsInBedFiles"); sys.exit(0)
            #|endfor
            
        #|endfor
        
        if debug_groupCTCFregionsInBedFiles:
            print ("^^^^^^^------------- groupCTCFregionsInBedFiles:listing")
        #
            
        if True: # debug_groupCTCFregionsInBedFiles:
            
            #########################################################
            print ("list of missing ctcf heatmaps: --------------------vvvvv")
            for mm in self.ctcf_missing:
                print ("%s" % mm)
            #
            
            print ("^^^^^-------------------- :list of missing ctcf heatmaps")
            print ("This may mean that they have no CTCFs at all....")
            #########################################################
            #print ("groupCTCFregionsInBedFiles: stop at 1"); sys.exit(0)
            
        #
        
    #
    
#    



def assemble_eheatMaps(bhm, tmpdir, debug = False):
    
    # bhm    -> class BuildHeatMaps
    # tmpdir -> class str
    # debug  -> class bool
    
    if not (type(bhm).__name__ == "BuildHeatMaps" or type(bhm).__name__ == "StitchHeatMaps"):
        # verify that the entry is the proper object
        print ("ERROR(assemble_eheatMaps): improper object call %s" % bhm.__name__)
        sys.exit(1)
    #
    
    
    show_make_Map = False
    # assemble heat maps in earnest.
    
    for kdmn in range(0, len(bhm.domain)):
        dmn_k = bhm.domain[kdmn]
        chrN = dmn_k.chrN
        gmin = dmn_k.gmin
        gmax = dmn_k.gmax
        
        if debug:
            print ("")
            print ("cluster: %6d  %5s  %10d   %10d" % (kdmn, chrN, gmin, gmax))
        #
        
        bhm.domain[kdmn].make_Map(show_make_Map)
        
        if debug:
            print (len(bhm.domain[kdmn].gmap[0]))
        #
        
        dmnflnm = "%s/%s_%d_%d.eheat" % (tmpdir, chrN, gmin, gmax)
        if debug: 
            print ("writing: vvvvvvv")
            print (dmnflnm)
            print (dmn_k.showDomain())
            print ("next is write_ExtHeatMap")
            print ("^^^^^^^^^^^^^^^^")
        #
        
        bhm.domain[kdmn].mtools.write_ExtHeatMap(dmnflnm, dmn_k, "0.0")                
        #print ("build_CTCF_map: stop at 2: "); sys.exit(0)
    #|endfor
    
    print ("finished assemble_eheatMaps")
    
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
                    if debug_main: print ("1", m)
                    k = m - 1 
                    break
                elif ft.check_ext(cl[m], EXTS):
                    if debug_main: print ("2", m)
                    k = m
                    inbedflnm += [cl[m]]
                else:
                    if debug_main: print ("3", m)
                    k = m 
                    break
                #
                
            #
            
            if debug_main: print ("inbedflnm: ", inbedflnm)
            
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
    #
    
    # now do the processing ... 
    dt = Data()
    
    # read in the bed file(s)
    for f in inbedflnm:
        if os.path.exists(f):
            dt.readfile(f)
        else:
            print ("ERROR: cannot open %s; the file missing" % f)
            sys.exit(1)
        #
        
    #|endfor
    
    dt.ordered_keylist()
    #print ("stop at main 1"); sys.exit(0)
    
    # ordered_keylist: make an ordered list in terms of chromosome
    # number and region. This was before I learned about "from
    # compilations import OrderedDict", but anyway.
    dt.display_Data(False)
    #print ("stop at main 2"); sys.exit(0)
    
    debug_fnctns = False # True # 
    bhm = BuildHeatMaps(dt)
    bhm.groupCTCFregionsInBedFiles(inbedflnm, debug_fnctns)
    tmpdir = "ctcftmp"
    #sys.exit(0)
    build_tmpDir(tmpdir)
    assemble_eheatMaps(bhm, tmpdir, debug_fnctns)
    print ("Done")
    
#

if __name__ == '__main__':
    # running the program
    main(sys.argv)

#
