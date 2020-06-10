#!/usr/bin/env python3

"""@@@

Main Program:  stitch_heatmaps.py 
               (originally called anal_CTCF.py (ANALyze CTCF))

Classes:       Anchor
               Domain
               BuildHeatMaps

Author:        Wayne Dawson
creation date: mostly 2016 and up to March 2017.
last update:   200311 some cosmetic upgrades
version:       0

Purpose:

Build eheat files out of CCDs stored in the *.bed files and the
respective heatmaps.


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
   > aseemble_heatmaps_and_CCDs.py -ff loops.CTCF.annotated.bed -dG

   where the file "loops.CTCF.annotated.bed" contains a list of data
   formatted according to Przemek. This has the format

   file format1 (loops.CTCF.annotated.withAandB.bed):

version: 1
# CTCF analysis from ChIA-PET data and epigenetic info
# active
# open 
# compartment A 
# compartment B

chromatin_tags:
chr      bgn     end    lCTCF   rCTCF  PETcnt  cmplx1  cmplx2  cmplx3   active          open            A       B 
chromatin_data:
chr1	838908	912011	R	L	11	0	4	6	0.0686291944243	0.426918183932	1.0	0.0
chr1	838908	922335	R	R	5	1	6	8	0.0601364066789	0.397329401752	1.0	0.0
chr1	838908	1000167	R	L	7	3	15	17	0.0910398799447	0.5333717808	1.0	0.0
chr1	918286	969271	R	R	5	0	4	4	0.119015396685	0.543963910954	1.0	0.0



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
from FileTools     import getHeadExt
from CUtils        import String


from ChromatinData import Data
from ChromatinData import make_file_heading

from BasicTools    import initialize_matrix
from HeatMapTools  import HeatMapTools

from assemble_heatmaps_and_CCDs import Domain
from assemble_heatmaps_and_CCDs import Anchor
from assemble_heatmaps_and_CCDs import BuildHeatMaps
from assemble_heatmaps_and_CCDs import build_tmpDir
from assemble_heatmaps_and_CCDs import assemble_eheatMaps


PROGRAM = "stitch_heatmaps.py"
EXTS    = ["bed"]
string = String(10)


"""@

loop_file.bed -- file containing references to heat maps for specific
                 regions of the chromatin. These heatmaps are stitched
                 together to generate the desired extended heatmap.


template_file.bed -- This is a reference file that defines the desired
                     regions of interest and helps explain how a given
                     region should be stitched together.

out_dir -- name of the output director


"""
# tag for printing
USAGE = ["\n\nUSAGE: %s -ff [loop_file].bed  [[loop_file2].bed ... ] -rr [template_file].bed [-out_dir]" % PROGRAM,
         "",
         "loop_file.bed -- file containing references to heat maps for specific",
         "                 regions of the chromatin. These heatmaps are stitched",
         "                 together to generate the desired extended heatmap.",
         "",
         "",
         "template_file.bed -- This is a reference file that defines the desired",
         "                     regions of interest and helps explain how a given",
         "                     region should be stitched together.",
         "",
         "out_dir -- name of the output director" ]
#


def usage():
    for uk in USAGE:
        print (uk)
    #|endfor
    
#




class StitchHeatMaps(BuildHeatMaps):
    def __init__(self, ctcfdata, refdt):
        BuildHeatMaps.__init__(self, ctcfdata)   # inherit InputSettings
        self.refData      = refdt
        self.ref_cdata    = refdt.cdata
        self.refkeylist   = refdt.keylist
        self.refDatatype  = refdt.datatype
        self.ref_missing  = []
        self.ctcf_missing = []
    #
    
    def stitchCTCFregions(self, inbedflnms, inrefflnm, debug = False):
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
        
        The reference or template file is also read and if no heat
        maps from the ctcf data can be found that match or fit in the
        region, then these missing regions are stored in ref_missing.
        
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
        debug_stitchCTCFregions = debug
        self.ctcf_missing = []
        
        if debug_stitchCTCFregions: 
            print ("stitchCTCFregions()")
            print (".... combining CCD info from following bed files:")
            for v in inbedflnms:
                print ("         <- %s" % v)
            #|endfor
            
            print ("reference AB regions file:")
            print ("         <- %s" % inrefflnm)
            #sys.exit(0)
        #
        
        
        ctcfkeylist = self.ctcfkeylist
        if len(self.ctcfkeylist) == 0:
            # not sure that it really works, but it would seem that if
            # there is some data loaded, then it can actually process
            # it at this step too. I think the program should just
            # stop if it doesn't have a set ctcfkeylist, but we will see.
            print ("stitchCTCFregions: setting ctcfkeylist")
            print ("length of the dataset: ", len(self.ctcf_cdata))
            ctcfkeylist = self.ctcfData.ordered_keylist()
        #
        
        print ("stitchCTCFregions -- length of ctcf keylist: ", len(self.ctcfkeylist))
        if debug_stitchCTCFregions: 
            print ("stitchCTCFregions:listing ctcf keys -------------vvvvvvv")
            print ("ctcfkeylist:")  # this is ordered
            print (ctcfkeylist)
            print ("ctcf-cdata.keys():")  # this is just the dictionary data
            print (self.ctcf_cdata.keys())
        
            print ("ctcfData.chrmsm_grp.keys():")  # This gives a list of the chromosome used
            print (self.ctcfData.chrmsm_grp.keys())
            print ("stitchCTCFregions:listing ctcf keys -------------^^^^^^^")
        #
        
        #print ("stop at shm 1"); sys.exit(0)
        
        
        refkeylist = self.refkeylist
        if len(self.refkeylist) == 0:
            # not sure that it really works, but it would seem that if
            # there is some data loaded, then it can actually process
            # it at this step too. I think the program should just
            # stop if it doesn't have a set ctcfkeylist, but we will see.
            print ("stitchCTCFregions: setting ctcfkeylist")
            print ("length of the dataset: ", len(self.ref_cdata))
            refkeylist = self.refData.ordered_keylist()
        #
        
        print ("stitchCTCFregions -- length of ctcf keylist:     ", len(self.ctcfkeylist))
        print ("stitchCTCFregions -- length of template keylist: ", len(self.refkeylist))
        if debug_stitchCTCFregions: 
            print ("stitchCTCFregions:listing AB keys -------------vvvvvvv")
            print ("refkeylist:")  # this is ordered
            print (refkeylist)
            print ("ab-cdata.keys():")  # this is just the dictionary data
            print (self.ref_cdata.keys())
            
            print ("data.chrmsm_grp.keys():")  # This gives a list of the chromosome used
            print (self.refData.chrmsm_grp.keys())
            print ("stitchCTCFregions:listing -------------vvvvvvv")
            print ("stitchCTCFregions:listing AB keys -------------^^^^^^^")
        #
        
        #print ("stop at shm 2"); sys.exit(0)
        
        domain_fixed = True
        for chrN in self.refData.chrmsm_grp.keys():
            for absegm_k in range(0, len(self.refData.chrmsm_grp[chrN].chrsegment)-1):
                
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                # domain is unassigned
                ab_pos = self.refData.chrmsm_grp[chrN].chrsegment[absegm_k]
                dmn_k = Domain(chrN, domain_fixed)
                dmn_k.gmin = ab_pos[0]
                dmn_k.gmax = ab_pos[1]
                dmn_k.length = ab_pos[1] - ab_pos[0]
                # Note that <class Domain> -> length is not used and
                # it is not presently clear to me if my intention was
                # to calculate the bp length of the bead length!
                
                gmin = dmn_k.gmin
                gmax = dmn_k.gmax
                
                if debug_stitchCTCFregions:
                    print ("dmn_k:\n%s" % dmn_k)
                    #sys.exit(0)
                #
                
                flag_found_CCD = False
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                for segm_k in range(0, len(self.ctcfData.chrmsm_grp[chrN].chrsegment)-1):
                    
                    """@
                    
                    Here, we really want to know the particular
                    segment that is in the bed file for the CTCF
                    data. we scan through the various entries in the
                    bed file, verify that the file actually exists,
                    and group it with any other files that overlap
                    with the AB file.  """
                    
                    chrpos = self.ctcfData.chrmsm_grp[chrN].chrsegment[segm_k]
                    
                    v = self.ctcfData.makekey([chrN, chrpos[0], chrpos[1]])
                    nPET   = self.ctcf_cdata[v].nPET
                    if debug_stitchCTCFregions: 
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
                    ctcfaa = Anchor(chrN, chrpos[0], chrpos[1], nPET, ss)
                    flhd = make_file_heading(chrN, chrpos[0], chrpos[1], self.res)
                    
                    # #######################################################
                    # verify that the files exist or not. If they are
                    # ctcf_missing, then at least make a record of it.
                    # #######################################################
                    flnm = flhd + ".heat"
                    if debug_stitchCTCFregions: 
                        print ("evaluating: ", os.path.exists(flnm), flnm)
                    #
                    
                    if not os.path.exists(flnm):
                        
                        # if the file is missing, there is nothing we
                        # can do, we cannot use this data to stitch
                        # togeether an AB region.
                        
                        self.ctcf_missing += [ flhd ]
                        continue
                    #
                    
                    if debug_stitchCTCFregions:
                        print ("current ctcfaa:")
                        print (ctcfaa)
                    #
                    
                    # now we determine if this file we just looked up
                    # fits within the range of AB requested.
                    
                    
                    
                    
                    if ((gmin <= chrpos[0]) and (chrpos[1] <= gmax)):
                        
                        """@
                        
                        Only one case is allowed here, that the heat
                        map exist within the boudaries of the refData
                        point.
                        
                        gmin <= chrpos_i < chrpos_j <= gmax --
                        
                        chrpos_ij is subsummed within both the current
                        domain boundaries """

                        flag_found_CCD = True
                        dmn_k.addAnchor(ctcfaa)
                        if debug_stitchCTCFregions:
                            print ("4 adding \"%s\" to dmn_k" % ctcfaa)
                        #
                        
                    #
                    
                    
                #|endfor
                if flag_found_CCD:
                    # only save if at least one heat map was found in
                    # the subregion
                    self.domain += [dmn_k]
                else:
                    flhd = make_file_heading(chrN, gmin, gmax, self.res)
                    self.ref_missing += [ flhd ]
                #
                
            #|endfor
            
            print ("\nResults from grouping CTCFs:")
            for kdmn in range(0, len(self.domain)):
                print ("cluster: %6d   %10d   %10d" \
                       % (kdmn, self.domain[kdmn].gmin, self.domain[kdmn].gmax))
                print ("chr        bgn         end    ctcfs      #PET ---- cplxty ----      open        active         A          B")
                print (self.domain[kdmn].showDomain())
                #print ("stop at 3 in stitchCTCFregions"); sys.exit(0)
            #|endfor
            
        #|endfor
        
        if debug_stitchCTCFregions:
            
            print ("^^^^^^^------------- stitchCTCFregions:listing")
            
            
            #########################################################
            print ("list of missing ctcf heatmaps: --------------------vvvvv")
            for mm in self.ctcf_missing:
                print ("%s" % mm)
            #
            
            print ("^^^^^-------------------- :list of missing ctcf heatmaps")
            print ("This may mean that they have no CTCFs at all....")
            print ("\n")
            print ("list of empty AB regions: --------------------vvvvv")
            for mm in self.ref_missing:
                print ("%s" % mm)
            #
            
            print ("^^^^^-------------------- :list of empty AB regions")
            
            #########################################################
            #print ("stitchCTCFregions: stop at 1"); sys.exit(0)
            
        #
        
    #
    
#


def disp_missing_files(refdt, ref_missing):
    
    sdata  = "# missing files from %s template\n" % refdt.inflnm
    sdata += "#  file name                                        chr       bgn         end         ctcf1  ctcf2    PET      C1         C2         C3        len/5kb   active    open       A      B\n"
    for mflnm in ref_missing:
        v = mflnm.split('_')
        bgn = int(v[1])
        end = int(v[2])
        sp = v[0] + '_' + string.hollerith(bgn) + '_' + string.hollerith(end)
        k = 0
        for seg in refdt.cdata[sp].segments:
            length = int(float(seg.end - seg.bgn)/5000. + 0.5)
            sv = ''
            if "active" in refdt.cdata[sp].state:
                sv += "%8.4f " % refdt.cdata[sp].state["active"]
            else:
                sv += "         "
            #
            
            if "open" in refdt.cdata[sp].state:
                sv += "%8.4f " % refdt.cdata[sp].state["open"]
            else:
                sv += "         "
            #
            
            if "A" in refdt.cdata[sp].state:
                sv += "%8.4f " % refdt.cdata[sp].state["A"]
            else:
                sv += "         "
            #
            
            if "B" in refdt.cdata[sp].state:
                sv += "%8.4f " % refdt.cdata[sp].state["B"]
            else:
                sv += "         "
            #
            
            if k == 0:
                sdata += ("%-50s  %s | %5d   %s\n" \
                          % (mflnm, \
                             refdt.cdata[sp].disp_data(), \
                             length, sv))
            else:
                sdata += ' '*50
                sdata += ("  %s | %5d   %s\n" % (refdt.cdata[sp].disp_data(), length, sv))
            #
        #|endfor
    #
    
    return sdata

#

def test1():
    # read in the reference AB bed file
    inrefflnm = "CCDs.annotated.w_active_open.bed"
    refdt = Data()
    if os.path.exists(inrefflnm):
        refdt.readfile(inrefflnm)
    else:
        print ("ERROR: cannot open %s; the file missing" % inrefflnm)
        sys.exit(1)
    #
    
    refdt.ordered_keylist()
    # ordered_keylist: make an ordered list in terms of chromosome
    # number and region. This was before I learned about "from
    # compilations import OrderedDict", but anyway.
    
    refdt.display_Data(False)
    
    ref_missing = [ "chr12_56906756_58303170_res5kb.eheat" ]
    
    sdata = disp_missing_files(refdt, ref_missing)
    print (sdata)
    
#


def main(cl):
    debug_main = True # False # 
    if len(cl) < 2:
        print ("ERROR: too few arguments")
        usage()
        sys.exit(1)
    #
    
    k = 1
    inbedflnm = []
    inrefflnm = ''
    out_dir   = 'reftmp' 
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
                elif ft.check_ext(cl[m], EXTS):
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
            if ft.check_ext(cl[k], EXTS):
                if debug_main:
                    print ("2", k)
                #
                
                inrefflnm = cl[k]
            else:
                print ("ERROR: reference file must be of type")
                if debug_main:
                    print ("3", m)
                #
                
                k = m 
                break
            #
            
            if debug_main:
                print ("inrefflnm: ", inrefflnm)
            #
            
        elif cl[k] == "-out_dir":
            k += 1
            out_dir = cl[k]
            
            if debug_main:
                print ("output directory: ", out_dir)
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

    print ("output directory:  %s" % out_dir)
    
    #print ("stop at main 1"); sys.exit(0)
    
    # now do the processing ...
    
    
    # read in the CTCF bed file(s) 
    ctcfdt = Data()
    for f in inbedflnm:
        if os.path.exists(f):
            ctcfdt.readfile(f)
        else:
            print ("ERROR: cannot open %s; the file missing" % f)
            sys.exit(1)
        #
        
    #|endfor
    
    ctcfdt.ordered_keylist()
    # ordered_keylist: make an ordered list in terms of chromosome
    # number and region. This was before I learned about "from
    # compilations import OrderedDict", but anyway.
    
    ctcfdt.display_Data(False)
    
    #print ("stop at main 2a"); sys.exit(0)
    
    # read in the reference AB bed file
    refdt = Data()
    if os.path.exists(inrefflnm):
        refdt.readfile(inrefflnm)
    else:
        print ("ERROR: cannot open %s; the file missing" % inrefflnm)
        sys.exit(1)
    #
    
    refdt.ordered_keylist()
    # ordered_keylist: make an ordered list in terms of chromosome
    # number and region. This was before I learned about "from
    # compilations import OrderedDict", but anyway.
    
    refdt.display_Data(False)
    
    #print ("stop at main 2b"); sys.exit(0)
    
    
    
    debug_fnctns = False # True # 
    bhm = StitchHeatMaps(ctcfdt, refdt)
    bhm.stitchCTCFregions(inbedflnm, inrefflnm, debug_fnctns)
    tmpdir = out_dir
    build_tmpDir(tmpdir)
    assemble_eheatMaps(bhm, tmpdir, debug_fnctns)

    refhd, ext = getHeadExt(inrefflnm)
    missing_reffls = refhd + "_missing_refFiles.dat"
    fp = open(missing_reffls, 'w')
    fp.write(disp_missing_files(refdt, bhm.ref_missing))
    fp.close()
    
    print ("Done")
    
    
#



if __name__ == '__main__':
    # running the program
    main(sys.argv)
    #test1()

#

