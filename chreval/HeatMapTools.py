#!/usr/bin/env python3

"""@@@

Main Module:    HeatMapTools.py 

Classes:        HeatMapData 
                GenerateHeatMapTools 
                HeatMapTools
                CSVmaps

Functions:      convert_to_basic_heatmap 
                (tools for analyzing ctcf weights)
                get_energydist_histogram 
                disp_energydist_histogram 

Author:        Wayne Dawson
creation date: 2016
last update:   200416 various updates and minor adjustments
version:       0.1

Purpose:

This is an important tool to be used with the chromatin analysis
programs. It originally was embedded in chreval. I have introduced at
as a separate module because it is used by many objects in this code.

Recently, I've also had to add some functions to directly analyze CTCF
weights. These are currently under the category "Functions" because
they are under development and are not yet systematic.

"""

from math import exp
from math import log

import sys
from copy import deepcopy
from FileTools  import FileTools
from FileTools  import getHeadExt
from Chromosome import Chromosome      # not used
from Chromosome import Chromatin       # not used
from Chromosome import Segment         # not used
from Chromosome import allowed_tags    # not used
from Chromosome import assign_bed_tags # this is used


####################################################################
#########################  HEATMAP TOOLS  ##########################
####################################################################

PROGRAM = "HeatMapTools.py"


def usage():
    print ("USAGE: %s -option args" % PROGRAM)
    print ("        option       arguments")
    print ("       -basic_hm     file.clust    -- read a cluster file")
    print ("       -fcsv         <file>.csv    -- format csv file")
    print ("       -o            <file>.heat   -- output from csv input")
    print ("       -opcsv_weight \"lin\"       -- linear rescaling")
    print ("                     \"exp\"       -- exponential rescaling")
    print ("       -opcsv_scale  N             -- maximum of the data set")
    print ("       -histogram <file>           -- histogram of heatmap <file>")
     
#



# This was introduced to handle different versions of heatmaps and
# types of content. Presumably, the program that creates the map would
# explain what the contents in the file are.

class HeatMapData(object):
    def __init__(self, version):
        self.version     = version
        self.resolution  = 5000
        self.position    =    0 # unset
        self.length      =    0
        self.heatmap     = []
        self.clusters    = []
        self.g_bgn         = 1e100
        self.g_end         = -1 
        self.set_HeatMapData = False
        self.fileformat = "none"
    #
    
    def set_heatmap(self, hmap):
        # should check the inputs; e.g., is it an NxN matrix?
        self.length = len(hmap[0])
        self.heatmap = hmap
        self.set_HeatMapData = True
    #
    
    def set_clusters(self, clusters):
        # should check the inputs; e.g., is it consistent with cluster
        # data? Does it have the right number of fields?, etc.
        self.clusters = clusters
    #
    
    def set_position(self, pos):
        # should check the inputs; e.g., is it a positive number?
        self.position = pos
    #
    
    def set_resolution(self, res):
        # should check the inputs; e.g., is it a positive number?
        self.resolution = res
    #
#


class GenerateHeatMapTools(object):
    # I constructed this here to provide the general starting
    # parameters for HeatMapTools.
    def __init__(self):
        self.source = "GenerateHeatMapTools"
        # PET cluster weight
        self.add_PET_wt          = False
        self.add_PET_wt_to_edges = False
        self.PETwt               = 100.0
        
        # weights for selecting out PET clusters 
        self.CTCF_scale          = 100.0
        self.ctcf_tthresh        = 0.20*self.CTCF_scale
        self.ctcf_cthresh        = 0.40*self.CTCF_scale
        
        # This started when I was confronted with the somewhat strange
        # data I got from Nenski that had counts in almost every bin
        # and very large numbers.
        
        self.pssbl_ctcf          = {}
        self.edge_ctcf           = {}
        self.from_Nenski         = False
        
        # this re-scales the data by some fraction
        self.rescale_wt          = 1.0
        
        # other handling matters
        self.allowed_extns       = ["heat", "data"]
    #
    
#




# HeatMapTools() is used to define the matrix points
class HeatMapTools(object):
    def __init__(self, cl = GenerateHeatMapTools()):
        self.debug_HeatMapTools  = False
        
        self.N                   = -1
        
        # PET cluster weight
        self.add_PET_wt          = cl.add_PET_wt
        self.add_PET_wt_to_edges = cl.add_PET_wt_to_edges
        self.PETwt               = cl.PETwt
        
        # weights for selecting out PET clusters 
        self.CTCF_scale          = cl.CTCF_scale
        self.ctcf_tthresh        = 0.20*self.CTCF_scale
        self.ctcf_cthresh        = 0.40*self.CTCF_scale
        
        # This started when I was confronted with the somewhat strange
        # data I got from Nenski that had counts in almost every bin
        # and very large numbers.
        
        self.pssbl_ctcf          = {}
        self.edge_ctcf           = {}
        
        # special considerations
        self.from_Nenski         = cl.from_Nenski
        
        
        
        # this re-scales the data by some fraction
        self.rescale_wt          = cl.rescale_wt
        
        self.allowed_extns       = cl.allowed_extns
        
        # general information
        self.hm_max = -1000.0
        self.hm_min = 1000.0
        self.wt_range = -1.0
        self.histogram = {}
        
        
        # filtering/rescaling functions
        self.rescale_wt   = 1.0
        
        
        self.flag_use_1_m_exp = False
        # This is related to problems I had with the Nencki data.
        # I though to take the integral of the exponential dye off
        
        # f(x) = (1 - exp(-b(x-1)))
        
        # but then I figured I should used the derivative
        
        # f(x) = a exp(-b(x-1)) + c.
        
        # Presently, f(x) is hard wired for using the exponential
        # weight because it makes the most physical sense. The hardest
        # thing to estimate in the first case is the exponent b.
        
        
    #
    
    def set_Nenski(self):
        self.from_Nenski = True
    #
    
    def set_CTCF_weights(self, ref_scale):
        self.ctcf_tthresh = 0.20*ref_scale
        self.ctcf_cthresh = 0.40*ref_scale
    #
    
    
    """@@@@@
    
    # ####################################################################
    # it might work to put this next function and its partner in
    # HeatMapTools. It seems like most outputs will be matrix data, so
    # it makes sense to generate similar style outputs and to have one
    # program that reads them and processes them accordingly.
    # ####################################################################
    # vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    
    This reads in the heatmap data and converts it into something that
    can be used for calculating the enthalpy with this data.

    """
    
    def read_heat(self, flnm):
        """
        read in the matrix data
        """
        
        gmtrx    = self.read_MatrixFile(flnm, self.allowed_extns, self.rescale_wt)
        self.N   = gmtrx.length
        hv       = gmtrx.heatmap
        clusters = gmtrx.clusters
        
        """@@@@
        
        
        self.pssbl_ctcf = self.pssbl_ctcf # potential CTCFs
        self.edge_ctcf  = self.edge_ctcf  # the edge CTCF
        
        ####  Distinguishing between singletons and PET clusters  ####
        
        From here, we have to decide about PET clusters. Presumably,
        when the matrix elements are something like 1 to 10, we can
        assume that the interactions are probably singletons. The CTCF
        clusters are mainly distinguished in having a much greater
        interaction frequency, basically numbers greater than 10.
        
        With respect to CTCF clusters, these are then divided into
        convergent structures ('c': the most common and therefore
        strongly bound structures) and tandem right or tandem left
        structures ('t': encountered less frequently and not as
        strongly bound).  The tandem CTCF PET clusters are assumed to
        range from about 10 to 40. The convergent CTCF PET clusters
        should have numbers greater than 40.
        
        Presently, the program only distinguishes the PET cluster by
        the magnitude of the information. Since the tandem right and
        left structures have the same (or similar) frequency, these
        are assigned the label 't' in this section for "tandem". If
        information can be obtained elsewhere about the PET cluster,
        this can be assigned later.  In general, this must be deduced
        from the sequence alone.
        
        Therefore, presently, all I can do is add a second set of
        labels that specify the character of the
        information. Presently, this is 's' for "singleton", 'c' for
        "convergent" and 't' for "tandem" CTCF PET clusters.
        
        # location of the CTCF clusters 
        
        """
        cv = [] # nothing present
        
        # here, we can search the file using the rules above.
        
        if self.debug_HeatMapTools:
            print (self.disp_fmatrix(hv, "enthalpy: "))
        #

        """@@@@
        
        #####  READJUSTMENT FOR PET CLUSTER WEIGHTS  #####
        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        
        Here, I presume that the end point it tied with CTCF (PET
        clusters). Since the general intensity of these interactions
        is 0 to 10, this is just a guess on the weight from the CTCF
        binding sites, but it comes out only to 4 kcal/mol, so setting
        it to 100 does not seem that horrible an idea in my opinion.
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        160914wkd: In my opinion, this should not be used anymore, but
        I leave it here anyway for the moment. At any rate, the reason
        that was initially put here was because there seemed to be no
        PET data, but now there generally is, so this section is
        unnecessary.
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        """
        
        if self.add_PET_wt:
            
            if self.add_PET_wt_to_edges:
                
                """@@@@@
                
                160602wkd: Initially, I thought that the location of
                the edge for the PET cluster was always at (0,N-1).
                However, it seems my understanding that this was not
                the correct. There can be some overhang in assigning
                the 1 kbp range where the PET cluster could be in bin
                (0,N-2) or bin (0,N-1). Possible other issues
                exist. Nevertheless, because it adds some additional
                options for debugging and other things, I decided to
                leave the option here.
                
                """
                
                hv[0  ][self.N-1] = self.PETwt
                hv[self.N-1][0  ] = self.PETwt
                cv = {(0, self.N-1) : self.PETwt }
                print ("Note: Additional weight added to PET cluster")
                print ("      borders at (%d,%d)" % (0, self.N-1))
                
            else:
                """@@@@
                
                160602wkd: This is probably a more accurate
                assumption, one which allows for the possibility of
                the edge being located at either (0,N-2) or (0,N-1),
                and allows further free play than that.
                
                """
                
                hv, cv = self.add_boundaryWt(self.N, hv, self.PETwt)
            # 
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        # #####  READJUSTMENT FOR PET CLUSTER WEIGHTS  ##### 
        
        if self.debug_HeatMapTools:
            print (self.disp_fmatrix(hv, "enthalpy: "))
        #
        
        return hv, cv, self.N
    #
    
    
    
    # outputs a string containing the contents of a matrix of integer
    # variables (here it matters what you put in!!!).
    def disp_imatrix(self, v,
                     name="matrix",
                     flag_x_on_dgnl = True,
                     flag_no_header = False):
        
        N = len(v)
        s = name + ": %d\n" % N
        for j in range(0,N):
            sj = ''
            
            for i in range(0,N):
                if i == j:
                    if flag_x_on_dgnl:
                        sj += "    x "
                    else:
                        sj += " %4d " % v[j][i]
                    #
                    
                else:
                    sj += " %4d " % v[j][i]
                #
                
            #|endfor
            
            s += sj + '\n'
        #|endfor
        
        return s
    #
    
    
    
    # outputs a string containing the contents of a matrix of integer
    # variables (here it matters what you put in!!!).
    def make_heatmap(self, v, flag_no_header = False):
        N = len(v)
        s = ''
        if not flag_no_header:
            s += "%d\n" % N
        #
        
        for j in range(0,N):
            sj = ''
            if not len(v[j]) == N:
                print ("ERROR(HeatMapTools.make_heatmap): undiscernable matrix dimesions.")
                print ("     matrix size %d, for column %d: %d" % (N, j, len(v[j])))
                sys.exit(1)
            #
            
            for i in range(0,N):
                if i == j:
                    sj += "0\t"
                else:
                    sj += "%d\t" % v[j][i]
                #
                
            #|endfor
            
            s += sj + '\n'
        #|endfor
        
        return s
    #
    
    
    # Outputs a string containing the contents of a matrix of float
    # variables (here it matters what you put in!!!).
    def disp_fmatrix(self, mtrx, name="matrix",
                     flag_x_on_dgnl = True,
                     flag_no_header = False):
        
        N = len(mtrx)
        s = ''
        if not flag_no_header:
            s += name + "  %d\n" % N
        #
        
        for j in range(0,N):
            for i in range(0,N):
                if i == j:
                    if flag_x_on_dgnl:
                        s += "        x "
                    else:
                        s += " %8.4g " % mtrx[j][i]
                    #
                    
                else:
                    s += " %8.4g " % mtrx[j][i]
                #
                
            #|endfor
            
            s += '\n'
        #|endfor
        
        return s
    #
    
    
    # This simply reads the heatmap data (raw and unfiltered). In the
    # diagonal region, it deletes anything in the diagonal regions.
    def read_heatmap(self,
                     flnm,
                     EXTS = ["heat", "eheat", "csv"],
                     PROGRAM = "read_MatrixFile",
                     flag_display = True):
        
        self.fileformat = "none"
        debug_read_heatmap = False # True # 
        flag_display = debug_read_heatmap
        flhd, ext = getHeadExt(flnm)
        self.fileformat = ext
        
        if debug_read_heatmap:
            print ("file name: %s" % flnm)
        #
        
        ft = FileTools()
        if not ft.check_ext(flnm, EXTS):
            print ("ERROR: %s only accepts files with extention '%s'." % (PROGRAM, EXTS))
            print ("       input file '%s' --> (extension(s) '%s')" % (flnm, ext))
            sys.exit(1)
        #
        
        try:
            fp = open(flnm, 'r')
        except IOError:
            print ("ERROR: cannot open file '%s'." % flnm)
            sys.exit(1)
        #
        
        lfp = fp.readline()
        fileInfoLine = lfp.strip().split(':')
        fp.close()
        
        # Determine what file type is being submitted.
        if len(fileInfoLine) > 1 and (ext == "eheat"):
            # means it is probably the extended format, unless by
            # chance this is the first line after the editing step.
            if fileInfoLine[0] == "file version":
                self.fileformat = "eheat"
            else:
                print ("ERROR: %s does not appear to be a recognizable file type" % flnm)
                print ("       first line: '%s'" % fileInfoLine)
                sys.exit(1)
            #
            
        elif ext == "heat" or ext == "clust" or ext == "cpif":
            if debug_read_heatmap:        
                print ("attempt to read basic format heatmap file")
            #
            
            try:
                if ext == "heat":
                    n = int(fileInfoLine[0])
                else:
                    s = fileInfoLine[0].strip().split()
                    n = int(s[1])
                #
                
            except(ValueError):
                print ("ERROR: first line of %s should contain an integer" % flnm)
                print ("       first line: ", fileInfoLine)
                sys.exit(1)
            #
            
        elif ext == "csv":
            if debug_read_heatmap:        
                print ("attempt to read a csv format heatmap file")
            #
            
            try:
                #print (fileInfoLine[0])
                line = fileInfoLine[0].strip().split(',')
                
                n = len(line)
                #print (n, line)
            except(ValueError):
                print ("ERROR: first line of %s does not contain recognizable data" % flnm)
                print ("       first line: ", fileInfoLine)
                sys.exit(1)
            #
            
        else:
            print ("ERROR: Unrecognized file format (%s)" % flnm)
            sys.exit(1)
        #
        
        gmtrx = None
        if self.fileformat == "eheat":
            if debug_read_heatmap:
                print ("going to read_extended_heatmap()")
            #
            
            gmtrx = self.read_extended_heatmap(flnm, PROGRAM, flag_display)
        elif self.fileformat == "heat" or ext == "clust" or ext == "cpif":
            if debug_read_heatmap:
                print ("going to read_basic_heatmap()")
            #
            
            gmtrx = self.read_basic_heatmap(flnm, PROGRAM, flag_display)
        elif self.fileformat == "csv":
            if debug_read_heatmap:
                print ("going to read_csv_heatmap()")
            #
            
            gmtrx = self.read_csv_heatmap(flnm, PROGRAM, flag_display)
        #
        
        self.N = gmtrx.length
        mtrx = gmtrx.heatmap
        clusters = gmtrx.clusters
        return gmtrx # mtrx, clusters, N
    #
    
    
    def read_basic_heatmap(self, flnm, PROGRAM = "read_heatmap", flag_display = False):
        """@
        
        This is the format that I've encountered from Przemek and
        Michal, so I make this a special format different from my own
        extended heatmap format. This simply reads the heatmap data
        (raw and unfiltered). In the diagonal region, it deletes
        anything in the diagonal regions.
        
        """
        
        flhd, ext = getHeadExt(flnm)
        
        gmtrx = HeatMapData("PrzemekS")
        
        if flag_display:
            print ("file name: %s" % flnm)
        #
        
        # this program is called by read_heatmap, so there SHOULD NOT
        # be any need to do the usual checking.
        fp = open(flnm, 'r')
        lfp = fp.readlines()
        fp.close()
        
        start = 1
        # read in the number of elements in the matrix
        try:
            s = lfp[0].strip().split()
            # print (s)
            N = int(s[len(s)-1]) # read the last field from line 0
            if len(s) > 2:
                nj = len(lfp)
                ni = len(s)
                # print (ni, nj)
                if ni == nj:
                    print ("heatmap: using Michal K format")
                    N = ni
                    start = 0
                else:
                    print ("ERROR: unrecognized format for heatmap file %s" % flnm)
                    sys.exit(1)
                #
                
            #
                    
            
        except ValueError:
            print ("ERROR: first line of %s cannot be read properly. " % flnm)
            print ("       Format should be '\"some text\"    integer'; i.e., '%s %d'. ")
            print ("       Please check the file and confirm.")
            sys.exit(1)
        #
        
        if flag_display: 
            print ("size of matrix: %d" % N)
        #
        
        gmtrx.length  = N
        gmtrx.heatmap = self.read_matrix(start, N, lfp)
        
        # purge unnecessary memory
        # ########################
        del lfp
        del fp
        # ########################
        
        
        return gmtrx
    #
    
    def read_extended_heatmap(self, flnm, PROGRAM = "read_heatmap", flag_display = False):
        """@ 
        
        This reads heatmap files with the extended information; with
        extension "eheat". Therefore, not only does the it read the
        file, it reads auxiliary information such has clusters,
        resolution, etc..
        
        """
        xsflag_display = True
        
        flhd, ext = getHeadExt(flnm)
        
        
        if flag_display:
            print ("file name: %s" % flnm)
        #
        
        # this program is called by read_heatmap, so there SHOULD NOT
        # be any need to do the usual checking.
        fp = open(flnm, 'r')
        lfp = fp.readlines()
        fp.close()
        
        # get version number
        s = lfp[0].strip().split(':')
        version =  s[1].strip()
        print ("reading version: %s" % version)
        gmtrx = HeatMapData(version)
        start = 1
        N = -1
        if version == "0.0": # presently the only version
            try:
                # note
                # s = lfp[1].strip().split(':')[1]
                # ss = s.strip().split()[0])
                
                # map resolution:
                gmtrx.resolution = int(lfp[1].strip().split(':')[1].strip().split()[0])
                # start position:
                gmtrx.position   = int(lfp[2].strip().split(':')[1].strip().split()[0])
                # chain length:
                gmtrx.length     = int(lfp[3].strip().split(':')[1].strip().split()[0])
                
                if flag_display:
                    print ("resolution: %10d" % gmtrx.resolution)
                    print ("position:   %10d" % gmtrx.position)
                    print ("length:     %10d" % gmtrx.length)
                #
                
            except ValueError:
                print ("ERROR: first lines of %s are not formatted properly. " % flnm)
                for k in range(0, 5):
                    print ("line %2d:  %s" % (k, lfp[k].strip()))
                #
                print ("       Please check the file and confirm.")
                sys.exit(1)
            #
            
            # now read the heat map
            start = 5
            gmtrx.heatmap = self.read_matrix(start, gmtrx.length, lfp)
            print ("loaded matrix")
            #print (gmtrx.heatmap)
            
            k = start + gmtrx.length
            if not lfp[k].strip() == "//":
                print ("ERROR: %s format has some incompatibilities" % flnm)
                sys.exit(1)
            #
            
            k += 1
            
            # read "genome segments" part 
            try:
                # begining and end of the segment
                s = lfp[k].strip().split(':')[1]
                gmtrx.g_bgn = int(s.strip().split()[0])
                gmtrx.g_end = int(s.strip().split()[1])
            except ValueError:
                print ("ERROR: %s, cluster line after heatmap data cannot be read properly." % flnm)
                for k in range(5+gmtrx.length, len(lfp)):
                    print ("line %2d:  %s" % (k, lfp[k].strip()))
                #
                print ("       Please check the file and confirm.")
                sys.exit(1)
            #
            
            k += 1
            
            print ("loading in CTCF loop info:")
            # read taglist
            taglist = lfp[k].strip().split(':')[1].strip().split()
            #print (taglist)
            
            
            print ("chr        bgn         end       lCTCF  rCTCF    PETcnt   cmplx1     cmplx2     cmplx3        &bgn      &end")
            k += 1
            while not lfp[k].strip() == "//":
                sbedfl = lfp[k].strip().split()
                
                #print ("k = %2d, sbedfl: %s" % (k, sbedfl))

                # ignore the information at the beginning and end of this list
                if  sbedfl[1] == "---------------------------------------------------------------------------------":
                    k += 1
                elif sbedfl[3] == "---------------------------------------------------------------------------------":
                    k += 1
                    
                else:
                    chrmtn = assign_bed_tags(taglist, sbedfl)
                    gmtrx.clusters += [chrmtn]
                    print ("%s" % (chrmtn.disp_data()))
                    #print ("k = %2d, %s" % (k, chrmtn.disp_data()))
                    k += 1
                    
                #
            #|endwhile
            
        #
        
        
        
        
        # purge unnecessary memory
        # ########################
        del lfp
        del fp
        # ########################
        
        
        return gmtrx
    #
    
    
    def read_csv_heatmap(self, flnm, PROGRAM = "read_csv_heatmap", flag_display = False):
        """ This is a recent format for reading ctcf weights. """
        flhd, ext = getHeadExt(flnm)
        
        gmtrx = HeatMapData("MichalL")
        
        if flag_display:
            print ("file name: %s" % flnm)
        #
        
        # this program is called by read_heatmap, so there SHOULD NOT
        # be any need to do the usual checking.
        fp = open(flnm, 'r')
        lfp = fp.readlines()
        fp.close()
        
        start = 0
        # read in the number of elements in the matrix
        try:
            s = lfp[0].strip().split(',')
            # print ("xx: ", s)
            N = len(s) # read the last field from line 0
            # print (N)
            if len(s) > 2:
                nj = len(lfp)
                ni = len(s)
                # print (ni, nj)
                if ni == nj:
                    print ("heatmap: using Michal Lazniewski format")
                    N = ni
                    start = 0
                else:
                    print ("ERROR: unrecognized format for heatmap file %s" % flnm)
                    sys.exit(1)
                #
                
            #
            
            
        except ValueError:
            print ("       Please check the file and confirm.")
            sys.exit(1)
        #
        
        if flag_display: 
            print ("size of matrix: %d" % N)
        #
        
        gmtrx.length  = N
        gmtrx.heatmap = self.read_matrix(start, N, lfp)
        
        # purge unnecessary memory
        # ########################
        del lfp
        del fp
        # ########################
        
        
        return gmtrx
    #
    
    
    # this does the actual reading of the heatmap matrix
    def read_matrix(self, start, N, lfp):
        debug_read_martix = False # True # 
        if debug_read_martix:
            print ("entered read_matrix():")
        #
        
        s = ''
        mtrx = []
        ij_min = 100000.0
        for k in range(start,N+start):
            mtrxj = []
            if self.fileformat == "csv":
                s = lfp[k].strip().split(',')
            else:
                s = lfp[k].strip().split()
            #
            
            
            if not len(s) == N:
                # probably should make sure that the dimensions are
                # exactly equal, but anyway, they should at least be
                # such that you take some part of the matrix!
                print ("ERROR: something wrong with the dimensions of the matrix!")
                sys.exit(1)
            #
            
            j = k - start
            for i in range(0,N):
                try:
                    w = float(s[i])
                    
                    if not (w >= 0.0 or ext == "clust"):
                        # just so stupid stuff doesn't get in.
                        print ("ERROR: Encountered a negative value in the chromatin ")
                        print ("       input data. ")
                        print ("       This cannot be heat map data for chromatin!")
                        print ("       => position (i,j) = (%d,%d), value = " % (i, j), w)
                        print ("       -- all data in the heat map must be positive")
                        print ("       Please verify or correct the file '%s'." % flnm)
                        sys.exit(1)
                    #
                    
                    if w < ij_min:
                        ij_min = w
                    #
                    
                    mtrxj += [w]
                    
                except ValueError:
                    if s[i] == 'x':
                        mtrxj += [0.0]
                    else:
                        print ("ERROR: matrix contains unrecognized terms")
                        print ("       from ", s)
                        print ("       cannot recognize '%s'" % s[i])
                        sys.exit(1)
                    #
                    
                #
                
            #|endfor
            
            mtrx += [mtrxj]
            
        #|endfor
        
        # this is a patch to fix some older data
        if ij_min < 0.0 and ext == "clust":
            shift = - ij_min
            print ("unshifting free energy data by ", shift)
            for j in range(0, N):
                for i in range(0, N):
                    if not mtrx[i][j] == 0.0:
                        mtrx[i][j] += shift
                        mtrx[j][i] += shift
                    #
                    
                #|endfor
                
            #|endfor
            
        #
        
        if debug_read_martix:
            print ("finished read_matrix()")
        #
        
        return mtrx
    #
    
    
    def disp_matrix(self, mtrx, N):
        s = ''
        for j in range(0, N):
            s_j = ''
            for i in range(0, N):
                s_j += "%8.4g\t" % mtrx[i][j]
            #|endfor
            
            s += s_j + '\n'
            
        #|endfor
        
        return s
    #
    
    
    
    # This simply reads the heatmap data (raw and unfiltered). In the
    # diagonal region, it deletes anything in the diagonal regions.
    def write_BscHeatMap(self, flnm, mtrx, N, fformat = "Michal_K"):
        # presently only the most simplest format of Michal Kadlof
        flhd, ext = getHeadExt(flnm)
        if self.debug_HeatMapTools:
            print ("output file name: %s" % flnm)
        #
        
        try:
            fp = open(flnm, 'w')
        except IOError:
            print ("ERROR: cannot open file '%s'." % flnm)
            sys.exit(1)
        #
        
        if fformat == "Michal_K":
            s = self.disp_matrix(mtrx, N)
            fp.write(s)
            fp.close()
        elif fformat == "v0.0":
            fp.write("%d\n" % N) # should really be first
            s = self.disp_matrix(mtrx, N)
            fp.write(s)
            fp.close()
            
        elif fformat == "cpif":
            fp.write("cpif  %d\n" % N) # should really be first
            s = self.disp_matrix(mtrx, N)
            fp.write(s)
            fp.close()
            
        elif fformat == "clust":
            fp.write("clust  %d\n" % N) # should really be first
            s = self.disp_matrix(mtrx, N)
            fp.write(s)
            fp.close()
            
        elif fformat == "heat":
            fp.write("%d\n" % N) # should really be first
            s = self.disp_matrix(mtrx, N)
            fp.write(s)
            fp.close()
            
        else:
            fp.write("%d\n" % N) # should really be first
            s = self.disp_matrix(mtrx, N)
            fp.write(s)
            fp.close()
        #
        
    #
    
    
    # This simply reads the heatmap data (raw and unfiltered). In the
    # diagonal region, it deletes anything in the diagonal regions.
    def write_ExtHeatMap(self, flnm, domain, version = "0.0"):
        # presently only the most simplest format of Michal Kadlof
        flhd, ext = getHeadExt(flnm)
        self.debug_HeatMapTools = True
        if self.debug_HeatMapTools:
            print ("output file name: %s" % flnm)
        #
        
        try:
            fp = open(flnm, 'w')
        except IOError:
            print ("ERROR: cannot open file '%s'." % flnm)
            sys.exit(1)
        #
        
        if version == "0.0":
            N = len(domain.gmap[0])
            v_bgn = int( 0.5 + (float(domain.gmin) / domain.resolution ))           
            fp.write("file version: %s\n" % version)
            fp.write("map resolution: %10d [bp]\n" % domain.resolution)
            fp.write("start position: %10d [sgmnts]\n" % v_bgn)
            fp.write("chain length:   %10d [sgmnts]\n" % N) # should really be first
            fp.write("heatmap:\n")
            s = self.disp_matrix(domain.gmap, N)
            fp.write(s)
            s =  "//\n"
            s += "genome segment:  %10d   %10d\n" % (domain.gmin, domain.gmax)
            s += ":chr      bgn      end    lCTCF rCTCF  PETcnt  cmplx1 cmplx2 cmplx3   open    active       A        B   vi   vj\n"
            s += domain.showDomain(False) # false = don't show header
            # assemble_heatmap_and_CCDs.Domain.showDomain() prints out
            # the contents of the particular line in the bed file.
            s += "//\n"
            fp.write(s)
            fp.close()
            
        else:
            print ("Sorry, no format for version '%s' is available presently" % version)
        #
        
    #
    
    
    
    """@@@@@
    
    
    # #############################################################
    # #############################################################
    
    161005wkd: I know this is a bit silly to do every operation
    separately; therefore force multiple N^2 steps rather than do
    these steps all at once.
    
    Unfortunately, we just have so many different types of data, that
    I found it annoying to have to constantly be writing in and
    checking some new routine for each type of new data. Sometimes,
    the placement of these things (particularly the histograms) need
    to be moved around when looking at some different issue.
    Moreover, make_heatmap needs some of the followign functions in
    this current arrangement, but other functions are not used. So it
    seems better to ask the program to obtain the relevant material,
    not do all sorts of additional strange things that we don't
    want. The modularity will allow me to use what I actually need.
    
    Therefore, I finally decided to modularize this section so that
    these functions can be called in a flexible way. Presently, they
    are still located in this routine, but I think that eventually I
    will locate them as needed in the programs that call this
    function.
    
    # #############################################################
    # #############################################################
    
    
    This reads things like the heatmap data, which always involves
    elements where i != j (no diagonal elements).
    
    """
    
    def read_MatrixFile(self, flnm, EXTS = ["heat", "eheat"],
                        rescale_wt = 1.0,
                        PROGRAM = "read_MatrixFile",
                        flag_display = True,
                        bgn_shift = 0):
        
        #gmtrx = HeatMapData(version)
        gmtrx = self.read_heatmap(flnm, EXTS, PROGRAM, flag_display)
        N = gmtrx.length
        clusters = gmtrx.clusters
        # scaling the data according to the problem and data type and
        # obtain the minimum and maximum of the heatmap data.
        
        if self.from_Nenski:
            if self.flag_use_1_m_exp:
                gmtrx.heatmap, self.hm_min, self.hm_max, self.wt_range \
                    = self.filter_HiC_using_1_m_exp(gmtrx.heatmap, N)
            else:
                gmtrx.heatmap, self.hm_min, self.hm_max, self.wt_range \
                    = self.filter_HiC_using_exp(gmtrx.heatmap, N)
            #
            
        #
        
        # full rescale relative to a self.rescale_wt
        self.rescale_wt = rescale_wt
        gmtrx.heatmap = self.rescale_mtrx(gmtrx.heatmap, N, self.rescale_wt)
        
        # obtain the final minimum and maximum weights
        self.hm_min, self.hm_max, self.wt_range = find_minmax(gmtrx.heatmap, N, True)
        # print (self.hm_min, self.hm_max, self.wt_range)
        
        # display a histogram of the data
        flag_use_minmax_hist = False
        flag_use_avg_hist    = False
        if flag_use_minmax_hist:
            self.histogram = self.make_histogram_maxmin(gmtrx.heatmap, N, "histogram_minmax.dat")
        #
        
        if flag_use_avg_hist:
            self.histogram = self.make_histogram_avg(gmtrx.heatmap, N, "histogram_avg.dat")
        #
        
        # find the CTCF points in the data
        if len(clusters) > 0:
            print ("used read in PET cluster information: %s" % flnm)
            self.set_CTCF_points(gmtrx, True) # flag_display
        else:
            print ("estimate the PET cluster information: %s" % flnm)
            self.estimate_CTCF_points(gmtrx.heatmap, N, flag_display, bgn_shift)
        #
        
        return gmtrx
    #
    
    
    # This reads things like the heatmap data, which always involves
    # elements where i != j (no diagonal elements).
    def read_MatrixFile_wt(self, flnm,
                           EXTS = ["heat", "eheat"],
                           PROGRAM = "read_MatrixFile"):
        
        gmtrx = self.read_heatmap(flnm, EXTS, PROGRAM)
        N     = gmtrx.length
        
        
        if self.from_Nenski:
            if self.flag_use_1_m_exp:
                gmtrx.heatmap, hm_min, hm_max, wt_range \
                    = self.filter_HiC_using_1_m_exp(gmtrx.heatmap, N)
            else:
                gmtrx.heatmap, hm_min, hm_max, wt_range \
                    = self.filter_HiC_using_exp(gmtrx.heatmap, N)
            #
            
        #
        
        
        # reweight the data to logarithmic weights
        gmtrx.heatmap, self.hm_min, self.hm_max, self.wt_range \
            = self.reweight_ln(gmtrx.heatmap, N)
        
        print ("hm_(min,max): ", self.hm_min, self.hm_max)
        
        return gmtrx
    #
    
    # find the CTCF points in the data
    def estimate_CTCF_points(self, mtrx, N, flag_display = True, bgn_shift = 0):
        irh = N; jrh = 0
        # identifies the largest region supported passing threshold
        # PET count. Maybe the purpose is not as important as it once
        # seemed, but I wanted to know the largest region encompassing
        # suffiently large PET counts.
        
        rhlist = {}
        # marks CTCF-like elements found directly _in_ the heat map.
        
        
        mx = 0.0;  i_mx = 0;  j_mx = 0
        # marks the larges PET count weight, wherever it is found
        
        for j in range(0,N):
            for i in range(0,N):
                w = mtrx[i][j]
                
                if w >= self.ctcf_tthresh:
                    # largest CTCF span
                    
                    if i < j:
                        rhlist.update({(i, j): w})
                        # this will help identify the maximum
                        # domain size for CTCF islands
                        if i < irh:
                            irh = i
                        #
                        
                        if j > jrh:
                            jrh = j
                        #
                        
                    #
                    
                #
                
                if w > mx: # largest PET count weight
                    if i < j:
                        mx = w
                        i_mx = i
                        j_mx = j
                    #
                    
                #
                
            #|endfor
            
        #|endfor
        
        if flag_display:
            print ("maximum matrix element at (%d,%d)[value= %8.3f]" % (i_mx, j_mx, mx))
        #
        
        # Filter out the domain boundary CTCF since that will
        # automatically be used for CTCF island formation. Moreover,
        # if it is the only one, there are no other potential CTCF
        # interactions from which to evaluate or build any islands.
        self.pssbl_ctcf = {}
        self.edge_ctcf = {}   # largest PET count span
        if len(rhlist) > 0:
            for rh in rhlist.keys():
                if not (rh[0] == irh and rh[1] == jrh):
                    # it is important to treat the border different
                    self.pssbl_ctcf.update({(rh[0],rh[1]) : rhlist[rh]})
                else:
                    self.edge_ctcf.update({(rh[0],rh[1]) : rhlist[rh]})
                #
                
            #|endfor
        #
        
        print (self.print_ctcf_dict("edge_ctcf ", self.edge_ctcf, bgn_shift))
        print (self.print_ctcf_dict("pssbl_ctcf", self.pssbl_ctcf, bgn_shift))
        # sys.exit(0)
        return i_mx, j_mx, mx
    #
    
    def print_ctcf_dict(self, title, cdict, bgn_shift = 0):
        s = "%12s  {" % title
        clist = list(cdict.keys())
        
        """@
        
        https://stackoverflow.com/questions/16819222/how-to-return-dictionary-keys-as-a-list-in-python
        
        to produce an indexed list of keys, you need to do
        list(dictionary.keys())
        
        """
        
        #print (s)
        #print ("clist", clist)
        if len(clist) > 1:
            for k in range(0, len(clist)-1):
                ck = clist[k]
                # print ("ck: ", ck, ck[0])
                s += "(%d, %d) : %d, " % (ck[0]+bgn_shift, ck[1]+bgn_shift, cdict[ck])
            #|endfor
            
            ck = clist[len(clist)-1]
            s += "(%d, %d) : %d}" % (ck[0]+bgn_shift, ck[1]+bgn_shift, cdict[ck])
            
        elif len(clist) == 1:
            ck = clist[0]
            s += "(%d, %d) : %d}" % (ck[0]+bgn_shift, ck[1]+bgn_shift, cdict[ck])
            
        else:
            s += " }"
        #
        
        #print (s)
        #sys.exit(0)
        return s
    #
    
    # find the CTCF points in the data
    def set_CTCF_points(self, gmtrx, flag_display = True):
        clusters = gmtrx.clusters
        mtrx = gmtrx.heatmap
        N = gmtrx.length
        
        i_mx = 0; j_mx = 0; mx    = 0.0
        irh  = N; jrh  = 0
        rhlist = {}
        
        # Filter out the domain boundary CTCF since that will
        # automatically be used for CTCF island formation. Moreover,
        # if it is the only one, there are no other potential CTCF
        # interactions from which to evaluate or build any islands.
        self.pssbl_ctcf = {}
        self.edge_ctcf = {}
        if len(clusters) > 0:
            k_mx   = -1
            for k in range(0, len(clusters)):
                vi = int(clusters[k].vi); vj = int(clusters[k].vj);
                nPET = clusters[k].nPET
                
                rhlist.update({(vi, vj): nPET})
                if nPET > mx:
                    mx = nPET
                    i_mx = int(vi); j_mx = int(vj); mx = nPET
                    k_mx   = k
                #
                
                # this will help identify the maximum
                # domain size for CTCF islands
                if vi < vj:
                    # this will help identify the maximum
                    # domain size for CTCF islands
                    if vi < irh:
                        irh = int(vi)
                    #
                    
                    if vj > jrh:
                        jrh = int(vj)
                    #
                    
                #
                
            #|endfor
            
            if flag_display:
                print ("maximum matrix element at  (%5d,%5d)[weight=   %8.1f]" % (i_mx, j_mx, mx))
                print ("total CTCF domain coverage (%5d,%5d)[length= %8d]" % (irh, jrh, (jrh-irh+1)))
            #
            
            # sometimes the boundary  element is missing
            #if not rhlist.has_key((irh,jrh)):
            if not ((irh,jrh) in rhlist):
                self.edge_ctcf.update({(irh,jrh): 5})
            #
            
            for k in range(0, len(clusters)):
                vi = int(clusters[k].vi); vj = int(clusters[k].vj); nPET = clusters[k].nPET
                if (vi == irh and vi == jrh):
                    # if k == k_mx: 
                    self.edge_ctcf.update({ (vi, vj) : nPET})
                else:
                    # it is important to treat the border different
                    self.pssbl_ctcf.update({ (vi, vj) : nPET})
                #
                
            #|endfor
            
        #
        
        print ("edge_ctcf:  ", self.edge_ctcf)
        print ("pssbl_ctcf: ", self.pssbl_ctcf)
        # sys.exit(0)
        return i_mx, j_mx, mx
    #
    
    
    
    
    # display a histogram of the matrix data: here it is the average
    # for a given |j-i|
    def make_histogram_avg(self, mtrx, N, flnm):
        histogram = {}
        for j in range(0,N):
            for i in range(0,N):
                w = mtrx[i][j]
                if i < j:
                    if (j-i) in histogram:
                        histogram[j-i] += [w]
                    
                    else:
                        histogram.update( { (j-i) : [w] })
                    #
                    
                #
                
            #|endfor
            
        #|endfor
        
        fp = open(flnm, 'w')
        
        for k in histogram.keys():
            wf = 0.0
            for w in histogram[k]:
                wf += w
            #
            
            wf /= float(len(histogram[k]))
            fp.write("%4d   %12.2f\n" % (k, wf))
            
        #|endfor
        
        fp.close()
        
        return histogram
    #
    
    # display a histogram of the matrix data: here it is the minimum
    # and maximum for a given |j-i|
    def make_histogram_maxmin(self, mtrx, N, flnm):
        histogram = {}
        for j in range(0,N):
            for i in range(0,N):
                w = mtrx[i][j]
                if i < j:
                    jmi = j - i
                    if jmi in histogram:
                        if histogram[jmi][0] > w:
                            histogram[jmi][0] = w
                        #
                        
                        if histogram[jmi][1] < w:
                            histogram[jmi][1] = w
                        #
                    
                    else:
                        histogram.update( { (jmi) : [w,w] })
                    #
                    
                #
                
            #|endfor
            
        #|endfor
        
        fp = open(flnm, 'w')
        for k in histogram.keys():
            w_min = histogram[k][0]
            w_max = histogram[k][1]
            fp.write("%4d   %12.2f   %12.2f\n" % (k, w_min, w_max))
        #|endfor
        
        fp.close()
        
        return histogram, w_min, w_max
    #
    
    
    
    """
    Scaling Hi-C data using f(|j-i|) = 1 - exp[-b(|j-i|-1)].  This is
    used to remove the homopolymer noise from the Hi-C data to make it
    possible to process the data using chreval and related
    programs. This should only be used with LINEAR DATA.  It should
    NOT be used with ChIA-PET data.

    """
    def filter_HiC_using_1_m_exp(self, mtrx, N):
        """@
        
        now settle on the data type
        
        Here, hm_max, hm_min, wt_range and mtrx are treated as local
        variables because there is the possibility that one does not
        wish to destroy the raw input data from the heatmap.
        
        """
        hm_max = -1000.0
        hm_min = 1000.0
        for j in range(0,N):
            for i in range(0,N):
                w = mtrx[i][j]*(1.0 - exp(-0.0035*(float(abs(j-i))-1.0)))
                w = int(w + 0.5)
                mtrx[i][j] = w
                
                if w > hm_max:
                    hm_max = w
                #
                
                if w < hm_min:
                    hm_min = w
                #
                
            #|endfor
            
        #|endfor
        
        wt_range = hm_max
        # this is LINEAR data; i.e., it ranges between 0 to hm_max
        return mtrx, hm_min, hm_max, wt_range
    #
    
    """
    Scaling Hi-C data using 1/{Nexp[-b(|j-i|-1)] + c}.  This is used
    to remove the homopolymer noise from the Hi-C data to make it
    possible to process the data using chreval and related
    programs. This should only be used with LINEAR DATA. It should NOT
    be used with ChIA-PET data.
    
    """
    def filter_HiC_using_exp(self, mtrx, N):
        """@
        
        now settle on the data type
        
        Here, hm_max, hm_min, wt_range and mtrx are treated as local
        variables because there is the possibility that one does not
        wish to destroy the raw input data from the heatmap.
        
        """
        
        hm_max = -1000.0
        hm_min =  1000.0
        for j in range(0,N):
            for i in range(0,N):
                w = mtrx[i][j]
                if abs(j-i) == 1:
                    # this position is rather special and particularly
                    # noisy. It is just the nearest neighbor
                    # interactions, so it is meaningless as well.
                    w = 0.0
                else:
                    w = w/(4451*exp(-0.374*float(abs(j-i))-1.0) + 55.1 )
                    w = int(w + 0.5)
                #
                
                mtrx[i][j] = w
                
                if w > hm_max:
                    hm_max = w
                #
                
                if w < hm_min:
                    hm_min = w
                #
                
            #|endfor
            
        #|endfor
        
        wt_range = hm_max 
        # this is LINEAR data; i.e., it ranges between 0 to hm_max
        return mtrx, hm_min, hm_max, wt_range
    #
    
    
    
    # full rescale relative to a self.PET_range
    def rescale_mtrx(self, mtrx, N, rescale_wt):
        # this guarantees that there will be no problems on whatever
        # data is rescaled.
        for j in range(0,N):
            for i in range(0,N):
                w = mtrx[i][j]
                w = int(w*rescale_wt+0.5)
                mtrx[i][j] = w
            #|endfor
            
        #|endfor
        
        return mtrx
    #
    
    
    def reweight_ln(self, mtrx, N):
        """@
        
        I think this could be done by hm_max = log(hm_max) and
        hm_min = log(hm_max). Anyway, it also has to be evaluated
        using find_minmax(), so either way, this step must be done
        somewhere.
        
        """
        hm_max = -1000.0
        hm_min = 1000.0
        v = 0.0
        for j in range(0,N):
            for i in range(0,N):
                w = mtrx[i][j]
                if not w == 0.0:
                    if w == 1:
                        v = 0.3
                    else:
                        v = log(w)
                    #
                    
                #
                
                if v > hm_max:
                    hm_max = v
                #
                
                if v < hm_min:
                    hm_min = v
                #
                
            #|endfor
            
        #|endfor
        
        wt_range = hm_max - hm_min
        if hm_max == hm_min:
            # this can happen if all the values in the matrix are 1,
            # for example.
            wt_range = 1.0
            print ("reseting wt_range = %.2f" % wt_range)
        #
        
        # reweight the data
        if hm_min < 0.0:
            # Data is from a *.clust file [min,max]. Here, min is
            # usually negative, max could also be negative.
            for j in range(0,N):
                for i in range(0,N):
                    w = mtrx[i][j]
                    if w <= 0.0:
                        w = hm_min / wt_range
                    else:
                        w = log(w) / wt_range
                    #
                    
                    mtrx[i][j] = w
                #|endfor
                
            #|endfor
            
        else:
            # Data is from a *.heat file [0, max]. Essentially
            # integers typically ranging between 0 and say 100, though
            # max can be larger than 100
            for j in range(0,N):
                for i in range(0,N):
                    w = mtrx[i][j] 
                    if w == 0.0:
                        w = 0.0 / wt_range
                    elif w == 1.0:
                        w = 0.3 / wt_range
                    else:
                        w = log(w) / wt_range
                    #
                    
                    mtrx[i][j] = w
                    
                #|endfor
                
            #|endfor
            
        #
        
        return mtrx, hm_min, hm_max, wt_range
    #
    
    
    def add_boundaryWt(self, N, mtrx, boundaryWt):
        i_min = 0; j_max = N-1
        
        if mtrx[i_min][j_max] == 0.0:
            if self.debug_HeatMapTools:
                print ("adjusting borders (%d,%d)" % (i_min, j_max))
            #
            
            j_max = 0
            i_min = N
            for j in range(N/2, N):
                for i in range(N/2, -1, -1):
                    if mtrx[i][j] > 0.0:
                        if self.debug_HeatMapTools:
                            print ("ij = ", i, j)
                        #
                        
                        if i <= i_min and j >= j_max:
                            i_min = i
                            j_max = j
                        #
                        
                    #
                    
                #|endfor
                
            #|endfor
        
        else:
            print ("no need to adjust the location of the border")
        #
        
        mtrx[i_min][j_max] = boundaryWt
        mtrx[j_max][i_min] = boundaryWt
        print ("Note: Additional weight added to PET cluster, borders at (%d,%d)" \
            % (i_min, j_max))
        flag_check = False #True
        if flag_check:
            sys.exit(0)
        #
        
        return mtrx, {(i_min, j_max): boundaryWt }
    #
    
    def normalize_matrix(self, mtrx, optimum = "neg"):
        # optimum = neg || pos: free energy is negative, scores are
        # often postive
        
        if not (optimum == "neg" or optimum == "pos"):
            print ("ERROR: option 'optimum' requires either a 'pos' or 'neg' (default)")
            print ("       Don't know what you want!")
            sys.exit(1)
        #
        
        N = len(mtrx)
        
        
        # read in the matrix data
        m_max = -1e10
        m_min =  1e10
        for j in range(0,N):
            for i in range(0,N):
                if not (i == j):
                    if mtrx[i][j] < m_min:
                        m_min = mtrx[i][j]
                    #
                    
                    if mtrx[i][j] > m_max:
                        m_max = mtrx[i][j]
                    #
                    
                else:
                    if not (mtrx[i][j] == 0.0):
                        """                        
                        I have to think about this, but I think we
                        should not allow digonal elements at least in
                        this construction.  Maybe it needs to be
                        renamed, but anyway, we the purpose of this
                        tool is to analyze matrices and these matrices
                        should all be ones that do not have diagonal
                        elements.
                        
                        """
                        print ("ERROR: diagonal element (%d,%d) is not zero!" % (i,j))
                        sys.exit(1)
                    #
                    
                #
                
            #|endfor
            
        #|endfor
        
        wt = m_max - m_min
        
        if optimum == "neg":
            wt = - wt
        #
        
        for j in range(0,N):
            for i in range(0,N):
                mtrx[i][j] /= wt
            #|endfor
            
        #|endfor
        
        return mtrx
    #

# end class HeatMapTools



def convert_to_basic_heatmap(flnm):
    """
    removes any artifacts that are in some heatmaps that were
    generated for visualization and produces the most crude form --
    the heatmap and nothing else.
    
    """
    mtools = HeatMapTools()
    gmtrx = mtools.read_heatmap(flnm, 
                                ["heat","eheat","clust"], "HeatMapTools")
    N        = gmtrx.length
    mtrx     = gmtrx.heatmap
    clusters = gmtrx.clusters
    
    # construct a new file name from the old one
    flhd, ext = getHeadExt(flnm)
    flnm_rev = flhd + "_new." + ext
    print ("new file name: ", flnm_rev)
    
    mtools.write_BscHeatMap(flnm_rev, mtrx, N, "Michal_K")
#


        


def get_energydist_histogram(hm, span = 10.0):
    N = len(hm)
    histogram = {}

    hm_max = -1.0
    for j in range(0,N):
        for i in range(0,j):
            if not hm[i][j] == 0:
                k = int(hm[i][j]/span)
                if hm[i][j] > hm_max:
                    hm_max = hm[i][j]
                #
                
                if k in histogram:
                    histogram[k][0] += 1
                    histogram[k][1] += hm[i][j]
                else:
                    histogram.update({k : [1, hm[i][j]] })
                #
                
            #
            
        #|endfor
        
    #|endfor

    hist_range = int(hm_max/span) + 1
    
    for k in range(0, hist_range):
        if k in histogram:
            histogram[k][1] /= float(histogram[k][0])
        #
        
    #|endfor
    
    print (disp_energydist_histogram(histogram, hist_range, span))
    return histogram, span
#


def disp_energydist_histogram(hist, hist_range, span = 10.0):
    s = "energy distribution histogram:\n"
    s += "   energy   counts     average\n"
    for k in range(0, hist_range):
        if k in hist:
            s +=  "%8.2f    %4d     %8.2f\n" \
                % (span*(0.5 + k), hist[k][0], hist[k][1])
        else:
            s += "%8.2f       0\n" % (span*(0.5 + k))
        #
        
    #|endfor
    
    return s
#


"""# 

This can be used anywhere after the matrix is built to scans through
the matrix elements and find the minimum and maximum value.

"""
def find_minmax(mtrx, N, flag_linear):
    if not N == len(mtrx):
        print ("ERROR: matrix length not matching with proposed size")
        print ("matrix length(%d)  vs N(%d)" % (len(mtrx), N))
        sys.exit(1)
    #
    
    """@
    
    now settle on the data type
    
    Here, hm_max, hm_min, wt_range and mtrx are treated as local
    variables because there is the possibility that one does not wish
    to destroy the raw input data from the heatmap.

    """
    hm_max = -1000.0
    hm_min =  1000.0
    for j in range(0,N):
        for i in range(0,N):
            w = mtrx[i][j]
            
            if w > hm_max:
                hm_max = w
            #
            
            if w < hm_min:
                hm_min = w
            #
            
        #|endfor
        
    #|endfor
    
    if flag_linear:
        # In THESE problems, the linear data is in the form of
        # ChIA-PET or Hi-C data and is ALWAYS expressed from 0 to
        # hm_max
        wt_range = hm_max
    else:
        # Presently, I use this when I rescale the linear data to
        # logarithmic data.
        wt_range = hm_max - hm_min
    #
    
    return hm_min, hm_max, wt_range
#


class CSVmaps(HeatMapTools):
    def __init__(self, method = "exp", scale = 50.0):
        HeatMapTools.__init__(self) # inherit HeatMapTools
        self.factor     = 300.0
        self.resolution = 5. # kb
        self.scale      = scale
        self.cutoff     = 1.0
        self.method     = method
        if self.method == "lin":
            self.cutoff = 1.0 
        #
        
        self.hm_max     = -100000.0
        self.hm_min     =  100000.0
        self.N          = -1
    #
    
    def process_csv_file(self, inflnm, outflnm, factor = 300.0):
        #a = HeatMapTools()
        gmtrx = self.read_heatmap(inflnm)
        self.N = gmtrx.length
        
        self.factor = factor
        
        if 0:
            #print (a.disp_matrix(gmtrx.heatmap, N))
            mn, mx, rng = find_minmax(gmtrx.heatmap, self.N, True)
            print ("min, max, range: ", mn, mx, rng)
            histogram1, w_min, w_max \
                = self.make_histogram_maxmin(gmtrx.heatmap,
                                             self.N,
                                             "histogram_minmax.dat")
            histogram2 \
                =  self.make_histogram_avg(gmtrx.heatmap,
                                           self.N,
                                           "histogram_avg.dat")
            histogram3, span \
                = get_energydist_histogram(gmtrx.heatmap, 5.0)
        #
        
        new_hm = []
        if self.method == "lin":
            new_hm = self.slice_and_dice(gmtrx.heatmap)
        else:
            new_hm = self.reweight_and_scale(gmtrx.heatmap)
        #
        
        get_energydist_histogram(new_hm, 2.0)
        self.write_BscHeatMap(outflnm, new_hm, self.N, "v0.0")
        
    #
    
    def reweight_and_scale(self, hm):
        # first simply purge everything less than the cutoff.
        self.N = len(hm)
        mn, mx, rng = find_minmax(hm, self.N, True)
        self.factor = mx / log(50.0)
        
        self.hm_max = -100000.0
        self.hm_min =  100000.0
        
        hmx = deepcopy(hm)
        for j in range(0,self.N):
            for i in range(0,self.N):
                hmx[i][j] = exp(hm[i][j]/self.factor)
                
                if hmx[i][j] > self.hm_max:
                    self.hm_max = hmx[i][j]
                #
                
                if hmx[i][j] < self.hm_min:
                    self.hm_min = hmx[i][j]
                #
                
            #|endfor
            
        #|endfor
        
        print ("after exp reweighting:")
        print ("hm_max: ", self.hm_max)
        print ("hm_min: ", self.hm_min)
        print (mx, self.factor)
        #sys.exit(0)
        
        # now rescale the data to a maximum value scale
        wt = self.scale/self.hm_max
        print ("wt:     ", wt)
        print ("cutoff: ", self.cutoff)
        print ("scale:  ", self.scale)
        print ("factor: ", self.factor)
        
        for j in range(0,self.N):
            for i in range(0,self.N):
                hmx[i][j] *= wt
                hmx[i][j] = float(int( (10.0*hmx[i][j]+0.5)/10.0))
                if hmx[i][j] < self.cutoff:
                    # print ("cut(%2d,%2d): %10.4g, wt(%.2f)" % (i, j, hmx[i][j], wt))
                    hmx[i][j] = 0.0
                #
                
            #|endfor
            
        #|endfor
        
        get_energydist_histogram(hmx, 2.0)
        #sys.exit(0)
        return hmx
    #
    
    
    
    def slice_and_dice(self, hm):
        
        # first simply purge everything less than the cutoff.
        self.N = len(hm)
        mn, mx, rng = find_minmax(hm, self.N, True)
        
        
        self.hm_min =  1000000.0
        self.hm_max = -1000000.0
        hmx = deepcopy(hm)
        for j in range(0,self.N):
            for i in range(0,self.N):
                
                # hmx[i][j] = hmx[i][j]**2
                if hmx[i][j] > self.hm_max:
                    self.hm_max = hmx[i][j]
                #
                
                if hmx[i][j] < self.hm_min:
                    self.hm_min = hmx[i][j]
                #
                
            #|endfor
            
        #|endfor
        
        # now rescale the data to a maximum value scale
        
        wt = self.scale/self.hm_max
        print ("wt = ", wt)
        for j in range(0,self.N):
            for i in range(0,self.N):
                hmx[i][j] *= wt
                hmx[i][j] = float(int( (10.0*hmx[i][j]+0.5)/10.0))
                if hmx[i][j] < self.cutoff:
                    hmx[i][j] = 0.0
                #
                
            #|endfor
            
        #|endfor
        
        return hmx
    #
    
#


def main(cl):
    n = len(cl)
    if  n == 1:
        usage()
        sys.exit(0)
    #
    inflnm = ''
    in_option = ''
    outflnm = 'test.heat'
    out_option = False
    csv_weight = "exp"
    csv_scale  = 50.0
    op_hist    = False
    k = 1
    while k < len(cl):
        arg = cl[k]
        if arg == '-basic_hm':
            in_option = "convert_to_basic_hm"
            k += 1
            if k < n:
                inflnm = cl[k]
                print (inflnm)
            #
        
        elif arg == '-fcsv':
            in_option = "format_csv_heatmap_file"
            k += 1
            if k < n:
                inflnm = cl[k]
                print (inflnm)
            #
            
        elif arg == '-opcsv_weight':
            k += 1
            if k < n:
                csv_weight = cl[k]
                print (csv_weight)
            #
            
        elif arg == '-opcsv_scale':
            k += 1
            if k < n:
                try:
                    csv_scale = float(cl[k])
                except:
                    print ("ERROR(HeatMapTools) opcsv_scale requires a float variable")
                    print ("                    (%s) was entered" % cl[k])
                #
                
                print (csv_scale)
            #
            
        elif arg == '-o':
            out_option = True
            k += 1
            if k < n:
                outflnm = cl[k]
                print (outflnm)
            #
            
        elif arg == '-histogram':
            op_hist = True
            in_option = "histogram"
            k += 1
            if k < n:
                inflnm = cl[k]
                print (inflnm)
            #
            
        elif arg == '-h' or arg == "--help":
            usage()
            sys.exit(0)
        else:
            print ("unrecognized argument at %d" % k)
            print (cl)
            usage()
            sys.exit(1)
        #
        
        k += 1
    #
    
    if in_option == "convert_to_basic_hm":
        convert_to_basic_heatmap(inflnm)
    
    elif in_option == "format_csv_heatmap_file":
        b = CSVmaps(csv_weight, csv_scale)
        b.process_csv_file(inflnm, outflnm, 100.0)
        
    elif op_hist:
        # construct a new file name from the old one
        mtools = HeatMapTools()
        
        gmtrx = mtools.read_heatmap(inflnm,["heat","eheat","clust", "csv"], "main")
        N        = gmtrx.length
        mtrx     = gmtrx.heatmap
        clusters = gmtrx.clusters
        get_energydist_histogram(gmtrx.heatmap, 2.0)
        mn, mx, rng = find_minmax(gmtrx.heatmap, N, True)
        print ("min(%12.2f), max(%12.2f), range(%12.2f)" % (mn, mx, rng))
    #
    
#


if __name__ == '__main__':
    main(sys.argv)
#
