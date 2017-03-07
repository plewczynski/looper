#!/usr/bin/env python

# This is an important tool to be used with the chromatin analysis
# programs. I have introduced at as a separate module because it is
# used within the suit of code.



import math
import sys
from FileTools import FileTools
from Chromosome import Chromosome
from Chromosome import Chromatin
from Chromosome import Segment
from Chromosome import allowed_tags
from Chromosome import assign_bed_tags


####################################################################
##########################  MATRIX TOOLS ###########################
####################################################################

class HeatMap:
    def __init__(self, version):
        self.version     = version
        self.resolution  = 5000
        self.position    =    0 # unset
        self.length      =    0
        self.heatmap     = []
        self.clusters    = []
        self.g_bgn         = 1e100
        self.g_end         = -1 
        self.set_HeatMap = False
    #
    
    def set_heatmap(self, hmap):
        # should check the inputs; e.g., is it an NxN matrix?
        self.length = len(hmap[0])
        self.heatmap = hmap
        self.set_HeatMap = True
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
    

# MatrixTools() is used to define the matrix points
class MatrixTools:
    def __init__(self):
        self.DEBUG  = False
        # special considerations
        self.from_Nenski = False
        
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
        
        
        # weights for selecting out PET clusters 
        self.set_CTCF_weights(100.0)
        self.pssbl_ctcf   = {}
        self.edge_ctcf    = {}
    #
    
    def set_Nenski(self):
        self.from_Nenski = True
    #
    
    def set_CTCF_weights(self, ref_scale):
        self.ctcf_tthresh = 0.20*ref_scale
        self.ctcf_cthresh = 0.40*ref_scale
    #
    
    def initialize_matrix(self, Var, N, w):
        # examples:
        
        # initialize_matrix(v, N, 0.0)
        # initialize_matrix(v, N, 1)
        
        # v = the matrix variable name
        
        # N = the size of the matrix (N x N).
        
        # w = the value to be used.
        
        # Because the python program doesn't care, the type can be
        # either float or integer. It can also be used to assign
        # vectors such as M-loops: cc = [(0,0)...]. Because this
        # simply acts as a place holder, any sort of object can be
        # inserted, as long as it is the same for all matrix elements.
        
        # However, this also mean that you have to be very careful
        # what you enter into the function, because python does not
        # care what you insert into this thing (integers, floats,
        # vectors, objects). If you don't pay attention to what goes
        # in, garbage is likely to follow ... python doesn't care. Of
        # course, this sort of thing could also be done with C++.
        
        
        Var = []
        for j in range(0, N):
            vj = []
            for i in range(0, N):
                vj += [ w ]
            #
            Var += [vj]
        return Var
    #
    
    
    
    # outputs a string containing the contents of a matrix of integer
    # variables (here it matters what you put in!!!).
    def disp_imatrix(self, v, name="matrix", flag_x_on_dgnl = True, flag_no_header = False):
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
            #
            s += sj + '\n'
        #
        return s
    #
    
    
    
    # outputs a string containing the contents of a matrix of integer
    # variables (here it matters what you put in!!!).
    def make_heatmap(self, v, flag_no_header = False):
        N = len(v)
        s = ''
        if not flag_no_header:
            s += "%d\n" % N
        for j in range(0,N):
            sj = ''
            if not len(v[j]) == N:
                print "ERROR(MatrixTools.make_heatmap): undiscernable matrix dimesions."
                print "     matrix size %d, for column %d: %d" % (N, j, len(v[j]))
                sys.exit(1)
            #
            
            for i in range(0,N):
                if i == j:
                    sj += "0\t"
                else:
                    sj += "%d\t" % v[j][i]
                #
            #
            s += sj + '\n'
        #
        return s
    #
    
    
    # Outputs a string containing the contents of a matrix of float
    # variables (here it matters what you put in!!!).
    def disp_fmatrix(self, mtrx, name="matrix", flag_x_on_dgnl = True, flag_no_header = False):
        N = len(mtrx)
        s = ''
        if not flag_no_header:
            s += name + "  %d\n" % N
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
            #
            s += '\n'
        #
        return s
    #
    
    
    # This simply reads the heatmap data (raw and unfiltered). In the
    # diagonal region, it deletes anything in the diagonal regions.
    def read_heatmap(self,
                     flnm,
                     EXTS = ["heat", "eheat"],
                     PROGRAM = "read_MatrixFile",
                     flag_display = True):
        
        debug_read_heatmap = flag_display
        flag_extended_heatmap = False
        
        splf = flnm.split('.')
        ext = splf[len(splf)-1]
        
        
        if debug_read_heatmap:
            print "file name: %s" % flnm
        #
        
        ft = FileTools()
        if not ft.check_ext(flnm, EXTS):
            print "ERROR: %s only accepts files with extention '%s'." % (PROGRAM, EXTS)
            print "       input file '%s' --> (extension(s) '%s')" % (flnm, ext)
            sys.exit(1)
        #
        
        try:
            fp = open(flnm, 'r')
        except IOError:
            print "ERROR: cannot open file '%s'." % flnm
            sys.exit(1)
        #
        lfp = fp.readline()
        fileInfoLine = lfp.strip().split(':')
        fp.close()
        
        
        # Determine what file type is being submitted.
        if debug_read_heatmap:        
            print "fileInfoLine: ", fileInfoLine
        #
        if len(fileInfoLine) > 1:
            # means it is probably the extended format, unless by
            # chance this is the first line after the editing step.
            if fileInfoLine[0] == "file version":
                flag_extended_heatmap = True
            else:
                print "ERROR: %s does not appear to be a recognizable file type" % flnm
                print "       first line: '%s'" % fileInfoLine
                sys.exit(1)
            #
        #
        
        gmtrx = None
        if flag_extended_heatmap:
            if debug_read_heatmap:
                print "going to read_extended_heatmap()"
            #
            gmtrx = self.read_extended_heatmap(flnm, PROGRAM, flag_display)
        else:
            if debug_read_heatmap:
                print "going to read_basic_heatmap()"
            #
            gmtrx = self.read_basic_heatmap(flnm, PROGRAM, flag_display)
        #
        N = gmtrx.length
        mtrx = gmtrx.heatmap
        clusters = gmtrx.clusters
        return gmtrx # mtrx, clusters, N
    #
    
    
    
    # This is the format that I've encountered from Przemek and
    # Michal, so I make this a special format different from my own
    # extended heatmap format. This simply reads the heatmap data (raw
    # and unfiltered). In the diagonal region, it deletes anything in
    # the diagonal regions.
    def read_basic_heatmap(self, flnm, PROGRAM = "read_heatmap", flag_display = False):
        splf = flnm.split('.')
        ext = splf[len(splf)-1]
        
        gmtrx = HeatMap("Przemek")

        if flag_display:
            print "file name: %s" % flnm
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
            # print s
            N = int(s[len(s)-1]) # read the last field from line 0
            if len(s) > 2:
                nj = len(lfp)
                ni = len(s)
                # print ni, nj
                if ni == nj:
                    print "heatmap: using Michal K format"
                    N = ni
                    start = 0
                else:
                    print "ERROR: unrecognized format for heatmap file %s" % flnm
                    sys.exit(1)
                #
            #
                    
            
        except ValueError:
            print "ERROR: first line of %s cannot be read properly. " % flnm
            print "       Format should be '\"some text\"    integer'; i.e., '%s %d'. "
            print "       Please check the file and confirm."
            sys.exit(1)
        #
        if flag_display: 
            print "size of matrix: %d" % N
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
    
    # This reads heatmap files with the extended information; with
    # extension "eheat". Therefore, not only does the it read the
    # file, it reads auxiliary information such has clusters,
    # resolution, etc..
    def read_extended_heatmap(self, flnm, PROGRAM = "read_heatmap", flag_display = False):
        splf = flnm.split('.')
        ext = splf[len(splf)-1]
        
        
        if flag_display:
            print "file name: %s" % flnm
        #
        
        # this program is called by read_heatmap, so there SHOULD NOT
        # be any need to do the usual checking.
        fp = open(flnm, 'r')
        lfp = fp.readlines()
        fp.close()
        
        s = lfp[0].strip().split(':')
        version =  s[1].strip()
        print "reading version: %s" % version
        gmtrx = HeatMap(version)
        start = 1
        N = -1
        if version == "0.0":
            try:
                # note
                # s = lfp[1].strip().split(':')[1]
                # ss = s.strip().split()[0])
                
                gmtrx.resolution = int(lfp[1].strip().split(':')[1].strip().split()[0])
                gmtrx.position   = int(lfp[2].strip().split(':')[1].strip().split()[0])
                gmtrx.length     = int(lfp[3].strip().split(':')[1].strip().split()[0])
                if flag_display:
                    print "resolution: %10d" % gmtrx.resolution
                    print "position:   %10d" % gmtrx.position
                    print "length:     %10d" % gmtrx.length
                #
            except ValueError:
                print "ERROR: first lines of %s are not formatted properly. " % flnm
                for k in range(0, 5):
                    print "line %2d:  %s" % (k, lfp[k].strip())
                #
                print "       Please check the file and confirm."
                sys.exit(1)
            #
            start = 5
            gmtrx.heatmap = self.read_matrix(start, gmtrx.length, lfp)
            print "loaded matrix"
            #print gmtrx.heatmap
            
            k = start + gmtrx.length
            if not lfp[k].strip() == "//":
                print "ERROR: %s format has some incompatibilities" % flnm
                sys.exit(1)
            #
            k += 1
            try:
                s = lfp[k].strip().split(':')[1]
                gmtrx.g_bgn = int(s.strip().split()[0])
                gmtrx.g_end = int(s.strip().split()[1])
            except ValueError:
                print "ERROR: %s, cluster line after heatmap data cannot be read properly." % flnm
                for k in range(5+gmtrx.length, len(lfp)):
                    print "line %2d:  %s" % (k, lfp[k].strip())
                #
                print "       Please check the file and confirm."
                sys.exit(1)
            #
            k += 1
            taglist = lfp[k].strip().split(':')[1].strip().split()
            print taglist
            k += 1
            while not lfp[k].strip() == "//":
                sbedfl = lfp[k].strip().split()
                chrmtn = assign_bed_tags(taglist, sbedfl)
                gmtrx.clusters += [chrmtn]
                print chrmtn.disp_data()
                k += 1
            #
        #
        
        
        
        
        # purge unnecessary memory
        # ########################
        del lfp
        del fp
        # ########################
        
        
        return gmtrx
    #
    
    # this does the actual reading of the heatmap matrix
    def read_matrix(self, start, N, lfp):
        debug_read_martix = False
        if debug_read_martix:
            print "entered read_matrix():"
        #
        
        mtrx = []
        ij_min = 100000.0
        for k in range(start,N+start):
            mtrxj = []
            s = lfp[k].strip().split()
            if not len(s) == N:
                # probably should make sure that the dimensions are
                # exactly equal, but anyway, they should at least be
                # such that you take some part of the matrix!
                print "ERROR: something wrong with the dimensions of the matrix!"
                sys.exit(1)
            #
            j = k - start
            for i in range(0,N):
                try:
                    w = float(s[i])
                    
                    if not (w >= 0.0 or ext == "clust"):
                        # just so stupid stuff doesn't get in.
                        print "ERROR: Encountered a negative value in the chromatin "
                        print "       input data. "
                        print "       This cannot be heat map data for chromatin!"
                        print "       => position (i,j) = (%d,%d), value = " % (i, j), w
                        print "       -- all data in the heat map must be positive"
                        print "       Please verify or correct the file '%s'." % flnm
                        sys.exit(1)
                    #
                    if w < ij_min:
                        ij_min = w
                    #
                    
                    mtrxj += [w]
                    #
                except ValueError:
                    if s[i] == 'x':
                        mtrxj += [0.0]
                    else:
                        print "ERROR: matrix contains unrecognized terms"
                        print "       from ", s
                        print "       cannot recognize '%s'" % s[i]
                        sys.exit(1)
                    #
                #
            #
            
            mtrx += [mtrxj]
        #
        
        # this is a patch to fix some older data
        if ij_min < 0.0 and ext == "clust":
            shift = - ij_min
            print "unshifting free energy data by ", shift
            for j in range(0, N):
                for i in range(0, N):
                    if not mtrx[i][j] == 0.0:
                        mtrx[i][j] += shift
                        mtrx[j][i] += shift
                    #
                #
            #
        #
        if debug_read_martix:
            print "finished read_matrix()"
        #
        return mtrx
    #
    
    
    def disp_matrix(self, mtrx, N):
        s = ''
        for j in range(0, N):
            s_j = ''
            for i in range(0, N):
                s_j += "%8.4g\t" % mtrx[i][j]
            s += s_j + '\n'
        return s
    #
    
    
    
    
    # This simply reads the heatmap data (raw and unfiltered). In the
    # diagonal region, it deletes anything in the diagonal regions.
    def write_BscHeatMap(self, flnm, mtrx, N, fformat = "Michal_K"):
        # presently only the most simplest format of Michal Kadlof
        splf = flnm.split('.')
        ext = splf[len(splf)-1]
        if self.DEBUG:
            print "output file name: %s" % flnm
        #
        
        try:
            fp = open(flnm, 'w')
        except IOError:
            print "ERROR: cannot open file '%s'." % flnm
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
        splf = flnm.split('.')
        ext = splf[len(splf)-1]
        if self.DEBUG:
            print "output file name: %s" % flnm
        #
        
        try:
            fp = open(flnm, 'w')
        except IOError:
            print "ERROR: cannot open file '%s'." % flnm
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
            s += domain.showDomain()
            s += "//\n"
            fp.write(s)
            fp.close()
        else:
            print "Sorry, no format for version '%s' is available presently" % version
        #
        
    #
    
    
    
    
    # #############################################################
    # #############################################################
    
    # 161005wkd: I know this is a bit silly to do every operation
    # separately; therefore force multiple N^2 steps rather than do
    # these steps all at once.
    
    # Unfortunately, we just have so many different types of data,
    # that I found it annoying to have to constantly be writing in and
    # checking some new routine for each type of new data. Sometimes,
    # the placement of these things (particularly the histograms) need
    # to be moved around when looking at some different issue.
    # Moreover, make_heatmap needs some of the followign functions in
    # this current arrangement, but other functions are not used. So
    # it seems better to ask the program to obtain the relevant
    # material, not do all sorts of additional strange things that we
    # don't want. The modularity will allow me to use what I actually
    # need.
    
    # Therefore, I finally decided to modularize this section so that
    # these functions can be called in a flexible way. Presently, they
    # are still located in this routine, but I think that eventually I
    # will locate them as needed in the programs that call this
    # function.
    
    # #############################################################
    # #############################################################
    
    
    # This reads things like the heatmap data, which always involves
    # elements where i != j (no diagonal elements).
    def read_MatrixFile(self, flnm, EXTS = ["heat", "eheat"], rescale_wt = 1.0, PROGRAM = "read_MatrixFile", flag_display = True, bgn_shift = 0):
        #
        #gmtrx = HeatMap(version)
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
                
        # full rescale relative to a self.rescale_wt
        self.rescale_wt = rescale_wt
        gmtrx.heatmap = self.rescale_mtrx(gmtrx.heatmap, N, self.rescale_wt)
        #
        # obtain the final minimum and maximum weights
        self.hm_min, self.hm_max, self.wt_range = self.find_minmax(gmtrx.heatmap, N, True)
        # print self.hm_min, self.hm_max, self.wt_range
        
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
            print "used read in PET cluster information: %s" % flnm
            self.set_CTCF_points(gmtrx, True) # flag_display
        else:
            print "estimate the PET cluster information: %s" % flnm
            self.estimate_CTCF_points(gmtrx.heatmap, N, flag_display, bgn_shift)
        #
        
        return gmtrx
    #
    
    
    # This reads things like the heatmap data, which always involves
    # elements where i != j (no diagonal elements).
    def read_MatrixFile_wt(self, flnm, EXTS = ["heat", "eheat"], PROGRAM = "read_MatrixFile"):
        #
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
        
        print "hm_(min,max): ", self.hm_min, self.hm_max
        
        return gmtrx
    #
    
    # find the CTCF points in the data
    def estimate_CTCF_points(self, mtrx, N, flag_display = True, bgn_shift = 0):
        irh = N; jrh = 0
        rhlist = {}
        mx = 0.0
        i_mx = 0
        j_mx = 0
        for j in range(0,N):
            for i in range(0,N):
                w = mtrx[i][j]
                #
                if w >= self.ctcf_tthresh:
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
                if w > mx:
                    if i < j:
                        mx = w
                        i_mx = i
                        j_mx = j
                    #
                #
            #
        #
        if flag_display:
            print "maximum matrix element at (%d,%d)[value= %8.3f]" % (i_mx, j_mx, mx)
        #
        
        # Filter out the domain boundary CTCF since that will
        # automatically be used for CTCF island formation. Moreover,
        # if it is the only one, there are no other potential CTCF
        # interactions from which to evaluate or build any islands.
        self.pssbl_ctcf = {}
        self.edge_ctcf = {}
        if len(rhlist) > 0:
            for rh in rhlist.keys():
                if not (rh[0] == irh and rh[1] == jrh):
                    # it is important to treat the border different
                    self.pssbl_ctcf.update({(rh[0],rh[1]) : rhlist[rh]})
                else:
                    self.edge_ctcf.update({(rh[0],rh[1]) : rhlist[rh]})
                #
            #
        #
        print self.print_ctcf_dict("edge_ctcf ", self.edge_ctcf, bgn_shift)
        print self.print_ctcf_dict("pssbl_ctcf", self.pssbl_ctcf, bgn_shift) 
        # sys.exit(0)
        return i_mx, j_mx, mx
    #
    
    def print_ctcf_dict(self, title, cdict, bgn_shift = 0):
        s = "%12s  {" % title
        clist = cdict.keys()
        if len(clist) > 1:
            for k in range(0, len(clist)-1):
                ck = clist[k]
                # print "ck: ", ck, ck[0]
                s += "(%d, %d) : %d, " % (ck[0]+bgn_shift, ck[1]+bgn_shift, cdict[ck])
            ck = clist[len(clist)-1]
            s += "(%d, %d) : %d}" % (ck[0]+bgn_shift, ck[1]+bgn_shift, cdict[ck])
        elif len(clist) == 1:
            ck = clist[0]
            s += "(%d, %d) : %d}" % (ck[0]+bgn_shift, ck[1]+bgn_shift, cdict[ck])
        else:
            s += " }"
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
                vi = int(clusters[k].vi); vj = int(clusters[k].vj); nPET = clusters[k].nPET 
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
            #
            
            if flag_display:
                print "maximum matrix element at  (%5d,%5d)[weight=   %8.1f]" % (i_mx, j_mx, mx)
                print "total CTCF domain coverage (%5d,%5d)[length= %8d]" % (irh, jrh, (jrh-irh+1))
            #
            
            # sometimes the boundary  element is missing
            if not rhlist.has_key((irh,jrh)):
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
            #
        #
        
        print "edge_ctcf:  ", self.edge_ctcf
        print "pssbl_ctcf: ", self.pssbl_ctcf
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
                    if histogram.has_key((j-i)):
                        histogram[j-i] += [w]
                    #
                    else:
                        histogram.update( { (j-i) : [w] })
                    #
                #
            #
        #
        fp = open(flnm, 'w')
        for k in histogram.keys():
            wf = 0.0
            for w in histogram[k]:
                wf += w
            wf /= float(len(histogram[k]))
            fp.write("%4d   %12.2f\n" % (k, wf))
        #
        fp.close()
        #
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
                    if histogram.has_key(jmi):
                        if histogram[jmi][0] > w:
                            histogram[jmi][0] = w
                        #
                        if histogram[jmi][1] < w:
                            histogram[jmi][1] = w
                        #
                    #
                    else:
                        histogram.update( { (jmi) : [w,w] })
                    #
                #
            #
        #
        fp = open(flnm, 'w')
        for k in histogram.keys():
            w_min = histogram[k][0]
            w_max = histogram[k][1]
            fp.write("%4d   %12.2f   %12.2f\n" % (k, w_min, w_max))
        #
        fp.close()
        #
        return histogram, w_min, w_max
    #
    
    
    # This can be used anywhere after the matrix is built to scans
    # through the matrix elements and find the minimum and maximum
    # value.
    def find_minmax(self, mtrx, N, flag_linear):
        if not N == len(mtrx):
            print "ERROR: matrix length not matching with proposed size"
            print "matrix length(%d)  vs N(%d)" % (len(mtrx), N)
            sys.exit(1)
        #
        # now settle on the data type
        
        # Here, hm_max, hm_min, wt_range and mtrx are treated as local
        # variables because there is the possibility that one does not
        # wish to destroy the raw input data from the heatmap.
        
        hm_max = -1000.0
        hm_min = 1000.0
        for j in range(0,N):
            for i in range(0,N):
                w = mtrx[i][j]
                #
                if w > hm_max:
                    hm_max = w
                if w < hm_min:
                    hm_min = w
                #
            #
        #
        if flag_linear:
            # In THESE problems, the linear data is in the form of
            # ChIA-PET or Hi-C data and is ALWAYS expressed from 0 to
            # hm_max
            wt_range = hm_max
        else:
            # Presently, I use this when I rescale the linear data to
            # logarithmic data.
            wt_range = hm_max - hm_min
        return hm_min, hm_max, wt_range
    #
    
    # Scaling Hi-C data using f(|j-i|) = 1 - exp[-b(|j-i|-1)].  This
    # is used to remove the homopolymer noise from the Hi-C data to
    # make it possible to process the data using chreval and related
    # programs. This should only be used with LINEAR DATA.  It should
    # NOT be used with ChIA-PET data.
    def filter_HiC_using_1_m_exp(self, mtrx, N):
        # now settle on the data type
        
        # Here, hm_max, hm_min, wt_range and mtrx are treated as local
        # variables because there is the possibility that one does not
        # wish to destroy the raw input data from the heatmap.
        hm_max = -1000.0
        hm_min = 1000.0
        for j in range(0,N):
            for i in range(0,N):
                w = mtrx[i][j]*(1.0 - math.exp(-0.0035*(float(abs(j-i))-1.0)))
                w = int(w + 0.5)
                mtrx[i][j] = w
                #
                if w > hm_max:
                    hm_max = w
                if w < hm_min:
                    hm_min = w
                #
            #
        #
        wt_range = hm_max
        # this is LINEAR data; i.e., it ranges between 0 to hm_max
        return mtrx, hm_min, hm_max, wt_range
    #
    
    # Scaling Hi-C data using 1/{Nexp[-b(|j-i|-1)] + c}.  This is used
    # to remove the homopolymer noise from the Hi-C data to make it
    # possible to process the data using chreval and related
    # programs. This should only be used with LINEAR DATA. It should
    # NOT be used with ChIA-PET data.
    def filter_HiC_using_exp(self, mtrx, N):
        # now settle on the data type
        
        # Here, hm_max, hm_min, wt_range and mtrx are treated as local
        # variables because there is the possibility that one does not
        # wish to destroy the raw input data from the heatmap.
        hm_max = -1000.0
        hm_min = 1000.0
        for j in range(0,N):
            for i in range(0,N):
                w = mtrx[i][j]
                if abs(j-i) == 1:
                    # this position is rather special and particularly
                    # noisy. It is just the nearest neighbor
                    # interactions, so it is meaningless as well.
                    w = 0.0
                else:
                    w = w/(4451*math.exp(-0.374*float(abs(j-i))-1.0) + 55.1 )
                    w = int(w + 0.5)
                #
                mtrx[i][j] = w
                #
                if w > hm_max:
                    hm_max = w
                if w < hm_min:
                    hm_min = w
                #
            #
        #
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
            #
        #
        return mtrx
    #
    
    
    def reweight_ln(self, mtrx, N):
        # I think this could be done by hm_max = math.log(hm_max) and
        # hm_min = math.log(hm_max). Anyway, it also has to be
        # evaluated using find_minmax(), so either way, this step must
        # be done somewhere.
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
                        v = math.log(w)
                    #
                if v > hm_max:
                    hm_max = v
                if v < hm_min:
                    hm_min = v
                #
                
            #
        #
        wt_range = hm_max - hm_min
        if hm_max == hm_min:
            # this can happen if all the values in the matrix are 1,
            # for example.
            wt_range = 1.0
            print "reseting wt_range = %.2f" % wt_range
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
                        w = math.log(w) / wt_range
                    mtrx[i][j] = w
                #
            #
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
                        w = math.log(w) / wt_range
                    mtrx[i][j] = w
                #
            #
        #
        return mtrx, hm_min, hm_max, wt_range
    #
    
    
    def add_boundaryWt(self, N, mtrx, boundaryWt):
        i_min = 0; j_max = N-1
        #
        if mtrx[i_min][j_max] == 0.0:
            if self.DEBUG:
                print "adjusting borders (%d,%d)" % (i_min, j_max)
            #
            j_max = 0
            i_min = N
            for j in range(N/2, N):
                for i in range(N/2, -1, -1):
                    if mtrx[i][j] > 0.0:
                        if self.DEBUG:
                            print "ij = ", i, j
                        #
                        if i <= i_min and j >= j_max:
                            i_min = i
                            j_max = j
                        #
                    #
                #
            #
        #
        else:
            print "no need to adjust the location of the border"
        
        mtrx[i_min][j_max] = boundaryWt
        mtrx[j_max][i_min] = boundaryWt
        print "Note: Additional weight added to PET cluster, borders at (%d,%d)" % (i_min, j_max)
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
            print "ERROR: option 'optimum' requires either a 'pos' or 'neg' (default)"
            print "       Don't know what you want!"
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
                #
                else:
                    if not (mtrx[i][j] == 0.0):
                        # I have to think about this, but I think we
                        # should not allow digonal elements at least
                        # in this construction.  Maybe it needs to be
                        # renamed, but anyway, we the purpose of this
                        # tool is to analyze matrices and these
                        # matrices should all be ones that do not have
                        # diagonal elements.
                        print "ERROR: diagonal element (%d,%d) is not zero!" % (i,j)
                        sys.exit(1)
                    #
                #
            #
        #
        wt = m_max - m_min
        
        if optimum == "neg":
            wt = - wt
        #
        
        for j in range(0,N):
            for i in range(0,N):
                mtrx[i][j] /= wt
            #
        #
        return mtrx
    #

# end class MatrixTools

def usage():
    print "USAGE: MatrixTools.py -option args"
    print "        option     arguments"
    print "       -basic_hm   file.clust"
#

def convert_to_basic_heatmap(flnm):
    """removes any artifacts that are in some heatmaps that were generated for visualization and produces the most crude form -- the heatmap and nothing else."""
    mtools = MatrixTools()
    gmtrx = mtools.read_heatmap(flnm, ["heat","eheat","clust"], "MatrixTools")
    N        = gmtrx.length
    mtrx     = gmtrx.heatmap
    clusters = gmtrx.clusters

    # construct a new file name from the old one
    splf = flnm.split('.')
    ext = splf[len(splf)-1]
    flhd = flnm[0:len(flnm)-len(ext) - 1]
    flnm_rev = flhd + "_new." + ext
    print "new file name: ", flnm_rev
    
    mtools.write_BscHeatMap(flnm_rev, mtrx, N, "Michal_K")
#


def main(cl):
    n = len(cl)
    if  n == 1:
        usage()
        sys.exit(0)
    #
    flnm = ''
    option = ''
    k = 1
    while k < len(cl):
        arg = cl[k]
        if arg == '-basic_hm':
            option = "convert_to_basic_hm"
            k += 1
            if k < n:
                flnm = cl[k]
                print flnm
            #
        else:
            print "unrecognized argument at %d" % k
            print cl
            sys.exit(1)
        #
        k += 1
    #
    if option == "convert_to_basic_hm":
        convert_to_basic_heatmap(flnm)
    #
#


if __name__ == '__main__':
    main(sys.argv)

