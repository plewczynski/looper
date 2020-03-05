#!/usr/bin/env python

"""@@@

Main Program:  anal_loops.py (ANALyze LOOPS)

Classes:       Anal

Author:        Wayne Dawson
creation date: mostly 2016 and up to March 2017.
last update:   180709
version:       0

anal_loops.py (ANALyze LOOPS)
    
    Anal is an analysis tool that runs batch calculations with
    Chreval.
    
    The input file is a list of loops formatted according to Przemek's
    approach (with extension 'bed') Teresa's approach (currently with
    extension 'txt').
    
    command line example:
    > anal_loops.py -ff loops.CTCF.annotated.bed -dG
    
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


    In general, bed type files should only be one instance.

    An alternative file input would involve the following command line

    command line example:
    > anal_loops.py -ff wayn_active_inactive.txt wayn_compartments.txt wayn_open.txt -dG

    where the various list of "txt" files (wayn_active_inactive.txt,
    etc.) are formatted according to Teresa.


 used as a identifier                 |used in the program ----------------------------------
 region of chromosome                 |specific location in region             state   length
 chr   bgn             end            |chr     bgn             end              
 chr1	10436116	10556545	chr1	10437400	10438000	active	600
 chr1	10436116	10556545	chr1	10438000	10438466	active	466
 chr1	10436116	10556545	chr1	10438558	10440200	active	1642
        

    These files often have different data related to open, active,
    inactive, compartment A or B, etc., so the results typically
    require more than one file.  Therefore, multiple file entries are
    allowed for both file types; however, in general, for "bed" files,
    there should only be one input file.

    The chromatin loop is specified by the region chr1 10436116 to
    10556545 and the zones listed above as "active" are bands between
    (10437400-10438000), (10438000-10438466) and
    (10438558-10440200). The percentage of active region is the sum of
    the lengths of the small segments divide by the length of the
    loop.

    There are several ways to weight the data.

    1. (default) The simplest is to just take the first is simply the
    first Boltzmann probability term (no argument). This assumes that
    a compact structure will have a very large Boltzmann probability
    compared to any other possible structure.

    basic comparision of stability
        
    Measure stability of a structure based on the magnitude of the
    first Boltzmann probability term. Presumably, the terms that have
    a very high Boltzmann probability also have the highest
    stability. However, it should be remembered that as the size of
    the loop becomes large, it is more likely that you will have
    multiple structures that are very similar in structure and free
    energy. As a result, I started searching for ways to combine
    similar structures together using various rules such as hamming
    distance, word similarity, variation in delta(TdS) and variation
    in delta(dG). Hence, even the first structure is shared between
    two or three or even a dozen different structures in the ensemble.
    
    Therefore, I have (presently) set up a variety of options for how
    to analyze the similarity of the structures in the ensemble.

    2. (option: -sim) This uses calculation of a similarity ratio of
    two sequences. It is only the pattern that is looked at.

    3. (-ham) Therefore, a second approach is to discriminate
    according to the Hamming distance of the structures.

    4. (option -TdS) Hamming distance has some problems of its
    own. For example, when you observe sliding of a chain, the Hamming
    distance simply measures the observable edits.

    ...(((....)))...
    ....(((...)))...
    Example 1

    ...(((....)))...
    ....(((....)))..
    Example 2


    Consider Example 1. From viewpoint of Hamming distance, the number
    of edits in Example 1 are two. However, in fact, three contacts
    have actually moved. Hamming distance may efficiently measure
    edits, but it does not measure change. Entropy (CLE) or free
    energy are probably better measure of this change.

    In example 2, the Hamming distance is 4, yet actually six
    position-changes have occurred. On the other hand, in terms of its
    cross linked entropy (CLE), there is no energy change at all
    because (given all interactions are identical) the number of
    contacts and their arrangement remain the same, only the postion
    of the contacts has shifted through a process of sliding.

    The cross linking entropy would observe some change in Example 1,
    but no change in Example 2.  In terms of energy change, it is not
    all that much and involves minor movements of the chain. It was my
    thought that simple sliding like this should not be penalized all
    that much in an ensemble of structures. It also seems like the
    Hamming distance is not particular accurate because it really is
    sliding of the chain of 3 residues in Example 1 or even 6 residues
    in Example 2. Therefore, I also introduced an entropic weight for
    measuring this difference.

    5. (option -dG) Finally, the free energy could greatly differ in
    Example 1 and Example 2, even though the CLE does no change at
    all, as in Example 2.  So I reasoned that the weight on the
    structures should be a function of the free energy (which includes
    the CLE contribution). So even though the ClE is unchanged in
    Example 2, the enthalpy could be quite different, and this would
    discriminate the flow of structures.  So how many equivalent
    structures (to within some set free energy) are observed is
    crucial.

"""

import sys
import os
import chreval
from GetOpts import GetOpts

from ChromatinData import Data
from ChromatinData import make_file_heading

from Functions import hamming_bin
from Functions import hamming_str
from Functions import similar

#
PROGRAM = "anal_loops.py"

def usage():
    print "USAGE: %s -ff file.bed file.bed ... " % PROGRAM
#



class Anal:
    def __init__(self, data):
        self.debug              = False
        self.data               = data
        self.oflnm              = ''
        self.iflnms             = ''
        self.cdata              = data.cdata
        self.keylist            = data.keylist
        self.datatype           = data.datatype
        self.res                = "res5kb"
        self.range_TddS         = 2.0
        self.range_ddG          = 2.0
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
    
    def disp_header(self, cl):
        n = len(self.data.flhd)
        if n > 1:
            header = '# files:\n'
        else:
            header = '# file:\n'
        for fk in range(0, len(self.data.flhd)):
            header += "#   %s.%s\n" % (self.data.flhd[fk], self.data.ext[fk])
        #
        header += '#\n'
        header += "# interactions:\n"
        for ntrctn in self.data.interaction:
            header += "#   %s\n" % ntrctn
        #

        header += "# comments:\n"
        if len(self.data.description) > 0:
            for dscrbk in self.data.description:
                if len(dscrbk) > 0:
                    for s in dscrbk:
                        header += "#   %s\n" % s
                    #
                #
                if n > 1:
                    header += "#   ----\n"
                #
            #
        #
        header += "# ----\n"
        
        # thermodynamic parameters
        # entropy parameters from the chreval
        header += "# thermodynamic parameters:\n"
        header += "#   entropy:\n"
        header += "#     xi       = %.3g\n" % cl.xi
        header += "#     seg_len  = %.3g\n" % cl.seg_len # [bps]
        header += "#     lmbd     = %.3g\n" % cl.lmbd    # [nt]
        header += "#     gmm      = %.3g\n" % cl.gmm     # [dimensionless]
        header += "#     T        = %.3g\n" % cl.T          # in Kelvin
        # constants: weights for the enthalpy terms 
        header += "#   enthalpy:\n"
        header += "#     febase   = %.3g\n" % cl.dHbase
        header += "#     feshift  = %.3g\n" % cl.dHshift
        header += "# ----\n"
        
        return header
        
        
    
    def anal_loops(self, cl):
        print "anal_loops()"
        keylist = self.keylist
        if len(self.keylist) == 0:
            print "setting keylist"
            print "length of the dataset: ", len(self.cdata)
            keylist = self.data.ordered_keylist()
        # print keylist
        # print self.cdata.keys()
        # print self.data.chrmsm_grp.keys()
        # sys.exit(0)
        self.oflnm              = cl.f_output
        self.iflnms             = cl.f_activity
        
        
        full_cmplxtydict = { 'cmplx1' : False,
                             'cmplx2' : False,
                             'cmplx3' : False }
        full_cmplxtylist = ['cmplx1', 'cmplx2', 'cmplx3']
        cmplxty = []
        
        # The size of "self.data.taglist" should be at least one --
        # from at least one input file. Otherwise, the program would
        # never get to this point!
        for k in range(0, len(self.data.taglist)):
            for tlk in self.data.taglist[k]:
                if full_cmplxtydict.has_key(tlk): 
                    full_cmplxtydict[tlk] = True
                #
            #
        #
        for fc in full_cmplxtylist:
            if full_cmplxtydict[fc]:
                cmplxty += [fc]
            #
        #
        
        
        full_epidict = { 'open'      : False,
                         'active'    : False,
                         'A'         : False,
                         'B'         : False,
                         'repressed' : False }
        
        full_epilist = ['open', 'active', 'A', 'B', 'repressed']
        epilist = []
        # self.data.taglist is the list that is used
        for k in range(0, len(self.data.taglist)):
            for tlk in self.data.taglist[k]:
                if full_epidict.has_key(tlk): 
                    full_epidict[tlk] = True
                #
            #
        #
        #
        for el in full_epilist:
            if full_epidict[el]:
                epilist += [el]
            #
        #
        
        p_active    = 0.0
        p_open      = 0.0
        p_repressed = 0.0
        p_A         = 0.0
        p_B         = 0.0
        missing     = []
        rPET_wt     = []
        wt_PETo     = cl.PETwt
        
        
        flag_basic        = cl.basic       # the first Boltzmann prob
        flag_similar      = cl.similarity  # prob wtd by similarity
        flag_hamming_bin  = False          # not used!!!!
        flag_hamming      = cl.hamming     # prob wtd by hamming distance
        flag_TdS          = cl.TdS         # prob wtd by CLE 
        self.range_TddS   = cl.TdS_range   # sets either default or set value
        flag_ddG          = cl.ddG         # prob wtd by d(dG)
        self.range_ddG    = cl.ddG_range   # sets either default or set value
        flag_allwts       = cl.allwts      # show all weights (default)
        
        
        header = self.disp_header(cl)
        
        
        # begin by setting up save file with header information
        rflnm = self.oflnm  # the output file name
        
        dlabel = "# chr     begin        end      ctcf1   ctcf2     nPET      "
        for cc in cmplxty:
            dlabel += "%s     " % cc
        #
        plabel = "len    dG_0      TdS_0     dGbar     TdSbar    p_max     "
        if flag_allwts:
            print "all:"
            plabel += "p_sim    cmplx  p_ham    cmplx  p_dTdS   cmplx  p_ddG    cmplx  "
        elif flag_similar:
            print "similar:"
            plabel = "len    dG_0      TdS_0     p_sim    cmplx  "
        elif flag_hamming:  # based on string edits
            print "Hamming:"
            plabel = "len    dG_0      TdS_0     p_ham    cmplx  "
        elif flag_TdS:
            print "dTdS:"
            plabel = "len    dG_0      TdS_0     p_dTdS   cmplx  "
        elif flag_ddG:
            print "ddG:"
            plabel = "len    dG_0      TdS_0     p_ddG    cmplx  "
        else:
            print "basic:"
        #
        fp = open(rflnm, 'w')
        edatfld = ''
        if self.datatype == 'type1':
            edatfld += "active    open   repressed      A        B"
        else:
            for kxs in epilist:
                edatfld += "%8s  " % kxs
            #
        #
        fp.write("%s" % header)
        fp.write("%s%s%s\n" % (dlabel, plabel, edatfld))
        fp.close()
        
        
        for ky in keylist:
            # this now shows that I can extract whatever datapoint I
            # desire from the set of experimental data. Therefore, I
            # can read in all the data and extract it according to my
            # needs.
            
            name = self.cdata[ky].name
            bgn  = self.cdata[ky].bgn
            end  = self.cdata[ky].end
                        
            flhd = make_file_heading(name, bgn, end, self.res)
            
            # #######################################################
            # verify that the files exist or not. If they are missing
            # then skip any further analysis.
            # #######################################################
            
            # print "evaluating: ", os.path.exists(cl.f_heatmap[0]), cl.f_heatmap[0]
            if not os.path.exists(flhd + ".heat"):
                missing += [ flhd ]
                continue
            #
            
            
            
            # ##################################################
            # ########   compute the FE of the heatmap  ########
            # ########      and epigenetic factors      ########
            # ##################################################
            
            cl.PETwt = wt_PETo
            if self.datatype == 'type1': 
                p_active, p_open, p_repressed, p_A, p_B = self.data.get_type1_prob(ky)
            elif self.datatype == 'type2':
                p_active, p_open, p_repressed, p_A, p_B = self.data.get_type2_prob(ky)
                cl.PETwt = self.get_PET_wt(ky, wt_PETo)
            #
            self.cdata[ky].state["active"]    = p_active
            self.cdata[ky].state["open"]      = p_open
            self.cdata[ky].state["repressed"] = p_repressed
            self.cdata[ky].state["A"]         = p_A
            self.cdata[ky].state["B"]         = p_B
            rPET_wt += [ [ cl.PETwt, flhd ] ] 
            
            # configuration settings
            cl.f_heatmap     = [flhd + ".heat"]
            cl.allowed_extns = cl.EXTS["chreval.py"]
            # have to override the extension settings when calling
            # programs such as chreval.py that involve heatmaps.
            
            
            print "\n\nfile name: %s.heat" % flhd
            manager = chreval.Manager()
            manager.runCalculations(cl)
            # manager.printResults()
            
            length = manager.N
            dt     = chreval.DispThreads(manager.calc.N)
            pr     = manager.trace.lt[0].p
            pr_h   = pr
            ham_include   = 0
            pr_s   = pr
            sim_include   = 0
            pr_TdS = pr
            dTdSp_include = 0
            pr_ddG = pr 
            ddGp_include  = 0
            TdS_0  = manager.trace.lt[0].TdS
            dG_0   = manager.trace.lt[0].dG
            s_0    = dt.makeLThreadDotBracket_1b(manager.trace.lt[0], 0)
            # True option produces only the structure string
            
            if flag_allwts:
                print "start-->[%4d]: %s   %8.3g" \
                    % (0, s_0,              pr)
                sim_include += 1
                ham_include += 1
                dTdSp_include += 1
                ddGp_include += 1
            elif flag_hamming_bin:
                print "start-->[%4d]: %s   %8.4f   %3d   %8.3g" \
                    % (0, s_0, 1.0, 0,      pr)
                ham_include += 1
            elif flag_hamming:  # based on string edits
                print "start-->[%4d]: %s   %8.4f   %3d   %8.3g" \
                    % (0, s_0, 1.0, 0,      pr)
                ham_include += 1
            elif flag_TdS:
                print "start-->[%4d]: %s   %8.4f   %3d   %8.3f   %8.3g" \
                    % (0, s_0, 1.0, 0, 0.0, pr)
                dTdSp_include += 1
            elif flag_ddG:
                print "start-->[%4d]: %s   %8.4f   %3d   %8.3f   %8.3g" \
                    % (0, s_0, 1.0, 0, 0.0, pr)
                ddGp_include += 1
            elif flag_similar:
                print "start-->[%4d]: %s   %8.4f   %8.3g" \
                    % (0, s_0, 1.0,         pr)
                sim_include += 1
            elif flag_basic:
                print "start-->[%4d]: %s   %8.3g" \
                    % (0, s_0,              pr)
            #
            
            for cnt in range(1, len(manager.trace.lt)):
                s_cnt = dt.makeLThreadDotBracket_1b(manager.trace.lt[cnt], True)
                p_cnt = manager.trace.lt[cnt].p
                if flag_allwts:
                    if self.debug:
                        print "        [%4d]: %s   %8.3g" \
                            % (cnt, s_cnt, p_cnt)
                    #
                    if similar(s_0, s_cnt) > 0.95:
                        if self.debug:
                            print "included[%4d]: %s   %8.4f   %8.3g" \
                                % (cnt, s_cnt, similar(s_0,s_cnt), p_cnt)
                        #
                        pr_s += manager.trace.lt[cnt].p
                        sim_include += 1
                    #
                    if hamming_str(s_0, s_cnt) < 5:
                        if self.debug:
                            print "included[%4d]: %s   %8.4f   %3d   %8.3g" \
                                % (cnt, s_cnt, similar(s_0,s_cnt), hamming_str(s_0,s_cnt), p_cnt)
                        #
                        pr_h += manager.trace.lt[cnt].p
                        ham_include += 1
                    #
                    TdS_cnt = manager.trace.lt[cnt].TdS
                    delta_TdS = TdS_0 - TdS_cnt
                    if -self.range_TddS < delta_TdS and delta_TdS < self.range_TddS:
                        if self.debug:
                            print "included[%4d]: %s   %8.4f   %3d   %8.3f   %8.3g" \
                                % (cnt, s_cnt, similar(s_0,s_cnt), hamming_str(s_0,s_cnt), delta_TdS, p_cnt)
                        #
                        pr_TdS += manager.trace.lt[cnt].p
                        dTdSp_include += 1
                    #
                    dG_cnt = manager.trace.lt[cnt].dG
                    delta_dG = dG_cnt - dG_0
                    if delta_dG < self.range_ddG:
                        if self.debug:
                            print "included[%4d]: %s   %8.4f   %3d   %8.3f   %8.3g" \
                                % (cnt, s_cnt, similar(s_0,s_cnt), hamming_str(s_0,s_cnt), delta_dG, p_cnt)
                        #
                        pr_ddG += manager.trace.lt[cnt].p
                        ddGp_include += 1
                    #
                elif flag_hamming_bin:
                    if hamming_bin(s_0, s_cnt) < 5:
                        if self.debug:
                            print "included[%4d]: %s   %8.4f   %3d   %8.3g" \
                                % (cnt, s_cnt, similar(s_0,s_cnt), hamming_bin(s_0,s_cnt), p_cnt)
                        #
                        pr_h += manager.trace.lt[cnt].p
                        ham_include += 1
                    else:
                        if self.debug:
                            print "        [%4d]: %s   %8.4f   %3d   %8.3g" \
                                % (cnt, s_cnt, similar(s_0,s_cnt), hamming_bin(s_0,s_cnt), p_cnt)
                        #
                    #
                elif flag_hamming: # based on string edits
                    if hamming_str(s_0, s_cnt) < 5:
                        if self.debug:
                            print "included[%4d]: %s   %8.4f   %3d   %8.3g" \
                                % (cnt, s_cnt, similar(s_0,s_cnt), hamming_str(s_0,s_cnt), p_cnt)
                        #
                        pr_h += manager.trace.lt[cnt].p
                        ham_include += 1
                    else:
                        if self.debug:
                            print "        [%4d]: %s   %8.4f   %3d   %8.3g" \
                                % (cnt, s_cnt, similar(s_0,s_cnt), hamming_str(s_0,s_cnt), p_cnt)
                        #
                    #
                elif flag_TdS:
                    TdS_cnt = manager.trace.lt[cnt].TdS
                    delta_TdS = TdS_0 - TdS_cnt
                    if -self.range_TddS < delta_TdS and delta_TdS < self.range_TddS:
                        if self.debug:
                            print "included[%4d]: %s   %8.4f   %3d   %8.3f   %8.3g" \
                                % (cnt, s_cnt, similar(s_0,s_cnt), hamming_str(s_0,s_cnt), delta_TdS, p_cnt)
                        #
                        pr_TdS += manager.trace.lt[cnt].p
                        dTdSp_include += 1
                    else:
                        if self.debug:
                            print "        [%4d]: %s   %8.4f   %3d   %8.3f   %8.3g" \
                                % (cnt, s_cnt, similar(s_0,s_cnt), hamming_str(s_0,s_cnt), delta_TdS, p_cnt)
                        #
                    #
                elif flag_ddG:
                    dG_cnt = manager.trace.lt[cnt].dG
                    delta_dG = dG_cnt - dG_0
                    if delta_dG < self.range_ddG:
                        if self.debug:
                            print "included[%4d]: %s   %8.4f   %3d   %8.3f   %8.3g" \
                                % (cnt, s_cnt, similar(s_0,s_cnt), hamming_str(s_0,s_cnt), delta_dG, p_cnt)
                        #
                        pr_ddG += manager.trace.lt[cnt].p
                        ddGp_include += 1
                    else:
                        if self.debug:
                            print "        [%4d]: %s   %8.4f   %3d   %8.3f   %8.3g" \
                                % (cnt, s_cnt, similar(s_0,s_cnt), hamming_str(s_0,s_cnt), delta_dG, p_cnt)
                        #
                    #
                elif flag_similar:
                    if similar(s_0, s_cnt) > 0.95:
                        if self.debug:
                            print "included[%4d]: %s   %8.4f   %8.3g" \
                                % (cnt, s_cnt, similar(s_0,s_cnt), p_cnt)
                        #
                        pr_s += manager.trace.lt[cnt].p
                        sim_include += 1
                    else:
                        if self.debug:
                            print "        [%4d]: %s   %8.4f   %8.3g" \
                                % (cnt, s_cnt, similar(s_0,s_cnt), p_cnt)
                        #
                    #
                elif flag_basic:
                    # just the first boltzmann probability 
                    if self.debug:
                        print "        [%4d]: %s   %8.3g" \
                            % (cnt, s_cnt, p_cnt)
                    #
                #
                else:
                    print "anal_loops(): Sorry, something is really fucked up"
                    sys.exit(1)
                #
            #
            
            # self.datatype == 'type1': // self.datatype == 'type2':
            # store the information in a temporary buffer
            
            dGbar = dG_0 / float(length)
            TdSbar = TdS_0 / float(length)
            output =  "%s  %5d  %8.2f  %8.2f  %8.3f  %8.3f   " \
                      % (self.cdata[ky].disp_data(), length, dG_0, TdS_0, dGbar, TdSbar)
            #
            if flag_allwts:
                output +=  "%8.5f  " % pr
                output +=  "%8.5f  %4d  " % (pr_s, sim_include)
                output +=  "%8.5f  %4d  " % (pr_h, ham_include)
                output +=  "%8.5f  %4d  " % (pr_TdS, dTdSp_include)
                output +=  "%8.5f  %4d  " % (pr_ddG, ddGp_include)
            elif flag_similar:
                output +=  "%8.5f  %4d  " % (pr_s, sim_include)
            elif flag_hamming:  # based on string edits
                output +=  "%8.5f  %4d  " % (pr_h, ham_include)
            elif flag_TdS:
                output +=  "%8.5f  %4d  " % (pr_TdS, dTdSp_include)
            elif flag_ddG:
                output +=  "%8.5f  %4d  " % (pr_ddG, ddGp_include)
            else:
                # boltzmann probabity of the first structure
                output +=  "%8.5f  " % pr
            #
            
            
            if self.datatype == 'type1':
                output += "  %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n" \
                          % (p_active, p_open, p_repressed, p_A, p_B)
            else:
                for kxs in epilist:
                    output += "  %8.5f" % self.cdata[ky].state[kxs]
                #
                output += "\n" 
            #
            
            fp = open(rflnm, 'a')
            fp.write(output)
            fp.close()
            
        #

        if not os.path.exists("loops_missing.txt"):
            fp = open("loops_missing.txt", 'w')
            s = "missing files:\n"
            for m in missing:
                s += m + '\n'
            fp.write(s)
            fp.close()
        #
        
        if not os.path.exists("loops_weights.txt"):
            fp = open("loops_weights.txt", 'w')
            fp.write("file                                          PET wts: \n")
            # 
            for r in rPET_wt:
                fp.write("%40s   %8.3f\n" % (r[1], r[0]))
            fp.close()
        #
        print "finished calculations"
        print "see results:           'loops_result.dat'"
        print "    missing files:     'loops_missing.txt'"
        print "    chromatin weights: 'loops_weights.txt'"
        print "DONE"
    #
    
#    

def main(cl):
    
    CL = GetOpts(PROGRAM)
    
    if len(cl) < 2:
        emsg = "ERROR: too few arguments"
        CL.error(emsg)
    #
    
    dt = Data()
    iflnms = CL.f_activity
    oflnm  = CL.f_output 
    print iflnms
    
    
    for f in iflnms:
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
    an.anal_loops(CL)
    
#

if __name__ == '__main__':
    # running the program
    main(sys.argv)
#

