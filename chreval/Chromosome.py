#!/usr/bin/python

"""@@@

Main Module:   Chromosome.py

Classes:       Segment
               Chromatin
               Chromosome

Author:        Wayne Dawson
creation date: parts 2016 (in chreval), made into a separate object 170426
last update:   180709
version:       0

Purpose:

this contains a number of important data structures that are used in
ChromatinData and a few other programs.


Comments:

Well, one thing is that it is kind of annoying that I have so many
different data formats that are different yet very similar. I would
like to streamline this, but presently, I think the philosophy is "if
it works, don't fix it!".

"""

import sys
import os
from CUtils import String

PROGRAM = "Chromosome.py"
string = String(10)


allowed_tags = ['chr', 'bgn', 'end',
                'lCTCF', 'rCTCF',
                'PETcnt', 'cmplx1', 'cmplx2', 'cmplx3',
                'open', 'active', 'A', 'B',
                'repressed',
                'vi', 'vj']

# tag definitions:
# 'chr'        chromosome index
# 'bgn'        beginning locus on the chromosome
# 'end'        ending locus on the chromosome
# 'lCTCF'      orientation of CTCF on the 5' side, sense strand and when applicable
# 'rCTCF'      orientation of CTCF on the 3' side, sense strand and when applicable
# 'PETcnt'     number of PET clusters counted for this loop
# 'cmplx1'     complexity of the CTCF anchors
# 'cmplx2'     complexity of the RNA polymerase II anchors
# 'cmplx3'     complexity of both CTCF or RNA polymerase II anchors
# 'open'       degree of openness of the chromatin in the locus
# 'active'     degree of activeness of the chromatin in the locus
# 'A'          is it compartment A?
# 'B'          is it compartment B?
# 'repressed'  For old *.txt files, opposite of active
# 'vi'         starting locus of the chromatin relative to the boundaries of the domain
# 'vj'         ending locus of the chromatin relative to the boundaries of the domain


def assign_bed_tags(taglist, sbedfl):
    name = ''    # name of the chromosome
    bgn = -1     # beginning locus on the chromosome
    end = -1     # ending locus on the chromosome
    ctcf1 = ''   # right side CTCF
    ctcf2 = ''   # right side CTCF
    nPET = 0     # number of PET counts
    cmplx = []   # complexity
    popen  = 0.0 # open
    pactv  = 0.0 # active
    cmprtA = 0.0 # compartment A
    cmprtB = 0.0 # compartment B
    vi     = -1
    vj     = -1
    cmplx = []
    # print "taglist: ", taglist
    # print "sbedfl:  ", sbedfl
    for tlk in range(0, len(taglist)):
        try:
            if taglist[tlk] == 'chr':
                name  = sbedfl[tlk]
            elif taglist[tlk] == 'bgn':
                bgn = int(sbedfl[tlk])
            elif taglist[tlk] == 'end':
                end = int(sbedfl[tlk])
            elif taglist[tlk] == 'lCTCF':
                ctcf1 = sbedfl[tlk]
            elif taglist[tlk] == 'rCTCF':
                ctcf2 = sbedfl[tlk]
            elif taglist[tlk] == 'PETcnt':
                nPET = int(sbedfl[tlk])
            elif taglist[tlk] == 'cmplx1':
                cmplx += [int(sbedfl[tlk])]
            elif taglist[tlk] == 'cmplx2':
                cmplx += [int(sbedfl[tlk])]
            elif taglist[tlk] == 'cmplx3':
                cmplx += [int(sbedfl[tlk])]
            elif taglist[tlk] == 'open':
                popen = float(sbedfl[tlk])
            elif taglist[tlk] == 'active':
                pactv = float(sbedfl[tlk])
            elif taglist[tlk] == 'A':
                cmprtA = float(sbedfl[tlk])
            elif taglist[tlk] == 'B':
                cmprtB = float(sbedfl[tlk])
            elif taglist[tlk] == 'vi':
                vi = float(sbedfl[tlk])
            elif taglist[tlk] == 'vj':
                vj = float(sbedfl[tlk])
            else:
                print "ERROR: undefined tag. "
                print "       allowed tags are the following:"
                for aa in allowed_tags:
                    print aa
                sys.exit(1)
            #
        except (ValueError):
            print "ERROR: incompatable data found for entry"
            print sbedfl
            sys.exit(1)
        #
    #
    chrmtn = Chromatin(name, bgn, end, ctcf1, ctcf2, nPET, cmplx)
    chrmtn.vi = vi; chrmtn.vj = vj
    chrmtn.add_segments(bgn, end, "-")
    chrmtn.state["active"] = pactv
    chrmtn.state["open"]   = popen
    chrmtn.state["A"]      = cmprtA
    chrmtn.state["B"]      = cmprtB
    
    
    return chrmtn
#



class Segment:
    def __init__(self, bgn, end, state):
        self.bgn    = bgn
        self.end    = end
        self.state  = state
        # Here, "state" this is used by the older *.txt data files
        # where a certain segment was either defined as "active" or
        # "inactive". So for this type of data, we have to store the
        # information as one of these types, not a specific
        # number. Based upon these types, we could then determine how
        # to the relative active state. For the *.bed files, "state"
        # is not used.
    #
    
    def disp(self):
        return "%10d  %10d  %10d  %s" % (self.bgn, self.end, (self.end - self.bgn), self.state)
    #
    
    def get_length(self):
        length = self.end - self.bgn
        return length
    #
#


# Serves as a container for Stores data on a particular chromosome, or
# part of the chromosome.
class Chromatin:
    """stores records of chromatin segments of a specific type"""
    def __init__(self, name, bgn, end, ctcf1 = '*', ctcf2 = '*', nPET = 0, cmplx = []):
        
        # For the "*.bed" files, this is the most important data
        # structure, Segment is mainly used to maintain compatibility
        # between the different data formats. In *.bed files, the
        # contents of Segment are exactly the same as self.bgn and
        # self.end and state is not used.
        
        self.name  = name       # name of the chromosome
        self.bgn   = int(bgn)   # beginning position
        self.end   = int(end)   # ending position
        self.ctcf1 = ctcf1      # CTCF direction at the beginning point
        self.ctcf2 = ctcf2      # CTCF direction at the ending point
        #  (vi,vj): The relative position in resolution units. This
        #  data available when data is read from an eheat file.
        self.vi    = -1         
        self.vj    = -1
        
        # number of PET counts
        if string.isFloat(nPET):
            self.nPET = float(nPET)
        else:
            self.nPET  = int(nPET)
        #
        
        self.cmplx = cmplx      # complexity
        self.state = { "active"    : 0.0,\
                       "open"      : 0.0,\
                       "repressed" : 0.0,\
                       "A"         : 0.0,\
                       "B"         : 0.0 }
        #
        self.segments = [] # mainly for type1 data
    #
    
    
    def add_segments(self, bgn, end, state):
        # print "add_segments"
        self.segments += [Segment(bgn, end, state)]
    #
    
    def disp_head(self):
        return "%5s  %10d  %10d" % (self.name, self.bgn, self.end)
    #
    
    def disp_data(self, fill_region = 0):
        s = ''
        s += "%5s  %10d  %10d      %s      %s" \
             % (self.name, self.bgn, self.end, self.ctcf1, self.ctcf2)
        s += "       p_%s" % string.field(self.nPET, 3)
        k = 1
        for c in self.cmplx:
            s += "    c%d_%s" % (k, string.field(c, 4))
            k += 1
        #
        
        # in case we have to add additional spaces for non existing
        # data and fit it within the same table as some other
        # data. This normally shouldn't be necessary, but there are a
        # variety of data types, and they could be drawn from a
        # variety of sources and mixed together for some reason.
        if fill_region > 0:
            for i in range(0, fill_region):
                s += "  c-_0     "
            #
        #
        
        if self.vi >= 0 and self.vj > 0:
            s += "  %8d  %8d" % (self.vi, self.vj)
        #
        
        return s
    #
#
        
class Chromosome:
    """stores information on the length and contents of each chromosome unit"""
    # this is the larger structure that contains the Chromatin sectors
    def __init__(self, name, bgn, end):
        self.name = name
        self.chrsegment = [(int(bgn), int(end))]
    #
    def add_chrsegment(self, bgn, end):
        self.chrsegment += [(int(bgn), int(end))]
    #
    def show_chr_fragments(self):
        s =  self.name + '\n'
        for chrsx in self.chrsegment:
            s += "%10d   %10d\n" % (chrsx[0], chrsx[1])
        return s
    #
#
