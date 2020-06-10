#!/usr/bin/env python3

"""@@@

Main Program:  MolSystem.py 

Classes:       MolSystem
               

Author:        Wayne Dawson
creation date: 200205
last update:   200211 (upgraded to python3), 200205 
version:       0


Purpose:

This keeps track of whether the system is "RNA" or "Chromatin" and
what sort of free energy function is being used. 

This doesn't currently store the parameter set itself, but it allows
the poosibility of knowing what that set is supposed to be. That is
useful when moving from sequence to structure direction.


"""


import sys
import os

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

from Constants import rpr2num
from Constants import lpr2num

from Seq import Seq

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



sysDefLabels = { "Chromatin" : 'c',
                 "RNA"       : 'A'}

dParamSet = { "Chromatin" : {  "genheat"   : 1,
                               "setheat"   : 2 },
                "RNA"       : { "Turner"   : 1,
                                "ViS"      : 2,
                                "gMatrix"  : 3,
                                "ParFile"  : 4 } }




dJobTypes   = {  "testing"    : 0,
                 "prediction" : 1,
                 "evaluation" : 2,
                 "generator"  : 3 }

dJobTasks   = {
    "testing"    : "debugging or checking some particular module",
    "prediction" : "call a structure prediction program",
    "evaluation" : "energy evaluation of a specified structure (not a prediction)",
    "generator"  : "generate a contact map or heat map"
}

dProgList = {"cantata"           : "evaluation",
             "sarabande"         : "prediction",
             "sonata"            : "prediction",
             "chreval"           : "prediction",
             "testing"           : "testing",
             "HeartoftheSunrise" : "generator", 
             "generator"         : "generator" }

dEXTS = { "RNA" :       ["par",  "gMtrx" ],
          "Chromatin" : ["heat", "eheat"] }


def explain_JobTypes():
    s  = ("-------------------------------------------------------------\n")
    s += ("Types of jobs and their respective tasks:\n")
    s += ("vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n")
    for v in dJobTasks.keys():
        s += ("%-10s   %-s\n" % (v, dJobTasks[v]))
    #
    return s
#

def check_system(molsys):
    flag_system_known = True
    if not molsys in sysDefLabels:
        print ("ERROR: unknown system \"%s\"" % molsys)
        
        s = ''.join([("%s " % vk) for vk in list(sysDefLabels)])
        print ("allowed types ", s)
        
        sys.exit(1)
        flag_system_known = False
    #

    return flag_system_known
#


def check_JobType(jobtype):
    flag_prot_type_known = True
    if not jobtype in dJobTypes:
        print ("ERROR: unknown calculation type \"%s\"" % jobtype)
        
        s = ''.join([("%s " % vk) for vk in list(dJobTypes)])
        print ("allowed types: ", s)
        
        sys.exit(1)
        flag_prot_type_known = False
    #

    return flag_prot_type_known
#

def check_paramFl(jobtype):
    flag_prot_type_known = True
    if not jobtype in dJobTypes:
        print ("ERROR: unknown calculation type \"%s\"" % jobtype)
        
        s = ''.join([("%s " % vk) for vk in list(dJobTypes)])
        print ("allowed types: ", s)
        
        sys.exit(1)
        flag_prot_type_known = False
    #

    return flag_prot_type_known
#




"""@

Generate a G C U type sequence from an input structure. Useful when
you need some sort of sequence to calculate a free energy but you
don't necessarily know what to put there. You give it a structure and
it assigns a G to '(', a 'C' to ')' and a 'U' to '.'.

"""
def genRNASeq(ss_seq):
    debug = False # True # 
    rnaseq = ''
    
    
    # input is a structural sequence that satisfies inputs to a
    # program like VARNA. This is a very dumb tool, if we have any of
    # the lpr2num keys, it wil generate a 'G' and if it has any of the
    # rpr2num keys, it will generate a 'C'. Everything else is a 'U'
    
    if len(ss_seq) > 0:
        lss = list(ss_seq)
        lsq = ['x' for k in range(len(lss))] # range starts at 0 anyway
        # print (lsq)
        
        for k in range(0, len(lss)):
            if   lss[k] in lpr2num:
                lsq[k] = 'G'
            elif lss[k] in rpr2num:
                lsq[k] = 'C'
            elif lss[k] == '.':
                lsq[k] = 'u'
            else:
                print ("ERROR; unrecognized symbol \"%s\" in " % lss[k])
                print ("       structure sequence \"%s\"" % ss_seq)
                sys.exit(1)
            #
        #
        
        rnaseq = ''.join(lsq)
        
        
    #
    if debug:
        print ("ss_seq: \"%s\"" % ss_seq)
        print ("rnaseq: \"%s\"" % rnaseq)
    #
    
    return rnaseq
#


"""@

Generate a x y c type sequence from an input structure. Useful when
you need some sort of sequence to calculate a free energy but you
don't necessarily know what to put there. You give it a structure and
it assigns a 'x' to '(', a 'y' to ')' and a 'c' to '.'.

"""
def genChrSeq(ss_seq):
    debug = False # True # 
    chrseq = ''
    
    # input is a structural sequence that satisfies inputs to a
    # program like VARNA. This is a very dumb tool, if we have any of
    # the lpr2num keys, it wil generate a 'G' and if it has any of the
    # rpr2num keys, it will generate a 'C'. Everything else is a 'U'
    
    if len(ss_seq) > 0:
        lss = list(ss_seq)
        lsq = ['c' for k in range(len(lss))] # range starts at 0 anyway
        # print (lsq)
        
        for k in range(0, len(lss)):
            if   lss[k] in lpr2num:
                if lss[k] == '{':
                    lsq[k] = 'W'
                else:    
                    lsq[k] = 'x'
                #
                
            elif lss[k] in rpr2num:
                if lss[k] == '}':
                    lsq[k] = 'Z'
                else:    
                    lsq[k] = 'y'
                #
                
            elif lss[k] == '.':
                lsq[k] = 'c'
            else:
                if lss[k] == '|':
                    lsq[k] = 'I'
                else:
                    print ("ERROR; unrecognized symbol \"%s\" in " % lss[k])
                    print ("       structure sequence \"%s\"" % ss_seq)
                    sys.exit(1)
                #
                
            #
            
        #|endfor
        
        chrseq = ''.join(lsq)
        
        
    #
    if debug:
        print ("ss_seq: \"%s\"" % ss_seq)
        print ("chrseq: \"%s\"" % chrseq)
    #
    
    return chrseq
#



class MolSystem(object):
    """@
    
    The advantage of this thing is that it manages all this stuff in a
    variety of ways. If you just use the system name without anything,
    it forces the default system paramTypes, if you specify a
    paramType, then it automatically sets the system too.
    
    So, passing this around to the various modules, it should actually
    get there with all the requred information all bundled
    properly. At least, that is the current idea.
    
    """
    
    def __init__(self):
        self.system   = "none"
        """@
        
        system: "RNA", or "Chromatin"
        
        200329: In general, system, jobtype, program, paramType have
        some kind of relation to each other. For example, if the
        jobtype is "evaluation", then the program might be "cantata",
        which also means that it is RNA and there are various free
        energy functions Turner, ViS, ParamFl, and gMatrix. However,
        the exact relation between them is not so easy to
        establish. Perhaps in the future it will be. The point is that
        we have a compact wiring that we pass through the system to
        tell all the submodules what we are doing, and from there, we
        can decide what we want to do with the information.
        
        I realize that setting up InputSettings and so forth becomes
        all the harder with this, but as we pass this module that is
        ensconced in the underbelly of this program, we have some way
        to be sure that all these modules know what the original
        calling program was and (presumably) what the intent was from
        the start. This is helpful.

        """
        
        self.jobtype = "unassigned"
        """@ 
        
        jobtype options: "prediction", "evaluation", "generator", or
        "testing"
        
        jobtype is very important, particularly in telling other
        modules that this is a prediction rather than some contrived
        test or an evaluation. I haven't firmly decided what is the
        difference between "testing" and "evaluation", but maybe there
        could be some difference.

        """
        
        self.program = "undefined"
        """@
        
        program is not important presently. Mainly useful in telling
        something about what ultimately called the program. So
        presently, its default (undefined) is also ok.
        
        """
        
        self.parFlnm      = "none"
        """@
        
        parameter configuration file: 
        
        RNA: extension is gMtrx or par
        
        Chr: extension is heat or eheat
        
        none means Turner rules or default Chromatin rules
        
        """
        
        self.paramType    = "none"
        # RNA free energy parameter options
        self.useTurner    = False
        self.useViS       = False
        self.useParamFl   = False
        self.usegMatrix   = False
        
        # chromatin free energy parameter options
        self.useChromatin = False
        
        
        # sequence/structure information 
        self.mseq         = ""    # class str
        self.mSeq         = None  # class Seq
        self.mstr         = ""    # class str
        
    #
    
    
    
    def set_system(self, s):
        check_system(s)
        self.system = s
        if s == "RNA":
            self.set_useTurner()
        elif s == "Chromatin":
            self.set_useChromatin()
        #
        
    #
    
    def set_program(self, prgm):
        # this is probably not a particularly important variable, but
        # maybe it might eventually be important for
        # something. Anyway, it is maybe nice to know from whence it
        # came.
        self.program = prgm
    #
    
    def set_JobType(self, ptype):
        jobtype = ptype
        if ptype in dProgList:
            jobtype = dProgList[ptype]
        #
        
        check_JobType(jobtype)        
        self.jobtype = jobtype # types are "prediction", "evaluation"
        
    #
    
    def set_mstr(self, struct):
        self.mstr = struct
    #
    
    def set_mseq(self, seq):
        self.mseq = seq
        self.mSeq = Seq(self.mseq, self.system)
    #
    
    
    def set_ParamType(self, parSet, parFlnm = "none"):
        flag_set = True
        if   parSet == "Chromatin":
            self.set_useChromatin()
            self.paramType = "genheat"
            # the default uses pseudo heat values from dot bracket
            # information
            
        elif parSet == "RNA":
            self.set_useTurner()
            self.paramType = "Turner"
            # the default Turner rules
            
        else:
            if self.system == "Chromatin":
                
                if parSet == "genheat":
                    #print ("xx genheat")
                    self.set_useChromatin()
                    self.paramType = "genheat"
                    # uses pseudo heat values from dot bracket information
                    
                elif parSet == "setheat":
                    #print ("yy setheat")
                    self.set_useParamFl("Chromatin", parFlnm)
                    self.set_useChromatin()
                    self.paramType = "setheat"
                    # information obtained from an experimental heat map
                    
                else:
                    flag_set = False
                #
                
            elif self.system == "RNA":
                
                
                if parSet == "Turner":
                    self.set_useTurner()
                    self.paramType = "Turner"
                    # the default di-nucleotide Turner rules
                    
                elif parSet == "ViS": 
                    self.set_useViS()
                    self.paramType = "ViS"
                    # special setting of di-nucleotide par file
                    
                elif parSet == "gMatrix":
                    print ("x2")
                    self.set_usegMatrix("RNA", parFlnm)
                    self.paramType = "gMatrix"
                    # A new format that allows mono/di/tri-nucleotide bps
                    
                elif parSet == "ParFile":
                    print ("y2")
                    self.set_useParamFl("RNA", parFlnm)
                    self.paramType = "ParFile"
                    # alternative di-nucleotide par file
                else:
                    flag_set = False
                #
                
            else:
                flag_set = False
            #
            
        #
        
        if not flag_set:
            print ("ERROR: unrecognized parameter set: %s" % parSet)
            s = "allowed param types for %s: " % self.system
            if self.system in dParamSet:
                for v in dParamSet[self.system].keys():
                    s += "%s " % v
                #|endfor
            else:
                s += "unknown"
            #
            
            print (s)
            
        #
        
    #
    
    
    
    def set_useTurner(self):
        self.system       = "RNA"
        self.useTurner    = True
        self.useViS       = False
        self.useParamFl   = False
        self.usegMatrix   = False
        self.useChromatin = False
    #
    
    def set_useViS(self):
        import RNARefPath
        
        self.system       = "RNA"
        RNA_file_path     = os.path.dirname(RNARefPath.__file__)
        self.parFlnm      = RNA_file_path + "/ViSparams.par"
        self.useTurner    = False
        self.useViS       = True
        self.useParamFl   = False
        self.usegMatrix   = False
        self.useChromatin = False
    #
    
    
    
    def set_useParamFl(self, system, parFlnm):
        # here we have to specify the system and the parameter file
        # name, otherwise, the program doesn't know what to do.
        
        self.system       = system
        self.parFlnm      = parFlnm
        if not os.path.isfile(self.parFlnm):
            print ("ERROR: cannot find %s parameter file (%s)." \
                   % (self.system, self.parFlnm))
            sys.exit(1)
        #
        
        self.useTurner    = False
        self.useViS       = False
        self.useParamFl   = True
        self.usegMatrix   = False
        self.useChromatin = False
    #
    
    
    
    def set_usegMatrix(self, system, parFlnm):
        # here we have to specify the system and the parameter file
        # name, otherwise, the program doesn't know what to do.
        
        self.system       = system
        self.parFlnm      = parFlnm
        if not os.path.isfile(self.parFlnm):
            print ("ERROR: cannot find %s parameter file (%s)." \
                   % (self.system, self.parFlnm))
            sys.exit(1)
        #
        
        self.useTurner    = False
        self.useViS       = False
        self.useParamFl   = False
        self.usegMatrix   = True
        self.useChromatin = False
    #
    
    
    
    def set_useChromatin(self):
        self.system       = "Chromatin"
        self.useTurner    = False
        self.useViS       = False
        self.useParamFl   = False
        self.usegMatrix   = False
        self.useChromatin = True
    #
    
    
    def get_system(self):
        return self.system
    #
    
    def show_FEdataSet(self):
        s = ''
        if self.useChromatin:
            s = "energy function from Chromatin"
        elif self.useTurner:
            s = "Energy function dinucleotide from Turner"
        elif self.useViS:
            s = "Energy function dinucleotide from ViS"
        elif self.usegMatrix:
            s = "Energy function gMatrix parameter from %s" % self.parFlnm
        elif self.useParamFl:
            s = "Energy function from par file from %s" % self.parFlnm
        else:
            s = "no energy function available"
            print (s)
            sys.exit(1)
        #
        
        return s
    #
    
    def __str__(self):
        
        s = "%s  %s  %s %s\n" % (self.system, self.paramType, self.parFlnm, self.jobtype)
        if len(self.mseq) > 0:
            s += "seq: %s\n" % self.mseq
        else:
            s += "seq; -\n"
        #
        
        if len(self.mstr) > 0:
            s += "str: %s" % self.mstr
        else:
            s += "str: -"
        #
            
        return s
    #
    
    def __repr__(self):
        return self.__str__()
    #
    
#


def test0():
    
    ss_str = "(((...)))"
    
    
    print ("generic MolSystem:")
    a = MolSystem()
    print (a)
    #print ("stop at 1a"); sys.exit(0)
    
    print ("Chromatin:")
    chrseq = genChrSeq(ss_str)
    print ("chrseq: ", chrseq)
    a.set_system("Chromatin")
    a.set_mseq(chrseq)
    a.set_mstr(ss_str)
    a.set_useChromatin()
    print (a)
    print ("=======")
    #print ("stop at 2a"); sys.exit(0)
    
    a.set_JobType("evaluation")
    a.set_program("MolSystem::test0") # "program" is not so fussy presently
    a.set_ParamType("genheat")
    print (a)
    print ("=======")
    #print ("stop at 3a"); sys.exit(0)
    
    a.set_JobType("prediction")
    a.set_program("chreval") # "program" is not so fussy presently
    # have to specify your own directory structure here
    a.set_ParamType("setheat", "/home/dawson/python/chromatin/tests/test.heat")
    print (a)
    print ("=======")
    #print ("stop at 4a"); sys.exit(0)
    
    print ("RNA:")
    rnaseq = genRNASeq(ss_str)
    print ("rnaseq: ", rnaseq)
    
    b = MolSystem()
    b.set_system("RNA")
    b.set_mseq(rnaseq)
    b.set_JobType("sarabande")
    b.set_program("MolSystem::test0") # "program" is not so fussy presently
    print (b)
    print ("=======")
    #print ("stop at 1b"); sys.exit(0)
    
    b.set_mstr(ss_str)
    b.set_JobType("cantata")
    b.set_program("cantata") # "program" is not so fussy presently
    print (b)
    print ("=======")
    #print ("stop at 2b"); sys.exit(0)
    
    b.set_ParamType("ViS")
    print (b)
    print ("=======")
    #print ("stop at 3b"); sys.exit(0)
    
    b.set_ParamType("fake")
    print (b)
    print ("=======")
    #print ("stop at 4b"); sys.exit(0)
    gmtrxflnm = "/home/dawson/python/RNA/test3s_1mm+1mmp1+2mm+3mm_3nt-w5x5.gMtrx"
    b.set_ParamType("gMatrix", gmtrxflnm)
    print (b)
    print ("=======")
    #print ("stop at 5b"); sys.exit(0)
    
    b.set_ParamType("ParFile", "ViSparams.par")
    print (b)
    print ("=======")
    #print ("stop at 6b"); sys.exit(0)
    
    print (explain_JobTypes())
    
#


if __name__ == '__main__':
    # running the program
    test0()
#    
