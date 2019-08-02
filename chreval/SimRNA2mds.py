#!/usr/bin/env python

"""@@@

Main Module:   SimRNA2mds.py 

Classes:       Rename 
               Atom 
               DataPDB

Functions:     doSimRNA2mds

Author:        Wayne Dawson
creation date: parts 2016, made into a separate object 170314
last update:   170314
version:       0

Purpose:

After the SimRNA simulation is finished, a file containing the coarse
grained structure of RNA is obtained in the form of a pdb file. For
chromatin, this file has too many residues and the RNA names for these
beads makes no sense. This tool is meant to convert the data to a
special format that is popular in mds type simulations, thus the name
SimRNA2mds.

"""

import sys
from os import path
from FileTools import FileTools

PROGRAM = "SimRNA2mds.py"
EXTS = ["pdb"]

# REF = {"P" : [] , "C4'" : [], "N1" : [], "C2" : [], "C4" : [], "N9" : [], "C6" : [] }
REF = {"P" : None }


def usage():
    print "USAGE: %s pdbfile" % PROGRAM
#

# ATOM      1  B   BEA A   1    1506.7671833.735  62.005  0.00  0.00           C
# ATOM      1  P     A A   1      41.995   5.331  84.010  1.00  0.00


class Rename:
    def __init__(self, aname, rname, atom):
        self.aname = aname
        self.rname = rname
        self.atom  = atom
    #
#

class Atom:
    def __init__(self):
        self.ndx   = -1
        self.aname  = ''
        self.rname = ''
        self.chain = ''
        self.rnumb = -1
        self.x     = None
        self.y     = None
        self.z     = None
        self.bfct  = 0.0
        self.var   = 0.0
        self.atom  = ''
    #
    def set_ATOM(self, ndx, aname, rname, chain, rnumb, x, y, z, bfct = 0.0, var = 0.0, atom = 'X'):
        self.ndx   = ndx
        self.aname = aname
        self.rname = rname
        self.chain = chain
        self.rnumb = rnumb
        self.x     = x
        self.y     = y
        self.z     = z
        self.bfct  = bfct
        self.var   = var
        self.atom  = atom
    #
    def disp_ATOM(self):
        #    ATOM      1  B   BEA A   1    1506.7671833.735  62.005  0.00  0.00           C
        s = "ATOM  %5d %4s %s %s%4s    %8.3f%8.3f%8.3f %5.2f %5.2f           %s" \
            % (self.ndx, \
               self.aname, \
               self.rname, \
               self.chain, \
               self.rnumb, \
               self.x, \
               self.y, \
               self.z, \
               self.bfct, \
               self.var, \
               self.atom)
        return s
        #
    #
#

class DataPDB:
    def __init__(self):
        self.atoms = []
        self.heteroatoms = []
        self.connect = []
        self.ref = REF
        self.ref['P'] = Rename('B  ', 'BEA', 'C')
    #

    def readpdb(self, flnm):
        try:
            fp = open(flnm, 'r')
        except IOError:
            print "ERROR: cannot open %s" % flnm
            sys.exit(1)
        #
        
        lfp = fp.readlines()
        b = 0.0; v =0.0; a = 'X'
        count = 0
        acnt = 0
        for lfpk in lfp:
            lbl = lfpk[0:6]
            # print lbl
            if lbl == "ATOM  ":
                anm = lfpk[12:16]
                rnm = lfpk[17:20]
                chn = lfpk[21:22]
                rnb = lfpk[22:26]
                # print ndx, anm, rnm, rnb
                if self.ref.has_key(anm.strip()):
                    acnt += 1
                    ndx = acnt
                    try:
                        x = float(lfpk[30:38])
                        y = float(lfpk[38:46])
                        z = float(lfpk[46:54])
                    except ValueError:
                        print "ERROR: problems reading coordinate data in pdb file %s" % flnm
                        print "       line number %d"
                        print "REMARK                       >_???.??|_???.??|_???.??|<"
                        print lfpk.strip()
                        sys.exit(1)
                    try:
                        b = float(lfpk[54:60])
                    except ValueError:
                        b = 0.0
                    try:
                        v = float(lfpk[60:66])
                    except ValueError:
                        v = 0.0
                    a = lfpk[72:80]
                    # print ndx, anm, rnm, rnb, len(a)
                    if len(a) == 0:
                        a = anm.strip()
                        atom = ''
                        for ak in a:
                            if not ak.isdigit() and ak.isalpha():
                                atom += ak
                        a = atom
                    #
                    r = Atom()
                    anm = self.ref['P'].aname
                    rnm = self.ref['P'].rname
                    a   = self.ref['P'].atom
                    
                    r.set_ATOM(ndx, anm, rnm, chn, rnb, x, y, z, b, v, a)
                    self.atoms += [ r ]
                #
            elif lbl == "CONECT":
                #print lfpk.strip()
                self.connect += [ lbl.strip() ]
            elif lbl == "HETATM":
                #print lfpk.strip()
                self.heteroatoms += [ lbl.strip() ]
                
            #
                                
            #
        #
    #
    def show_ATOM(self):
        s = ''
        for xyzk in self.atoms:
            s += xyzk.disp_ATOM() + '\n'
        return s
    def show_CONECT(self):
        #CONECT    1    2
        #CONECT    2    1    3
        #CONECT    3    2    4
        # .....
        #CONECT  168  167  169
        #CONECT  169  168  170
        #CONECT  170  169  171
        #CONECT  171  170
        s = ''
        ix = 1
        s = "CONECT%5d%5d\n" % (ix, ix+1)
        for i in range(1, len(self.atoms)-1):
            s += "CONECT%5d%5d%5d\n" % (i+1, i, i+2)
        ix = len(self.atoms)
        s += "CONECT%5d%5d\n" % (ix, ix-1)
        return s
    #
#
    


def doSimRNA2mds(cl):

    if len(cl) == 1:
        usage()
        sys.exit(0)
    #
    if cl[1] == '-h' or cl[1] == '--help':
        usage()
        sys.exit(0)
    #
    pdbflnm_in = cl[1]
    if not path.isfile(pdbflnm_in):
        print "ERROR: cannot find %s" % pdbflnm_in
        usage()
        sys.exit(1)
    #
    
    ft = FileTools()
    ft.check_ext(pdbflnm_in, EXTS, PROGRAM)
    pdbflnm_out = ft.flhd + "_v.pdb"
    if len(cl) == 3:
        pdbflnm_out = cl[2]
        ft.check_ext(pdbflnm_out, EXTS, PROGRAM)
    #
    
    dp = DataPDB()
    dp.readpdb(pdbflnm_in)
    
    s = dp.show_ATOM()
    fp = open(pdbflnm_out, "w")
    fp.write(s)
    fp.write("TER\n")
    s = dp.show_CONECT()
    fp.write(s)
    fp.write("END\n")
    fp.close()
#

if __name__ == '__main__':
    # running the program
    doSimRNA2mds(sys.argv)


