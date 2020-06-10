# chreval

# README #

This README would normally document whatever steps are necessary to
get your application up and running.

### What is this repository for? ###

This repository contains the following executable programs.


chreval.py (CHRomatin EVALuation program for heatmaps)

+ this is the actually dynamic programming algorithm part of the
program. It calculates the optimal and suboptmal structures of
chromatin the Boltzmann distribution of the principal structures.


analyze_loops.py

+ this is a driver for carrying out analysis of the bed files with
corresponding heat maps. For a simple example, see the directory
tests/test_analyze_loops in this distribution and, for complete
calculations of AB compartments and CCDs using this tool, see the
directory results_in_manuscript/ in this distribution.


make_heatmap.py

+ generates a heatmap of the specified input heatmap


my_generation.py

+ generates a heatmap file based on user specified contacts (and
weights if desired). For information on the input format, run

> my_generation.py -hExFile
> my_generation.py -hExSeq



SimRNA_make_polyA.py

+ generates a poly(A) sequence of a specified length. These sequences
are used to help generate 3D structures from the *.simres files
generated while running chreval.py.



SimRNA2mds.py

+ converts the SimRNA3.21 generated pdb output files into a single
bead model (for building 3D representations of chromatin
.... presently, the coordinates are _not_ rescaled).
 


_Chreval_ is the main driver for calculating the free energy of
observed heatmaps.  Chreval is a dynamic programming method for
determining the most probably chromatin structure arrangement as well
as the distribution of chromatin structure arrangements as a function
of free energy of various motifs found in the structure.

_Analyze_loops_ is for handling a group of heatmap files with the
properly. The format is Przemek's bed format. The program calls
Chreval. You should make sure that the files listed in the bed file
also exist in your directory. Please see the directories tests and
results_in_manuscript for examples of using this tool.

_Make_heatmap_ is a tool to generate heatmaps either from *.heat files
or the output from Chreval *.clust

_My_generation_ is a program to generate heatmaps. Entries must be
given in dot bracket format. An example is provided in the command line
by using the flag -hExFile. 


* how to run Chreval?

To run, the minimum information you need is a heat map data file.
This file contains experimental data from a source (particularly
ChIA-PET but possibly Hi-C). This data (particularly ChIA-PET) is
highly correlated with positions of the CTCF binding proteins and
cohesin. The CTCF sites largely correspond to regions that form loops
and typically regulate the chromatin.

The heat maps consists of a symmetric matrix where the indices (i,j),
corresponding to row and column positions, indicate points on the
chromatin chain where segment i and j interact with each other. For
simplicity, we assume i < j, and look only at the upper triangle. A
reflection of matrix is found across the diagonal. The intensity of a
given point (i,j) is proportional to the frequency that this
particular interaction was encountered in the experimental data, so
small numbers mean only weak interactions, and large numbers mean
strong interactions.

Once you obtain a heat make, or have created one using
my_generation.py, you can run this program with the following command
line example:

> chreval.py -f myexample.heat

Specifically, in the directory "tests", the following command can be
used to calculate the example heatmap:

> cd tests
> chreval.py -f chr10_64313472_64921344_res5kb.heat

or, another example

> cd tests/test_analyze_loops/eheat_files
> chreval.py -f chr1_1890973_2316695.eheat


The file chr10_64313472_64921344_res5kb.heat has the labeling
information "chrN_x_y_res5kb.heat", where N is the chromosome number,
x is the starting position and y is the ending position. Further, the
"res5kb" indicates the resolution of the grid that was generated. If
the grid size is smaller or larger, an option can be used to change
this grid size from the default (presently 5 kb). The file must end
with the extension heat or the newer form "eheat" (extended heatmap
file).




The output directory from chreval in the above example will be
chr10_64313472_64921344_res5kb; i.e., the directory name
"chrN_x_y_res5kb". This directory contains (in separate files) the top
layer of suboptimal structures within some specifable energy range
from the mimumum free energy (default is 10 kcal/mol) or a fractional
percentage of the free energy. These files have the extension "DBN"
and can be read by the 3rd party program VARNA. Additionally, Chevral
is set up to provide heatmap, pairing info, and restraint files for
SimRNA 3D calculations. Two additional files are
chrN_x_y_res5kb_BDwt.clust that contains a matrix with the Boltzmann
probabilities for different interactions and
chrN_x_y_res5kb_summary.txt that contains a shorthand list of the
secondary structures.

There are a variety of additional options. Please run

> chreval.py -h 

to obtain additional information on additional command line options.


* How to run analyze_loops?

Examples of using analyze_loops.py are provided in the directories
"tests/test_anal_loops" and "results_in_manuscript" included in this
distribution.

> cd tests/test_analyze_loops
> analyze_loops.py -ff test_loops.CTCF.withAandB.annotated.bed

The program anal_loops.py will look up files in the same directory as
the *.bed file and try to compute the free energy using the object
Chreval.

For more information on how to run the program, please run

> analyze_loops.py -h


* How to run make_heatmap.py? Other programs

Here are some additional command line arguments of the other programs
in this set.


> make_heatmap.py chr1_1890973_2316695.eheat + extended heatmap files
+ makes a 2D heat map of the file. contain more detailed information. See the directory
"tests/test_analyze_loops/eheat_files"

> make_heatmap.py chr10_64313472_64921344_res5kb.heat
+ makes a heat map of the file. 


* How to run my_generation.py?

Here is an example.

> my_generation.py -seq ".ABCDE.((..)).abcde" "{.................}"

Note that it doesn't matter that one of the structures overlaps the
other one.


* How to run the SimRNA packages to obtain 3D structures from Chreval outputs?

You must download the executable version of SimRNA from the following
website

http://genesilico.pl/software/stand-alone/simrna

+ Copy the "data" directory and config.dat file to a separate directory where you want to build the 3D structure.

+ Copy the relevant simres file to that same directory. 

in data, change all the values in histograms3D_3.list to 0, except for the last four that have 0.1 in as the

> cat histograms3D_3.list
./data/AA3.hist    0.0
./data/AC3.hist    0.0

...

./data/A_3_exvol.hist 0.1
...

in the config.dat file, change the parameter ETA_THETA_WEIGHT from 0.4 to 0.0

> cat config.dat
...

ETA_THETA_WEIGHT 0.00


Build a sequence of the proper length (N) for poly(A)

> SimRNA_make_polyA.py N  > myseq.seq


Now run a replica exchange Monte Carlo simulation using SimRNA

SimRNA3 -s myseq.seq -r myseq.simres -c config.dat -E 10 -o myseq >& myseq.log &

To generate sequences, use the following

> cat myseq_x_??.trafl >> myseq_x.trafl
> clustering myseq_x.trafl 0.1 15.0

this generates the files

myseq_x_thrs15.00A_clust01.trafl


to convert the first trajectory in this file to a pdb representation:

> SimRNA_trafl2pdbs myseq_x_01-000001.pdb myseq_x.trafl 1

this generates the file
myseq_x_thrs15.00A_clust01-000001.pdb

Now we convert that to a single bead representation

SimRNA2mds.py myseq_x_thrs15.00A_clust01-000001.pdb

which generates a file

myseq_x_thrs15.00A_clust01-000001_v.pdb

this pdb file can be viewed in a recognizable for using chimera, vmd,
or pymol


----

> SimRNA_make_polyA.py 10
aaaaaaaaaa

> SimRNA2mds.py SimRNA_pdboutput.pdb

+ generates a PDB file formatted in a way so that (currently) the
phosphate atom is treated as the main binding interaction. The
particular atom can be changed in the program, and even more than one
atom displayed, if actually desired.

----




Version 1.0
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up

The package can run as is, but some applications may require the
installation of the following python packages, if not already
installed.

matplotlib, numpy, random and argparse


* Configuration

Standard python 2.7

* Dependencies

matplot, numpy, random and argparse; otherwise, none

* Database configuration

None

* How to run tests

some examples are provided in the directory test

* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

to consult about bugs and other issues of the code, please contact 

Wayne Dawson

Laboratory of Functional and Structural Genomeics

Center of New Technologies

University of Warsaw

Banacha 2C, 02-098 Warsaw

email w.dawson@cent.uw.edu.pl

* Repo owner or admin

Wayne Dawson
* Other community or team contact
