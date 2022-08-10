   BWOLO


Bwolo is a fast and memory efficient tool for approximate string matching over
the DNA alphabet. It is specifically designed for short patterns (< 50nt) over
large texts (10^8-10^9) with a high error rate (7%-15%).

Bwolo is implemented in C++ and freely distributed under the GNU lesser general
public license (LGPL).


WEBSITE
=======

The latest version of Bwolo can be downloaded from http://bioinfo.lifl.fr/bwolo


TL:DR (Using Bwolo in a hurry)
==============================

Before building you need :
  - SeqAn 1.4 in your building path (SeqAn.de)
  - Zlib installed
  - POSIX realtime extension (certainly ever installed)
  - G++ compatible with C++11

To build Fasta2Fmi and Bwolo, just use the command “make”.

Index the text :
fasta2Fmi -f [TextFasta] -i [OutputIndex]

Search one or several patterns :
bwolo -i [TextIndex] -p [pattern] -s “-[score]”
or
bwolo -i [TextIndex] -f [PatternMultiFasta] -s “-[score]”

Each output line is written as
[OccurenceBeginPosition]:[OccurenceEndPosition]:[distance]:[OccurenceSequence]:[PatternSequence]


BUILDING FROM SOURCE
====================

Before building
---------------

Download and install the library SeqAn (seqan.de) at
"http://packages.seqan.de/" or at "https://github.com/seqan/seqan".
Make sure it is in your building path.

Bwolo and fasta2Fmi need at least the version 1.4 of the library.
They have been tested with the version 1.4 and 1.4.1.

Other requirements :
  - ZLib library if it is not yet installed.
  - POSIX realtime extension, but it is certainly installed.
  - G++, C++11 or at least C++0x compatible (GCC>=4.3).

Building fasta2Fmi and Bwolo
----------------------------

Download the source archive from http://bioinfo.lifl.fr/bwolo. Decompress the
archive. Run the GNU 'make' command with no argument. If you use a version of
G++ between 4.3 and 4.6 run 'make STD=c++0x'.

You can now copy the binaries found in the respective folders to a folder in
your PATH variable.


HOW TO LAUNCH BWOLO
===================

Index the text with fasta2Fmi
-----------------------------

Before all, you need to index your text with fasta2Fmi. The text must be in
FASTA format (the NCBI FASTA specification is available at this address
http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml). Note that if your FASTA file
contains several sequences, only the first one will be indexed.

fasta2Fmi uses the SeqAn implementation of the FM-Index and its algorithm to
build it.  The index will be stored into several files (for wavelet tree,
samples of the suffix array and original text).
The options are :

  -f --fasta [required][required_argument:text] : Fasta file of the text.
  -i --index [required][required_argument:text] : Output prefix for files of the index.
  -v --verbose [optional] : Verbose mode.

__Example__

fasta2Fmi -f myText.fasta -i ./indexFolder/myText

This will index the first sequence in the file "myText.fasta" and save the
created index into the "./indexFolder/"  folder in many files.
All file names begin with "myText".

Search your pattern(s)
----------------------

You can search a single pattern by using the -p option or a collection of
patterns available in a multiFASTA file using the -f option.

For the text, Bwolo needs not the original FASTA file. It only needs the files
generated previously by fasta2Fmi. You should indicate the prefix of the file
name of the index as an argument of the -i option.

The maximal number of allowed errors is given in the score option. The score is
a negative value. The default value is "-3". It means that the search algorithm
will return all occurrences containing up to 3 errors. An error can be either a
deletion, an insertion or a substitution.

Bwolo returns the results in the standard output. Each line represents an
occurrence, in this format:

[OccurenceBeginPosition]:[OccurenceEndPosition]:[distance]:[OccurenceSequence]:[PatternSequence]

Occurrences are sorted by pattern (in the order of the fasta), then by start
position.

The options are :

  -p --pattern [optional*][required_argument:text] : pattern to search.
  -f --fasta [optional*][required_argument:text] : fasta file containing the patterns to search.
  -i --index [required][required_argument:text] : Output file for the index.
  -s --score [optional][requiered_argument:int] : minimal alignment score. Default : -3
  -v --verbose [optional] : Verbose mode.
  -h --help print this message

* : requires one of them

__Examples__

bwolo -s -2 -f myPatterns.fasta -i ./indexFolder/myText

It will search all the patterns contained in the FASTA file "myPatterns.fasta"
in the indexed text in "./indexFolder/myText".

bwolo -s -2 -p "TGACAGAAGAGAGTGAGCAC" -i ./indexFolder/myText

It will search the pattern "TGACAGAAGAGAGTGAGCAC" in the indexed text in
"./indexFolder/myText".


ADDITIONAL INFORMATION
======================

Changing the alphabet
---------------------

You can edit the two files MirTypeConfig.h in the folders
"./bwolo/include/" and "./fasta2Fmi/src/" before compiling  to change
the alphabet used by modifying  these lines :

typedef Dna BwoloTIndexAlphabet;
typedef Dna BwoloTPatternAlphabet;

You can replace “Dna” values by “Dna5”, Rna”, “Rna5”, “Iupac”
“AminoAcid” or “char”. This is at your own risk.


KNOWN ISSUES
============

Some overlapping occurrences may be returned. In this case, they are 2 different
alignments of the same occurrence. They differ by more than k nucleotides,
thus they are considered different.

CONTACT
=======

If you have problems or questions, please contact us at bwolo@univ-lille1.fr


AUTHORS
=======

Christophe Vroland
