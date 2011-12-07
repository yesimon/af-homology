Alignment Free Homology Detection
=================================
These are a bunch of scripts that help detect sequence homology
without the use of alignment. These mainly align a given training
set of regions from hg18 to danRer5 identified using conserved
flanks even though the middle region cannot be aligned.

Requirements
------------
 * Python 2.7, or 2.6 + argparse
 * Perl
 * Bash
 * Biopython
 * numpy, scipy, scikits.statsmodels, matplotlib
 * twoBitToFa (UCSC Kent)
 * faRc (UCSC Kent)

Pipeline
--------
To get the 2bit genomes for hg18 and danRer5 (over 1GB):
    
    ./get_data.sh

Obtain the formatted file used for most utilities here:

    ./extract_sequences.sh

Perform the alignment for D2z and HexMCD (takes a long time):

    cat hg18.toDanRer5.seqs.txt | ./af.py d2z
    cat hg18.toDanRer5.seqs.txt | ./af.py hexmcd

These produce pickled files, which can be converted to .dat as below:

    ./pkl_to_dat d2z.pkl > d2z.dat

For the alignments of HexDiff, HexYMF:

    cat hg18.toDanRer5.seqs.txt | ./hexdiff.pl region 1.0  > hexdiff.dat
    cat hg18.toDanRer5.seqs.txt | hexymf/hexymf.pl > hexymf.dat

Use scoring to generate ranked score reports using either the
overlap or ranked_peaks algorithm:

    ./scoring.py -f d2z.dat ranked_peaks > d2z_overlap

Generate cdf plots:

    ./scoring_cdf.py -f d2z_overlap hexmcd_overlap

You can also plot an individual cne's impulse function:

    ./plot.py -f d2z.dat <cne_name>
