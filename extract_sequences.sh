#!/usr/bin/env bash

# Get a sorted list of all genes.
GENES=$(sort -k4,4 hg18.bejscHumanCNE.toDanRer5Region.txt | cut -f4)

# Get a list of zebrafish negative strand sequences - faRc used for reverse complement.
NEG_SEQS=$(twoBitToFa -seqList=<(sort -k4,4 hg18.bejscHumanCNE.toDanRer5Region.txt | cut -f8 | grep ,- | ./negToPos.py | cut -f2) -noMask danRer5.2bit stdout | faRc stdin stdout -keepCase | ./fastaToTxt.py)

# Get a list of zebrafish positive strand sequences.
POS_SEQS=$(twoBitToFa -seqList=<(sort -k4,4 hg18.bejscHumanCNE.toDanRer5Region.txt | cut -f8 | grep ,+ | ./negToPos.py | cut -f2) -noMask danRer5.2bit stdout | ./fastaToTxt.py)

# All zebrafish sequences.
DANRER_SEQS=$(cat <(paste <(sort -k4,4 hg18.bejscHumanCNE.toDanRer5Region.txt | grep ,- | cut -f4) <(echo "$NEG_SEQS")) <(paste <(sort -k4,4 hg18.bejscHumanCNE.toDanRer5Region.txt | grep ,+ | cut -f4) <(echo "$POS_SEQS")) | sort -k1,1)

# All human sequences.
HG_SEQS=$(paste <(echo "$GENES") <(twoBitToFa -seqList=<(cut -f1,2,3 hg18.bejscHumanCNE.toDanRer5Region.txt | ./fieldsToCoord.py) -noMask hg18.2bit stdout | ./fastaToTxt.py) <(twoBitToFa -seqList=<(cut -f5,6,7 hg18.bejscHumanCNE.toDanRer5Region.txt | ./fieldsToCoord.py) -noMask hg18.2bit stdout | ./fastaToTxt.py) | sort -k1,1)

# Join human and zebrafish sequences.
SEQS=$(join -t $'\t' <(echo "$HG_SEQS") <(echo "$DANRER_SEQS"))

DANRER_CORRECT_COORDS=$(join -t $'\t' <(cat hg18.bejscHumanCNE.toDanRer5Region.txt | awk '{print $1"+"$2"+"$3"\t"$4}' | sort) <(cat hg18.danRer5.aligningRegions.txt | awk '{print $1"+"$2"+"$3"\t"$4}' | sort) | cut -f2,3 | sort -k1,1)

# Format is
# Name hg18chr hg18start hg18end danRer5Coord hg18Seq hg18Seq_long danRer5_seq danRer5_correct_coords
join -t $'\t' <(join -t $'\t' <(sort -k4,4 hg18.bejscHumanCNE.toDanRer5Region.txt | awk '{print $4"\t"$1"\t"$2"\t"$3"\t"$8}') <(echo "$SEQS")) <(echo "$DANRER_CORRECT_COORDS") > hg18.toDanRer5.seqs.txt


