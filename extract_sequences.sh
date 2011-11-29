#!/usr/bin/env bash

# Get a sorted list of all genes.
GENES=$(sort -k4,4 hg18.bejscHumanCNE.toDanRer5Region.txt | cut -f4)

# Get a list of zebrafish negative strand sequences - faRc used for reverse complement.
NEG_SEQS=$(twoBitToFa -seqList=<(sort -k4,4 hg18.bejscHumanCNE.toDanRer5Region.txt | cut -f8 | grep ,- | ./negToPos.py | cut -f2) -noMask danRer5.2bit stdout | faRc stdin stdout -keepCase | ./fastaToTxt.py)

# Get a list of zebrafish positive strand sequences.
POS_SEQS=$(twoBitToFa -seqList=<(sort -k4,4 hg18.bejscHumanCNE.toDanRer5Region.txt | cut -f8 | grep ,+ | ./negToPos.py | cut -f2) -noMask danRer5.2bit stdout | ./fastaToTxt.py)

# All zebrafish sequences.
DANRER_SEQS=$(cat <(paste <(sort -k4,4 hg18.bejscHumanCNE.toDanRer5Region.txt | grep ,- | cut -f4) <(echo "$NEG_SEQS")) <(paste <(sort -k4,4 hg18.bejscHumanCNE.toDanRer5Region.txt | grep ,+ | cut -f4) <(echo "$POS_SEQS")))

# All human sequences.
HG_SEQS=$(paste <(echo "$GENES") <(twoBitToFa -seqList=<(cut -f1,2,3 hg18.bejscHumanCNE.toDanRer5Region.txt | ./fieldsToCoord.py) -noMask hg18.2bit stdout | ./fastaToTxt.py) <(twoBitToFa -seqList=<(cut -f5,6,7 hg18.bejscHumanCNE.toDanRer5Region.txt | ./fieldsToCoord.py) -noMask hg18.2bit stdout | ./fastaToTxt.py))

# Join human and zebrafish sequences.
SEQS=$(join -t $'\t' <(echo "$HG_SEQS" | sort -k1,1) <(echo "$DANRER_SEQS" | sort -k4,4))

echo "$SEQS" > hg18.toDanRer5.seqs.txt
