#!/usr/bin/env bash

# All zebrafish sequences.
DANRER_SEQS=$(twoBitToFa -seqList=<(cut -f8 hg18.bejscHumanCNE.toDanRer5Region.txt | \
    awk 'sub("..$", "")') -noMask danRer5.2bit stdout | ./fastaToTxt.py)

# All human sequences.
HG_SEQS=$(twoBitToFa -seqList=<(cut -f1,2,3 hg18.bejscHumanCNE.toDanRer5Region.txt | \
    ./fieldsToCoord.py) -noMask hg18.2bit stdout | ./fastaToTxt.py)
# All extended human sequences.
HG_SEQS_EXTENDED=$(twoBitToFa -seqList=<(cut -f5,6,7 hg18.bejscHumanCNE.toDanRer5Region.txt | \
    ./fieldsToCoord.py) -noMask hg18.2bit stdout | ./fastaToTxt.py)

ALL_SEQS=$(paste \
    <(cat hg18.bejscHumanCNE.toDanRer5Region.txt | awk '{print $4"\t"$1"\t"$2"\t"$3"\t"$8}') \
    <(echo "$HG_SEQS") \
    <(echo "$HG_SEQS_EXTENDED") \
    <(echo "$DANRER_SEQS"))

# Get the correct danrer coordinates.
DANRER_CORRECT_COORDS=$(join -t $'\t' <(cat hg18.bejscHumanCNE.toDanRer5Region.txt | awk '{print $1"+"$2"+"$3"\t"$4}' | sort) <(cat hg18.danRer5.aligningRegions.txt | awk '{print $1"+"$2"+"$3"\t"$4}' | sort) | sort -k2,2 | cut -f3)

# Format is
# Name hg18chr hg18start hg18end danRer5Coord hg18Seq hg18Seq_long danRer5_seq danRer5_correct_coords
paste <(echo "$ALL_SEQS" | sort) <(echo "$DANRER_CORRECT_COORDS") > hg18.toDanRer5.seqs.txt
