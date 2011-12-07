#!/usr/bin/env bash

# extract the human dna
cat hg18.bejscHumanCNE.toDanRer5Region.txt | awk '{print $1":"$2"-"$3}' > tmp.dat
twoBitToFa hg18.2bit tmp.fa -seqList=tmp.dat -noMask
cat tmp.fa | perl faToSeq.pl > first.txt
rm tmp.dat tmp.fa

# extract the extended human dna
cat hg18.bejscHumanCNE.toDanRer5Region.txt | awk '{print $5":"$6"-"$7}' > tmp.dat
twoBitToFa hg18.2bit tmp.fa -seqList=tmp.dat -noMask
cat tmp.fa | perl faToSeq.pl > second.txt
rm tmp.dat tmp.fa

# extract extended zebrafish dna
cat hg18.bejscHumanCNE.toDanRer5Region.txt | awk '{print $8}' | cut -d, -f1 > tmp.dat
twoBitToFa danRer5.2bit tmp.fa -seqList=tmp.dat -noMask
cat tmp.fa | perl faToSeq.pl > third.txt
rm tmp.dat tmp.fa

# switch over to our seq format (not totally sure why orginal format isn't fine)
cat hg18.bejscHumanCNE.toDanRer5Region.txt | awk '{print $4"\t"$1"\t"$2"\t"$3"\t"$8}' > tmp.txt

# and finally put it all together
paste tmp.txt first.txt second.txt third.txt > hg18.toDanRer5.seqs.txt
rm tmp.txt first.txt second.txt third.txt
