#!/usr/bin/env bash

cat hg18.bejscHumanCNE.toDanRer5Region.txt | sort -k1,1 -k2,2n | cut -f4,8 | cut -d, -f1 > tmp1

cat hg18.danRer5.aligningRegions.txt | sort -k1,1 -k2,2n > tmp2

paste tmp1 tmp2 > ans.txt
rm tmp1 tmp2

