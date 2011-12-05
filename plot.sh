#!/usr/bin/env bash

cat $1 | grep $2 | cut -f2 > tmp.dat
perl ./plot.pl $2
rm tmp.dat

