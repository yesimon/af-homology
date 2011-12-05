#!bin/bash

# this script downloads all the necessary data files. Typically these files
# are too large to be checked into version control. Be sure to add any
# new data files here, or else you might break other developer's code.

wget -nc ftp://hgdownload.cse.ucsc.edu/goldenPath/hg18/bigZips/hg18.2bit
wget -nc http://hgdownload.cse.ucsc.edu/gbdb/danRer5/danRer5.2bit

