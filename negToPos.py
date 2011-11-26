#!/usr/bin/env python
import re
import sys
from util import read_fields

line_tups = read_fields(f=open('danRer5.lengths.txt'))
zebrafish_lengths = dict([(tup[0], int(tup[1].replace(',',''))) for tup in line_tups])

coord_regex = re.compile(r'(?P<chrom>\w+):(?P<start>\d+)-(?P<end>\d+),(?P<dir>[+-])')
line_tups = read_fields()
coords = [coord_regex.match(l[0]).groupdict() for l in line_tups]
for coord in coords:
    if coord['dir'] == '-':
        length, start, end = zebrafish_lengths[coord['chrom']], int(coord['start']), int(coord['end'])
        coord['start'] = length - end
        coord['end'] = length - start
    sys.stdout.write("{chrom}:{start}-{end}\n".format(**coord))
