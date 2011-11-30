#!/usr/bin/env python
import re
import sys
from util import read_fields

import argparse
parser = argparse.ArgumentParser(description='Append fields with positive strand coordinates.')
parser.add_argument("-f", type=int, default=0, help='Field number of the strand coordinates.')
OPTS = parser.parse_args()
line_tups = read_fields(f=open('danRer5.lengths.txt'))
zebrafish_lengths = dict([(tup[0], int(tup[1].replace(',',''))) for tup in line_tups])

coord_regex = re.compile(r'(?P<chrom>\w+):(?P<start>\d+)-(?P<end>\d+),(?P<dir>[+-])')
line_tups = read_fields()
for l in line_tups:
    coord = coord_regex.match(l[OPTS.f-1]).groupdict()
    if coord['dir'] == '-':
        length, start, end = zebrafish_lengths[coord['chrom']], int(coord['start']), int(coord['end'])
        coord['start'] = length - end
        coord['end'] = length - start
    sys.stdout.write('\t'.join(l + ['{chrom}:{start}-{end}\n'.format(**coord)]))
