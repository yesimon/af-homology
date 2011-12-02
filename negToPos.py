#!/usr/bin/env python
import sys
from util import read_fields, parse_coords, COORD_FORMATS

line_tups = read_fields(f=open('danRer5.lengths.txt'))
zebrafish_lengths = dict([(tup[0], int(tup[1].replace(',',''))) for tup in line_tups])

def forward_strand_zebrafish(coord):
    """Compute forward strand coordinates of zebrafish. Returns coord dict."""
    if coord['dir'] == '-':
        length, start, end = zebrafish_lengths[coord['chrom']], coord['start'], coord['end']
        coord['start'] = length - end
        coord['end'] = length - start
    return coord

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Append fields with positive strand coordinates.')
    parser.add_argument("-f", type=int, default=1, help='Field number of the strand coordinates.')
    OPTS = parser.parse_args()
    line_tups = read_fields()
    for l in line_tups:
        coord = forward_strand_zebrafish(parse_coords(l[OPTS.f-1]))
        sys.stdout.write('\t'.join(l + [COORD_FORMATS['forward'].format(**coord)]) + '\n')


if __name__ == '__main__':
    main()
