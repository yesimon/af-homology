#!/usr/bin/env python
import sys
from d2z import add_d2z_arguments, D2z
from hexmcd import add_hexmcd_arguments, HexMCD
from util import read_fields, index_coords, parse_coords, overlap_coords, COORD_FORMATS
from negToPos import forward_strand_zebrafish

def row_search(OPTS, m, line_tups):
    for l in line_tups:
        a, b, c, = l[OPTS.a-1], l[OPTS.b-1], forward_strand_zebrafish(parse_coords(l[OPTS.c-1]))
        correct_coord = parse_coords(l[OPTS.valid-1])
        m.fit([a])
        hits = m.scan([b], n=None, reverse_complement=True)[0]
        coords = [(index_coords(c, index, l=m.l), score) for index, score in hits]
        # hit_index = 0
        # valid_coord = None
        # hit = None
        # for coord in coords:
        #     if overlap_coords(coord, correct_coord):
        #         hit = hit_index
        #         valid_coord = c
        #         break
        #     hit_index += 1
        # if hit == None:
        #     continue
        for coord, score in coords:
            stats = '\t'.join([str(score)])
            coord_str = COORD_FORMATS['forward'].format(**coord)
            sys.stdout.write('\t'.join([l[0], coord_str, stats]) + '\n')

def main():
    import argparse
    parser = argparse.ArgumentParser('Harness for alignment free homology.', add_help=False)
    parser.add_argument('-a', type=int, default=6, help='Field number of A (training) sequences.')
    parser.add_argument('-b', type=int, default=8, help='Field number of B (test) sequences.')
    parser.add_argument('-c', type=int, default=5, help='Field number of test coordinates.')
    parser.add_argument('--valid', type=int, default=9, help='Field number of valid test coordinates.')
    parser.add_argument('-l', type=int, default=None, help='Length of scanning window. Defaults to the average of training sequences.')

    subparsers = parser.add_subparsers(help='Model algorithm to use.', dest='model')
    d2z_parser = subparsers.add_parser('d2z', help='D2z scoring metric.')
    d2z_parser = add_d2z_arguments(d2z_parser)
    hexmcd_parser = subparsers.add_parser('hexmcd', help='HexMCD algorithm.')
    hexmcd_parser = add_hexmcd_arguments(hexmcd_parser)
    # Add more parsers here.

    OPTS = parser.parse_args()
    line_tups = read_fields()

    if OPTS.model == 'd2z':
        m = D2z(**vars(OPTS))
    if OPTS.model == 'hexmcd':
        m = HexMCD(bg_list=[l[6] for l in line_tups], **vars(OPTS))
    row_search(OPTS, m, line_tups)
    return
    a_seqs = [l[OPTS.a-1] for l in line_tups]
    m.fit(a_seqs)
    for l in line_tups:
        b, c, correct_coord = l[OPTS.b-1], parse_coords(l[OPTS.c-1]), parse_coords(l[OPTS.valid-1])
        hits = m.scan([b], n=None, reverse_complement=True)[0]
        coords = [forward_strand_zebrafish(index_coords(c, index, l=m.l)) for index, score in hits]
        hit_index = 0
        valid_coord = None
        hit = None
        for coord in coords:
            if overlap_coords(coord, correct_coord):
                hit = hit_index
                valid_coord = c
                break
            hit_index += 1
        if hit == None:
            continue
        stats = '\t'.join([str(hit), str(len(b))])
        coord_str = COORD_FORMATS['dense'].format(**valid_coord)
        sys.stdout.write('\t'.join([l[0], coord_str, stats]) + '\n')

if __name__ == '__main__':
    main()
