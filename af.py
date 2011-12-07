#!/usr/bin/env python
import sys
import cPickle as pickle
from d2z import add_d2z_arguments, D2z
from hexmcd import add_hexmcd_arguments, HexMCD
from util import read_fields, parse_coords, progress

def row_search(OPTS, m, line_tups):
    """Writes the scan results to a pickle file.

    The pickled object is a  dictionary with cne names as keys and list as value.
    Each list contains tuples of (coord, index) where coord is a coord dict
    and index is int index starting from 0.
    """
    rows = {}
    for i, l in enumerate(line_tups):
        name, a, b, c, = l[0], l[OPTS.a-1], l[OPTS.b-1], parse_coords(l[OPTS.c-1])
        m.fit([a])
        hits = m.scan([b], n=None, reverse_complement=True, sort=False)[0]
        scores = [score for index, score in hits]
        rows[name] = scores
        progress(50, (float(i)+1)/len(line_tups)*100, pre="Processing genes")
    pkl_filename = OPTS.model + '.pkl'
    sys.stdout.write("Writing pickle file %s. Do not exit!\n" % pkl_filename)
    pkl = open(OPTS.model + '.pkl', 'wb')
    pickle.dump(rows, pkl, pickle.HIGHEST_PROTOCOL)

def main():
    import argparse
    parser = argparse.ArgumentParser('Harness for alignment free homology.', add_help=False,
                                     epilog='')
    parser.add_argument('--input', type=argparse.FileType('rb'), default=sys.stdin,
                        help='Input file e.g. hg18.toDanRer5.seqs.txt.')
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
    line_tups = read_fields(f=OPTS.input)

    if OPTS.model == 'd2z':
        m = D2z(**vars(OPTS))
    if OPTS.model == 'hexmcd':
        m = HexMCD(bg_list=[l[6] for l in line_tups], smoothing='ones', **vars(OPTS))
    row_search(OPTS, m, line_tups)


if __name__ == '__main__':
    main()
