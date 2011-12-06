#!/usr/bin/env python
from collections import defaultdict
import sys
import matplotlib.pyplot as plt
from util import read_fields, parse_coords
from negToPos import forward_strand_zebrafish

def parse_dat(line_tups):
    cne_dict = defaultdict(list)
    for cne, score in line_tups:
        cne_dict[cne].append(score)
    return cne_dict

def plot_cne(name, scores, valid):
    plt.plot(scores)
    p = plt.axvspan(valid[0], valid[1], facecolor='r', alpha=0.4)
    plt.xlabel('alignment')
    plt.ylabel('score')
    plt.title('%s' % name)
    plt.show()

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help='Input file e.g. d2z.dat.')
    parser.add_argument('--extra_data', type=argparse.FileType('r'), default='hg18.toDanRer5.seqs.txt',
                        help='Extra data file e.g. hg18.toDanRer5.seqs.txt.')
    parser.add_argument('--valid', type=int, default=9, help='Field number of valid test coordinates.')
    OPTS = parser.parse_args()
    cne_dict = parse_dat(read_fields(f=OPTS.infile))
    line_tups = read_fields(f=OPTS.extra_data)
    cne_valids = {}
    for l in line_tups:
        danrer_co, valid_co = parse_coords(l[4]), forward_strand_zebrafish(parse_coords(l[OPTS.valid-1]))
        valid_indices = (valid_co['start'] - danrer_co['start'], valid_co['end'] - danrer_co['start'])
        cne_valids[l[0]] = valid_indices
    one = cne_dict.keys()[0]
    plot_cne(one, cne_dict[one], cne_valids[one])


if __name__ == '__main__':
    main()
