#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt
from util import read_fields, parse_coords, parse_dat
from peakdetect import peakdetect, smooth

def plot_cne(name, scores, valid):
    lookahead = len(scores) / 50
    scores = np.array([float(x) for x in scores])
    signal = smooth(scores, window_len=40, window='bartlett')
    maxima = peakdetect(signal, look5ahead=lookahead)[0]
    m_x = np.array([m[0] for m in maxima])
    m_y = np.array([m[1] for m in maxima])
    plt.plot(range(len(signal)), signal, 'k', m_x, m_y, 'bo')
    plt.axvspan(valid[0], valid[1], facecolor='r', alpha=0.4)
    plt.xlabel('alignment')
    plt.ylabel('score')
    plt.title(name)
    fig = plt.gcf()
    fig.set_size_inches(8, 5)
    plt.savefig(name + '.png', dpi=140)

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('identifier', help='Identifier name e.g. cne.100899.FST.')
    parser.add_argument('-f', '--file', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help='Input file e.g. d2z.dat.')
    parser.add_argument('--extra_data', type=argparse.FileType('r'), default='hg18.toDanRer5.seqs.txt',
                        help='Extra data file e.g. hg18.toDanRer5.seqs.txt.')
    parser.add_argument('--valid', type=int, default=9, help='Field number of valid test coordinates.')
    OPTS = parser.parse_args()
    cne_dict = parse_dat(read_fields(f=OPTS.file))
    line_tups = read_fields(f=OPTS.extra_data)
    cne_valids = {}
    for l in line_tups:
        danrer_co, valid_co = parse_coords(l[4]), parse_coords(l[OPTS.valid-1])
        valid_indices = (valid_co['start'] - danrer_co['start'], valid_co['end'] - danrer_co['start'])
        cne_valids[l[0]] = valid_indices
    plot_cne(OPTS.identifier, cne_dict[OPTS.identifier], cne_valids[OPTS.identifier])


if __name__ == '__main__':
    main()
