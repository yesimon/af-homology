#!/usr/bin/env python
import argparse
import sys
from collections import defaultdict
from itertools import product
import numpy as np
from util import AFModel

def generate_hex_map():
    A = 'ACTG'
    seqs = [''.join(seq) for seq in product(A, repeat=6)]
    return dict([(seq, i) for i, seq in enumerate(seqs)])

def generate_hex_fuzzy(hex_map):
    A = set(['A', 'C', 'T', 'G'])
    hex_fuzzy = defaultdict(list)
    for seq in hex_map.iterkeys():
        for i in range(len(seq)):
            hex_fuzzy[seq].extend([seq[:i] + l + seq[i+1:] for l in A if l != seq[i]])
    return hex_fuzzy

HEX_MAP = generate_hex_map()
HEX_FUZZY = generate_hex_fuzzy(HEX_MAP)

def transition_matrix(seqs, w=0.25, k=6):
    """Generate numpy transition matrix."""
    tm = np.zeros([len(HEX_MAP), len(HEX_MAP)])
    for seq in seqs:
        for i in range(len(seq)-k-1):
            s, t = seq[i:i+k], seq[i+1:i+k+1]
            tm[HEX_MAP[s], HEX_MAP[t]] += 1
            for s_prime in HEX_FUZZY[s]:
                tm[HEX_MAP[s_prime], HEX_MAP[t]] += w
    return tm

class HexMCD(AFModel):
    def __init__(self, bg_file=None, bg_list=None):
        """Initialize background markov model from bg."""
        assert bg_file or bg_list
        self.w = 0.25
        self.k = 6
        if bg_file:
            self.bg_tm = None
        if bg_list:
            self.bg_tm = transition_matrix(bg_list, w=self.w, k=self.k)

    def fit(self, X):
        self.tm = transition_matrix(X, w=self.w, k=self.k)
        return self

    def prob(self, seq):
        for i in range(len(seq)-self.k-1):
            s, t = seq[i:i+self.k], seq[i+1:i+self.k+1]
            #self.tm[HEX_MAP[s], HEX_MAP[t]] / 

    def score(self, X):
        for x in X:



def add_hexmcd_arguments(parser, main=False):
    parser.add_argument('--bg', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    if not main:
        return parser
    return parser
