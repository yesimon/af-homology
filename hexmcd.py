#!/usr/bin/env python
import argparse
import math
import sys
from collections import defaultdict
from itertools import product
import numpy as np
from util import AFModel, shuffle_string

def generate_hex_map():
    hex_map = defaultdict(float)
    A = 'ACTG'
    seqs = [''.join(seq) for seq in product(A, repeat=5)]
    hex_map.update(dict([(seq, i) for i, seq in enumerate(seqs)]))
    return hex_map

def generate_hex_fuzzy(hex_map):
    A = set(['A', 'C', 'T', 'G'])
    hex_fuzzy = defaultdict(list)
    for seq in hex_map.iterkeys():
        for i in range(len(seq)):
            hex_fuzzy[seq].extend([seq[:i] + l + seq[i+1:] for l in A if l != seq[i]])
    return hex_fuzzy

HEX_MAP = generate_hex_map()
HEX_FUZZY = generate_hex_fuzzy(HEX_MAP)

def transition_matrix(seqs, w=0.25, k=5, smoothing=None):
    """Generate numpy transition matrix."""
    if smoothing == 'ones':
        tm = np.ones([len(HEX_MAP), len(HEX_MAP)])
    else:
        tm = np.zeros([len(HEX_MAP), len(HEX_MAP)])
    for seq in seqs:
        for i in range(len(seq)-k-1):
            s, t = seq[i:i+k], seq[i+1:i+k+1]
            tm[HEX_MAP[s], HEX_MAP[t]] += 1
            for s_prime in HEX_FUZZY[s]:
                tm[HEX_MAP[s_prime], HEX_MAP[t]] += w
    return tm

class HexMCD(AFModel):
    def __init__(self, bg_file=None, bg_list=None, smoothing=None, *args, **kwargs):
        """Initialize background markov model from bg."""
        assert bg_file or bg_list
        super(HexMCD, self).__init__(*args, **kwargs)
        self.w = 0.25
        self.k = 5
        self.smoothing = smoothing
        if bg_file:
            raise Exception
        if bg_list:
            shuffled = [shuffle_string(s) for s in bg_list]
            self.bg_tm = transition_matrix(shuffled, w=self.w, k=self.k,
                                           smoothing=self.smoothing)
        self.bg_tm_rowsum = np.sum(self.bg_tm, axis=1)

    def fit(self, X):
        super(HexMCD, self).fit(X)
        self.tm = transition_matrix(X, w=self.w, k=self.k, smoothing=self.smoothing)
        self.tm_rowsum = np.sum(self.tm, axis=1)
        return self

    def prob(self, seq):
        S = 0
        for i in range(len(seq)-self.k-1):
            s, t = seq[i:i+self.k], seq[i+1:i+self.k+1]
            a_plus = self.tm[HEX_MAP[s], HEX_MAP[t]] / self.tm_rowsum[HEX_MAP[s]]
            a_minus = self.bg_tm[HEX_MAP[s], HEX_MAP[t]] / self.bg_tm_rowsum[HEX_MAP[s]]
            if a_minus and a_plus: S += math.log(a_plus/a_minus)
        return S

    def score(self, X):
        return [self.prob(x) for x in X]

def add_hexmcd_arguments(parser, main=False):
    parser.add_argument('--bg', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    if not main:
        return parser
    raise Exception
