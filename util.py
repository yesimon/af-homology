#!/usr/bin/env python
import re
import sys
from operator import itemgetter

dense_regex = re.compile(r'(?P<chrom>\w+):(?P<start>\d+)-(?P<end>\d+),(?P<dir>[+-])')
forward_regex = re.compile(r'(?P<chrom>\w+):(?P<start>\d+)-(?P<end>\d+)')
bed_regex = re.compile(r'(?P<chrom>\w+)\t(?P<start>\d+)\t(?P<end>\d+)\t(?P<name>\S+)\t(?P<score>\S+)\w(?P<dir>[+-])')
bed_short_regex = re.compile(r'(?P<chrom>\w+)\t(?P<start>\d+)\t(?P<end>\d+)')

COORD_FORMATS = {
    'forward': '{chrom}:{start}-{end}',
    'dense': '{chrom}:{start}-{end},{dir}',
    'bed_long': '{chrom}\t{start}\t{end}\t{name}\t{score}\t{dir}',
    'bed': '{chrom}\t{start}\t{end}',
}

def set_missing(d, k, v):
    if not k in d: d[k] = v
    return d

def read_fields(f=sys.stdin, sep=None):
    return [l.strip().split(sep) for l in f.read().splitlines() if l]

def parse_coords(co):
    """Given a strand coordinate like 'chrX:1241-1888,+' or bed format return a
    dictionary of relevant coordinates. Raises exception on error.
    """
    if not isinstance(co, basestring):
        co = '\t'.join(co[:6])
    mo = None
    if not mo: mo = dense_regex.match(co)
    if not mo: mo = forward_regex.match(co)
    if not mo: mo = bed_regex.match(co)
    if not mo: mo = bed_short_regex.search(co)
    if not mo: raise Exception
    return set_missing(mo.groupdict(), 'dir', '+')

def index_coords(co, index, l=None):
    """Returns modified coordinates with forward index and window length l
    accounting for strand direction."""
    if co['dir'] == '+':
        co['start'] += index
        co['end'] += index
    elif co['dir'] == '-':
        co['start'] -= index
        co['end'] -= index
    if l: co['end'] = co['start'] + l
    return co

class AFModel(object):
    def __init__(self, l=None):
        """Initialize the model.

        Parameters:
          l (int): Sliding window length.
        """
        self.l = l

    def fit(self, X):
        """Fit the model against a list X of sequence strings.

        Parameters:
          X (list): List of sequence strings.

        Returns:
          (instance) self
        """
        return self

    def score(self, X):
        """Score a list X of sequence strings directly against the model.

        Parameters:
          X (list): List of sequence strings.

        Returns:
          (list) Float scores corresponding to input strings.
        """
        return

    def scan(self, X, n=20):
        """Scan a list X of sequence strings with a sliding window.

        Parameters:
          X (list): List of sequence strings.
          n (int): Number of top hits.

        Returns:
          (list) List of (list) of n best (index, score) tuples. If n == 1,
          return list of tuples directly.
        """
        results = []
        for x in X:
            hits = list(enumerate(self.score([x[i:i+self.l] for i in range(len(x)-self.l)])))
            hits.sort(key=itemgetter(1))
            results.append(hits[:n])
        return results
