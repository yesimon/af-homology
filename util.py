#!/usr/bin/env python
from collections import defaultdict
import math
import re
import sys
from operator import itemgetter
from random import shuffle
from Bio.Seq import Seq

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

def shuffle_string(seq):
    seq_list = list(seq)
    shuffle(seq_list)
    return ''.join(seq_list)

def read_fields(f=sys.stdin, sep=None):
    return [l.strip().split(sep) for l in f.read().splitlines() if l]

def parse_dat(line_tups):
    cne_dict = defaultdict(list)
    for cne, score in line_tups:
        cne_dict[cne].append(score)
    return cne_dict

def parse_extra_data(line_tups):
    extra = {}
    for l in line_tups:
      extra[l[0]] = {
        'name': l[0],
        'dr_co': parse_coords(l[4]),
        'dr_valid_co': parse_coords(l[8]),
        'dr_seq': l[7],
        'hg_co': parse_coords('\t'.join(l[1:4])),
        'hg_seq': l[5],
      }
    return extra

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
    d = mo.groupdict()
    d['start'] = int(d['start'])
    d['end'] = int(d['end'])
    return d

def index_coords(co, index, l=None):
    """Returns modified coordinates with forward index and window length l.
    Works for native + and - coordinates. Do not use with - coordinates
    coverted to + beause they will be backwards."""
    c = co.copy()
    c['start'] += index
    c['end'] += index
    if l: c['end'] = c['start'] + l
    return c

def overlap_coords(c1, c2):
    """Consider whether c1 overlaps c2."""
    within_c2 = lambda c: c2['start'] < c < c2['end']
    return within_c2(c1['start']) or within_c2(c1['end'])

def progress(width, percent, pre=None):
    # width defines bar width
    # percent defines current percentage
    marks = math.floor(width * (percent / 100.0))
    spaces = math.floor(width - marks)
    loader = '[' + ('=' * int(marks)) + (' ' * int(spaces)) + ']'
    if pre: loader = pre + ' ' + loader
    sys.stdout.write("%s %d%%\t\t\r" % (loader, percent))
    if percent >= 100:
        sys.stdout.write("\n")
    sys.stdout.flush()

class AFModel(object):
    def __init__(self, l=None, *args, **kwargs):
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
        self.l = self.l or int(float(len(''.join(X)))/len(X))
        return self

    def score(self, X):
        """Score a list X of sequence strings directly against the model.

        Parameters:
          X (list): List of sequence strings.

        Returns:
          (list) Float scores corresponding to input strings.
        """
        return

    def scan(self, X, n=20, reverse_complement=True, sort=True):
        """Scan a list X of sequence strings with a sliding window.

        Parameters:
          X (list): List of sequence strings.
          n (int): Number of hits. None to return all hits.
          reverse_complement (bool): Also scan reverse complement of sequence.
          sort (bool): If true, sort the sequences by highest score.

        Returns:
          (list) List of (list) of n best (index, score) tuples. If n == 1,
          return list of tuples directly.
        """
        results = []
        for x in X:
            hits = list(enumerate(self.score([x[i:i+self.l] for i in range(len(x)-self.l+1)])))
            x_len = len(x)
            x_rc = Seq(x).reverse_complement().tostring()
            hits_rc = list(enumerate(self.score([x_rc[i:i+self.l] for i in range(len(x)-self.l+1)])))
            hits_rc = [(x_len - i - self.l, score) for i, score in hits_rc]
            hits_d = dict(hits)
            hits_rc_d = dict(hits_rc)
            for i, score in hits_d.iteritems():
                if hits_rc_d[i] > score:
                    hits_d[i] = hits_rc_d[i]
            hits = list(hits_d.items())
            if sort: hits.sort(key=itemgetter(1), reverse=True)
            if n == 1: results.append(hits[0])
            elif n == None: results.append(hits)
            else: results.append(hits[:n])
        return results
