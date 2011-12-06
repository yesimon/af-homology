#!/usr/bin/env python
from collections import defaultdict
from operator import itemgetter
import sys
import matplotlib.pyplot as plt
import numpy as np
from util import read_fields, parse_coords, parse_dat, overlap_coords, index_coords
from peakdetect import peakdetect, smooth

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

def overlap(cne_dict, extra):
  window_len = sum([len(d['hg_seq']) for d in extra.itervalues()])/float(len(extra))
  results = defaultdict(dict)
  for cne, scores in cne_dict.iteritems():
    scores = sorted(enumerate(scores), key=itemgetter(1), reverse=True)
    for hit_rank, (i, score) in enumerate(scores):
      if overlap_coords(index_coords(extra[cne]['dr_co'], i, l=window_len),
                        extra[cne]['dr_valid_co']):
        results[cne]['rank'] = hit_rank + 1
        results[cne]['sequence_length'] = len(scores)
        break
  return results

def ranked_peaks(cne_dict, extra):
  def dist(cne, i):
    v_avg = (extra[cne]['dr_valid_co']['start'] + extra[cne]['dr_valid_co']['end']) / 2.0
    v_i = v_avg - extra[cne]['dr_co']['start']
    return abs(i - v_i)
  results = defaultdict(dict)
  for cne, scores in cne_dict.iteritems():
    scores = np.array([float(x) for x in scores])
    signal = smooth(scores, window_len=40, window='bartlett')
    maxima = peakdetect(signal, lookahead=500)[0]
    maxima.sort(key=itemgetter(1), reverse=True)
    if not maxima:
      print cne
      continue
    hit_rank = np.argmin(np.array([dist(cne, i) for i, score in maxima]))
    results[cne]['rank'] = hit_rank + 1
    results[cne]['n_peaks'] = len(maxima)
  return results

def main():
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument('--infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                      help='Input file e.g. d2z.dat.')
  parser.add_argument('--extra_data', type=argparse.FileType('r'), default='hg18.toDanRer5.seqs.txt',
                    help='Extra data file e.g. hg18.toDanRer5.seqs.txt.')
  parser.add_argument('scoring', choices=['ranked_peaks', 'overlap'])
  OPTS = parser.parse_args()
  line_tups = read_fields(f=OPTS.infile)
  cne_dict = parse_dat(line_tups)
  extra = parse_extra_data(read_fields(f=OPTS.extra_data))
  if OPTS.scoring == 'ranked_peaks':
    results = ranked_peaks(cne_dict, extra)
    for cne, result in sorted(results.iteritems()):
      sys.stdout.write('\t'.join([cne, str(result['rank']),
                                  str(result['n_peaks'])]) + '\n')
  elif OPTS.scoring == 'overlap':
    results = overlap(cne_dict, extra)
    for cne, result in sorted(results.iteritems()):
      sys.stdout.write('\t'.join([cne, str(result['rank']),
                                  str(result['sequence_length'])]) + '\n')


if __name__ == '__main__':
  main()
