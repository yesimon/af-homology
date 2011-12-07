#!/usr/bin/env python
import cPickle as pickle
import sys

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', type=argparse.FileType('rb'), default=sys.stdin)
    OPTS = parser.parse_args()
    cne_dict = pickle.load(OPTS.infile)

    for cne, results in cne_dict.iteritems():
        for score in results:
            sys.stdout.write('%s\t%s\n' % (cne, score))

if __name__ == '__main__':
    main()
