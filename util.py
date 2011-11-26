#!/usr/bin/env python
import sys

def read_fields(f=sys.stdin, sep=None):
    return [l.strip().split(sep) for l in f.read().split('\n') if l]
