#!/usr/bin/env python
import sys
from util import read_fields

line_tups = read_fields(f=sys.stdin)
for l in line_tups:
    sys.stdout.write('%s:%s-%s\n' %(l[0], l[1], l[2]))
