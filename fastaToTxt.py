#!/usr/bin/env python
import sys
from Bio import SeqIO

for seq_record in SeqIO.parse(sys.stdin, "fasta"):
    sys.stdout.write('%s\n' % seq_record.seq)
