#!/usr/bin/env python
"""
D2z reference here: http://bioinformatics.oxfordjournals.org/content/23/13/i249.full
"""
import sys
from collections import Counter
from util import read_fields

def k_word_counts(seq, k):
    return Counter([seq[i:i+k] for i in range(len(seq)-k)])

def letter_prob(seq):
    return dict([(letter, float(count)/len(seq)) for letter, count in Counter(seq).iteritems()])

def d2(A, B, k):
    N_A, N_B = k_word_counts(A, k), k_word_counts(B, k)
    return sum([N_A[a]*N_B[a] for a in N_A if a in N_B])

def g(f_a, f_b, x, y):
    return sum([f_a[a]**x * f_b[a]**y for a in f_a if a in f_b])

def E_d2(A, B, k, f_a, f_b):
    return (len(A)-k+1) * (len(B)-k+1) * g(f_a, f_b, 1, 1) ** k

def V_d2(A, B, k, f_a, f_b):
    nbar1, nbar2 = len(A) - k + 1, len(B) - k + 1
    qbar1, qbar2 = len(A) - 2 * k + 2, len(B) - 2 * k + 2
    p2 = sum([f_a[a] * f_b[a] for a in f_a if a in f_b])
    p31 = sum([f_a[a] * f_b[a] * f_a[a] for a in f_a if a in f_b])
    p32 = sum([f_a[a] * f_b[a] * f_b[a] for a in f_a if a in f_b])
    # Crabgrass with l=0 (complete overlap)
    pow1 = (pow(p32, k) - pow(p2, 2*k)) * nbar1 * qbar2 * (qbar2 - 1)
    pow2 = (pow(p31, k) - pow(p2, 2*k)) * nbar2 * qbar1 * (qbar1 - 1)
    variance = pow1 + pow2
    # Crabgrass with l>0
    for l in range(1, k):
        variance += 2*(nbar1-l)*qbar2*(qbar2-1)*(pow(p2,2*l)*pow(p32,k-l)-pow(p2,2*k))+2*(nbar2-l)*qbar1*(qbar1-1)* \
                    (pow(p2,2*l)*pow(p31,k-l)-pow(p2,2*k))
    # Accordian main diagonal
    variance += nbar1*nbar2*(pow(p2,k)-pow(p2,2*k)) # l=0 term
    for l in range(1, k):
        variance += 2*(nbar1-l)*(nbar2-l)*(pow(p2,k+l)-pow(p2,2*k)) # l>0 terms
    return variance

def d2z(A, B, k):
    f_a, f_b = letter_prob(A), letter_prob(B)
    return (d2(A, B, k) - E_d2(A, B, k, f_a, f_b))/(V_d2(A, B, k, f_a, f_b) ** 0.5)

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Compute d2z scores.')
    parser.add_argument("-a", type=int, default=1, help='Field number of A (training) sequences.')
    parser.add_argument("-b", type=int, default=2, help='Field number of B (test) sequences.')
    parser.add_argument("-k", type=int, default=4, help='k-mer lengths.')
    parser.add_argument("-l", type=int, default=None, help='Length of scanning window. Defaults to the average of training sequences.')
    OPTS = parser.parse_args()
    line_tups = read_fields()
    a_seqs = [l[OPTS.a-1] for l in line_tups]
    a_seq = ''.join(a_seqs)
    window = float(len(a_seq))/len(a_seqs)
    for l in line_tups:
        b_seq = l[OPTS.b-1]
        sys.stdout.write('%s\n' % d2z(a_seq, b_seq, OPTS.k))


if __name__ == '__main__':
    main()
