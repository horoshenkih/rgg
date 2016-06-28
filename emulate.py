#!/usr/bin/env python

import sys
from argparse import ArgumentParser
from collections import Counter

import matplotlib.pyplot as plt

from graph import HypRG

def main():
    parser = ArgumentParser()
    parser.add_argument('n', type=int, help='number of vertices')
    parser.add_argument('--alpha', type=float, help='alpha', default=1.)
    parser.add_argument('-C', type=float, help='C', default=0.)
    parser.add_argument('-f', help='outfile')

    args = parser.parse_args()

    if args.f:
        out_f = open(args.f, 'w')
    else:
        out_f = sys.stdout

    n = args.n
    alpha = args.alpha
    C = args.C
    g = HypRG(n,alpha=alpha,C=C)

    for e in g.edges():
        e_fmt = []
        for v in e:
            e_fmt.append("{0:.5f},{1:.5f}".format(*v))
        out_f.write(' '.join(e_fmt) + '\n')

if __name__ == '__main__':
    main()
