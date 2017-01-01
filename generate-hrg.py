#!/usr/bin/python

import sys
from argparse import ArgumentParser
from collections import Counter

import matplotlib.pyplot as plt

from lib.graph import HypRG

def main():
    parser = ArgumentParser(description='Emulate hyperbolic random graph', epilog='Returns only edges (without isolated vertices!)')
    parser.add_argument('n', type=int, help='number of vertices')
    parser.add_argument('--alpha', type=float, help='alpha', default=1.)
    parser.add_argument('-C', type=float, help='C', default=0.)
    parser.add_argument('-f', help='outfile')
    parser.add_argument('-s', '--seed', help='random seed', type=int)

    args = parser.parse_args()

    if args.f:
        out_f = open(args.f, 'w')
    else:
        out_f = sys.stdout

    n = args.n
    alpha = args.alpha
    C = args.C
    seed = 0 if args.seed is None else args.seed
    g = HypRG(n, alpha=alpha, C=C, seed=seed)

    for e in g.edges():
        e_fmt = []
        for v in e:
            e_fmt.append("{0:.3f},{1:.3f}".format(*v))
        out_f.write(' '.join(e_fmt) + '\n')

if __name__ == '__main__':
    main()
