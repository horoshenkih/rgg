from argparse import ArgumentParser
from collections import Counter

import matplotlib.pyplot as plt

from graph import HypRG

def main():
    parser = ArgumentParser()
    parser.add_argument('n', type=int, help='number of vertices')
    parser.add_argument('--alpha', type=float, help='alpha', default=1.)
    parser.add_argument('-C', type=float, help='C', default=0.)
    parser.add_argument('--pg', action='store_true', help='plot graph')
    parser.add_argument('--pd', action='store_true', help='plot degree distribution')

    # parse options
    args = parser.parse_args()
    n = args.n
    alpha = args.alpha

    # construct graph
    C = args.C
    g = HypRG(n,alpha=alpha,C=C)
    v = g.vertices()
    e = g.edges()

    print "number of vertices: {}".format(len(v))
    print "number of edges: {}".format(len(e))

    # determine plots
    plots = [
        ['graph', args.pg],
        ['degrees', args.pd],
    ]
    active_plots = [pl for pl in plots if pl[1] is True]
    plot2idx = {pl[0]: pl_i+1 for pl_i, pl in enumerate(active_plots)}
    n_plots = len(plot2idx)

    # draw graph
    if args.pg:
        plot_idx = plot2idx['graph']
        # vertices
        r, phi = zip(*v)
        ax = plt.subplot(1, n_plots, plot_idx, projection='polar')
        ax.plot(phi, r, marker='o', linewidth=0)

        # edges
        for (v1, v2) in e:
            r1, phi1 = v1
            r2, phi2 = v2
            ax.plot((phi1, phi2), (r1, r2), color='g')

    # analyze degree distribution
    deg = g.degrees()
    n_isolated_vertices = len([d for d in deg if d == 0])
    deg_dist = [i for i in Counter(deg).iteritems() if i[0] > 0]
    degs, counts = zip(*deg_dist)
    freqs = [float(c) / n for c in counts]

    print "n isolated vertices: {}".format(n_isolated_vertices)
    components = g.components()
    print "n components: {}".format(len(components))
    print "n nontrivial components: {}".format(len(components) - n_isolated_vertices)
    print "largest compoment size: {}".format(len(components[0]))

    # plot degree distribution
    if args.pd:
        plot_idx = plot2idx['degrees']
        ax_deg = plt.subplot(1, n_plots, plot_idx)
        ax_deg.loglog(degs, freqs, marker='o', linewidth=0)
        ax_deg.grid(True)

    if n_plots > 0:
        plt.show()

if __name__ == '__main__':
    main()
