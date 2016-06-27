from argparse import ArgumentParser
from collections import Counter

import matplotlib.pyplot as plt

from graph import HypRG

def main():
    parser = ArgumentParser()
    parser.add_argument('-n', type=int, help='number of vertices', default=100)
    parser.add_argument('--alpha', type=float, help='alpha', default=1.)
    parser.add_argument('-C', type=float, help='C', default=0.)

    args = parser.parse_args()
    n = args.n
    alpha = args.alpha
    g = HypRG(n,alpha=alpha)
    v = g.vertices()
    e = g.edges()

    # vertices
    print "number of vertices: {}".format(len(v))
    r, phi = zip(*v)
    ax = plt.subplot(121, projection='polar')
    ax.plot(phi, r, marker='o', linewidth=0)

    # edges
    print "number of edges: {}".format(len(e))
    for (v1, v2) in e:
        r1, phi1 = v1
        r2, phi2 = v2
        ax.plot((phi1, phi2), (r1, r2), color='g')

    deg = g.degrees()
    n_isolated_vertices = deg[0]
    deg_dist = [i for i in Counter(deg).iteritems() if i[0] > 0]
    degs, counts = zip(*deg_dist)
    freqs = [float(c) / n for c in counts]
    ax_deg = plt.subplot(122)
    ax_deg.set_xscale("log")
    ax_deg.set_yscale("log")
    ax_deg.errorbar(degs, freqs)

    print "size of giant component: {}".format(n - n_isolated_vertices)
    plt.show()

if __name__ == '__main__':
    main()
