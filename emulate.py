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
    C = args.C
    g = HypRG(n,alpha=alpha,C=C)
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
    n_isolated_vertices = len([d for d in deg if d == 0])
    deg_dist = [i for i in Counter(deg).iteritems() if i[0] > 0]
    degs, counts = zip(*deg_dist)
    freqs = [float(c) / n for c in counts]
    ax_deg = plt.subplot(122)
    ax_deg.loglog(degs, freqs, marker='o', linewidth=0)
    ax_deg.grid(True)

    print "n isolated vertices: {}".format(n_isolated_vertices)
    components = g.components()
    print "n components: {}".format(len(components))
    print "largest compoment size: {}".format(len(components[0]))
    plt.show()

if __name__ == '__main__':
    main()
