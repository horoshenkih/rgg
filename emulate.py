from argparse import ArgumentParser
from collections import Counter

import matplotlib.pyplot as plt

from graph import RGG, HypRG

def main():
    parser = ArgumentParser()
    parser.add_argument('-n', type=int, help='number of vertices', default=100)
    parser.add_argument('--alpha', type=float, help='alpha', default=1.)
    parser.add_argument('-C', type=float, help='C', default=0.)
    if False:
        rgg = RGG(
            t=100,
            #cluster_attachment='uniform'
        )
        print rgg.clusters
        print sorted(rgg.degrees())

    args = parser.parse_args()
    g = HypRG(args.n,alpha=args.alpha)
    v = g.vertices()
    e = g.edges()
    # vertices
    print(len(v))
    r, phi = zip(*v)
    ax = plt.subplot(111, projection='polar')
    ax.plot(phi, r, marker='o', linewidth=0)
    # edges
    print(len(e))
    for (v1, v2) in e:
        r1, phi1 = v1
        r2, phi2 = v2
        ax.plot((phi1, phi2), (r1, r2), color='g')

    plt.show()
    #deg = g.degrees()
    #print Counter(deg)

if __name__ == '__main__':
    main()
