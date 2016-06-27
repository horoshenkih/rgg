from graph import RGG, HypRG

if False:
    rgg = RGG(
        t=100,
        #cluster_attachment='uniform'
    )
    print rgg.clusters
    print sorted(rgg.degrees())

g = HypRG(10000)
print(len(g.V))
