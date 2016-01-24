from graph import RGG

rgg = RGG(
    t=100,
    #cluster_attachment='uniform'
)
print rgg.clusters
print sorted(rgg.degrees())
