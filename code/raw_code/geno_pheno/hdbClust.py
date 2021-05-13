#!/usr/bin/env python3
#hdbClust.py [distance matrix] [output file (optional)]

import numpy, hdbscan, sys

distance_matrix = numpy.loadtxt(open(sys.argv[1], "rb"), delimiter=",", skiprows=0)

clusterer = hdbscan.HDBSCAN(metric='precomputed', min_cluster_size=20)
clusterer.fit(distance_matrix)
print(clusterer.labels_)
