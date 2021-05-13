#!/usr/bin/env python3

#blah.py[directory/ containing mouse cluster .bed files] [peaks from dataset for each cluster]

import sys
from scipy.stats import hypergeom
from glob import glob
from os.path import basename

bedGlob = glob(sys.argv[1] + "*.bed")
sizeMap = {} #Sizes of each mouse cluster
for fileName in bedGlob:
    name = basename(fileName)[:-4]
    bedFile = open(fileName, "r")
    lenFile = len([True for l in bedFile.readlines() if len(l.strip()) > 0])
    bedFile.close()
    sizeMap[name] = lenFile

nMousePeaks = sum([sizeMap[x] for x in sizeMap]) #All peaks in mouse clusters

peakCountMap = {} #Peaks from the neuron type in each cluster
nPeaksInClusters = 0 # All peaks from neuron type in a mouse cluster

inFile = open(sys.argv[2], "r")
for line in inFile:
    tokens = line.strip().split()
    if len(tokens) == 2 and tokens[1] in sizeMap:
        nPeaksInClusters += int(tokens[0])
        peakCountMap[tokens[1]] = int(tokens[0])
inFile.close()

for cls in sizeMap:
    if cls not in peakCountMap:
        peakCountMap[cls] = 0

print("Reading peak counts from", sys.argv[2])
print("Found", len(sizeMap), "clusters containing", nMousePeaks, "enhancers with mouse orthologs")
print("Found", nPeaksInClusters, "peaks in a mouse-active cluster")
print("Cluster\tp-value\tPeaks in Cluster\tCluster Size")
for cls in sizeMap:
    prb = hypergeom.sf(peakCountMap[cls], nMousePeaks, nPeaksInClusters, sizeMap[cls])
    print(cls + "\t" + str(prb) + "\t" + str(peakCountMap[cls]) + "\t" + str(sizeMap[cls]))
print()

        
