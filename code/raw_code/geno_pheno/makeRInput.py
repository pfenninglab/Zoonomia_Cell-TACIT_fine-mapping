#!/usr/bin/env python3
#makeRInput.py [predictions file] [clusters file] [output file (optional)]

import sys

clusterFile = open(sys.argv[2], "r")
clusterDict = {}
i = 1
for line in clusterFile:
    tokens = line.strip().split()
    for t in tokens:
        clusterDict[t] = i
    i += 1
clusterFile.close()
predFile = open(sys.argv[1], "r")
outFile = sys.out
if len(sys.argv > 3):
    outFile = open(sys.argv[3], "w")
    
for line in predFile:
    tokens = line.strip().split()
    if tokens[0] in clusterDict:
        outFile.write(str(clusterDict[tokens[0]]) + "\t")
        outFile.write("\t".join(tokens[1:]) + "\n")
predFile.close()
outFile.close()
