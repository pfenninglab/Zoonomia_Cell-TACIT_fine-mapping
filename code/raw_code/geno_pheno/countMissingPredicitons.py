#!/usr/bin/env python3
#foo.py [Predictions] [List of species]

import sys

columnOrderFile = open(sys.argv[2], "r")
columnOrder = [x.strip() for x in columnOrderFile if len(x.strip()) > 0]
columnOrderFile.close()

missingCount = {}
for sp in columnOrder:
    missingCount[sp] = 0

querySums = {}
queryAffectedCounts = {}
queryHeavyCounts = {}
queryBadCounts = {}

predFile = open(sys.argv[1], "r")

for line in predFile:
    te = line.strip().split("\t")
    tokens = te[1:]
    suffix = te[0].split("_")[1]
    count = 0
    for i in range(len(tokens)):
        if tokens[i] == "-1" and suffix in ["mo", "ra"]:
            missingCount[columnOrder[i]] += 1
            count += 1
    if suffix in querySums:
        querySums[suffix] += count
    else:
        querySums[suffix] = count
    if count > 0:
        if suffix in queryAffectedCounts:
            queryAffectedCounts[suffix] += 1
        else:
            queryAffectedCounts[suffix] = 1
    if count > 21:
        if suffix in queryHeavyCounts:
            queryHeavyCounts[suffix] += 1
        else:
            queryHeavyCounts[suffix] = 1
    if count > 109:
        if suffix in queryBadCounts:
            queryBadCounts[suffix] += 1
        else:
            queryBadCounts[suffix] = 1
predFile.close()
            
for sp in sorted(columnOrder):
    print(sp + "\t" + str(missingCount[sp]))

for s in queryAffectedCounts:
    print(s + "\t" + str(querySums[s]) + "\t" + str(queryAffectedCounts[s]) + "\t" + str(queryHeavyCounts[s]) + "\t" + str(queryBadCounts[s]), file=sys.stderr)