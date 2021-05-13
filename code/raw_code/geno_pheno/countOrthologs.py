#!/usr/bin/env python3

#countOrthologs.py [enhancer predictions] [species order] [species of interest]

import sys

specFile = open(sys.argv[3], "r")
species = [x.strip().replace("_", " ") for x in specFile.readlines() if not x.isspace()]
specFile.close()

orthologCounts = {}
for s in species:
    orthologCounts[s] = 0 

orderFile = open(sys.argv[2], "r")
order = [x.strip() for x in orderFile.readlines() if not x.isspace()]
orderFile.close

predFile = open(sys.argv[1])
for line in predFile:
    tokens = line.strip().split()[1:]
    if len(tokens) == 0:
        continue
    for i in range(len(tokens)):
        spec = order[i]
        if spec in species and tokens[i] != "-1":
             orthologCounts[spec] += 1
predFile.close()

for spec in species:
    print(spec + "\t" + str(orthologCounts[spec]))
