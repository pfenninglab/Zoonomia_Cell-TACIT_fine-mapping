#!/usr/bin/env python3
#assignUniqueNames.py [input file] [output file (optional)]

import sys

inFile = open(sys.argv[1], "r")
outFile = sys.stdout
if len(sys.argv) > 2:
    outFile = open(sys.argv[2], "w")
    
i = 1
for line in inFile:
    name = line.split()[0];
    newName = "peak"+str(i)
    outFile.write(line.replace(name, newName, 1))
    i += 1
inFile.close()
outFile.close()
