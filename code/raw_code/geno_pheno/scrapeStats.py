#!/usr/bin/env python3

import sys, urllib.request

inFile = open(sys.argv[1], "r")
inFile.readline() #Header
for line in inFile:
    tokens = line.split("\t")
    species = tokens[1]
    genomeFile = tokens[15]
    statFile = genomeFile.split()[1].replace("_genomic.fna.gz", "_assembly_stats.txt")
    cn50 = "-"
    cl50 = "-"
    if statFile[:3] == "ftp":
        r = urllib.request.urlopen(statFile)
        stats = r.read().decode()
        lines = stats.split("\n")
        for line in lines:
            if len(line) > 0 and line[0] != "#":
                tokens = line.split()
                if tokens[4] == "contig-N50":
                    cn50 = tokens[5]
                elif tokens[4] == "contig-L50":
                    cl50 = tokens[5]
    print(species + "\t" + cn50 + "\t" + cl50)
inFile.close()