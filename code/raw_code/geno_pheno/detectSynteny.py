#!/usr/bin/env python3

#foo.py [Instruction file] [list of species]
#Instruction file is of the form species \t peak file \t list of ortholog .gz files
#Precondition: peak files are 1) sorted within each chromosome and has contiguous chromosomes
#Precondition: both peak & ortholog files have peak names in column 4 under the same naming scheme
#Note: The orthologs of the first and last peaks that have an ortholog
#      on each chromosome will never be rejected
#Rejected 3119459 out of 19090876 orthologs


import sys, gzip
def findSubstring(s, l):
    matches = []
    for x in l:
        if s in x:
            matches.append(x)
    if len(matches) == 0:
        print("cannot file file for", s, file=sys.stderr)
        return ""
    else: 
        return min(matches, key=len) #This solves species/subspecies problem


columnOrderFile = open(sys.argv[2], "r")
columnOrder = [x.strip() for x in columnOrderFile if len(x.strip()) > 0]
columnOrderFile.close()

print("Original Species\tOrtholog Species\tPeak Name")

insFile = open(sys.argv[1], "r")
rCount = 0
count = 0
for line in insFile:
    ins = line.strip().split("\t")
    if len(ins) != 3:
        print(tokens2, file=sys.stderr)    
        continue
    
    baseSpecies = ins[0]
    orthFileListFile = open(ins[2], "r")
    orthFiles = [x.strip() for x in orthFileListFile if len(x.strip()) > 0]
    orthFileListFile.close()
    
    for otherSpecies in columnOrder:
        if otherSpecies == baseSpecies:
            continue
            
        orthFileName = findSubstring(otherSpecies.replace(" ", "_"), orthFiles)
        #print("For", otherSpecies, "using", orthFileName, file=sys.stderr)
        orthChrDict = {}
        orthFile = gzip.open(orthFileName, "rb")
        orthFileLines = orthFile.read().decode().split("\n")
        orthFile.close()
        for line in orthFileLines:
            tokens = line.strip().split("\t")
            if len(tokens) >= 4:
                orthChrDict[tokens[3]] = tokens[0]
        
        bedFile = open(ins[1])
        oldChr = ""
        olderOrth = [] #Peak, chr
        oldOrth = [] #Peak, chr
        for line in bedFile:
            tokens = line.strip().split("\t")
            if len(tokens) < 4:
                continue
            chr = tokens[0]
            peak = tokens[3]
            if peak not in orthChrDict:
                continue
            orthChr = orthChrDict[peak]
            count += 1
            if chr == oldChr:
                if len(olderOrth) > 0 and len(oldOrth) > 0 and olderOrth[1] != oldOrth[1] and oldOrth[1] != orthChr:
                    print(baseSpecies + "\t" + otherSpecies + "\t" + peak)
                    rCount += 1
                olderOrth = oldOrth
                oldOrth = [peak, orthChr]
            else:
                oldChr = chr
                olderOrth = []
                oldOrth = []
        bedFile.close()
print("Rejected", rCount, "out of", count, "orthologs", file=sys.stderr)