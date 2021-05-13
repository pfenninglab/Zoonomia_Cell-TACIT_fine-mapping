#!/usr/bin/env python3
#assignUniqueNames.py [input file] [output file (optional)]


inFile = open("..\\data\\SomeEnhancers.txt", "r")
species = "Macaque"
outFile = open("..\\data\\" + species + "LookupTable.txt", "w")

i = 1
for line in inFile:
    name = line.split()[0];
    newName = "peak"+str(i)
    if species == "Mouse" and 1 <=i and i <= 12006:
        outFile.write(name + "\t" + newName + "\n")
    elif species == "Human" and 12007 <=i and i <= 13984:
        outFile.write(name + "\t" + newName + "\n")
    elif species == "Macaque" and 13985 <=i and i <= 22521:
        outFile.write(name + "\t" + newName + "\n")
    i += 1

inFile.close()
outFile.close()
    
