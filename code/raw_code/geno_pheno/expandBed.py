#! /usr/bin/env python3
# ... [bed file] [radius] > [new bed file]

import sys

inFile = open(sys.argv[1], "r")
radius = int(sys.argv[2])
for line in inFile:
	tokens = line.strip().split("\t")
	if len(tokens) >= 3:
		new1 = str(max(int(tokens[1]) - radius, 0))
		new2 = str(int(tokens[2]) + radius) #Possible problem if too large
		print(tokens[0] + "\t" + new1 + "\t" + new2 + "\t" + "\t".join(tokens[3:]))
	else:
		print(line, end="")
inFile.close()
