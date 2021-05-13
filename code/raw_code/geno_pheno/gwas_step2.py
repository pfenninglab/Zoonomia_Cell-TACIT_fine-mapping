#!/usr/bin/env python

#As argument, give the background bed file, which is assumed to be named foo_ProjectName.bed
#In particular, this constrains the project name to not include an underscore.
#Additionally, the .bed directory cannot be root

from glob import glob
from os.path import abspath, split, basename
import sys


bgAbsName = abspath(sys.argv[1])
bedPath, bgFileName = split(bgAbsName)
bedGlob = bedPath + "/*.bed"
bgName, te = bgFileName.rsplit("_", 1)
prName = te[:-4] #remove .bed
annDir = bedPath + "/annotations/" + prName + "/"

with open(bedPath + "/" + prName + ".ldcts", "w") as f:
    for file in glob(bedGlob):
        name = basename(str(file))
        if name != bgFileName:
            f.write(name[:-4] + " " + annDir + prName + "_" + name[:-3] + "," + annDir + bgName + "_" + prName + ".\n")


