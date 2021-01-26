import pyBigWig
bw = pyBigWig.open("/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/phyloP/200m_scoresPhyloP_20201221.bigWig")
for line in open("foo.bed"):
    cols = line.strip().split()
    vals = bw.values(cols[0], int(cols[1]), int(cols[2]))
    # Do something with the values...
bw.close()