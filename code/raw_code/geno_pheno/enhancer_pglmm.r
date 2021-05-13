#Script to preform PGLMM analysis on enhancer activity predictions. This script is for Windows/RStudio or command-line use, depedning on which lines are (un)commented.

library(ape)
library(nlme)
#source("D:\\Documents\\Daniel's folder\\CMU\\Pfenning\\custom_rer.r")
source("/home/dschaffe/enhancer-clustering/scripts/custom_rer.r")

#Positional Arguments: tree, predictions, species of prediction columns, traits, output file template, output file number and inital line in preds file, step size in preds file, number of trials, seed to use, trait column (pos), trait column (neg)
#E.g. if putput file given is /path/to/foo.csv and number is i, output is in /path/to/foo_ri.csv
#Output file # must be the same as the line to start in the input file (laziness)
#Template must end in .csv
#What trait to use is currently still hardcoded
# given last three arguments i,j,k,
#does k shiffles for enhancer i, i+j, i+2j, ... until end of file is reached


args <- commandArgs()
#print(args[14])
seed = as.integer(args[14])
set.seed(seed)
#print(seed)
#print(args)

#ATTN: Change method to read.tree for newick tree and read.nexus for nexus tree
tree <- read.nexus(file = args[6])
#tree <- read.nexus(file = "200mammals_x_lessgc40_75consensus.tree")

#For reading Mammal diet spreadsheet
#traits = read.csv(file= args[9])
#traits = read.csv("MammalDIET_v1.0.csv")
#trait.col = args[15]
#trait.col = "Vertebrate"
#trait.all = traits[[trait.col]]
#iv = which(trait.all != 3)
#trait = replace(trait.all[iv], which(trait.all[iv] == 2), 1) #Currently 1 & 2 -> 1, 0 -> 0
#genera = traits$Genus[iv]
#species = traits$Species[iv]
#trait.species = paste(genera, species, sep="_") #Alert! duplicates. Not sure if matters.

#For reading Meyer lab data
traits = read.csv(file= args[9])
#traits = read.csv("meyer_diet_2-14-21.csv")
trait.col = args[15]
#trait.col = "Vert..vs.herb"
trait.all = traits[[trait.col]]
iv = which(trait.all %in% c(0,1))
trait = as.numeric(as.character(trait.all[iv]))
species.spaces = traits$Species.Name[iv]
trait.species = sub(" ", "_", species.spaces)
names(trait) = trait.species

#Possible TODO - replace reading with reading one line at a time. Might reduce memory usage by each run
preds = read.csv(file = args[7], header = F, sep = "\t")
#preds = read.csv(file = "liver\\MoRaMaCoPiEnhancersNRSummit.txt", header = F, sep = "\t")
names(preds)[1] = "Enhancer"
te = read.csv(file = args[8], header=F)
#te = read.csv(file = "BoreoeutheriaSpeciesNamesNew.txt", header=F)
pred.species = sub(" ", "_", te$V1)
names(preds)[2:(length(pred.species)+1)] = pred.species

common.species = intersect(intersect(pred.species, tree$tip.label), trait.species)
te = which(trait.species %in% common.species)
tree.common = keep.tip(tree, common.species)


########

row_init = as.integer(args[11])
row_step = as.integer(args[12])
num_shuffles = as.integer(args[13])
#row_init = 1
#row_step = 1
#num_shuffles = 1 #Should be 1 for non-random

max_iter = (nrow(preds)-row_init) %/% row_step
n = (max_iter + 1) * num_shuffles
t = 0.5 * length(common.species)

enh.names = character(n)
p.vals = double(n)
coeffs = double(n)
random = T
index = 1

#t.tree = 0
#t.fit = 0
ptm <- proc.time()
for (m in 0:max_iter) {
  e = row_init + m*row_step
  name = as.character(preds[e, 1])
  cur.preds = preds[e, 2:(length(pred.species)+1)]
  good.preds = cur.preds[which(cur.preds != -1)]
  int.species = intersect(names(good.preds), common.species)
  l = length(int.species)
  if (l >= t) {
    int.preds = good.preds[int.species]
    int.trait = trait[int.species]
    int.tree = keep.tip(tree.common, int.species)
    int.tree.di = multi2di(int.tree)
    if (random) {
      fg.species = names(int.trait[which(int.trait == 1)])
      fg.trait.tree = customForeground2Tree(fg.species, int.tree.di, clade="all")
      fg.node.count = sum(fg.trait.tree$edge.length == 1)
      bg.species = names(int.trait[which(int.trait == 0)])
    }
    for (f in 1:num_shuffles) {
      if (random) {
        root.species = sample(bg.species, 1)
        #ptm2 <- proc.time()
        fg.species.shuffled = customSimBinPhenoVec(int.tree.di, root.species, int.trait, fg.node.count)
        #t.tree = t.tree + proc.time() - ptm2
        int.trait = double(l)
        names(int.trait) = int.species
        int.trait[fg.species.shuffled] = 1
      }
      dat <- data.frame(Y=int.trait, X1=as.double(int.preds), row.names = int.species)
      #ptm3 <- proc.time()
      m <- binaryPGLMM(Y ~ X1, phy=int.tree.di, data=dat)
      #t.fit = t.fit + proc.time() - ptm3
      enh.names[index] = name
      p.vals[index] = m$B.pvalue[2]
      coeffs[index] = m$B[2]
      index = index + 1
    }
  }
}
proc.time() - ptm

dat = data.frame(Enhancer = enh.names[1:index-1], Pvalue = p.vals[1:index-1], Coeff = coeffs[1:index-1])
#write.csv(dat, "liver\\pglmm_results", row.names = FALSE)
write.csv(dat, sub(".csv", paste("_r", args[11], "_s", args[14], ".csv", sep=""),  args[10]), row.names = FALSE)


#Temoprary - report # of orthologs per species
#counts = double(length(pred.species))
#for (e in 1:nrow(preds)) {
#  cur.preds = preds[e, 2:(length(pred.species)+1)]
#  good.preds.indices = which(cur.preds != -1)
#  counts[good.preds.indices] = counts[good.preds.indices] + 1
#}
#dat = data.frame(Species=pred.species, Ortholog_Count=counts)
#write.csv(dat, "liver\\liver_ortholog_counts.csv", row.names=F)


