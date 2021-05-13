#Script to preform PGLMM analysis on enhancer activity predictions. This script is for Windows/RStudio or command-line use, depedning on which lines are (un)commented.

library(ape)
library(phylolm)
#source("D:\\Documents\\Daniel's folder\\CMU\\Pfenning\\fast_bin_perm.r")

#Positional Arguments: tree, predictions, species of prediction columns, traits, output file template, output file number and inital line in preds file, step size in preds file, number of trials, seed to use, trait column (pos), trait column (neg)
#E.g. if output file given is /path/to/foo.csv and number is i, output is in /path/to/foo_ri.csv
#Output file # must be the same as the line to start in the input file (laziness)
#Template must end in .csv
#What trait to use is currently still hardcoded
# given last three arguments i,j,k,
#does k shiffles for enhancer i, i+j, i+2j, ... until end of file is reached


args <- commandArgs()
#print(args[14])
source(paste(args[16], "/fast_bin_perm.r", sep=""))
seed = as.integer(args[14])
set.seed(seed)
#print(seed)
#print(args)

#ATTN: Change method to read.tree for newick tree and read.nexus for nexus tree
tree <- read.tree(file = args[6])
#tree <- read.tree(file = "Zoonomia_ChrX_lessGC40_241species_30Consensus.tree")

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
#traits = read.csv("meyer_diet_4-7-21.csv")
trait.col = args[15]
#trait.col = "True_carn_v_all"
trait.all = traits[[trait.col]]
iv = which(trait.all %in% c(0,1))
trait = as.numeric(as.character(trait.all[iv]))
species.spaces = traits$Species.Name[iv]
trait.species = gsub(" ", "_", species.spaces)
names(trait) = trait.species

#Possible TODO - replace reading with reading one line at a time. Might reduce memory usage by each run
preds = read.csv(file = args[7], header = F, sep = "\t")
#preds = read.csv(file = "liver\\MoRaMaCoPiEnhancersNRSummit.txt", header = F, sep = "\t")
names(preds)[1] = "Enhancer"
te = read.csv(file = args[8], header=F)
#te = read.csv(file = "BoreoeutheriaSpeciesNamesNew.txt", header=F)
pred.species = gsub(" ", "_", te$V1)
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
#num_shuffles = 4 #Should be 0 for non-random

random = T
if (num_shuffles == 0) {
  random=F
  num_shuffles = 1
}

max_iter = (nrow(preds)-row_init) %/% row_step
n = (max_iter + 1) * num_shuffles
t = 0.5 * length(common.species)

enh.names = character(n)
p.vals = double(n)
coeffs = double(n)
index = 1

#t.tree = 0
#t.fit = 0
options(warn = -1) #suppress warning from phyloglm that boundaries of parameters are reached
ptm <- proc.time()
for (m in 0:max_iter) {
  e = row_init + m*row_step
  name = as.character(preds[e, 1])
  cur.preds = preds[e, 2:(length(pred.species)+1)]
  good.preds = cur.preds[which(cur.preds != -1)]
  int.species = intersect(names(good.preds), common.species)
  l = length(int.species)
  int.trait = trait[int.species]
  sum.trait = sum(int.trait)
  if (l >= t && sum.trait > 0 && sum.trait < l) { #50%
    int.preds = good.preds[int.species]
    int.tree = keep.tip(tree.common, int.species)
    int.tree.di = multi2di(int.tree)
    if (random) {
      leafMap=makeLeafMap(int.tree.di)
      fg.species = names(int.trait[which(int.trait == 1)])
      bg.species = names(int.trait[which(int.trait == 0)])
      fg.leaf.count = length(fg.species)
      fg.internal.count = countInternal(int.tree.di, leafMap, fg.species)
      rate.matrix=ratematrix(int.tree.di, int.trait)
    }
    for (f in 1:num_shuffles) {
      if (random) {
        fg.species.shuffled = fastSimBinPhenoVec(int.tree.di, tips=fg.leaf.count, fg.internal.count, rm=rate.matrix, leafBitMaps=leafMap)
        int.trait = double(l)
        names(int.trait) = int.species
        int.trait[fg.species.shuffled] = 1
      }
      dat <- data.frame(Y=int.trait, X=as.double(int.preds), row.names = int.species)
      m <- phyloglm(Y ~ X, data = dat, phy=int.tree.di,  method = "logistic_MPLE")
      enh.names[index] = name
      m.coeff = summary(m)$coefficients
      p.vals[index] = m.coeff[8]
      coeffs[index] = m.coeff[2]
      index = index + 1
    }
  }
}
proc.time() - ptm
options(warn = 1)


dat = data.frame(Enhancer = enh.names[1:index-1], Pvalue = p.vals[1:index-1], Coeff = coeffs[1:index-1])
#write.csv(dat, "liver\\pglmm_results", row.names = FALSE)
write.csv(dat, sub(".csv", paste("_r", args[11], "_s", args[14], ".csv", sep=""),  args[10]), row.names = FALSE)


