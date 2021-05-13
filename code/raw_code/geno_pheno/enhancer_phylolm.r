#Script to preform PGLMM analysis on enhancer activity predictions. This script is for Windows/RStudio or command-line use, depedning on which lines are (un)commented.

library(ape)
library(geiger)
library(stringr)
library(phylolm)

#Copied from RERconverge
#'Generates a permulated continuous phenotype given an observed continuous phenotype and a phylogeny
#' @param namedvec A named numeric vector with phenotype values for each speices
#' @param treewithbranchlengths A rooted phylogenetic tree with the same species as namedvec and branch lengths representing average evolutionary rate.
#' @param rm A precomputed rate matrix
#' @return A vector with permulated phenotype values
#' @export
simpermvec=function (namedvec, treewithbranchlengths, rm=NULL) 
{
  if(is.null(rm)) {
    rm = ratematrix(treewithbranchlengths, namedvec)
  }
  sims = sim.char(treewithbranchlengths, rm, nsim = 1)
  nam = rownames(sims)
  s = as.data.frame(sims)
  simulatedvec = s[, 1]
  names(simulatedvec) = nam
  simsorted = sort(simulatedvec)
  realsorted = sort(namedvec)
  l = length(simsorted)
  c = 1
  while (c <= l) {
    simsorted[c] = realsorted[c]
    c = c + 1
  }
  simsorted
}


#Positional Arguments: 
#    tree, 
#    predictions, 
#    species of prediction columns, 
#    traits, 
#    output file template, 
#    output file number and inital line in preds file, 
#    step size in preds file, 
#    number of trials, 
#    seed to use, 
#    trait column
#E.g. if putput file given is /path/to/foo.csv and number is i, output is in /path/to/foo_ri.csv
#Output file # must be the same as the line to start in the input file (laziness)
#Template must end in .csv
#What trait to use is currently still hardcoded
# given last three arguments i,j,k,
#does k shiffles for enhancer i, i+j, i+2j, ... until end of file is reached


args <- commandArgs()
seed = as.integer(args[14])
set.seed(seed)
#print(args)
#print(seed)

#ATTN: Change method to read.tree for newick tree and read.nexus for nexus tree
tree <- read.tree(file = args[6])
#tree <- read.tree(file = "Zoonomia_ChrX_lessGC40_241species_30Consensus.tree")

##For reading panthiera
#traits = read.csv(file= args[9])
#traits = read.csv("panthiera.csv")
#trait.col = args[15]
#trait.col = "X9.1.GestationLen.d"
#trait.all = traits[[trait.col]]
#iv = which(trait.all > 0)
#trait = trait.all[iv]
#species.spaces = traits$Species[iv]
#trait.species = sub(" ", "_", species.spaces) #Alert! duplicates. Not sure if matters.

##For reading zoonomia sheet
traits = read.csv(file= args[9])
#traits = read.csv("teeling_longevity_2-12-21.csv")
trait.col = args[15]
#trait.col = "LQ"
trait.all = traits[[trait.col]]
iv = which(trait.all != "NA") 
trait = trait.all[iv]
species.lower = traits$species.binomial[iv]
trait.species = str_to_title(species.lower)
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
#num_shuffles = 1 #Should be 1 for non-random

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

options(warn = -1) #suppress warning from phylolm
ptm <- proc.time()
for (m in 0:max_iter) {
  e = row_init + m*row_step
  name = as.character(preds[e, 1])
  cur.preds = preds[e, 2:(length(pred.species)+1)]
  good.preds = cur.preds[which(cur.preds != -1)]
  int.species = intersect(names(good.preds), common.species)
  l = length(int.species)
  if (l >= t) { #50%
    int.preds = good.preds[int.species]
    int.trait = trait[int.species]
    int.tree = keep.tip(tree.common, int.species)
    int.tree.di = multi2di(int.tree)
    if (random) {
      rate.matrix=ratematrix(int.tree.di, int.trait)
      int.trait.real = int.trait
      names(int.trait.real) = int.species
    }
  
    for (f in 1:num_shuffles) {
      if (random) {
        int.trait = simpermvec(int.trait.real, int.tree.di, rm=rate.matrix)
      }
      dat <- data.frame(X = as.double(int.preds), Y = int.trait)
      m <- phylolm(Y ~ X, data = dat, phy=int.tree.di, model = "BM")
      res = as.data.frame(summary(m)$tTable)
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
#write.csv(dat, "liver\\pgls_lq_results.csv", row.names = FALSE)
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


