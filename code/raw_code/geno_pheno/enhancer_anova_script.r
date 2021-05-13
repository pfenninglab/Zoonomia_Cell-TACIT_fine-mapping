library(phytools)

#Positional Arguments: tree, predictions, species of prediction columns, traits, output file template, output file number and inital line in preds file, step size in preds file, number of sims, seed to use

args <- commandArgs()

print(args[14])
seed = as.integer(args[14])
set.seed(seed)
print(seed)
#print(args)
tree <- read.tree(file = args[6])

preds = read.csv(file = args[7], header = F, sep = "\t")
names(preds)[1] = "Enhancer"
te = read.csv(file = args[8], header=F)
pred.species = sub(" ", "_", te$V1)

traits = read.csv(file = args[9])
traitAll = traits$Animal
iv = which(traitAll != 3)
trait = replace(traitAll[iv], traitAll[iv]==2, 1) #Currently 1 & 2 -> 1, 0 -> 0
genera = traits$Genus[iv]
species = traits$Species[iv]
trait.species = paste(genera, species, sep="_") #Alert! duplicates. Not sure if matters.
unique.traits = trait[which(unique(trait.species) %in% trait.species)] #Filter to only species in Zoonomia
#####################

row_init = as.integer(args[11])
row_step = as.integer(args[12])
num_sims = as.integer(args[13])
max_iter = (nrow(preds)-row_init) %/% row_step

n = (max_iter + 1) * num_sims
t = 0.5 * length(pred.species)


#for (rn in 0:2) {
#rand_index = as.integer(args[6]) + (rn*333) 


enh.names = character(n)
p.vals = double(n)

index = 1
for (m in 0:max_iter) {
  e = row_init + m*row_step
  name = as.character(preds[e, 1])
  cur.preds = preds[e, 2:(length(pred.species)+1)]
  good.preds = cur.preds[which(cur.preds != -1)]
  good.species = pred.species[which(cur.preds != -1)]
  int.species = good.species[which(good.species %in% trait.species)]
  int.preds = as.double(good.preds[which(good.species %in% trait.species)])
  
  if (length(int.preds) >= t) { #50%
    int.trait = trait[match(int.species, trait.species)]
    int.tree = keep.tip(tree, int.species)
    int.preds.treeorder = int.preds[match(int.species, int.tree$tip.label)]
    int.trait.treeorder = int.trait[match(int.species, int.tree$tip.label)]
	res = suppressWarnings(phylANOVA(int.tree, int.trait.treeorder, int.preds.treeorder, nsim=num_sims))
	enh.names[index] = name
    p.vals[index] = res$Pf
	index = index + 1
  }
}

dat = data.frame(Enhancer = enh.names[1:index-1], Pvalue = p.vals[1:index-1])
write.csv(dat, sub(".csv", paste("_r", args[11], "_s", args[14], ".csv", sep=""),  args[10]), row.names = FALSE)

