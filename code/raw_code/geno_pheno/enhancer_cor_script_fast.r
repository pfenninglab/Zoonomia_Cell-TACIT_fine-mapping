library(ape)
library(nlme)

#Positional Arguments: tree, predictions, species of prediction columns, traits, output file template, output file number and inital line in preds file, step size in preds file, number of trials, seed to use
#E.g. if putput file given is /path/to/foo.csv and number is i, output is in /path/to/foo_ri.csv
#Output file # must be the same as the line to start in the input file (laziness)
#Template must end in .csv
#What trait to use is currently still hardcoded
# given last three arguments i,j,k,
#does k shiffles for enhancer i, i+j, i+2j, ... until end of file is reached


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

traits = read.csv(file= args[9])

########

traitAll.a = traits$X17.1.MaxLongevity.m
traitAll.b = traits$X5.1.AdultBodyMass.g
indices = intersect(which(traitAll.a != -999), which(traitAll.b != -999))
trait = traitAll.a[indices]  - (12 * 4.88* (traitAll.b[indices])^0.153) 
species.spaces <- traits$Species[indices]
trait.species = sub(" ", "_", species.spaces) #Alert! duplicates. Not sure if matters.

row_init = as.integer(args[11])
row_step = as.integer(args[12])
num_shuffles = as.integer(args[13])
max_iter = (nrow(preds)-row_init) %/% row_step

n = (max_iter + 1) * num_shuffles
t = 0.5 * length(pred.species)


#for (rn in 0:2) {
#rand_index = as.integer(args[6]) + (rn*333) 


enh.names = character(n)
p.vals = double(n)
coeffs = double(n)

index = 1
for (m in 0:max_iter) {
  e = row_init + m*row_step
  name = as.character(preds[e, 1])
  cur.preds = preds[e, 2:(length(pred.species)+1)]
  good.preds = cur.preds[which(cur.preds != -1)]
  good.species = pred.species[which(cur.preds != -1)]
  int.species = good.species[which(good.species %in% trait.species)]
  int.preds = good.preds[which(good.species %in% trait.species)]
  int.trait = trait[match(int.species, trait.species)]
  if (length(int.preds) >= t) { #50%
    for (f in 1:num_shuffles) {
	  shuffle.trait = sample(int.trait)
      dat <- data.frame(Species = int.species, X = as.double(int.preds), Y = shuffle.trait)
      m1 <- gls(Y ~ X, dat, correlation=corBrownian(1, tree, form = ~Species))
      res = as.data.frame(summary(m1)$tTable)
	  enh.names[index] = name
      p.vals[index] = res$p[2]
	  coeffs[index] = res$Value[2]
      index = index + 1
    }
  }
}

dat = data.frame(Enhancer = enh.names[1:index-1], Pvalue = p.vals[1:index-1], Coeff = coeffs[1:index-1])
write.csv(dat, sub(".csv", paste("_r", args[11], "_s", args[14], ".csv", sep=""),  args[10]), row.names = FALSE)



