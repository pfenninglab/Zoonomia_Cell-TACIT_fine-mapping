library(ape)
library(nlme)
#library(rcompanion)

#Positional Arguments: tree, predictions, species of prediction columns, traits, output file template, output file number
#E.g. if putput file given is /path/to/foo.csv and number is i, output is in /path/to/foo_ri.csv
#Template must end in .csv
#What trait to use is currently still hardcoded



args <- commandArgs() #For now, only some i; everything else hardcoded :-|

tree <- read.tree(file = args[6])

preds = read.csv(file = args[7], header = F, sep = "\t")
names(preds)[1] = "Enhancer"
te = read.csv(file = args[8], header=F)
pred.species = sub(" ", "_", te$V1)

traits = read.csv(file= args[9])

########

shuffle = T
threshold = 0.5

traitAll.a = traits$X17.1.MaxLongevity.m
traitAll.b = traits$X5.1.AdultBodyMass.g
indices = intersect(which(traitAll.a != -999), which(traitAll.b != -999))
trait = traitAll.a[indices]  - (12 * 4.88* (traitAll.b[indices])^0.153) 
species.spaces <- traits$Species[indices]
trait.species = sub(" ", "_", species.spaces) #Alert! duplicates. Not sure if matters.


n = nrow(preds)
t = threshold * length(pred.species)

#for (rn in 0:2) {
#rand_index = as.integer(args[6]) + (rn*333) 


names = character(n)
p.vals = double(n)
#prs = double(n) #Pseudo R^2
#like.p.vals = double(n) #Liklehood ratio test p-values

for (e in 1:n) {
  name = as.character(preds[e, 1])
  names[e] = name
  cur.preds = preds[e, 2:(length(pred.species)+1)]
  good.preds = cur.preds[which(cur.preds != -1)]
  good.species = pred.species[which(cur.preds != -1)]
  int.species = good.species[which(good.species %in% trait.species)]
  int.preds = good.preds[which(good.species %in% trait.species)]
  int.trait = trait[match(int.species, trait.species)]
  if (shuffle) {
    int.trait = sample(int.trait)
  }
  if (length(int.preds) >= t) { #50%
    dat <- data.frame(Species = int.species, X = as.double(int.preds), Y = int.trait)
    m1 <- gls(Y ~ X, dat, correlation=corBrownian(1, tree, form = ~Species))
    res = as.data.frame(summary(m1)$tTable)
    p.vals[e] = res$p[2]
    #te = nagelkerke(m1)
    #prs[e] = te$P[2] #Cox & Snell pseudo R-squared
    #like.p.vals[e] = te$L[4]
  } else {
    p.vals[e] = -1.0
  }
}

filt.names = names[which(p.vals != -1)]
filt.p.vals = p.vals[which(p.vals != -1)]
cor.p.vals = filt.p.vals * length(filt.p.vals)
#filt.prs = prs[which(p.vals != -1)]
#filt.like = like.p.vals[which(p.vals != -1)]
#cor.like = filt.like * length(filt.p.vals)

#dat = data.frame(Enhancer = filt.names, Pvalue = filt.p.vals, Corrected_Pvalue = cor.p.vals, Pseudo_R2 = filt.prs, Liklikood_Pvalue = filt.like, Corrected_Like_Pvalue = cor.like)
dat = data.frame(Enhancer = filt.names, Pvalue = filt.p.vals, Corrected_Pvalue = cor.p.vals)
write.csv(dat, sub(".csv", paste("_r", args[11], ".csv", sep=""),  args[10]), row.names = FALSE)



