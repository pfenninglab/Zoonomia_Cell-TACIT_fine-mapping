## conda activate r4
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(here)
library(data.table)
library(tidyverse)
library(ape)     #Phylogenetic tree processing
library(geiger)  #Brownian motion simulation
library(stringr) #Species name manipulation
library(phylolm) #Phylogeny-corrected correlation
library(optparse)

source(here('code/final_code/TACIT-main/fast_bin_perm.r'))

CODEDIR='code/raw_code/geno_pheno'
DATADIR='data/raw_data/geno_pheno'
PLOTDIR='figures/exploratory/geno_pheno'
dir.create(here(DATADIR, 'tables', 'permulations_out'), showWarnings = F, recursive = T) 

## command line options
option_list = list(
  make_option(c("-c", "--cell"), type="numeric", default=NULL, 
              help="index of cell type to run", metavar="character"),
  make_option(c("-t", "--trait"), type="numeric", default=NULL, 
              help="index of trait to run", metavar="character"),
  make_option(c("-p", "--perm"), type="numeric", default=NULL, 
              help="index of permulation set to run", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Read in tree and traits
tree = read.tree(file = here('data/raw_data/ldsc_caudate_zoonomia/traits_pgls', 
                             'Zoonomia_ChrX_lessGC40_241species_30Consensus.tree')) 
traits = read.csv(file= here('data/raw_data/ldsc_caudate_zoonomia/traits_pgls',
                             'traitAnnotations_caudate_phylolm.csv'))

## cell types to do associations over
celltypes = c('MSN_D1', 'MSN_D2', 'MSN_SN', 'Oligo', 'OPC')
nPerms = 10*10^seq(6)

# celltype = celltypes[1]
# idx_trait = names(traits)[-c(1:2)][2]
# numPermute = nPerms[1]

celltype = celltypes[opt$c]
idx_trait = names(traits)[-c(1:2)][opt$t]
numPermute = nPerms[opt$p]

print(paste0('Working on celltype: ',celltype))
print(paste0('Working on trait: ', idx_trait))
print(paste0('Working on permulations batch: ', numPermute))

## read in the predictions
pred_fn = here('data/raw_data/ldsc_caudate_zoonomia/rdas',
             paste0('Corces2020.',celltype,'.enhPeaks.avgCNN.predictions.rds'))
pred_df = readRDS(pred_fn)

## grab the trait, subset to common species
trait = traits %>% dplyr::select(all_of(c('species.binomial', idx_trait))) %>%
  filter(species.binomial %in% names(pred_df)) %>%
  deframe() %>% enframe() %>% filter(!is.na(value)) %>% deframe()
trait.species = names(trait)
isBinary = length(table(trait))==2

# drop species w/ too many NAs
pred_df2 = pred_df %>% dplyr::select( all_of(trait.species), name) %>%
  column_to_rownames('name')
pred_df2 = pred_df2[, colSums(is.na(pred_df2)) < .9 * nrow(pred_df2)]

# drop peaks w/ too many NAs
pred_df2 = pred_df2[rowSums(is.na(pred_df2)) < 0.25 * ncol(pred_df2), ]
pred.species = names(pred_df2)
peaks = setNames(seq(nrow(pred_df2)), rownames(pred_df2))
numNA = rowSums(!is.na(pred_df2))

## subset tree
common.species = intersect(intersect(pred.species, tree$tip.label), trait.species)
tree.common = keep.tip(tree, common.species)

## grab the original run, calculate # of permutations per peak
pgls_df = here(DATADIR, 'rdas', 'phylolm_out', 
               paste0(celltype, '.', idx_trait, '.phylolm.rds')) %>% readRDS()
pgls_df = pgls_df %>% mutate(nPerm = pmin(100 * 10 ^ round(-log10(p.value)), 1e7)) %>%
  filter(nPerm == numPermute)
pgls_df = pgls_df[sample(nrow(pgls_df)), ] # do permulations out of order
peaks = pgls_df %>% mutate(num = seq(n())) %>% 
  dplyr::select(peak, num) %>% deframe()

## final file to put this batch of permulations
save_fn = here(DATADIR, 'tables', 'permulations_out', 
               paste0(celltype, '.', idx_trait,'.permulations.', numPermute,'.txt.gz'))
perm_types = c( rep("NULL", 2), 'character', rep("NULL", 4))
perm_cols = c('peak')

if(file.exists(save_fn)){
  # check cols in progress
  peak_tbl = fread(save_fn, header= F, colClasses = perm_types, col.names = perm_cols) %>%
    deframe() %>% table()
  # see which peaks w/ enough permulations
  peak_tbl = peak_tbl[peak_tbl >= numPermute] 
  # only compute permulations on peaks that haven't been permulated
  peaks = peaks[!names(peaks) %in% names(peak_tbl)]
}

if(numPermute >= 1e6) { ## for the big ones, do in chunks to reduce RAM
  n = 1000
  numPermute = round(numPermute / n)
} else if(numPermute >= 1e5) { ## for the big ones, do in chunks to reduce RAM
  n = 100
  numPermute = round(numPermute / n)
} else {
  n = 1
}

## chunk by # of permutations, peakwise
for(idxPeak in peaks){
  print(paste('On peak', idxPeak, 'of', length(peaks), '.'))
  
  # get predictions for each species
  Y = pred_df2[idxPeak, !is.na(pred_df2[idxPeak,])] %>% t() %>% 
    data.frame() %>% rownames_to_column() %>% deframe()
  # subset to non-NA species
  species = Reduce('intersect', list(names(Y), common.species, names(trait)))
  pred = Y[species]; trait2 = trait[species]
  tree2 = keep.tip(tree.common, species) 
  tree.di = multi2di(tree2)
  rate.matrix=ratematrix(tree2, trait2)
  
  for(i in seq(n)){
    if(isBinary){
      leafMap=makeLeafMap(tree.di)
      fg.species = names(trait2[which(trait2 == 1)])
      bg.species = names(trait2[which(trait2 == 0)])
      fg.leaf.count = length(fg.species)
      fg.internal.count = countInternal(tree.di, leafMap, fg.species)
      
      ## permulate binary traits
      peak_perm_list = parallel::mclapply(seq(numPermute), function(x){ 
        fg.species.shuffled = fastSimBinPhenoVec(tree.di, tips=fg.leaf.count, 
                                                 fg.internal.count, rm=rate.matrix, 
                                                 leafBitMaps=leafMap)
        Y2 = setNames(double(length(species)), species)
        Y2[fg.species.shuffled] = 1
        dat = data.frame(X = pred[names(Y2)], Y = Y2)
        ## test permulated traits
        m <- phyloglm(Y ~ X, data = dat, phy=tree.di, method = "logistic_MPLE")
        ret = data.frame(summary(m)$coefficients)[2,]
        rownames(ret) = names(idxPeak)
        return(ret)
      }, mc.cores = min(parallel::detectCores(), 8))
    } else{
      ## permulate continuous traits
      peak_perm_list = parallel::mclapply(seq(numPermute), function(x){ 
        Y2 = simpermvec(trait2, tree.di, rm=rate.matrix)
        dat = data.frame(X = pred[names(Y2)], Y = Y2)
        ## test permulated traits
        m <- phylolm(Y ~ X, data = dat, phy=tree.di, model = "BM")
        ret = data.frame(summary(m)$coefficients)[2,]
        rownames(ret) = names(idxPeak)
        return(ret)
      }, mc.cores = min(parallel::detectCores(), 8))
    }
    
    ## export each peak's permulated tests and append to save file
    peak_perm_df  = peak_perm_list %>% rbindlist() %>% 
      mutate(celltype = celltype, trait = idx_trait, peak = names(peaks)[idxPeak]) %>%
      relocate(c(celltype, trait, peak), .before = everything()) %>%
      write_tsv(save_fn, append = T)
  }
}