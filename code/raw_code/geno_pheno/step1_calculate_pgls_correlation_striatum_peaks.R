## conda activate r4
library(ape)     #Phylogenetic tree processing
library(geiger)  #Brownian motion simulation
library(stringr) #Species name manipulation
library(phylolm) #Phylogeny-corrected correlation

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(here)
library(data.table)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)
library(tidyverse)
library(tidymodels)
library(BiocParallel)
library(swfdr)

options(MulticoreParam=quote(MulticoreParam(workers=8)))

CODEDIR='code/raw_code/geno_pheno'
DATADIR='data/raw_data/geno_pheno'
PLOTDIR='figures/exploratory/geno_pheno'
dir.create(here(DATADIR, 'rdas', 'phylolm_out'), showWarnings = F, recursive = T) 

# Read in tree and traits
tree = read.tree(file = here('data/raw_data/ldsc_caudate_zoonomia/traits_pgls', 
                             'Zoonomia_ChrX_lessGC40_241species_30Consensus.tree')) 
traits = read.csv(file= here('data/raw_data/ldsc_caudate_zoonomia/traits_pgls',
                             'traitAnnotations_caudate_phylolm.csv'))

## cell types to do associations over
celltypes = c('MSN_D1', 'MSN_D2', 'MSN_SN', 'Oligo', 'OPC')

idx_trait = "Brain.resid_eg" 
celltype = 'MSN_D1'
idxPeak = 1

for(celltype in celltypes){
  ## read in the predictions
  pred_fn = here('data/raw_data/ldsc_caudate_zoonomia/rdas',
                 paste0('Corces2020.',celltype,'.enhPeaks.avgCNN.predictions.rds'))
  pred_df = readRDS(pred_fn)
  
  for (idx_trait in names(traits)[-c(1:2)]) {
    save_fn = here(DATADIR, 'rdas', 'phylolm_out', 
                   paste0(celltype, '.', idx_trait, '.phylolm.rds'))
    print(paste0('Working on celltype: ',celltype, ', Trait: ', idx_trait))
    if(!file.exists(save_fn)){
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
      ## read in the predictions
      common.species = intersect(intersect(pred.species, tree$tip.label), trait.species)
      tree.common = keep.tip(tree, common.species)
      
      ## main PGLS run
      if(isBinary){
        pgls_df = parallel::mclapply(peaks, function(idxPeak){
          # get predictions for each species
          Y = pred_df2[idxPeak, !is.na(pred_df2[idxPeak,])]
          # subset to non-NA species
          species = Reduce('intersect', list(names(Y), common.species, names(trait)))
          pred = Y[species]; trait = trait[species]
          tree.di = keep.tip(tree.common, species) 
          tree.di = multi2di(tree.di)
          dat <- data.frame(X = as.double(pred), Y = trait, row.names = species)
          # run the trait-prediction association test
          m <- phyloglm(Y ~ X, data = dat, phy=tree.di,  method = "poisson_GEE")
          ret = data.frame(summary(m)$coefficients)[2,]
          rownames(ret) = names(idxPeak)
          return(ret)
        }, mc.cores = 8) %>% rbindlist(idcol = 'peak')
      } else {
        pgls_df = parallel::mclapply(peaks, function(idxPeak){
          # get predictions for each species
          Y = pred_df2[idxPeak, !is.na(pred_df2[idxPeak,])]
          # subset to non-NA species
          species = Reduce('intersect', list(names(Y), common.species, names(trait)))
          pred = Y[species]; trait = trait[species]
          tree.di = keep.tip(tree.common, species) 
          tree.di = multi2di(tree.di)
          dat <- data.frame(X = as.double(pred), Y = trait, row.names = species)
          # run the trait-prediction association test
          m <- phylolm(Y ~ X, data = dat, phy=tree.di, model = "BM")
          ret = data.frame(summary(m)$coefficients)[2,]
          rownames(ret) = names(idxPeak)
          return(ret)
        }, mc.cores = 8) %>% rbindlist(idcol = 'peak')
      }
      
      ## add the trait-cell type and per-peak QC numbers
      pgls_df = pgls_df %>% 
        mutate(celltype = celltype, trait = idx_trait, 
               numMappable = numNA[peak], FDR = p.adjust(p.value, 'fdr')) %>%
        relocate(c(celltype, trait, peak, numMappable), .before = everything())
      saveRDS(pgls_df, save_fn)
    }
  }
}
