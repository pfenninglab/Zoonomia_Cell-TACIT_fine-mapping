# to be run in the root github directory
LABEL='cnn_enhancer_ortholog'
DATADIR=file.path('/home/bnphan/projects/snATAC_cross_species_caudate/data/raw_data/',LABEL)
DATADIR=file.path('data/raw_data/',LABEL)

### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(rtracklayer))
library(tidyverse)
library(here)

source(here('code/raw_code/hal_scripts/narrowPeakFunctions.R'))
source(here('code/raw_code/hal_scripts/gen_enh_ortholog_sets.R'))

## chromosomal splits
testSet = c('chr1','chr2')
folds = list(fold1 = c('chr6', 'chr13', 'chr21'),
             fold2 = c('chr7', 'chr14', 'chr18'),
             fold3 = c('chr11','chr17', 'chrX'),
             fold4 = c('chr9', 'chr12'),
             fold5 = c('chr8', 'chr10'))
split = c('test', 'valid')

###########################################
# load in the human positive peak files
genome = 'hg38'
load(here(DATADIR, 'rdas', paste('caudate_positive_sequences', 'hg38','.rda', sep = '.')))
load(here(DATADIR, 'rdas', paste('caudate_nonEnh_sequences', 'hg38','.rda', sep = '.')))

### find cell type peaks not in other cell types non-reproducible peaks
humanCelltypeOnly_enhList = lapply(names(human_enhList), function(cell){
  excl_gr = GenomicRanges::reduce(unlist(GRangesList(human_pos_list)[names(human_peakList) %ni% cell]))
  oo = findOverlaps(query = human_enhList[[cell]], subject = excl_gr)
  keepIdx = which(!seq_along(human_enhList[[cell]]) %in% queryHits(oo))
  return(human_enhList[[cell]][keepIdx])
}) %>% GRangesList()
names(humanCelltypeOnly_enhList) = names(human_enhList)
lengths(humanCelltypeOnly_enhList)

humanCelltypeOnly_enhList_split = lapply(folds, function(fold){
  ret = lapply(humanCelltypeOnly_enhList, splitPeakSet, testSet = testSet, validSet = fold, useCol = 'col')
  return(ret)
})

### find non-celltype peaks not overlapping non-reproducible cell type ####
humanNonCelltype_enhList = lapply(names(human_enhList), function(cell){
  nonCelltype_gr = unlist(human_enhList[names(human_enhList) %ni% cell])
  oo = findOverlaps(query = nonCelltype_gr, subject = human_pos_list[[cell]])
  keepIdx = which(!seq_along(nonCelltype_gr) %in% unique(queryHits(oo)))
  return(nonCelltype_gr[keepIdx])
}) %>% GRangesList()
names(humanNonCelltype_enhList) = names(human_enhList)
lengths(humanNonCelltype_enhList)

humanNonCelltype_enhList_split = lapply(folds, function(fold){
  ret = lapply(humanNonCelltype_enhList, splitPeakSet, testSet = testSet, validSet = fold, useCol = 'col')
  return(ret)
})

# export postive sequences to summit-centered 501bp.fasta file
system(paste('mkdir -p',  file.path(DATADIR, 'fasta')))
for(cell in names(human_enhList)){
  for(fold in names(folds)){
    # write the positives
    pos.fasta_fn = file.path(DATADIR, 'fasta', 
                             paste(paste0(genome, 'CelltypeOnly'), cell, fold, split, 'positive.fa.gz', sep = '_'))
    posFasta = mapply(writeGRangesToFasta, gr = humanCelltypeOnly_enhList_split[[fold]][[cell]][split],  
                      file = pos.fasta_fn, genome = genome)
    neg.fasta_fn = file.path(DATADIR, 'fasta', 
                             paste(paste0(genome, 'NonCelltype'), cell, fold, split, 'negative.fa.gz', sep = '_'))
    negFasta = mapply(writeGRangesToFasta, gr = humanNonCelltype_enhList_split[[fold]][[cell]][split],  
                      file = neg.fasta_fn, genome = genome)
  }}

# save_fn = file.path(DATADIR, 'rdas', paste('caudate_sequences', paste0(genome, 'CelltypeOnlyNonCelltype'),'rda', sep = '.'))
# save(humanCelltypeOnly_enhList, humanNonCelltype_enhList, file = save_fn )


###########################################
# load in the mouse positive peak files
genome = 'mm10'
load(here(DATADIR, 'rdas', paste('caudate_positive_sequences', genome,'.rda', sep = '.')))
load(here(DATADIR, 'rdas', paste('caudate_nonEnh_sequences', genome,'.rda', sep = '.')))

### find cell type peaks not in other cell types non-reproducible peaks
mouseCelltypeOnly_enhList = lapply(names(mouse_enhList), function(cell){
  excl_gr = GenomicRanges::reduce(unlist(GRangesList(mouse_pos_list)[names(mouse_peakList) %ni% cell]))
  oo = findOverlaps(query = mouse_enhList[[cell]], subject = excl_gr)
  keepIdx = which(!seq_along(mouse_enhList[[cell]]) %in% queryHits(oo))
  return(mouse_enhList[[cell]][keepIdx])
}) %>% GRangesList()
names(mouseCelltypeOnly_enhList) = names(mouse_enhList)
lengths(mouseCelltypeOnly_enhList)

mouseCelltypeOnly_enhList_split = lapply(folds, function(fold){
  ret = lapply(mouseCelltypeOnly_enhList, splitPeakSet, testSet = testSet, validSet = fold, useCol = 'col')
  return(ret)
})

### find non-celltype peaks not overlapping non-reproducible cell type ####
mouseNonCelltype_enhList = lapply(names(mouse_enhList), function(cell){
  nonCelltype_gr = unlist(mouse_enhList[names(mouse_enhList) %ni% cell])
  oo = findOverlaps(query = nonCelltype_gr, subject = mouse_pos_list[[cell]])
  keepIdx = which(!seq_along(nonCelltype_gr) %in% unique(queryHits(oo)))
  return(nonCelltype_gr[keepIdx])
}) %>% GRangesList()
names(mouseNonCelltype_enhList) = names(mouse_enhList)
lengths(mouseNonCelltype_enhList)

mouseNonCelltype_enhList_split = lapply(folds, function(fold){
  ret = lapply(mouseNonCelltype_enhList, splitPeakSet, testSet = testSet, validSet = fold, useCol = 'col')
  return(ret)
})

# export postive sequences to summit-centered 501bp.fasta file
system(paste('mkdir -p',  file.path(DATADIR, 'fasta')))
for(cell in names(mouse_enhList)){
  for(fold in names(folds)){
    # write the positives
    pos.fasta_fn = file.path(DATADIR, 'fasta', 
                             paste(paste0(genome, 'CelltypeOnly'), cell, fold, split, 'positive.fa.gz', sep = '_'))
    posFasta = mapply(writeGRangesToFasta, gr = mouseCelltypeOnly_enhList_split[[fold]][[cell]][split],  
                      file = pos.fasta_fn, genome = genome)
    neg.fasta_fn = file.path(DATADIR, 'fasta', 
                             paste(paste0(genome, 'NonCelltype'), cell, fold, split, 'negative.fa.gz', sep = '_'))
    negFasta = mapply(writeGRangesToFasta, gr = mouseNonCelltype_enhList_split[[fold]][[cell]][split],  
                      file = neg.fasta_fn, genome = genome)
  }}

# save_fn = file.path(DATADIR, 'rdas', paste('caudate_sequences', paste0(genome, 'CelltypeOnlyNonCelltype'),'rda', sep = '.'))
# save(mouseCelltypeOnly_enhList, mouseNonCelltype_enhList, file = save_fn )






###########################################
# load in the monkey positive peak files
genome = 'rheMac10'
load(here(DATADIR, 'rdas', paste('caudate_positive_sequences', genome,'.rda', sep = '.')))
load(here(DATADIR, 'rdas', paste('caudate_nonEnh_sequences', genome,'.rda', sep = '.')))

### find cell type peaks not in other cell types non-reproducible peaks
rhesusCelltypeOnly_enhList = lapply(names(rhesus_enhList), function(cell){
  excl_gr = GenomicRanges::reduce(unlist(GRangesList(rhesus_pos_list)[names(rhesus_peakList) %ni% cell]))
  oo = findOverlaps(query = rhesus_enhList[[cell]], subject = excl_gr)
  keepIdx = which(!seq_along(rhesus_enhList[[cell]]) %in% queryHits(oo))
  return(rhesus_enhList[[cell]][keepIdx])
}) %>% GRangesList()
names(rhesusCelltypeOnly_enhList) = names(rhesus_enhList)
lengths(rhesusCelltypeOnly_enhList)

rhesusCelltypeOnly_enhList_split = lapply(folds, function(fold){
  ret = lapply(rhesusCelltypeOnly_enhList, splitPeakSet, testSet = testSet, validSet = fold, useCol = 'col')
  return(ret)
})

### find non-celltype peaks not overlapping non-reproducible cell type ####
rhesusNonCelltype_enhList = lapply(names(rhesus_enhList), function(cell){
  nonCelltype_gr = unlist(rhesus_enhList[names(rhesus_enhList) %ni% cell])
  oo = findOverlaps(query = nonCelltype_gr, subject = rhesus_pos_list[[cell]])
  keepIdx = which(!seq_along(nonCelltype_gr) %in% unique(queryHits(oo)))
  return(nonCelltype_gr[keepIdx])
}) %>% GRangesList()
names(rhesusNonCelltype_enhList) = names(rhesus_enhList)
lengths(rhesusNonCelltype_enhList)

rhesusNonCelltype_enhList_split = lapply(folds, function(fold){
  ret = lapply(rhesusNonCelltype_enhList, splitPeakSet, testSet = testSet, validSet = fold, useCol = 'col')
  return(ret)
})

# export postive sequences to summit-centered 501bp.fasta file
system(paste('mkdir -p',  file.path(DATADIR, 'fasta')))
for(cell in names(rhesus_enhList)){
  for(fold in names(folds)){
    # write the positives
    pos.fasta_fn = file.path(DATADIR, 'fasta', 
                             paste(paste0(genome, 'CelltypeOnly'), cell, fold, split, 'positive.fa.gz', sep = '_'))
    posFasta = mapply(writeGRangesToFasta, gr = rhesusCelltypeOnly_enhList_split[[fold]][[cell]][split],  
                      file = pos.fasta_fn, genome = genome)
    neg.fasta_fn = file.path(DATADIR, 'fasta', 
                             paste(paste0(genome, 'NonCelltype'), cell, fold, split, 'negative.fa.gz', sep = '_'))
    negFasta = mapply(writeGRangesToFasta, gr = rhesusNonCelltype_enhList_split[[fold]][[cell]][split],  
                      file = neg.fasta_fn, genome = genome)
  }}

# save_fn = file.path(DATADIR, 'rdas', paste('caudate_sequences', paste0(genome, 'CelltypeOnlyNonCelltype'),'rda', sep = '.'))
# save(rhesusCelltypeOnly_enhList, rhesusNonCelltype_enhList, file = save_fn )
