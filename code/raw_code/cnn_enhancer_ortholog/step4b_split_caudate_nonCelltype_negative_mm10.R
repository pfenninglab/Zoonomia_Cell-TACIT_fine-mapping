# to be run in the root github directory
LABEL='cnn_enhancer_ortholog'
setwd('code/raw_code/cnn_enhancer_ortholog')
PROJDIR=file.path('../../../data/raw_data',LABEL)

### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(rtracklayer))
suppressMessages(library(tidyverse))
suppressMessages(library(ArchR))
source('../hal_scripts/narrowPeakFunctions.R')
source('../hal_scripts/gen_enh_ortholog_sets.R')

## parameters for annotation and filtering
fromTSS = c(-20000,20000)
include = c('Distal.Intergenic','Intron')

## chromosomal splits
testSet = c('chr1','chr2')
folds = list(fold1 = c('chr6', 'chr13', 'chr21'),
             fold2 = c('chr7', 'chr14', 'chr18'),
             fold3 = c('chr11','chr17', 'chrX'),
             fold4 = c('chr9', 'chr12'), fold5 = c('chr8', 'chr10'))

genome = 'mm10'

##############################################
# get the positives, and nonreproducible peaks
load( file.path(PROJDIR, 'rdas', paste('caudate_positive_sequences', genome,'.rda', sep = '.')) )
ARCHDIR=file.path('../../../data/raw_data/mm10/BICCN_mouse_caudoputamen',
                  'ArchR_BICCN_CP_labeled','PeakCalls','ReplicateCalls')
nonRepPeak_list = lapply(names(mouse_peakList), function(cell){
  fn = list.files(path = ARCHDIR, pattern = cell, full.names = T)
  gr = lapply(fn, function(f) gr = readRDS(f)) %>%GRangesList() %>% unlist() 
  start(gr) = start(gr) -250; end(gr) = end(gr) + 250
  gr = GenomicRanges::reduce(gr)
  return(gr)
})
names(nonRepPeak_list) = names(mouse_peakList)

# cSNAIL pooled rep peaks
cSNAIL_celltypes = c('INT_Pvalb','MSN_D1','MSN_D2')
nonRepPeak_list[cSNAIL_celltypes] = 
  lapply(cSNAIL_celltypes, function(cell){
    fn = list.files(path =file.path('../../../data/raw_data/mm10', 'Mouse_cSNAIL_D1D2',
                      'peaks_nonreproducible'), pattern = cell, full.names = T)
    gr = import(fn)
  })



##########################################################
# gather nonCelltype peaks (enhancers in other cell types)
nonCelltype_peakList = lapply(names(mouse_enhList), function(name){
  # reproducible peaks in other cell types except target
  nonCelltype_gr = unlist(GRangesList(mouse_enhList[!names(mouse_enhList) %in% name]))
  # exclude peaks in non-reproducible of target cell type
  oo = findOverlaps(query = nonRepPeak_list[[name]], subject = nonCelltype_gr)
  ind = which(! seq_along(nonCelltype_gr) %in% subjectHits(oo))
  ret_gr = nonCelltype_gr[ind]
  return(ret_gr)
}) %>% GRangesList()
names(nonCelltype_peakList) = names(mouse_enhList)

lengths(nonCelltype_peakList) / 1000
#  Astro INT_Pvalb Microglia    MSN_D1    MSN_D2    MSN_SN     Oligo       OPC 
#  196.011  31.762   250.448    40.762    54.527   189.537  214.214     214.419

## ratio of positives to nonCelltype enhancers
lengths(mouse_enhList) / lengths(nonCelltype_peakList)
#     Astro INT_Pvalb Microglia    MSN_D1    MSN_D2    MSN_SN     Oligo       OPC 
# 0.09393861 1.81005604 0.02196863 1.85714636 1.50620793 0.05803616 0.06458028  0.03940416

# split non-enhancer orthologs
chainFile =file.path("/home/bnphan/resources/liftOver_chainz", paste0(genome,'ToHg38.over.chain'))
nonCelltype_peakList = mapply(getLiftedChr, p = nonCelltype_peakList, chainFile = chainFile) %>% GRangesList()
negativeSet = nonCelltype_peakList
negativeSplit_list = lapply(folds, function(fold){
  lapply(negativeSet, splitPeakSet, testSet = testSet, validSet = fold, useCol = 'col')
})



#############################################################
# export postive sequences to summit-centered 501bp.fasta file
split = names(negativeSplit_list[[1]][[1]])
system(paste('mkdir -p',  file.path(PROJDIR, 'fasta')))
for(cell in names(negativeSet)){
  for(fold in names(folds)){
    # write the negatives
    neg.fasta_fn = file.path(PROJDIR, 'fasta', 
                             paste(genome, cell, fold, split, 'nonCelltypeNeg.fa.gz', sep = '_'))
    negFasta = mapply(writeGRangesToFasta, gr = negativeSplit_list[[fold]][[cell]],  
                      file = neg.fasta_fn, genome = genome)
  }}

system(paste('mkdir -p',  file.path(PROJDIR, 'rdas')))
save_fn = file.path(PROJDIR, 'rdas', paste('caudate_nonCelltype_sequences', genome,'.rda', sep = '.'))
save( negativeSet, negativeSplit_list, file = save_fn )



