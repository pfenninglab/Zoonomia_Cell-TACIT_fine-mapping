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
             fold4 = c('chr9', 'chr12'),
             fold5 = c('chr8', 'chr10'))

genome = 'hg38'

############################################
# get the positives, and nonreproducible peaks
load( file.path(PROJDIR, 'rdas', paste('caudate_positive_sequences', genome,'.rda', sep = '.')) )

ARCHDIR=file.path('../../../data/raw_data/hg38/Corces_2020',
                  'ArchR_Corces2020_caudate_labeled','PeakCalls','ReplicateCalls')

nonRepPeak_list = lapply(names(human_peakList), function(cell){
  fn = list.files(path = ARCHDIR, pattern = cell, full.names = T)
  gr = lapply(fn, function(f) gr = readRDS(f)) %>%GRangesList() %>% unlist() 
  start(gr) = start(gr) -250; end(gr) = end(gr) + 250
  gr = GenomicRanges::reduce(gr)
  return(gr)
})
names(nonRepPeak_list) = names(human_enhList)


##########################################################
# gather nonCelltype peaks (enhancers in other cell types)
nonCelltype_peakList = lapply(names(human_enhList), function(name){
  # reproducible peaks in other cell types except target
  nonCelltype_gr = unlist(GRangesList(human_enhList[!names(human_enhList) %in% name]))
  # exclude peaks in non-reproducible of target cell type
  oo = findOverlaps(query = nonRepPeak_list[[name]], subject = nonCelltype_gr)
  ind = which(! seq_along(nonCelltype_gr) %in% subjectHits(oo))
  ret_gr = nonCelltype_gr[ind]
  return(ret_gr)
}) %>% GRangesList()
names(nonCelltype_peakList) = names(human_enhList)

lengths(nonCelltype_peakList) / 1000
#  Astro INT_Pvalb Microglia    MSN_D1    MSN_D2    MSN_SN     Oligo       OPC 
#  165.041 146.308 211.345 115.490 117.401 152.812 192.811 169.926

## ratio of positives to nonCelltype enhancers
lengths(human_enhList) / lengths(nonCelltype_peakList)
#     Astro INT_Pvalb Microglia    MSN_D1    MSN_D2    MSN_SN     Oligo       OPC 
# 0.2265679 0.2878175 0.1249521 0.4086847 0.3995026 0.2193872 0.1390636 0.1797724

# export postive sequences to summit-centered 501bp fasta file
negativeSet = summitCenter(nonCelltype_peakList, width = 501)
negativeSplit_list = lapply(folds, function(fold) {
  ret = lapply(negativeSet, splitPeakSet, testSet = testSet, validSet = fold)
  return(ret)
})


#############################################################
# export postive sequences to summit-centered 501bp fasta file
split = names(negativeSplit_list[[1]][[1]])
system(paste('mkdir -p',  file.path(PROJDIR, 'fasta')))
for(cell in names(negativeSet)){
  for(fold in names(folds)){
    # write the negatives
    neg_fasta_fn = file.path(PROJDIR, 'fasta', 
                             paste(genome, cell, fold, split, 'nonCelltypeNeg.fa', sep = '_'))
    negFasta = mapply(writeGRangesToFasta, gr = negativeSplit_list[[fold]][[cell]],  
                      file = neg_fasta_fn, genome = genome)
  }}

system(paste('mkdir -p',  file.path(PROJDIR, 'rdas')))
save_fn = file.path(PROJDIR, 'rdas', paste('caudate_nonCelltype_sequences', genome,'.rda', sep = '.'))
save(negativeSet, negativeSplit_list, file = save_fn )


