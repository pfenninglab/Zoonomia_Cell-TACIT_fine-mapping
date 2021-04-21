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

genome = 'rheMac10' # do everything in rheMac8 coordinates, export in rheMac10
chainFile1 =file.path("/home/bnphan/resources/liftOver_chainz", 
                      paste0('rheMac8ToRheMac10.over.chain'))

############################################
# get the positives, and nonreproducible peaks
load( file.path(PROJDIR, 'rdas', paste('caudate_positive_sequences', genome,'.rda', sep = '.')) )

ARCHDIR=file.path('../../../data/raw_data/rheMac10/Stauffer_caudate',
                  'ArchR_Stauffer_caudate_labeled','PeakCalls','ReplicateCalls')
nonRepPeak_list = lapply(names(rhesus_enhList), function(cell){
  fn = list.files(path = ARCHDIR, pattern = cell, full.names = T)
  gr = lapply(fn, function(f) gr = readRDS(f)) %>%GRangesList() %>% unlist() 
  start(gr) = start(gr) -250; end(gr) = end(gr) + 250
  gr = GenomicRanges::reduce(gr)
  return(gr)
})
names(nonRepPeak_list) = names(rhesus_enhList)


##########################################################
# gather nonCelltype peaks (enhancers in other cell types)
nonCelltype_peakList = lapply(names(rhesus_enhList), function(name){
  # reproducible peaks in other cell types except target
  nonCelltype_gr = unlist(GRangesList(rhesus_enhList[!names(rhesus_enhList) %in% name]))
  # exclude peaks in non-reproducible of target cell type
  oo = findOverlaps(query = nonRepPeak_list[[name]], subject = nonCelltype_gr)
  ind = which(! seq_along(nonCelltype_gr) %in% subjectHits(oo))
  ret_gr = nonCelltype_gr[ind]
  return(ret_gr)
}) %>% GRangesList()
names(nonCelltype_peakList) = names(rhesus_enhList)

lengths(nonCelltype_peakList) / 1000
#  Astro INT_Pvalb Microglia    MSN_D1    MSN_D2    MSN_SN     Oligo       OPC 
# 151.706 159.605 192.131  96.914  97.481 130.438 184.720 160.388

## ratio of positives to non-enhancer orthologs
lengths(rhesus_enhList) / lengths(nonCelltype_peakList)
#     Astro INT_Pvalb Microglia    MSN_D1    MSN_D2    MSN_SN     Oligo       OPC 
# 0.2723030 0.1035870 0.1721794 0.5232474 0.5206040 0.2267744 0.1231756 0.1817343

# split non-enhancer orthologs
chainFile =file.path("/home/bnphan/resources/liftOver_chainz", paste0(genome,'ToHg38.over.chain'))
nonCelltype_peakList = mapply(getLiftedChr, p = nonCelltype_peakList, chainFile = chainFile) %>% GRangesList()
negativeSet = summitCenter(nonCelltype_peakList, width = 501)
negativeSplit_list = lapply(folds, function(fold){
  lapply(negativeSet, splitPeakSet, testSet = testSet, validSet = fold, useCol = 'col')
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
save( negativeSet, negativeSplit_list, file = save_fn )




