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

######################
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

# combine w/ reproducible peaks
mouse_pos_list = mapply(function(gr1, gr2){
  gr = GRangesList(list(gr1, gr2)) %>% unlist() %>%
    GenomicRanges::reduce()
  return(gr)
}, gr1 = mouse_peakList, gr2 = nonRepPeak_list)

############################################
### get the human mapped to mouse peaks ####
human_hal2mm10_fn =file.path('../../../data/raw_data/','hg38') %>%
  list.files(path = ., full.names = T, recursive = T, 
             pattern = 'ToMus_musculus.HALPER.narrowPeak.gz')
names(human_hal2mm10_fn) = basename(human_hal2mm10_fn) %>% ss('\\.', 2)
human_hal2mm10_fn = human_hal2mm10_fn[names(mouse_pos_list)]
human_hal2mm10_peakList = lapply(human_hal2mm10_fn, rtracklayer::import) %>% GRangesList()

# annotate and keep intron, distal intergenic
# exclude mapped mouse peak overlapping human peak in same cell type
human_hal2mm10_peakList = annotatePeaks(human_hal2mm10_peakList, fromTSS = fromTSS, genome = genome)
human_hal2mm10_enhList = filterPeaks(human_hal2mm10_peakList, include = include)
human_nonEnhList = mapply(getNonEnhOrthPeaks, inPeaks = human_hal2mm10_enhList, 
                          excludePeaks = mouse_pos_list)

#############################################
### get the monkey mapped to mouse peaks ####
rhesus_hal2mm10_fn =file.path('../../../data/raw_data/','rheMac10') %>%
  list.files(path = ., full.names = T, recursive = T, 
             pattern = 'ToMus_musculus.HALPER.narrowPeak.gz')
names(rhesus_hal2mm10_fn) = basename(rhesus_hal2mm10_fn) %>% ss('\\.', 2)
rhesus_hal2mm10_fn = rhesus_hal2mm10_fn[names(mouse_pos_list)]
rhesus_hal2mm10_peakList = lapply(rhesus_hal2mm10_fn, rtracklayer::import) %>% GRangesList()

# annotate and keep intron, distal intergenic
# exclude mapped mouse peak overlapping human peak in same cell type
rhesus_hal2mm10_peakList = annotatePeaks(rhesus_hal2mm10_peakList, fromTSS = fromTSS, genome = genome)
rhesus_hal2mm10_enhList = filterPeaks(rhesus_hal2mm10_peakList, include = include)
rhesus_nonEnhList = mapply(getNonEnhOrthPeaks, inPeaks = rhesus_hal2mm10_enhList, 
                          excludePeaks = mouse_pos_list)


##############################################
# combine sets of non-enh orthologs and save #
# number of non-enhancer orthologs
nonEnh_peakList = GRangesList(mapply(c, rhesus_nonEnhList, human_nonEnhList))
all.equal(lengths(rhesus_nonEnhList) + lengths(human_nonEnhList), lengths(nonEnh_peakList))
lengths(nonEnh_peakList) / 1000
#  Astro INT_Pvalb Microglia    MSN_D1    MSN_D2    MSN_SN     Oligo       OPC 
# 30.240    18.920    21.636    29.357    29.951    26.325    19.015    26.243 

## ratio of positives to non-enhancer orthologs
lengths(mouse_enhList) / lengths(nonEnh_peakList)
#     Astro INT_Pvalb Microglia    MSN_D1    MSN_D2    MSN_SN     Oligo       OPC 
# 0.6088955 3.0386364 0.2542984 2.5786354 2.7421121 0.4178538 0.7275309 0.3219525 

# get the hg38 chromosome for splitting
chainFile =file.path("/home/bnphan/resources/liftOver_chainz", paste0(genome,'ToHg38.over.chain'))
nonEnh_peakList = mapply(getLiftedChr, p = nonEnh_peakList, chainFile = chainFile) %>% GRangesList()

# export postive sequences to summit-centered 501bp fasta file
negativeSet = summitCenter(nonEnh_peakList, width = 501)
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
                             paste(genome, cell, fold, split, 'nonEnhNeg.fa', sep = '_'))
    negFasta = mapply(writeGRangesToFasta, gr = negativeSplit_list[[fold]][[cell]],  
                      file = neg_fasta_fn, genome = genome)
  }}

system(paste('mkdir -p',  file.path(PROJDIR, 'rdas')))
save_fn = file.path(PROJDIR, 'rdas', paste('caudate_nonEnh_sequences', genome,'.rda', sep = '.'))
save( mouse_pos_list, negativeSet, negativeSplit_list, file = save_fn )



