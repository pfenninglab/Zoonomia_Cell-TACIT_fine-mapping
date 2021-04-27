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

######################
# get the positives, and nonreproducible peaks
load( file.path(PROJDIR, 'rdas', paste('caudate_positive_sequences', genome,'.rda', sep = '.')) )

ARCHDIR=file.path('../../../data/raw_data/rheMac10/Stauffer_caudate',
                  'ArchR_Stauffer_caudate_labeled','PeakCalls','ReplicateCalls')
nonRepPeak_list = lapply(names(rhesus_pos_list), function(cell){
  fn = list.files(path = ARCHDIR, pattern = cell, full.names = T)
  gr = lapply(fn, function(f) gr = readRDS(f)) %>%GRangesList() %>% unlist() 
  start(gr) = start(gr) -250; end(gr) = end(gr) + 250
  gr = GenomicRanges::reduce(gr)
  return(gr)
})
names(nonRepPeak_list) = names(rhesus_pos_list)

# combine w/ reproducible peaks
rhesus_pos_list = mapply(function(gr1, gr2){
  gr = GRangesList(list(gr1, gr2)) %>% unlist() %>%
    GenomicRanges::reduce()
  return(gr)
}, gr1 = rhesus_pos_list, gr2 = nonRepPeak_list)

############################################
### get the mouse mapped to monkey peaks ####
mouse_hal2rm8_fn =file.path('../../../data/raw_data/','mm10') %>%
  list.files(path = ., full.names = T, recursive = T, 
             pattern = 'ToMacaca_mulatta.HALPER.narrowPeak.gz')
mouse_hal2rm8_fn = mouse_hal2rm8_fn[ss(basename(mouse_hal2rm8_fn), '\\.') %in% 
                                c('BICCN_CP', 'Pfenning_Cpu')]
# exclude BICCN_CP.MSN_D1 and BICCN_CP.MSN_D2 & halper peaks for now
mouse_hal2rm8_fn = mouse_hal2rm8_fn[!grepl('BICCN_CP.MSN_D|BICCN_CP.INT', mouse_hal2rm8_fn)]
names(mouse_hal2rm8_fn) = basename(mouse_hal2rm8_fn) %>% ss('\\.', 2)
mouse_hal2rm8_fn = mouse_hal2rm8_fn[names(rhesus_enhList)]
mouse_hal2rm8_peakList = lapply(mouse_hal2rm8_fn, rtracklayer::import) %>% GRangesList()
mouse_hal2rm8_peakList = lapply(mouse_hal2rm8_peakList, convertHalChrName, 
                                chrOut = 'UCSC', species = 'Macaca_mulatta') %>% GRangesList()

# annotate and keep intron, distal intergenic in rheMac10
# exclude mapped mouse peak overlapping human peak in same cell type
mouse_hal2rm10_peakList = lapply(mouse_hal2rm8_peakList, liftOver_narrowPeak, 
                                chainFile = chainFile1) %>% GRangesList()
mouse_hal2rm10_peakList = annotatePeaks(mouse_hal2rm10_peakList, fromTSS = fromTSS, genome = genome)
mouse_hal2rm10_enhList = filterPeaks(mouse_hal2rm10_peakList, include = include)
mouse_nonEnhList = mapply(getNonEnhOrthPeaks, inPeaks = mouse_hal2rm10_enhList, 
                          excludePeaks = rhesus_pos_list)
lengths(mouse_hal2rm10_enhList) / lengths(mouse_hal2rm10_peakList)
lengths(mouse_nonEnhList) / lengths(mouse_hal2rm10_peakList)


#############################################
### get the human mapped to monkey peaks ####
human_hal2rm8_fn =file.path('../../../data/raw_data/','hg38') %>%
  list.files(path = ., full.names = T, recursive = T, 
             pattern = 'ToMacaca_mulatta.HALPER.narrowPeak.gz')
names(human_hal2rm8_fn) = basename(human_hal2rm8_fn) %>% ss('\\.', 2)
human_hal2rm8_fn = human_hal2rm8_fn[names(rhesus_enhList)]
human_hal2rm8_peakList = lapply(human_hal2rm8_fn, rtracklayer::import) %>% 
  lapply(convertHalChrName, chrOut = 'UCSC', species = 'Macaca_mulatta') %>% GRangesList()

# annotate and keep intron, distal intergenic
# exclude mapped mouse peak overlapping human peak in same cell type
human_hal2rm10_peakList = lapply(human_hal2rm8_peakList, liftOver_narrowPeak, 
                                 chainFile = chainFile1) %>% GRangesList()
human_hal2rm10_peakList = annotatePeaks(human_hal2rm10_peakList, fromTSS = fromTSS, genome = genome)
human_hal2rm10_enhList = filterPeaks(human_hal2rm10_peakList, include = include)
human_nonEnhList = mapply(getNonEnhOrthPeaks, inPeaks = human_hal2rm10_enhList, 
                          excludePeaks = rhesus_pos_list)
lengths(human_hal2rm10_enhList) / lengths(human_hal2rm10_peakList)
lengths(human_nonEnhList) / lengths(human_hal2rm10_peakList)


###################################
# number of non-enhancer orthologs
nonEnh_peakList = GRangesList(mapply(c, human_nonEnhList, mouse_nonEnhList))
lengths(nonEnh_peakList) / 1000
#  Astro INT_Pvalb Microglia    MSN_D1    MSN_D2    MSN_SN     Oligo       OPC 
# 15.905    34.485    10.650    34.660    35.357    13.789    15.710    11.861

## ratio of positives to non-enhancer orthologs
lengths(rhesus_enhList) / lengths(nonEnh_peakList)
#     Astro INT_Pvalb Microglia    MSN_D1    MSN_D2    MSN_SN     Oligo       OPC 
# 2.5972964 0.4794258 3.1061972 1.4630698 1.4353311 2.1451882 1.4483132 2.4574656

# split non-enhancer orthologs
chainFile =file.path("/home/bnphan/resources/liftOver_chainz", paste0(genome,'ToHg38.over.chain'))
nonEnh_peakList = mapply(getLiftedChr, p = nonEnh_peakList, chainFile = chainFile) %>% GRangesList()
negativeSet = summitCenter(nonEnh_peakList, width = 501)
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
                             paste(genome, cell, fold, split, 'nonEnhNeg.fa.gz', sep = '_'))
    negFasta = mapply(writeGRangesToFasta, gr = negativeSplit_list[[fold]][[cell]],  
                      file = neg.fasta_fn, genome = genome)
  }}

system(paste('mkdir -p',  file.path(PROJDIR, 'rdas')))
save_fn = file.path(PROJDIR, 'rdas', paste('caudate_nonEnh_sequences', genome,'.rda', sep = '.'))
save( rhesus_pos_list, negativeSet, negativeSplit_list, file = save_fn )




