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

######################
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
names(nonRepPeak_list) = names(human_peakList)
human_pos_list = mapply(function(gr1, gr2){
  gr = GRangesList(list(gr1, gr2)) %>% unlist() %>%
    GenomicRanges::reduce()
  return(gr)
}, gr1 = human_peakList, gr2 = nonRepPeak_list)


############################################
### get the mouse mapped to human peaks ####
mouse_hal2hg38_fn =file.path('../../../data/raw_data/','mm10') %>%
  list.files(path = ., full.names = T, recursive = T, 
             pattern = 'ToHomo_sapiens.HALPER.narrowPeak.gz')
mouse_hal2hg38_fn = mouse_hal2hg38_fn[ss(basename(mouse_hal2hg38_fn), '\\.') %in% 
                                        c('BICCN_CP', 'Pfenning_Cpu')]
# exclude BICCN_CP.MSN_D1 and BICCN_CP.MSN_D2 & halper peaks for now
mouse_hal2hg38_fn = mouse_hal2hg38_fn[!grepl('BICCN_CP.MSN_D|BICCN_CP.INT', mouse_hal2hg38_fn)]
names(mouse_hal2hg38_fn) = basename(mouse_hal2hg38_fn) %>% ss('\\.', 2)
mouse_hal2hg38_fn = mouse_hal2hg38_fn[names(human_peakList)]
mouse_hal2hg38_peakList = lapply(mouse_hal2hg38_fn, rtracklayer::import) %>% GRangesList()

# annotate and keep intron, distal intergenic
# exclude mapped mouse peak overlapping human peak in same cell type
mouse_hal2hg38_peakList = annotatePeaks(mouse_hal2hg38_peakList, fromTSS = fromTSS, genome = genome)
mouse_hal2hg38_enhList = filterPeaks(mouse_hal2hg38_peakList, include = include)
mouse_nonEnhList = mapply(getNonEnhOrthPeaks, inPeaks = mouse_hal2hg38_enhList, 
                          excludePeaks = human_pos_list)


#############################################
### get the monkey mapped to human peaks ####
rhesus_hal2hg38_fn =file.path('../../../data/raw_data/','rheMac10') %>%
  list.files(path = ., full.names = T, recursive = T, 
             pattern = 'ToHomo_sapiens.HALPER.narrowPeak.gz')
names(rhesus_hal2hg38_fn) = basename(rhesus_hal2hg38_fn) %>% ss('\\.', 2)
rhesus_hal2hg38_fn = rhesus_hal2hg38_fn[names(human_peakList)]
rhesus_hal2hg38_peakList = lapply(rhesus_hal2hg38_fn, rtracklayer::import) %>% GRangesList()

# annotate and keep intron, distal intergenic
# exclude mapped mouse peak overlapping human peak in same cell type
rhesus_hal2hg38_peakList = annotatePeaks(rhesus_hal2hg38_peakList, fromTSS = fromTSS, genome = genome)
rhesus_hal2hg38_enhList = filterPeaks(rhesus_hal2hg38_peakList, include = include)
rhesus_nonEnhList = mapply(getNonEnhOrthPeaks, inPeaks = rhesus_hal2hg38_enhList, 
                           excludePeaks = human_pos_list)


###################################
# number of non-enhancer orthologs
nonEnh_peakList = GRangesList(mapply(c, rhesus_nonEnhList, mouse_nonEnhList))
all.equal(lengths(rhesus_nonEnhList) + lengths(mouse_nonEnhList), lengths(nonEnh_peakList))
lengths(nonEnh_peakList) / 1000
#  Astro INT_Pvalb Microglia    MSN_D1    MSN_D2    MSN_SN     Oligo       OPC 
# 17.952    22.956    15.263    41.330    42.594    13.914    10.984    11.774 

## ratio of positives to non-enhancer orthologs
lengths(human_enhList) / lengths(nonEnh_peakList)
#    Astro INT_Pvalb Microglia    MSN_D1    MSN_D2    MSN_SN     Oligo       OPC 
# 2.082943  1.834379  1.730197  1.142003  1.101141  2.409444  2.441096  2.594530 

# export postive sequences to summit-centered 501bp fasta file
negativeSet = summitCenter(nonEnh_peakList, width = 501)
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
                             paste(genome, cell, fold, split, 'nonEnhNeg.fa', sep = '_'))
    negFasta = mapply(writeGRangesToFasta, gr = negativeSplit_list[[fold]][[cell]],  
                      file = neg_fasta_fn, genome = genome)
  }}

system(paste('mkdir -p',  file.path(PROJDIR, 'rdas')))
save_fn = file.path(PROJDIR, 'rdas', paste('caudate_nonEnh_sequences', genome,'.rda', sep = '.'))
save( human_pos_list, negativeSet, negativeSplit_list, file = save_fn )


