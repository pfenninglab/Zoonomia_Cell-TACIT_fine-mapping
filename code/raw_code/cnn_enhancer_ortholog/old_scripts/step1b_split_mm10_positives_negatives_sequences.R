# to be run in the root github directory
LABEL='cnn_enhancer_ortholog'
setwd('code/raw_code/cnn_enhancer_ortholog')
PROJDIR=file.path('../../../data/raw_data/',LABEL)

### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(rtracklayer))
library(tidyverse)
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

genome = 'mm10'

############################
### get the human peaks ####
mouse_peak_fn = file.path('../../../data/raw_data/',genome) %>%
  list.files(path = ., full.names = T, recursive = T, pattern = 'narrowPeak.gz')
mouse_peak_fn = mouse_peak_fn[!grepl('HALPER.narrowPeak.gz', mouse_peak_fn)]
# exclude BICCN_CP.MSN_D1 and BICCN_CP.MSN_D2 & halper peaks for now
mouse_peak_fn = mouse_peak_fn[!grepl('BICCN_CP.MSN_D|BICCN_CP.INT|Consensus', 
                                     mouse_peak_fn)]
mouse_peak_fn = mouse_peak_fn[ss(basename(mouse_peak_fn), '\\.') %in% 
                                c('BICCN_CP', 'Pfenning_Cpu')]
names(mouse_peak_fn) = basename(mouse_peak_fn) %>% ss('\\.', 2)
mouse_peak_fn = mouse_peak_fn[order(names(mouse_peak_fn))]
mouse_peakList= lapply(mouse_peak_fn, rtracklayer::import) %>% GRangesList()
mouse_peakList[c('INT_Pvalb','MSN_D1','MSN_D2')] = 
  lapply(mouse_peakList[c('INT_Pvalb','MSN_D1','MSN_D2')], function(gr){
    gr$name = paste0(seqnames(gr), ':',start(gr),'-',end(gr),':',gr$peak)
    return(gr)
  }) %>% GRangesList() 

# annotate peaks, filter out exons, promoters
mouse_peakList = annotatePeaks(mouse_peakList, fromTSS = fromTSS, genome = genome)
mouse_enhList = filterPeaks(mouse_peakList, include = include)
lengths(mouse_enhList) / lengths(mouse_peakList) ## 55% are noncoding

### get the mouse mapped to human peaks, used in splitting sets
mouse_hal2hg38_fn = mouse_peak_fn %>% 
  gsub(pattern = 'peak',replacement = 'halper') %>%
  gsub(pattern = '.narrowPeak.gz', replacement = '.Mus_musculusToHomo_sapiens.HALPER.narrowPeak.gz')
mouse_hal2hg38_peakList = lapply(mouse_hal2hg38_fn, rtracklayer::import) %>% GRangesList()

#tranfer hg38 chr to `col` column in the mouse_enhList
mouse_enhList = mapply(transferColumn, toPeaks = mouse_enhList, colOut = 'col',
                       fromPeaks = mouse_hal2hg38_peakList) %>% GRangesList()
mouse_enhList = mouse_enhList %>% lapply(as.data.frame) %>% 
  bind_rows(.id = 'cell_type') %>% arrange(seqnames, start, end, peak) %>% 
  fill(contains('col')) %>% group_by(cell_type) %>% group_split() %>%
  lapply(GRanges) %>% GRangesList()
names(mouse_enhList) = sapply(mouse_enhList, function(x) unique(x$cell_type))
table(hg38 = mouse_enhList[[2]]$col, mm10 = seqnames(mouse_enhList[[2]]))

############################################
### get the mouse mapped to human peaks ####
human_hal2mm10_fn =file.path('../../../data/raw_data/','hg38') %>%
  list.files(path = ., full.names = T, recursive = T, 
             pattern = 'ToMus_musculus.HALPER.narrowPeak.gz')
names(human_hal2mm10_fn) = basename(human_hal2mm10_fn) %>% ss('\\.', 2)
human_hal2mm10_fn = human_hal2mm10_fn[names(mouse_enhList)]
human_hal2mm10_peakList = lapply(human_hal2mm10_fn, rtracklayer::import) %>% GRangesList()

# annotate and keep intron, distal intergenic
# exclude mapped mouse peak overlapping human peak in same cell type
human_hal2mm10_peakList = annotatePeaks(human_hal2mm10_peakList, fromTSS = fromTSS, genome = genome)
human_hal2mm10_enhList = filterPeaks(human_hal2mm10_peakList, include = include)
human_nonEnhList = mapply(getNonEnhOrthPeaks, inPeaks = human_hal2mm10_enhList, 
                          excludePeaks = mouse_peakList)
human_nonEnhList = lapply(human_nonEnhList, function(gr){
  gr$col = ss(gr$name, ':',2)
  return(gr)
}) %>% GRangesList()

#############################################
### get the monkey mapped to human peaks ####
rhesus_hal2mm10_fn =file.path('../../../data/raw_data/','rheMac10') %>%
  list.files(path = ., full.names = T, recursive = T, 
             pattern = 'ToMus_musculus.HALPER.narrowPeak.gz')
names(rhesus_hal2mm10_fn) = basename(rhesus_hal2mm10_fn) %>% ss('\\.', 2)
rhesus_hal2mm10_fn = rhesus_hal2mm10_fn[names(mouse_enhList)]
rhesus_hal2mm10_peakList = lapply(rhesus_hal2mm10_fn, rtracklayer::import) %>% GRangesList()

# annotate and keep intron, distal intergenic
# exclude mapped mouse peak overlapping human peak in same cell type
rhesus_hal2mm10_peakList = annotatePeaks(rhesus_hal2mm10_peakList, fromTSS = fromTSS, genome = genome)
rhesus_hal2mm10_enhList = filterPeaks(rhesus_hal2mm10_peakList, include = include)
rhesus_nonEnhList = mapply(getNonEnhOrthPeaks, inPeaks = rhesus_hal2mm10_enhList, 
                          excludePeaks = mouse_peakList)

## read in rhesus peaks mapped to hg38 to get chr for split
rhesus_hal2hg38_fn = rhesus_hal2mm10_fn %>%
  gsub(pattern = 'ToMus_musculus.HALPER.narrowPeak.gz', 
       replacement = 'ToHomo_sapiens.HALPER.narrowPeak.gz')
rhesus_hal2hg3_peakList = lapply(rhesus_hal2hg38_fn, rtracklayer::import) %>% GRangesList()

#tranfer hg38 chr to `col` column in the rhesus_nonEnhList
rhesus_nonEnhList = mapply(transferColumn, toPeaks = rhesus_nonEnhList, colOut = 'col',
                       fromPeaks = rhesus_hal2hg3_peakList) %>% GRangesList()
table(hg38 = rhesus_nonEnhList[[1]]$col, mm10 = seqnames(rhesus_nonEnhList[[1]]))

#######################
# number of postives  #
lengths(mouse_enhList) / 1000
#  Astro INT_Pvalb Microglia    MSN_D1    MSN_D2    MSN_SN     Oligo       OPC 
# 18.413    57.491     5.502    75.701    82.129    11.000    13.834     8.449 

# combine sets of non-enh orthologs and save #
# number of non-enhancer orthologs
nonEnh_peakList = GRangesList(mapply(c, rhesus_nonEnhList, human_nonEnhList))
all.equal(lengths(rhesus_nonEnhList) + lengths(human_nonEnhList), lengths(nonEnh_peakList))
lengths(nonEnh_peakList) / 1000
#  Astro INT_Pvalb Microglia    MSN_D1    MSN_D2    MSN_SN     Oligo       OPC 
# 35.493    23.554    23.298    35.705    35.123    30.423    21.757    30.529 

## ratio of positives to non-enhancer orthologs
lengths(mouse_enhList) / lengths(nonEnh_peakList)
#     Astro INT_Pvalb Microglia    MSN_D1    MSN_D2    MSN_SN     Oligo       OPC 
# 0.5187784 2.4408168 0.2361576 2.1201792 2.3383253 0.3615686 0.6358413 0.2767533


###############################################
# split positives and negatives across folds
positiveSet = summitCenter(mouse_enhList, width = 501)
positiveSplit_list = lapply(folds, function(fold){
  ret = lapply(positiveSet, splitPeakSet, testSet = testSet, 
               validSet = fold, useCol = 'col')
  return(ret)
})

# export postive sequences to summit-centered 501bp fasta file
negativeSet = summitCenter(nonEnh_peakList, width = 501)
negativeSplit_list = lapply(folds, function(fold){
  ret = lapply(negativeSet, splitPeakSet, testSet = testSet, 
               validSet = fold, useCol = 'col')
  return(ret)
})


#############################################################
# export postive sequences to summit-centered 501bp fasta file
split = names(negativeSplit_list[[1]][[1]])
dir.create(file.path(PROJDIR, 'fasta'))
for(cell in names(positiveSet)){
  for(fold in names(folds)){
    # write the positives
    pos_fasta_fn = file.path(PROJDIR, 'fasta', 
                             paste(genome, cell, fold, split, 'positive.fa', sep = '_'))
    posFasta = mapply(writeGRangesToFasta, gr = positiveSplit_list[[fold]][[cell]],  
                      file = pos_fasta_fn, genome = genome)
    # write the negatives
    neg_fasta_fn = file.path(PROJDIR, 'fasta', 
                             paste(genome, cell, fold, split, 'negative.fa', sep = '_'))
    negFasta = mapply(writeGRangesToFasta, gr = negativeSplit_list[[fold]][[cell]],  
                      file = neg_fasta_fn, genome = genome)
  }}


system(paste('mkdir -p',  file.path(PROJDIR, 'rdas')))
save_fn = file.path(PROJDIR, 'rdas', paste('cnn_enhancer_non-enhancer_split',
                                           genome,'.rda', sep = '_'))
save(positiveSet, positiveSplit_list, 
     negativeSet, negativeSplit_list, file = save_fn )



