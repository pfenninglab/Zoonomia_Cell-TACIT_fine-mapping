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


############################
### get the human peaks ####
genome = 'hg38'
human_peak_fn = file.path('../../../data/raw_data/',genome,'Corces_2020', 'peak') %>%
  list.files(path = ., full.names = T, pattern = 'narrowPeak.gz')
human_peak_fn = human_peak_fn[ss(basename(human_peak_fn), '\\.') == 'Corces2020_caudate']
names(human_peak_fn) = basename(human_peak_fn) %>% ss('\\.', 2)
human_peak_fn = human_peak_fn[! names(human_peak_fn) %in% c('Consensus')]
human_peakList= lapply(human_peak_fn, rtracklayer::import) %>% GRangesList()

# annotate peaks, filter out exons, promoters
human_peakList = annotatePeaks(human_peakList, fromTSS = fromTSS, genome = genome)
human_enhList = filterPeaks(human_peakList, include = include)
human_enhList = mapply(transferColumn, toPeaks = human_enhList, colOut = 'col',
                       fromPeaks = human_enhList) %>% GRangesList()
lengths(human_enhList) / lengths(human_peakList) 

# summit-center and 
human_enhList = summitCenter(human_enhList, width = 501)
human_enhList_split = lapply(folds, function(fold){
  ret = lapply(human_enhList, splitPeakSet, testSet = testSet, validSet = fold, useCol = 'col')
  return(ret)
})

# export postive sequences to summit-centered 501bp.fasta file
split = names(human_enhList_split[[1]][[1]])
system(paste('mkdir -p',  file.path(PROJDIR, 'fasta')))
# write the positives w/o splits
pos.fasta_fn = file.path(PROJDIR, 'fasta', paste(genome, names(human_enhList), 'positive.fa.gz', sep = '_'))
posFasta = mapply(writeGRangesToFasta, gr = human_enhList, file = pos.fasta_fn, genome = genome)
for(cell in names(human_enhList)){
  for(fold in names(folds)){
    # write the positives
    pos.fasta_fn = file.path(PROJDIR, 'fasta', 
                             paste(genome, cell, fold, split, 'positive.fa.gz', sep = '_'))
    posFasta = mapply(writeGRangesToFasta, gr = human_enhList_split[[fold]][[cell]],  
                      file = pos.fasta_fn, genome = genome)
}}

system(paste('mkdir -p',  file.path(PROJDIR, 'rdas')))
save_fn = file.path(PROJDIR, 'rdas', paste('caudate_positive_sequences', genome,'.rda', sep = '.'))
save(human_peakList, human_enhList, file = save_fn )



############################
### get the mouse peaks ####
genome = 'mm10'
mouse_peak_fn = file.path('../../../data/raw_data/',genome) %>%
  list.files(path = ., full.names = T, recursive = T, pattern = 'narrowPeak.gz')
mouse_peak_fn = mouse_peak_fn[!grepl('HALPER.narrowPeak.gz|peaks_nonreproducible', mouse_peak_fn)]
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
lengths(mouse_enhList) / lengths(mouse_peakList)

### get the rough chromosome numbers for each peak
chainFile =file.path("/home/bnphan/resources/liftOver_chainz", paste0(genome,'ToHg38.over.chain'))
mouse_enhList = mapply(getLiftedChr, p = mouse_enhList, chainFile = chainFile) %>% GRangesList()

# summit-center  
mouse_enhList = summitCenter(mouse_enhList, width = 501)
mouse_enhList_split = lapply(folds, function(fold){
  ret = lapply(mouse_enhList, splitPeakSet, testSet = testSet, validSet = fold, useCol = 'col')
  return(ret)
})

# write the positives w/o splits
split = names(mouse_enhList_split[[1]][[1]])
pos.fasta_fn = file.path(PROJDIR, 'fasta', paste(genome, names(mouse_enhList), 'positive.fa.gz', sep = '_'))
posFasta = mapply(writeGRangesToFasta, gr = mouse_enhList, file = pos.fasta_fn, genome = genome)
for(cell in names(mouse_enhList)){
for(fold in names(folds)){
  pos.fasta_fn = file.path(PROJDIR, 'fasta', paste(genome, cell, fold, split, 'positive.fa.gz', sep = '_'))
  posFasta = mapply(writeGRangesToFasta, gr = mouse_enhList_split[[fold]][[cell]],  
                    file = pos.fasta_fn, genome = genome)
}}

system(paste('mkdir -p',  file.path(PROJDIR, 'rdas')))
save_fn = file.path(PROJDIR, 'rdas', paste('caudate_positive_sequences', genome,'.rda', sep = '.'))
save(mouse_peakList, mouse_enhList, file = save_fn )





#######################################################################
### get the macaque peaks, note narrowPeaks in rheMac8 coordinates ####
genome = 'rheMac10'
rhesus_fn =file.path('../../../data/raw_data/',genome) %>%
  list.files(path = ., full.names = T, recursive = T, pattern = 'narrowPeak.gz')
rhesus_fn = rhesus_fn[!grepl('HALPER.narrowPeak.gz', rhesus_fn)]
names(rhesus_fn) = basename(rhesus_fn) %>% ss('\\.', 2)
rhesus_fn = rhesus_fn[!grepl('Consensus', rhesus_fn)]
rhesus_fn = rhesus_fn[order(names(rhesus_fn))]
rhesus_peakList = lapply(rhesus_fn, rtracklayer::import) %>% GRangesList()
rhesus_peakList = lapply(rhesus_peakList, convertHalChrName, chrOut = 'UCSC', 
                         species = 'Macaca_mulatta') %>% GRangesList()

## lift rheMac8 coordinates to rheMac10
chainFile1 =file.path("/home/bnphan/resources/liftOver_chainz", paste0('rheMac8ToRheMac10.over.chain'))
rhesus_peakList = lapply(rhesus_peakList, liftOver_narrowPeak, chainFile = chainFile1) %>% GRangesList()

# annotate peaks, filter out exons, promoters
# need to convert genBank chr names to UCSC chr names
rhesus_peakList = annotatePeaks(rhesus_peakList, fromTSS = fromTSS, genome = genome)
rhesus_enhList = filterPeaks(rhesus_peakList, include = include)
lengths(rhesus_enhList) / lengths(rhesus_peakList)

## read in rhesus peaks mapped to hg38 to get chr for split
chainFile =file.path("/home/bnphan/resources/liftOver_chainz", paste0(genome,'ToHg38.over.chain'))
rhesus_enhList = mapply(getLiftedChr, p = rhesus_enhList, chainFile = chainFile) %>% GRangesList()

# summit-center  
rhesus_enhList = summitCenter(rhesus_enhList, width = 501)
library(BSgenome.Mmulatta.UCSC.rheMac10)
seqinfo(rhesus_enhList) = seqinfo(BSgenome.Mmulatta.UCSC.rheMac10)
rhesus_enhList = lapply(rhesus_enhList, trim) %>% GRangesList()
rhesus_enhList_split = lapply(folds, function(fold){
  ret = lapply(rhesus_enhList, splitPeakSet, testSet = testSet, validSet = fold, useCol = 'col')
  return(ret)
})

# write the positives w/o splits
split = names(rhesus_enhList_split[[1]][[1]])
pos.fasta_fn = file.path(PROJDIR, 'fasta', paste(genome, names(rhesus_enhList), 'positive.fa.gz', sep = '_'))
posFasta = mapply(writeGRangesToFasta, gr = rhesus_enhList, file = pos.fasta_fn, genome = genome)
for(cell in names(rhesus_enhList)){
  for(fold in names(folds)){
    pos.fasta_fn = file.path(PROJDIR, 'fasta', paste(genome, cell, fold, split, 'positive.fa.gz', sep = '_'))
    posFasta = mapply(writeGRangesToFasta, gr = rhesus_enhList_split[[fold]][[cell]],  
                      file = pos.fasta_fn, genome = genome)
  }}

system(paste('mkdir -p',  file.path(PROJDIR, 'rdas')))
save_fn = file.path(PROJDIR, 'rdas', paste('caudate_positive_sequences', genome,'.rda', sep = '.'))
save(rhesus_peakList, rhesus_enhList, file = save_fn )




