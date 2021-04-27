# to be run in the root github directory
LABEL='cnn_enhancer_ortholog'
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

### get the mouse mapped to human peaks ####
mouse_hal2hg38_fn =here('data/raw_data/','mm10') %>%
  list.files(path = ., full.names = T, recursive = T, 
             pattern = 'ToHomo_sapiens.HALPER.narrowPeak.gz')
mouse_hal2hg38_fn = mouse_hal2hg38_fn[ss(basename(mouse_hal2hg38_fn), '\\.') %in% 
                                        c('BICCN_CP', 'Pfenning_Cpu')]
# exclude BICCN_CP.MSN_D1 and BICCN_CP.MSN_D2 & halper peaks for now
mouse_hal2hg38_fn = mouse_hal2hg38_fn[!grepl('BICCN_CP.MSN_D|BICCN_CP.INT', mouse_hal2hg38_fn)]
names(mouse_hal2hg38_fn) = basename(mouse_hal2hg38_fn) %>% ss('\\.', 2)
mouse_hal2hg38_fn = mouse_hal2hg38_fn[names(human_peakList)]
mouse_hal2hg38_peakList = lapply(mouse_hal2hg38_fn, rtracklayer::import) %>% GRangesList()

### get the monkey mapped to human peaks ####
rhesus_hal2hg38_fn =here('data/raw_data/','rheMac10') %>%
  list.files(path = ., full.names = T, recursive = T, 
             pattern = 'ToHomo_sapiens.HALPER.narrowPeak.gz')
names(rhesus_hal2hg38_fn) = basename(rhesus_hal2hg38_fn) %>% ss('\\.', 2)
rhesus_hal2hg38_fn = rhesus_hal2hg38_fn[names(human_peakList)]
rhesus_hal2hg38_peakList = lapply(rhesus_hal2hg38_fn, rtracklayer::import) %>% GRangesList()

### find human peaks not overlapping mouse/monkey peaks ####
humanOnly_enhList = mapply(function(keep, excl1, excl2){
  excl_gr = c(excl1, excl2)
  oo = findOverlaps(query = keep, subject = excl_gr)
  keepIdx = which(!seq_along(keep) %in% queryHits(oo))
  return(keep[keepIdx])
}, keep = human_enhList, excl1 = mouse_hal2hg38_peakList, excl2 = rhesus_hal2hg38_peakList) %>%
  GRangesList()

humanOnly_enhList_split = lapply(folds, function(fold){
  ret = lapply(humanOnly_enhList, splitPeakSet, testSet = testSet, validSet = fold, useCol = 'col')
  return(ret)
})

# export postive sequences to summit-centered 501bp.fasta file
system(paste('mkdir -p',  file.path(PROJDIR, 'fasta')))
for(cell in names(human_enhList)){
  for(fold in names(folds)){
    # write the positives
    pos.fasta_fn = file.path(PROJDIR, 'fasta', 
                             paste(paste0(genome, 'Only'), cell, fold, split, 'positive.fa.gz', sep = '_'))
    posFasta = mapply(writeGRangesToFasta, gr = humanOnly_enhList_split[[fold]][[cell]][split],  
                      file = pos.fasta_fn, genome = genome)
  }}

save_fn = file.path(PROJDIR, 'rdas', paste('caudate_positive_sequences', paste0(genome, 'Only'),'rda', sep = '.'))
save(humanOnly_enhList, file = save_fn )


###########################################
# load in the mouse positive peak files
genome = 'mm10'
load(here(DATADIR, 'rdas', paste('caudate_positive_sequences', 'mm10','.rda', sep = '.')))

### get the human mapped to mouse peaks ####
human_hal2mm10_fn =file.path('data/raw_data/','hg38') %>%
  list.files(path = ., full.names = T, recursive = T, 
             pattern = 'ToMus_musculus.HALPER.narrowPeak.gz')
names(human_hal2mm10_fn) = basename(human_hal2mm10_fn) %>% ss('\\.', 2)
human_hal2mm10_fn = human_hal2mm10_fn[names(mouse_peakList)]
human_hal2mm10_peakList = lapply(human_hal2mm10_fn, rtracklayer::import) %>% GRangesList()

### get the monkey mapped to mouse peaks ####
rhesus_hal2mm10_fn =file.path('data/raw_data/','rheMac10') %>%
  list.files(path = ., full.names = T, recursive = T, 
             pattern = 'ToMus_musculus.HALPER.narrowPeak.gz')
names(rhesus_hal2mm10_fn) = basename(rhesus_hal2mm10_fn) %>% ss('\\.', 2)
rhesus_hal2mm10_fn = rhesus_hal2mm10_fn[names(mouse_peakList)]
rhesus_hal2mm10_peakList = lapply(rhesus_hal2mm10_fn, rtracklayer::import) %>% GRangesList()

### find mouse peaks not overlapping human/monkey peaks ####
mouseOnly_enhList = mapply(function(keep, excl1, excl2){
  excl_gr = c(excl1, excl2)
  oo = findOverlaps(query = keep, subject = excl_gr)
  keepIdx = which(!seq_along(keep) %in% queryHits(oo))
  return(keep[keepIdx])
}, keep = mouse_enhList, excl1 = human_hal2mm10_peakList, excl2 = rhesus_hal2mm10_peakList) %>%
  GRangesList()

mouseOnly_enhList_split = lapply(folds, function(fold){
  ret = lapply(mouseOnly_enhList, splitPeakSet, testSet = testSet, validSet = fold, useCol = 'col')
  return(ret)
})

# export postive sequences to summit-centered 501bp.fasta file
for(cell in names(mouseOnly_enhList)){
  for(fold in names(folds)){
    # write the positives
    pos.fasta_fn = file.path(PROJDIR, 'fasta', 
                             paste(paste0(genome, 'Only'), cell, fold, split, 'positive.fa.gz', sep = '_'))
    posFasta = mapply(writeGRangesToFasta, gr = mouseOnly_enhList_split[[fold]][[cell]][split],  
                      file = pos.fasta_fn, genome = genome)
  }}

save_fn = file.path(PROJDIR, 'rdas', paste('caudate_positive_sequences', paste0(genome, 'Only'),'rda', sep = '.'))
save(mouseOnly_enhList, file = save_fn )







###########################################
# load in the monkey positive peak files
genome = 'rheMac10'
chainFile1 =file.path("/home/bnphan/resources/liftOver_chainz", 
                      paste0('rheMac8ToRheMac10.over.chain'))

load(here(DATADIR, 'rdas', paste('caudate_positive_sequences', 'rheMac10','.rda', sep = '.')))

### get the mouse mapped to monkey peaks ####
mouse_hal2rm8_fn =file.path('data/raw_data/','mm10') %>%
  list.files(path = ., full.names = T, recursive = T, 
             pattern = 'ToMacaca_mulatta.HALPER.narrowPeak.gz')
mouse_hal2rm8_fn = mouse_hal2rm8_fn[ss(basename(mouse_hal2rm8_fn), '\\.') %in% 
                                      c('BICCN_CP', 'Pfenning_Cpu')]
mouse_hal2rm8_fn = mouse_hal2rm8_fn[!grepl('BICCN_CP.MSN_D|BICCN_CP.INT', mouse_hal2rm8_fn)]
names(mouse_hal2rm8_fn) = basename(mouse_hal2rm8_fn) %>% ss('\\.', 2)
mouse_hal2rm8_fn = mouse_hal2rm8_fn[names(rhesus_enhList)]
mouse_hal2rm8_peakList = lapply(mouse_hal2rm8_fn, rtracklayer::import) %>% GRangesList()
mouse_hal2rm8_peakList = lapply(mouse_hal2rm8_peakList, convertHalChrName, 
                                chrOut = 'UCSC', species = 'Macaca_mulatta') %>% GRangesList()
mouse_hal2rm10_peakList = lapply(mouse_hal2rm8_peakList, liftOver_narrowPeak, 
                                 chainFile = chainFile1) %>% GRangesList()


### get the human mapped to monkey peaks ####
human_hal2rm8_fn =file.path('data/raw_data/','hg38') %>%
  list.files(path = ., full.names = T, recursive = T, 
             pattern = 'ToMacaca_mulatta.HALPER.narrowPeak.gz')
names(human_hal2rm8_fn) = basename(human_hal2rm8_fn) %>% ss('\\.', 2)
human_hal2rm8_fn = human_hal2rm8_fn[names(rhesus_enhList)]
human_hal2rm8_peakList = lapply(human_hal2rm8_fn, rtracklayer::import) %>% 
  lapply(convertHalChrName, chrOut = 'UCSC', species = 'Macaca_mulatta') %>% GRangesList()
human_hal2rm10_peakList = lapply(human_hal2rm8_peakList, liftOver_narrowPeak, 
                                 chainFile = chainFile1) %>% GRangesList()


### find monkey peaks not overlapping human/mouse peaks ####
rhesusOnly_enhList = mapply(function(keep, excl1, excl2){
  excl_gr = c(excl1, excl2)
  oo = findOverlaps(query = keep, subject = excl_gr)
  keepIdx = which(!seq_along(keep) %in% queryHits(oo))
  return(keep[keepIdx])
}, keep = rhesus_enhList, excl1 = human_hal2rm10_peakList, excl2 = mouse_hal2rm10_peakList) %>%
  GRangesList()

rhesusOnly_enhList_split = lapply(folds, function(fold){
  ret = lapply(rhesusOnly_enhList, splitPeakSet, testSet = testSet, validSet = fold, useCol = 'col')
  return(ret)
})

# export postive sequences to summit-centered 501bp.fasta file
for(cell in names(rhesusOnly_enhList)){
  for(fold in names(folds)){
    # write the positives
    pos.fasta_fn = file.path(PROJDIR, 'fasta', 
                             paste(paste0(genome, 'Only'), cell, fold, split, 'positive.fa.gz', sep = '_'))
    posFasta = mapply(writeGRangesToFasta, gr = rhesusOnly_enhList_split[[fold]][[cell]][split],  
                      file = pos.fasta_fn, genome = genome)
  }}

save_fn = file.path(PROJDIR, 'rdas', paste('caudate_positive_sequences', paste0(genome, 'Only'),'rda', sep = '.'))
save(rhesusOnly_enhList, file = save_fn )




