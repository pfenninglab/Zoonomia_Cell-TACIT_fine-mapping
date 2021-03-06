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

genome = 'rheMac10' # do everything in rheMac8 coordinates, export in rheMac10



############################
### get the human peaks ####
rhesus_fn =file.path('../../../data/raw_data/','rheMac10') %>%
  list.files(path = ., full.names = T, recursive = T, pattern = 'narrowPeak.gz')
rhesus_fn = rhesus_fn[!grepl('HALPER.narrowPeak.gz', rhesus_fn)]
names(rhesus_fn) = basename(rhesus_fn) %>% ss('\\.', 2)
rhesus_fn = rhesus_fn[!grepl('Consensus', rhesus_fn)]
rhesus_fn = rhesus_fn[order(names(rhesus_fn))]
rhesus_peakList = lapply(rhesus_fn, rtracklayer::import) %>% GRangesList()
rhesus_peakList = lapply(rhesus_peakList, convertHalChrName, chrOut = 'UCSC', 
                         species = 'Macaca_mulatta') %>% GRangesList()

# annotate peaks, filter out exons, promoters
# need to convert genBank chr names to UCSC chr names
rhesus_peakList = annotatePeaks(rhesus_peakList, fromTSS = fromTSS, genome = genome)
rhesus_enhList = filterPeaks(rhesus_peakList, include = include)
lengths(rhesus_enhList) / lengths(rhesus_peakList) ## 55% are noncoding

## read in rhesus peaks mapped to hg38 to get chr for split
rhesus_hal2hg38_fn = rhesus_fn %>% gsub(pattern = 'peak',replacement = 'halper') %>%
  gsub(pattern = '.narrowPeak.gz', replacement = '.Macaca_mulattaToHomo_sapiens.HALPER.narrowPeak.gz')
rhesus_hal2hg3_peakList = lapply(rhesus_hal2hg38_fn, rtracklayer::import) %>% GRangesList()

#tranfer hg38 chr to `col` column in the rhesus_nonEnhList
rhesus_enhList = mapply(transferColumn, toPeaks = rhesus_enhList, colOut = 'col',
                           fromPeaks = rhesus_hal2hg3_peakList) %>% GRangesList()
table(hg38 = rhesus_enhList[[1]]$col, mm10 = seqnames(rhesus_enhList[[1]]))



############################################
### get the mouse mapped to human peaks ####
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

# annotate and keep intron, distal intergenic
# exclude mapped mouse peak overlapping human peak in same cell type
mouse_hal2rm8_peakList = annotatePeaks(mouse_hal2rm8_peakList, fromTSS = fromTSS, genome = genome)
mouse_hal2rm8_enhList = filterPeaks(mouse_hal2rm8_peakList, include = include)
mouse_nonEnhList = mapply(getNonEnhOrthPeaks, inPeaks = mouse_hal2rm8_enhList, 
                          excludePeaks = rhesus_peakList)

### get the mouse mapped to human peaks, used in splitting sets
mouse_hal2hg38_fn = mouse_hal2rm8_fn %>% 
  gsub(pattern = 'ToMacaca_mulatta.HALPER.narrowPeak.gz', 
       replacement = 'ToHomo_sapiens.HALPER.narrowPeak.gz')
mouse_hal2hg38_peakList = lapply(mouse_hal2hg38_fn, rtracklayer::import) %>% GRangesList()

#tranfer hg38 chr to `col` column in the mouse_nonEnhList
mouse_nonEnhList = mapply(transferColumn, toPeaks = mouse_nonEnhList, colOut = 'col',
                       fromPeaks = mouse_hal2hg38_peakList) %>% GRangesList()
mouse_nonEnhList = mouse_nonEnhList %>% lapply(as.data.frame) %>% 
  bind_rows(.id = 'cell_type') %>% arrange(seqnames, start, end, peak) %>% 
  fill(contains('col')) %>% group_by(cell_type) %>% group_split() %>%
  lapply(GRanges) %>% GRangesList()
names(mouse_nonEnhList) = sapply(mouse_nonEnhList, function(x) unique(x$cell_type))
table(hg38 = mouse_nonEnhList[[2]]$col, rm8 = seqnames(mouse_nonEnhList[[2]]))




#############################################
### get the monkey mapped to human peaks ####
human_hal2rm8_fn =file.path('../../../data/raw_data/','hg38') %>%
  list.files(path = ., full.names = T, recursive = T, 
             pattern = 'ToMacaca_mulatta.HALPER.narrowPeak.gz')
names(human_hal2rm8_fn) = basename(human_hal2rm8_fn) %>% ss('\\.', 2)
human_hal2rm8_fn = human_hal2rm8_fn[names(rhesus_enhList)]
human_hal2rm8_peakList = lapply(human_hal2rm8_fn, rtracklayer::import) %>% 
  lapply(convertHalChrName, chrOut = 'UCSC', species = 'Macaca_mulatta') %>% GRangesList()

# annotate and keep intron, distal intergenic
# exclude mapped mouse peak overlapping human peak in same cell type
human_hal2rm8_peakList = annotatePeaks(human_hal2rm8_peakList, fromTSS = fromTSS, genome = genome)
human_hal2rm8_enhList = filterPeaks(human_hal2rm8_peakList, include = include)
human_nonEnhList = mapply(getNonEnhOrthPeaks, inPeaks = human_hal2rm8_enhList, 
                          excludePeaks = rhesus_enhList)

### get the mouse mapped to human peaks, used in splitting sets
human_fn = human_hal2rm8_fn %>% 
  gsub(pattern = 'halper',replacement = 'peak') %>%
  gsub(pattern = '.Homo_sapiensToMacaca_mulatta.HALPER.narrowPeak.gz', replacement = '.narrowPeak.gz')
human_peakList = lapply(human_fn, rtracklayer::import) %>% GRangesList()

#tranfer hg38 chr to `col` column in the mouse_enhList
human_nonEnhList = mapply(transferColumn, toPeaks = human_nonEnhList, colOut = 'col',
                       fromPeaks = human_peakList) %>% GRangesList()
table(hg38 = human_nonEnhList[[2]]$col, rm8 = seqnames(human_nonEnhList[[2]]))



###########################
# number of postives 
lengths(rhesus_enhList) / 1000
#  Astro INT_Pvalb Microglia    MSN_D1    MSN_D2    MSN_SN     Oligo       OPC 
# 65.574    29.263    59.072    72.363    72.601    44.666    39.496    48.376 

# combine sets of non-enh orthologs and save #
# number of non-enhancer orthologs
nonEnh_peakList = GRangesList(mapply(c, human_nonEnhList, mouse_nonEnhList))
all.equal(lengths(human_nonEnhList) + lengths(mouse_nonEnhList), lengths(nonEnh_peakList))
lengths(nonEnh_peakList) / 1000
#  Astro INT_Pvalb Microglia    MSN_D1    MSN_D2    MSN_SN     Oligo       OPC 
# 35.103    75.464    23.679    63.727    64.621    35.555    32.652    29.842

## ratio of positives to non-enhancer orthologs
lengths(rhesus_enhList) / lengths(nonEnh_peakList)
#     Astro INT_Pvalb Microglia    MSN_D1    MSN_D2    MSN_SN     Oligo       OPC 
# 1.8680455 0.3877743 2.4946999 1.1355156 1.1234893 1.2562509 1.2096043 1.6210710 



##############################################
## liftOver rheMac8 coordinates to rheMac10 
## split positives and negatives across folds
chainFile =file.path("/home/bnphan/resources/liftOver_chainz",
                     'rheMac8ToRheMac10.over.chain')

# export postive sequences to summit-centered 501bp fasta file
positiveSet = summitCenter(rhesus_enhList, width = 501)
positiveSet_rheMac10 = lapply(positiveSet, liftOver_narrowPeak, chainFile = chainFile)
positiveSplit_list = lapply(folds, function(fold){
  ret = lapply(positiveSet_rheMac10, splitPeakSet, testSet = testSet, 
               validSet = fold, useCol = 'col')
  return(ret)
})

# export postive sequences to summit-centered 501bp fasta file
negativeSet = summitCenter(nonEnh_peakList, width = 501)
negativeSet_rheMac10 = lapply(negativeSet, liftOver_narrowPeak, chainFile = chainFile)
negativeSplit_list = lapply(folds, function(fold){
  ret = lapply(negativeSet_rheMac10, splitPeakSet, testSet = testSet, 
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

dir.create(file.path(PROJDIR, 'rdas'))
save_fn = file.path(PROJDIR, 'rdas', 
                    paste('cnn_enhancer_non-enhancer_split', genome,'.rda', sep = '_'))
save(positiveSet_rheMac10, positiveSplit_list, 
     negativeSet_rheMac10, negativeSplit_list, file = save_fn )


