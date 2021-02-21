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

genome = 'hg38'

############################
### get the human peaks ####
human_peak_fn = file.path('../../../data/raw_data/',genome,'Corces_2020', 'peak') %>%
  list.files(path = ., full.names = T, pattern = 'narrowPeak.gz')
human_peak_fn = human_peak_fn[ss(basename(human_peak_fn), '\\.') == 'Corces2020_caudate']
names(human_peak_fn) = basename(human_peak_fn) %>% ss('\\.', 2)
human_peak_fn = human_peak_fn[! names(human_peak_fn) %in% c('Consensus')]
human_peakList= lapply(human_peak_fn, rtracklayer::import) %>% GRangesList()

# annotate peaks, filter out exons, promoters
human_peakList = annotatePeaks(human_peakList, fromTSS = fromTSS, genome = genome)
human_enhList = filterPeaks(human_peakList, include = include)
lengths(human_enhList) / lengths(human_peakList) ## 25% are in intron/distal

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
mouse_hal2hg38_fn = mouse_hal2hg38_fn[names(human_enhList)]
mouse_hal2hg38_peakList = lapply(mouse_hal2hg38_fn, rtracklayer::import) %>% GRangesList()

# annotate and keep intron, distal intergenic
# exclude mapped mouse peak overlapping human peak in same cell type
mouse_hal2hg38_peakList = annotatePeaks(mouse_hal2hg38_peakList, fromTSS = fromTSS, genome = genome)
mouse_hal2hg38_enhList = filterPeaks(mouse_hal2hg38_peakList, include = include)
mouse_nonEnhList = mapply(getNonEnhOrthPeaks, inPeaks = mouse_hal2hg38_enhList, 
                          excludePeaks = human_peakList)

#############################################
### get the monkey mapped to human peaks ####
rhesus_hal2hg38_fn =file.path('../../../data/raw_data/','rheMac10') %>%
  list.files(path = ., full.names = T, recursive = T, 
             pattern = 'ToHomo_sapiens.HALPER.narrowPeak.gz')
names(rhesus_hal2hg38_fn) = basename(rhesus_hal2hg38_fn) %>% ss('\\.', 2)
rhesus_hal2hg38_fn = rhesus_hal2hg38_fn[names(human_enhList)]
rhesus_hal2hg38_peakList = lapply(rhesus_hal2hg38_fn, rtracklayer::import) %>% GRangesList()

# annotate and keep intron, distal intergenic
# exclude mapped mouse peak overlapping human peak in same cell type
rhesus_hal2hg38_peakList = annotatePeaks(rhesus_hal2hg38_peakList, fromTSS = fromTSS, genome = genome)
rhesus_hal2hg38_enhList = filterPeaks(rhesus_hal2hg38_peakList, include = include)
rhesus_nonEnhList = mapply(getNonEnhOrthPeaks, inPeaks = rhesus_hal2hg38_enhList, 
                          excludePeaks = human_peakList)

######################
# number of postives 
lengths(human_enhList) / 1000
#  Astro INT_Pvalb Microglia    MSN_D1    MSN_D2    MSN_SN     Oligo       OPC 
#   37.393    42.110    26.408    47.199    46.902    33.525    26.813    30.548 

# combine sets of non-enh orthologs and save #
# number of non-enhancer orthologs
nonEnh_peakList = GRangesList(mapply(c, rhesus_nonEnhList, mouse_nonEnhList))
all.equal(lengths(rhesus_nonEnhList) + lengths(mouse_nonEnhList), lengths(nonEnh_peakList))
lengths(nonEnh_peakList) / 1000
# Astro INT_Pvalb Microglia    MSN_D1    MSN_D2    MSN_SN     Oligo       OPC 
#    24.256    29.507    19.756    53.509    54.487    18.585    14.193    16.035

## ratio of positives to non-enhancer orthologs
lengths(human_enhList) / lengths(nonEnh_peakList)
# Astro INT_Pvalb Microglia    MSN_D1    MSN_D2    MSN_SN     Oligo       OPC 
#1.5415980 1.4271190 1.3367078 0.8820759 0.8607925 1.8038741 1.8891707 1.9050826 

###############################################
# split positives and negatives across folds
positiveSet = summitCenter(human_enhList, width = 501)
positiveSplit_list = lapply(folds, function(fold){
  ret = lapply(positiveSet, splitPeakSet, testSet = testSet, validSet = fold)
  return(ret)
})

# export postive sequences to summit-centered 501bp fasta file
negativeSet = summitCenter(nonEnh_peakList, width = 501)
negativeSplit_list = lapply(folds, function(fold){
  ret = lapply(negativeSet, splitPeakSet, testSet = testSet, validSet = fold)
  return(ret)
})

#############################################################
# export postive sequences to summit-centered 501bp fasta file
split = names(negativeSplit_list[[1]][[1]])
system(paste('mkdir -p',  file.path(PROJDIR, 'fasta')))
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
save_fn = file.path(PROJDIR, 'rdas', 
                    paste('cnn_enhancer_non-enhancer_split', genome,'.rda', sep = '_'))
save(positiveSet, positiveSplit_list, 
     negativeSet, negativeSplit_list, file = save_fn )


