# to be run in the root github directory
LABEL='meuleman_dhs_index'
setwd('code/raw_code/meuleman_dhs_index')
PROJDIR=file.path('../../../data/raw_data/',LABEL)

#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F)
suppressMessages(library(rtracklayer))
suppressMessages(library(tidyverse))
source('../hal_scripts/narrowPeakFunctions.R')

######################################
# read in the human hg38 cCRE regions
human_peak_fn = list.files(path = file.path(PROJDIR, 'peaks'), full.names = T,
                           pattern = 'narrowPeak.gz')
human_peak_fn = human_peak_fn[ss(basename(human_peak_fn), '\\.') == 'DHS_Index_and_Vocabulary']
names(human_peak_fn) = basename(human_peak_fn) %>% 
  ss('DHS_Index_and_Vocabulary.', 2) %>% ss('\\.narrowPeak', 1)
human_peakList= lapply(human_peak_fn, rtracklayer::import) %>% GRangesList()

######################################################
## find they hg38 cCREs that can be halper-ed to mm10
human_hal2mm10_fn = human_peak_fn %>% 
  gsub(pattern = '.narrowPeak.gz', 
       replacement = '.Homo_sapiensToMus_musculus.HALPER.narrowPeak.gz') %>% 
  gsub(pattern = 'peaks/', replacement = 'halper/')
all(file.exists(human_hal2mm10_fn))
human_hal2mm10_peakList = lapply(human_hal2mm10_fn, rtracklayer::import) %>% 
  GRangesList()

human_mappable_peakList = mapply(function(hum, hg38tomm10){
  return(hum[hum$name %in% hg38tomm10$name])
}, hum = human_peakList, hg38tomm10 = human_hal2mm10_peakList) %>% 
  GRangesList()

# see how many peaks are mappable
lengths(human_mappable_peakList) / lengths(human_peakList)

# write mappable peaks to narrowPeak files
out_narrowPeak_fn = human_peak_fn %>% 
  gsub(pattern = 'DHS_Index_and_Vocabulary', 
       replacement = 'DHS_Index_and_Vocabulary_mappedToMm10')
outList = mapply(write_GRangesToNarrowPeak,gr = human_mappable_peakList, 
                 file = out_narrowPeak_fn, genome = 'hg38')

