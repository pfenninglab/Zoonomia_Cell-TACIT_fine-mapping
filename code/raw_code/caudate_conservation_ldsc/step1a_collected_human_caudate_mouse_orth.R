# to be run in the root github directory
LABEL='caudate_conservation_ldsc'
setwd('code/raw_code/caudate_conservation_ldsc')
PROJDIR=file.path('../../../data/raw_data/',LABEL)

#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(rtracklayer))
library(tidyverse)
source('../hal_scripts/narrowPeakFunctions.R')

######################################
# read in the human hg38 cCRE regions
human_peak_fn = file.path('../../../data/raw_data/','hg38','Corces_2020', 'peak') %>%
  list.files(path = ., full.names = T, pattern = 'narrowPeak.gz')
human_peak_fn = human_peak_fn[ss(basename(human_peak_fn), '\\.') == 'Corces2020_caudate']
names(human_peak_fn) = basename(human_peak_fn) %>% ss('\\.', 2)
human_peakList= lapply(human_peak_fn, rtracklayer::import) %>% GRangesList()

######################################################
## find they hg38 cCREs that can be halper-ed to mm10
human_hal2mm10_fn = human_peak_fn %>% 
  gsub(pattern = '.narrowPeak.gz', 
       replacement = '.Homo_sapiensToMus_musculus.HALPER.narrowPeak.gz') %>% 
  gsub(pattern = 'peak', replacement = 'halper')
human_hal2mm10_peakList = lapply(human_hal2mm10_fn, rtracklayer::import) %>% 
  GRangesList()

human_mappable_peakList = mapply(function(hum, hg38tomm10){
  return(hum[hum$name %in% hg38tomm10$name])
}, hum = human_peakList, hg38tomm10 = human_hal2mm10_peakList) %>% 
  GRangesList()

# see how many peak are mappable
lengths(human_mappable_peakList) / lengths(human_peakList)

# write mappable peak to narrowPeak files
out_narrowPeak_fn = human_peak_fn %>% 
  gsub(pattern = 'Corces2020_caudate.', replacement = 'Corces2020_caudate_mappedToMm10.')
# outList = mapply(write_GRangesToNarrowPeak,gr = human_mappable_peakList,
#                  file = out_narrowPeak_fn, genome = 'hg38')

#########################################################
## find mouse peaks that can be halpered to hg38
mouse_peak_fn =file.path('../../../data/raw_data/','mm10') %>%
  list.files(path = ., full.names = T, pattern = 'narrowPeak.gz', recursive = T)
mouse_peak_fn = mouse_peak_fn[ss(basename(mouse_peak_fn), '\\.') %in% 
                                c('BICCN_CP', 'Pfenning_Cpu')]
# exclude BICCN_CP.MSN_D1 and BICCN_CP.MSN_D2 & halper peaks for now
mouse_peak_fn = mouse_peak_fn[!grepl('Mus_musculusTo', mouse_peak_fn)]
mouse_peak_fn = mouse_peak_fn[!grepl('BICCN_CP.MSN_D|BICCN_CP.INT', mouse_peak_fn)]
names(mouse_peak_fn) = basename(mouse_peak_fn) %>% ss('\\.', 2)
mouse_peak_fn = mouse_peak_fn[names(human_peak_fn)]
mouse_peakList = lapply(mouse_peak_fn, rtracklayer::import) %>% GRangesList()

# get the hg38 coordinates of mouse peak
mouse_hal2hg38_fn = mouse_peak_fn %>% 
  gsub(pattern = '.narrowPeak.gz', 
       replacement = '.Mus_musculusToHomo_sapiens.HALPER.narrowPeak.gz') %>% 
  gsub(pattern = 'peak', replacement = 'halper')
mouse_hal2hg38_peakList = lapply(mouse_hal2hg38_fn, rtracklayer::import) %>% 
  GRangesList()

# percent of mouse orthologs in hg38
lengths(mouse_hal2hg38_peakList) / lengths(mouse_peakList)

#########################################################
## find halpered mouse peak that can be halpered to hg38
ortholog_peakList = mapply(function(human, mm10toHum){
  # find halpered human peak in mm10 that overlaps mouse peak
  oo = findOverlaps(query = human, subject = mm10toHum)
  h2 = human[unique(queryHits(oo))]
  m2h2 = mm10toHum[unique(subjectHits(oo))]
  ret = GenomicRanges::reduce(c(h2, m2h2))
  return(ret)
}, human = human_peakList, mm10toHum = mouse_hal2hg38_peakList) %>% 
  GRangesList()

# percent of mouse cCRE mapped to human w/ human cCRE
lengths(ortholog_peakList) / lengths(mouse_hal2hg38_peakList)

# percent of human cCRE w/ mouse cCRE mapped to human
lengths(ortholog_peakList) / lengths(human_peakList)

#number of orthologs
lengths(ortholog_peakList)

# write mappable peak to narrowPeak files
out_narrowPeak3_fn = human_peak_fn %>% 
  gsub(pattern = 'Corces2020_caudate', replacement = 'Corces2020_caudate_hgMmOrth')
outList = mapply(write_GRangesToNarrowPeak,gr = ortholog_peakList,
                 file = out_narrowPeak3_fn, genome = 'hg38')

##################################################
## export granges list objects for making plots ##
human_mouse_orthologs = list(hgPeaks = human_peakList,
                               hg2Mm = human_mappable_peakList,
                               MmPeaks = mouse_peakList, 
                               Mm2Hg = mouse_hal2hg38_peakList,
                               hgMmOrth = ortholog_peakList)
save_fn = PROJDIR=file.path('../../../data/raw_data/','caudate_conservation_ldsc',
                            'rdas/human_mouse_orthologs_peakList.rda')
save(human_mouse_orthologs, file = save_fn, compress = T)
