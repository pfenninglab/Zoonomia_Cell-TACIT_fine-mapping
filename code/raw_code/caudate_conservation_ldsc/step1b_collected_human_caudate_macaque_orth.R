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
## find they hg38 cCREs that can be halper-ed to rheMac10
human_hal2rheMac10_fn = human_peak_fn %>% 
  gsub(pattern = '.narrowPeak.gz', 
       replacement = '.Homo_sapiensToMacaca_mulatta.HALPER.narrowPeak.gz') %>% 
  gsub(pattern = 'peak', replacement = 'halper')
human_hal2rheMac10_peakList = lapply(human_hal2rheMac10_fn, rtracklayer::import) %>% 
  GRangesList()

human_mappable_peakList = mapply(function(hum, hg38torheMac10){
  return(hum[hum$name %in% hg38torheMac10$name])
}, hum = human_peakList, hg38torheMac10 = human_hal2rheMac10_peakList) %>% 
  GRangesList()

# see how many peak are mappable
lengths(human_mappable_peakList) / lengths(human_peakList)

## write mappable peak to narrowPeak files
out_narrowPeak_fn = human_peak_fn %>%
  gsub(pattern = 'Corces2020_caudate.', replacement = 'Corces2020_caudate_mappedToRheMac10.')
# outList = mapply(write_GRangesToNarrowPeak,gr = human_mappable_peakList,
#                  file = out_narrowPeak_fn, genome = 'hg38')

#########################################################
## find rhesus peaks that can be halpered to hg38
rhesus_peak_fn =file.path('../../../data/raw_data/','rheMac10') %>%
  list.files(path = ., full.names = T, pattern = 'narrowPeak.gz', recursive = T)
rhesus_peak_fn = rhesus_peak_fn[ss(basename(rhesus_peak_fn), '\\.') %in% 
                                c('Stauffer_caudate', 'Pfenning_Cpu')]
rhesus_peak_fn = rhesus_peak_fn[!grepl('Macaca_mulattaTo', rhesus_peak_fn)]
names(rhesus_peak_fn) = basename(rhesus_peak_fn) %>% ss('\\.', 2)
rhesus_peak_fn = rhesus_peak_fn[names(human_peak_fn)]
rhesus_peakList = lapply(rhesus_peak_fn, rtracklayer::import) %>% GRangesList()

# get the hg38 coordinates of rhesus peak
rhesus_hal2hg38_fn = rhesus_peak_fn %>% 
  gsub(pattern = '.narrowPeak.gz', 
       replacement = '.Macaca_mulattaToHomo_sapiens.HALPER.narrowPeak.gz') %>% 
  gsub(pattern = 'peak', replacement = 'halper')
rhesus_hal2hg38_peakList = lapply(rhesus_hal2hg38_fn, rtracklayer::import) %>% 
  GRangesList()

# percent of rhesus orthologs in hg38
lengths(rhesus_hal2hg38_peakList) / lengths(rhesus_peakList)

#########################################################
## find halpered rhesus peak that can be halpered to hg38
ortholog_peakList = mapply(function(human, rheMac10toHum){
  # find halpered human peak in rheMac10 that overlaps rhesus peak
  oo = findOverlaps(query = human, subject = rheMac10toHum)
  h2 = human[unique(queryHits(oo))]
  m2h2 = rheMac10toHum[unique(subjectHits(oo))]
  ret = GenomicRanges::reduce(c(h2, m2h2))
  return(ret)
}, human = human_peakList, rheMac10toHum = rhesus_hal2hg38_peakList) %>% 
  GRangesList()

# percent of rhesus cCRE mapped to human w/ human cCRE
lengths(ortholog_peakList) / lengths(rhesus_hal2hg38_peakList)

# percent of human cCRE w/ rhesus cCRE mapped to human
lengths(ortholog_peakList) / lengths(human_peakList)

#number of orthologs
lengths(ortholog_peakList)

# write mappable peak to narrowPeak files
out_narrowPeak3_fn = human_peak_fn %>% 
  gsub(pattern = 'Corces2020_caudate', replacement = 'Corces2020_caudate_hgRmOrth')
# outList = mapply(write_GRangesToNarrowPeak,gr = ortholog_peakList, 
#                  file = out_narrowPeak3_fn, genome = 'hg38')

##################################################
## export granges list objects for making plots ##
human_macaque_orthologs = list(hgPeaks = human_peakList,
                               hg2Rm = human_mappable_peakList,
                               rmPeaks = rhesus_peakList, 
                               rm2Hg = rhesus_hal2hg38_peakList,
                               hgRmOrth = ortholog_peakList)
save_fn = PROJDIR=file.path('../../../data/raw_data/','caudate_conservation_ldsc',
                            'rdas/human_macaque_orthologs_peakList.rda')
save(human_macaque_orthologs, file = save_fn, compress = T)



