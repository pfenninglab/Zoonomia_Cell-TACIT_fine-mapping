# to be run in the root github directory
LABEL='mouse-human_peak_orthologs'
setwd('code/raw_code/mouse-human_peak_orthologs')
HUMDIR=file.path('../../../data/raw_data/hg38/Corces_2020')
MUSDIR1=file.path('../../../data/raw_data/mm10/BICCN_CP')
MUSDIR2=file.path('../../../data/raw_data/mm10/Mouse_cSNAIL_D1D2')

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
human_peak_fn = list.files(path = file.path(HUMDIR, 'peak'), full.names = T,
                           pattern = 'narrowPeak.gz')
names(human_peak_fn) = basename(human_peak_fn) %>% ss('\\.', 2)
human_peak_fn = 
  human_peak_fn[basename(human_peak_fn) %>% ss('\\.',1) == 'Corces2020_caudate']
human_peakList= lapply(human_peak_fn, rtracklayer::import) %>% GRangesList()

######################################################
## find they hg38 cCREs that can be halper-ed to mm10
human_hal2mm10_fn = human_peak_fn %>% 
  gsub(pattern = '.narrowPeak.gz', 
       replacement = '.Homo_sapiensToMus_musculus.HALPER.narrowPeak.gz') %>% 
  gsub(pattern = 'peak', replacement = 'halper')
human_hal2mm10_fn[!file.exists(human_hal2mm10_fn)]
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
  gsub(pattern = 'Corces2020_caudate', replacement = 'Corces2020_caudate_mappedToMm10')
outList = mapply(write_GRangesToNarrowPeak,gr = human_mappable_peakList, 
                 file = out_narrowPeak_fn, genome = 'hg38')

#########################################################
## find halpered mouse peaks that can be halpered to hg38
mouse_peak_fn = 
  c(list.files(path = file.path(MUSDIR1, 'peak'), full.names = T,pattern = 'narrowPeak.gz'),
    list.files(path = file.path(MUSDIR2, 'peak'), full.names = T,pattern = 'narrowPeak.gz'))
names(mouse_peak_fn) = basename(mouse_peak_fn) %>% ss('.narrowPeak.gz', 1)
mouse_peak_fn = 
  mouse_peak_fn[basename(mouse_peak_fn) %>% ss('\\.',1) %in% c('BICCN_CP','Pfenning_Cpu')]
mouse_peak_fn = mouse_peak_fn[
  !names(mouse_peak_fn) %in% c("BICCN_CP.MSN_D1","BICCN_CP.MSN_D2",
                                   "Pfenning_Cpu.bulk","Pfenning_Cpu.INT_SST")]
mouse_peak_fn = mouse_peak_fn[order(ss(names(mouse_peak_fn), '\\.', 2))]
mouse_peakList = lapply(mouse_peak_fn, rtracklayer::import) %>% GRangesList()

# get the hg38 coordinates of mouse peaks
mouse_hal2hg38_fn = mouse_peak_fn %>% gsub(pattern = '.narrowPeak.gz', 
       replacement = '.Mus_musculusToHomo_sapiens.HALPER.narrowPeak.gz') %>% 
  gsub(pattern = 'peak', replacement = 'halper')
mouse_hal2hg38_fn[!file.exists(mouse_hal2hg38_fn)]
mouse_hal2hg38_peakList = lapply(mouse_hal2hg38_fn, rtracklayer::import) %>% 
  GRangesList()

# percent of mouse orthologs in hg38
lengths(mouse_hal2hg38_peakList) / lengths(mouse_peakList)

# write mappable peaks to narrowPeak files
out_narrowPeak2_fn = mouse_peak_fn %>% 
  gsub(pattern = 'BICCN_CP', 
       replacement = 'BICCN_CP_mappedToHg38')
outList = mapply(write_GRangesToNarrowPeak,gr = mouse_hal2hg38_peakList, 
                 file = out_narrowPeak2_fn, genome = 'hg38')

#########################################################
## find halpered mouse peaks that can be halpered to hg38
mouse_ortholog_peakList = mapply(function(human, mm10toHum){
  start(mm10toHum) = start(mm10toHum) - 250
  end(mm10toHum) = end(mm10toHum) + 250
  # find halpered human peak in mm10 that overlaps mouse peak
  oo = findOverlaps(query = human, subject = mm10toHum)
  return(human[unique(queryHits(oo))])
}, human = human_peakList, mm10toHum = mouse_hal2hg38_peakList) %>% 
  GRangesList()

# percent of mouse cCRE mapped to human w/ human cCRE
lengths(mouse_ortholog_peakList) / lengths(mouse_hal2hg38_peakList)

# percent of human cCRE w/ mouse cCRE mapped to human
lengths(mouse_ortholog_peakList) / lengths(human_peakList)

# write mappable peaks to narrowPeak files
out_narrowPeak3_fn = human_peak_fn %>% 
  gsub(pattern = 'Corces2020_caudate', replacement = 'BICCN_CP_hasOrthologInHg38')
outList = mapply(write_GRangesToNarrowPeak,gr = mouse_ortholog_peakList, 
                 file = out_narrowPeak3_fn, genome = 'hg38')


#########################################################
## find halpered human peaks that can be halpered to mm10
human_ortholog_peakList = mapply(function(human, hg38tomm10, mouse){
  # find halpered human peak in mm10 that overlaps mouse peak
  start(hg38tomm10) = start(hg38tomm10) - 250
  end(hg38tomm10) = end(hg38tomm10) + 250
  oo = findOverlaps(query = hg38tomm10, subject = mouse)
  hg38tomm10 = hg38tomm10[unique(queryHits(oo))]
  return(human[human$name %in% hg38tomm10$name])
}, human = human_peakList, hg38tomm10 = human_hal2mm10_peakList, 
mouse = mouse_peakList) %>% GRangesList()

# percent of human cCRE mapped to mouse w/ mouse cCRE
lengths(human_ortholog_peakList) / lengths(human_hal2mm10_peakList)

# percent of mouse cCRE w/ human cCRE mapped to mouse
lengths(human_ortholog_peakList) / lengths(mouse_peakList)

# percent of human cCRE mapped to mouse w/ mouse cCRE mapped to human
lengths(human_ortholog_peakList) / lengths(mouse_ortholog_peakList)

# write mappable peaks to narrowPeak files
out_narrowPeak4_fn = human_peak_fn %>% 
  gsub(pattern = 'Corces2020_caudate', replacement = 'Corces2020_caudate_hasOrthologInMm10')
outList = mapply(write_GRangesToNarrowPeak,gr = human_ortholog_peakList, 
                 file = out_narrowPeak4_fn, genome = 'hg38')
