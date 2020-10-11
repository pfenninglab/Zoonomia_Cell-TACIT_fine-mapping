# to be run in the root github directory
LABEL='Stauffer_caudate'
setwd('code/raw_code/rheMac10_Stauffer_caudate')
PROJDIR=file.path('../../../data/raw_data/rheMac10',LABEL)

#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(rtracklayer))

source('../hal_scripts/narrowPeakFunctions.R')

#########################################################
### get the peaks from experiment and halpered peaks ####
narrowPeak_fn = list.files(file.path(PROJDIR,'peak'), full.names = T,
                           pattern ='.narrowPeak.gz')
names(narrowPeak_fn) = ss(basename(narrowPeak_fn), '\\.',2)
halper_narrowPeak_fn = list.files(file.path(PROJDIR,'halper'), full.names = T,
                                pattern ='HALPER.narrowPeak.gz')

peakList_fn = lapply(narrowPeak_fn, function(n){
  tmp = c(n, grep(ss(basename(n), '\\.',2), halper_narrowPeak_fn, value = T))
  val = grepl('ToHuman', tmp)
  if(any(val))
    # put the human coordinates first
    tmp = c(tmp[val], tmp[!val])
  names(tmp) = ss(basename(tmp), '\\.',2)
  return(tmp)
})

peakList = lapply(peakList_fn, function(ls){
  GRangesList(lapply(ls, import))
})

#####################################################
### get human peaks w/ mouse & macaque orthologs ####
orthologPeakRanges = lapply(peakList, keepOrthologs, idxReturn = 1)
ORTHDIR=file.path(PROJDIR,'..','..','cross_species_peak_orthologs')
ortholog_fn = file.path(ORTHDIR, 'peaks',paste0(LABEL,'_orthologPeakList.rds'))
saveRDS(orthologPeakRanges, file = ortholog_fn)

