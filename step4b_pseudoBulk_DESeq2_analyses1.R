# to be run in the root github directory
setwd('code/raw_code/cross_species_peak_orthologs')
PROJDIR=file.path('../../../data/raw_data/cross_species_peak_orthologs')

#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F, repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(rtracklayer)); suppressMessages(library(DESeq2))

source('../hal_scripts/narrowPeakFunctions.R')

########################################################################
## load in ranged summarized experiment file and create DESeq2 object ##
orthPeakRSE_fn = file.path(PROJDIR, 'rdas', 'orthologPeak_pseudoBulkRSE.rds')
orthPeakRSE = readRDS(file = orthPeakRSE_fn)

