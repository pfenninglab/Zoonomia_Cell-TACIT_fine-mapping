# to be run in the root github directory
setwd('code/raw_code/cross_species_peak_orthologs')
PROJDIR=file.path('../../../data/raw_data/cross_species_peak_orthologs')
TMPDIR=file.path('/scratch/bnphan/glmGamPoi')
system(paste('mkdir -p', TMPDIR))

#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F, repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(rtracklayer)); suppressMessages(library(glmGamPoi))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(GenomicRanges))
suppressMessages(library(HDF5Array))
suppressMessages(library(tidyverse))
suppressMessages(library(BiocParallel))
register(MulticoreParam(round(parallel::detectCores()/3)))

source('../hal_scripts/narrowPeakFunctions.R')

#############################
### load peak count data ####
orthPeakH5_fn = file.path(PROJDIR, 'rdas', 'orthologPeakH5')
if(!file.exists(orthPeakH5_fn)){
  matListRDS_fn = file.path(PROJDIR, 'rdas', 'multispeciesMergedMatrixList.rds')
  peakListRDS_fn = file.path(PROJDIR, 'rdas', 'multispeciesMergedPeakList.rds')
  
  matList = readRDS(file = matListRDS_fn)
  peakList = readRDS(file = peakListRDS_fn)
  
  #############################################################
  ### combine the data & make summarized experiment object ####
  counts <- do.call('cbind', lapply(matList, function(x) assays(x)[[1]]))
  rownames(counts) = names(peakList[["Consensus"]]) = peakList[["Consensus"]]$name
  rowData = peakList[["Consensus"]]
  keepCols = Reduce('intersect',lapply(matList, function(x) names(colData(x))))
  
  colData <- do.call('rbind', lapply(matList, function(x) colData(x)[,keepCols]))
  colData = as.data.frame(colData)
  colData$log10nFrags = log10(colData$nFrags)
  
  orthPeakRSE = SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts), rowRanges=rowData, colData=colData)
  saveHDF5SummarizedExperiment(orthPeakRSE, dir = orthPeakH5_fn)
} 


#######################################
## copy to tmp directory for fast IO ##
tmpOrthPeakH5_fn = file.path(TMPDIR, 'orthologPeakH5')
system(paste('cp -R', orthPeakH5_fn,tmpOrthPeakH5_fn))
setwd(TMPDIR)

orthPeakRSE = loadHDF5SummarizedExperiment(tmpOrthPeakH5_fn)

## keep only cell types in all species ##
tt = table(orthPeakRSE$Species, orthPeakRSE$Clusters2)
keepCellTypes = colnames(tt)[apply(tt, 2, prod) > 0]
orthPeakRSE = orthPeakRSE[, orthPeakRSE$Clusters2 %in% keepCellTypes]
keepCellTypes

## subset to only one  ##
keepSamples = c('14_1018.CAUD', 'CAUD_WS1H_STA682A131','CEMBA171214_4D','CEMBA180813_5E')
orthPeakRSE = orthPeakRSE[, orthPeakRSE$Sample %in% keepSamples]


## scale and center confounding variables
colData = as.data.frame(colData(orthPeakRSE))
nSample = ncol(orthPeakRSE)
colData = colData %>% mutate(
  log10nFrags.scaled = (log10nFrags - mean(log10nFrags))/sd(log10nFrags),
  TSSEnrichment.scaled = (TSSEnrichment - mean(TSSEnrichment))/sd(TSSEnrichment))


##################################################
## fit glmGamPoi model on open chromatin counts ##
fitGlmGamPoi_fn =file.path('/home/bnphan/projects/snATAC_cross_species_caudate',
                           'data/raw_data/cross_species_peak_orthologs/rdas', 
                           paste0('orthologPeakGlmGamPoiFit_N',nSample, '.rds'))

if(!file.exists(fitGlmGamPoi_fn)){
fit <- glm_gp(design = ~Species + Clusters2 + Species:Clusters2,
  data = orthPeakRSE, on_disk = TRUE, col_data = colData, verbose = TRUE)

saveRDS(fit, file = fitGlmGamPoi_fn, compress = T)
} else{
  fit = readRDS(fitGlmGamPoi_fn)
}

###############################################
## find species conserved D1 and D2 clusters ##
res <- test_de(fit, reduced_design = ~Species + Clusters2)


