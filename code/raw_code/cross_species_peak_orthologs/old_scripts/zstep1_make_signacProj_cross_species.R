# to be run in the root github directory
setwd('code/raw_code/cross_species_peak_orthologs')
PROJDIR=file.path('../../../data/raw_data/cross_species_peak_orthologs')

#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F, repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(Signac)); suppressMessages(library(Seurat))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(EnsDb.Hsapiens.v79))

source('../hal_scripts/narrowPeakFunctions.R')
source('../hal_scripts/seuratFunctions.R')

### load peak count data ####
matListRDS_fn = file.path(PROJDIR, 'rdas', 'multispeciesMergedMatrixList.rds')
peakListRDS_fn = file.path(PROJDIR, 'rdas', 'multispeciesMergedPeakList.rds')

matList = readRDS(file = matListRDS_fn)
peakList = readRDS(file = peakListRDS_fn)

#####################################################
### combine the data & make Signac/Seurat object ####
counts <- do.call('cbind', lapply(matList, function(x) assays(x)[[1]]))
rownames(counts) = GRangesToString(peakList[["Consensus"]], sep = c(":", "-"))
keepCols = Reduce('intersect',lapply(matList, function(x) names(colData(x))))
metadata <- do.call('rbind', lapply(matList, function(x) colData(x)[,keepCols]))
metadata = as.data.frame(metadata)

obj_seurat  = makeChomAssaySeurat(
  counts, metadata, EnsDb.Hsapiens.v79, 
  group.by = 'Sample', min.cells = 10, verbose = TRUE)

## save the seurat object
saveRDS_fn = file.path(PROJDIR, 'rdas', 'multispeciesMergedSeurat.rds')
system(paste('mkdir -p', dirname(saveRDS_fn)))
saveRDS(obj_seurat, file = saveRDS_fn)


