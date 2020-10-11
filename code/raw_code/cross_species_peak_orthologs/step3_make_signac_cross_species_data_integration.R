# to be run in the root github directory
setwd('code/raw_code/cross_species_peak_orthologs')
PROJDIR=file.path('../../../data/raw_data/cross_species_peak_orthologs')

#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F, repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(Signac)); suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(EnsDb.Hsapiens.v79))

source('../hal_scripts/narrowPeakFunctions.R')

### load peak count data ####
saveRDA_fn = file.path(PROJDIR, 'rdas', 'multispeciesMergedMatrixList.rda')
load(file = saveRDA_fn)

#####################################################
### combine the data & make Signac/Seurat object ####
counts <- do.call('cbind', lapply(matList, function(x) assays(x)[[1]]))
keepCols = Reduce('intersect',lapply(matList, function(x) names(colData(x))))
metadata <- do.call('rbind', lapply(matList, function(x) colData(x)[,keepCols]))
metadata = as.data.frame(metadata)

chrom_assay <- CreateChromatinAssay(
  counts = counts, ranges = peakList[["Consensus"]], sep = c(":", "-",'.'), 
  genome = 'hg38', fragments = NULL, min.cells = 0, min.features = 0)

obj_seurat <- CreateSeuratObject( counts = chrom_assay, assay = "peaks", 
                                  metadata = metadata)

### extract gene annotations from EnsDb ####
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v79)

# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(obj_seurat) <- annotations

#########################################
### perform TF-IDF on the peak counts ###
obj_seurat <- RunTFIDF(obj_seurat)
obj_seurat <- FindTopFeatures(obj_seurat, min.cutoff = 'q0')
obj_seurat <- RunSVD(obj_seurat)

## save the seurat object
saveRDS_fn = file.path(PROJDIR, 'rdas', 'multispeciesMergedSeurat.rds')
system(paste('mkdir -p', dirname(saveRDS_fn)))
saveRDS(obj_seurat, file = saveRDS_fn)
