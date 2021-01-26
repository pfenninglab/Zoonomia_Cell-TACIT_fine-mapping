# to be run in the root github directory
setwd('code/raw_code/cross_species_peak_orthologs')
PROJDIR=file.path('../../../data/raw_data/cross_species_peak_orthologs')

#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F, repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(rtracklayer)); suppressMessages(library(ArchR))
suppressMessages(library(biomaRt)); suppressMessages(library(Seurat))
suppressMessages(library(SummarizedExperiment))

source('../hal_scripts/narrowPeakFunctions.R')

#########################################
### read in Corces2020 ArchR project ####
addArchRGenome("hg38")
PROJDIR2='../../../data/raw_data/hg38/Corces_2020'
ARCHDIR2=file.path(PROJDIR2,paste0('ArchR_Corces2020_caudate_labeled'))
humanProj = loadArchRProject(ARCHDIR2, showLogo = F)

# get peak matrix into SummarizedExperiment object
humanGeneScoreMat = getMatrixFromProject(
  ArchRProj = humanProj, useMatrix = "GeneScoreMatrix",  binarize = FALSE)
colData(humanGeneScoreMat)$Species = 'hg38'

##############################
### read in ArchR project ####
GENOMEDIR='/home/bnphan/resources/genomes/rheMac10'
load(file.path(GENOMEDIR,'rheMac10_liftoff_GRCh38.p13_ArchR_annotations.rda'))
PROJDIR4='../../../data/raw_data/rheMac10/Stauffer_caudate'
ARCHDIR4=file.path(PROJDIR4,paste0('ArchR_Stauffer_caudate_labeled'))
macaqueProj = loadArchRProject(ARCHDIR4, showLogo = F)

# get peak matrix into SummarizedExperiment object
macaqueGeneScoreMat = getMatrixFromProject(
  ArchRProj = macaqueProj, useMatrix = "GeneScoreMatrix",  binarize = FALSE)
colData(macaqueGeneScoreMat)$Species = 'rheMac10'

#############################################
### read in mouse BICCN CP ArchR project ####
addArchRGenome("mm10")
PROJDIR3='../../../data/raw_data/mm10/BICCN_mouse_caudoputamen'
ARCHDIR3=file.path(PROJDIR3,paste0('ArchR_BICCN_CP_labeled'))
mouseProj = loadArchRProject(ARCHDIR3, showLogo = F)

# get peak matrix into SummarizedExperiment object
mouseGeneScoreMat = getMatrixFromProject(
  ArchRProj = mouseProj, useMatrix = "GeneScoreMatrix",  binarize = FALSE)
colData(mouseGeneScoreMat)$Species = 'mm10'

## use biomaRt to convert mouse gene names to human gene names ###
mouse = useMart(dataset = "mmusculus_gene_ensembl", "ensembl")
human = useMart(dataset = "hsapiens_gene_ensembl", "ensembl")

grep('sapien', listAttributes(mart=mouse)[,1], value = T)
grep('symbol', listAttributes(mart=mouse)[,1], value = T)

mouseGenes = rowData(mouseGeneScoreMat)$name
humanGenes = intersect(rowData(humanGeneScoreMat)$name, 
                       rowData(macaqueGeneScoreMat)$name)
bm = getLDS(attributes=c('mgi_symbol'), filters = 'mgi_symbol',
           values = mouseGenes, mart = mouse,
           attributesL=c('hgnc_symbol'), filtersL = 'hgnc_symbol',
           valuesL= humanGenes, martL = human)

######################################################
# rearrange geneScoreMatrices by shared gene orthologs
humanGeneScoreMat = humanGeneScoreMat[
  match(bm$HGNC.symbol, rowData(humanGeneScoreMat)$name), ]

macaqueGeneScoreMat = macaqueGeneScoreMat[
  match(bm$HGNC.symbol, rowData(macaqueGeneScoreMat)$name), ]
rowData(macaqueGeneScoreMat) = rowData(humanGeneScoreMat)

mouseGeneScoreMat = mouseGeneScoreMat[
  match(bm$MGI.symbol, rowData(mouseGeneScoreMat)$name), ]
rowData(mouseGeneScoreMat) = rowData(humanGeneScoreMat)

## save the gene score matrices
matListRDS_fn = file.path(PROJDIR, 'rdas', 'geneScoreMergedMatrixList.rds')

matList = list(humanGeneScoreMat, macaqueGeneScoreMat, mouseGeneScoreMat)
system(paste('mkdir -p', dirname(matListRDS_fn)))
saveRDS(matList, compress = T, file = matListRDS_fn)


library(future)
plan("multiprocess", gc = TRUE, workers = round(parallel::detectCores()/4))
options(future.globals.maxSize = 10 * 1024^3)
# plan('sequential')

##############################################################
### combine the data & integrate w/ Seurat reciprocal PCA ####
matList = readRDS(matListRDS_fn)

counts <- do.call('cbind', lapply(matList, function(x) assays(x)[[1]]))
rownames(counts) = rowData(matList[[1]])$name
keepCols = Reduce('intersect',lapply(matList, function(x) names(colData(x))))
metadata <- do.call('rbind', lapply(matList, function(x) colData(x)[,keepCols]))
metadata = as.data.frame(metadata)
metadata$log10nFrags = log10(metadata$nFrags)

# merge the gene score matrices
obj_seurat <- CreateSeuratObject( counts, project = "CrossSpeciesGeneScores",
  assay = "RNA", names.field = 1, names.delim = "_", meta.data = metadata)

# use CPM normalization, no log transform
obj_seurat <- NormalizeData(obj_seurat, normalization.method = 'RC', 
                            scale.factor = 1e6, verbose = TRUE)

# scale data using NB rather than linear
obj_seurat <- ScaleData(obj_seurat, verbose = TRUE, split.by = 'Sample',
                        do.scale = T, do.center = F, model.use = 'poisson',
                        vars.to.regress = c('log10nFrags','TSSEnrichment'))

# regular HVG, PCA, and such
obj_seurat <- FindVariableFeatures(obj_seurat, verbose = FALSE)
obj_seurat <- RunPCA(obj_seurat, verbose = FALSE)
obj_seurat <- RunUMAP(obj_seurat, dims = 1:30, verbose = FALSE)

## save the seurat object
geneScoreMergedRDS_fn = file.path(PROJDIR, 'rdas', 'geneScoreMergedSeurat.rds')
saveRDS(obj_seurat, file = geneScoreMergedRDS_fn)

#######################################################
# split up by sample for reciprocal PCA integration ##
obj_seurat.list <- SplitObject(obj_seurat, split.by = "Sample")
obj_seurat.list <- lapply(X = obj_seurat.list, FUN = function(x) {
  # x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = obj_seurat.list)
obj_seurat.list <- lapply(X = obj_seurat.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(
  # using one mouse and human as reference
  object.list = obj_seurat.list, reference = c(2,5), 
  reduction = "rpca",dims = 1:30)

obj_seurat.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
obj_seurat.integrated <- ScaleData(obj_seurat.integrated, verbose = FALSE, 
                                   vars.to.regress = c('nFrags'))
obj_seurat.integrated <- RunPCA(obj_seurat.integrated, verbose = FALSE)
obj_seurat.integrated <- RunUMAP(obj_seurat.integrated, dims = 1:30)

geneScoreRDS_fn = file.path(PROJDIR, 'rdas', 'geneScoreIntegratedSeurat.rds')
saveRDS(obj_seurat.integrated, file = geneScoreRDS_fn)


