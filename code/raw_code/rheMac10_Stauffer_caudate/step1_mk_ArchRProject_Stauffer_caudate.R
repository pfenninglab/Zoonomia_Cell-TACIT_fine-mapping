suppressMessages(library(ArchR))
library(tidyverse)
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

##################################
### load rheMac10 ArchR genome ###
GENOMEDIR='/home/bnphan/resources/genomes/rheMac10'
load(file.path(GENOMEDIR,'rheMac10_liftoff_GRCh38.p13_ArchR_annotations.rda'))

##################################
### set Arrow File parameters ####
addArchRThreads(threads = 15)

### find the pre-processed arrow files
PROJDIR='../../../data/raw_data/rheMac10/Stauffer_caudate'
ARCHDIR=file.path(PROJDIR,'ArchR_Stauffer_caudate')
ArrowFiles = list.files(path = file.path(PROJDIR,'arrow'),
                        pattern = '.arrow', full.names = TRUE)
print(ArrowFiles)

##################################
### make an ArchR Projects ######
proj = ArchRProject( 
  ArrowFiles = ArrowFiles, outputDirectory = ARCHDIR,
  geneAnnotation = geneAnnotation, #must be the custom rheMac10 version
  genomeAnnotation = genomeAnnotation, #must be the custom rheMac10 version
  copyArrows = TRUE, showLogo = FALSE)

proj$Subject = ss(basename(proj$Sample),'_', 3)
proj$Region = ss(basename(proj$Sample),'_', 1)

proj = filterDoublets( ArchRProj = proj,  cutEnrich = .5)

####################
# add iterative LSI
set.seed(1234)
# add iterative LSI
proj <- addIterativeLSI( 
  ArchRProj = proj, useMatrix = "TileMatrix", 
  name = "IterativeLSI",
  LSIMethod = 2, #"tf-logidf","log(tf-idf)", "logtf-logidf"
  iterations = 4, # increase this if noticing subtle batch effects
  scaleTo = 3000,
  selectionMethod = 'var',
  clusterParams = list( # See Seurat::FindClusters
    resolution = .2, # lower this if noticing subtle batch effects
    sampleCells = 10000,  n.start = 10), 
  varFeatures = 150000, # also can reduce this if noticing subtle batch effects
  dimsToUse = 1:30, force = TRUE)


# add harmony batch correction #
proj <- addHarmony( ArchRProj = proj, reducedDims = "IterativeLSI",
                    name = "Harmony", groupBy = "Sample",force = TRUE)

# add imputation
proj <- addImputeWeights(proj, reducedDims = "IterativeLSI")


# add umap
proj <- addUMAP( ArchRProj = proj, reducedDims = "IterativeLSI", 
                 name = "UMAP", nNeighbors = 30, minDist = 0.5, 
                 metric = "cosine", force = TRUE)

# add clusters
proj <- addClusters(input = proj, reducedDims = "IterativeLSI", 
                    method = "Seurat", name = "Clusters", 
                    resolution = .5, force = TRUE)
proj = saveArchRProject(ArchRProj = proj)

