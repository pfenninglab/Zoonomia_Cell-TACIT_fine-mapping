suppressMessages(library(ArchR))
library(tidyverse)
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

##################################
### set Arrow File parameters ####
addArchRThreads(threads = 24)
addArchRGenome("mm10")
PROJDIR='../../../data/raw_data/mm10/BICCN_mouse_caudoputamen'

### find the pre-processed arrow files
ARCHDIR=file.path(PROJDIR,'ArchR_BICCN_CP')
ArrowFiles = list.files(path = file.path(PROJDIR,'arrow'),
                        pattern = '.arrow', full.names = TRUE)
print(ArrowFiles)

##################################
### make an ArchR Projects ######
proj = ArchRProject( ArrowFiles = ArrowFiles,
                     outputDirectory = ARCHDIR,
                     copyArrows = TRUE, showLogo = FALSE)
proj = filterDoublets( ArchRProj = proj,  cutEnrich = .5)

# read in the metadata and filter
meta_fn =file.path('/projects/pfenninggroup/singleCell/BICCN_mouse_CATlas_snATAC-seq',
                   'tables', 'BICCN_mouse_CATlas__snATAC-seq_metadata.RDS') 
meta_df = readRDS(meta_fn)
meta_df = meta_df %>% filter(SubRegion %in% c('CP'))
meta_df$cellID = gsub('\\.','#', meta_df$cellID)

# rearrange cells
cellsKeep = intersect(proj$cellNames, meta_df$cellID)
proj = proj[cellsKeep,]
meta_df = meta_df[match(cellsKeep, meta_df$cellID),]
all.equal(meta_df$cellID, proj$cellNames)

for(i in seq(ncol(meta_df))){
  proj = addCellColData(
    ArchRProj = proj,
    data = meta_df[,i],
    cells = cellsKeep,
    name = names(meta_df)[i],
    force = T
  )
}



####################
# add iterative LSI
set.seed(1234)
# add iterative LSI
proj <- addIterativeLSI(
  ArchRProj = proj, useMatrix = "TileMatrix", 
  name = "IterativeLSI",
  LSIMethod = 2, #"tf-logidf", "log(tf-idf)", "logtf-logidf"
  iterations = 6, # increase this if noticing subtle batch effects
  scaleTo = 3000,
  selectionMethod = 'var',
  clusterParams = list( # See Seurat::FindClusters
    resolution = c(.1, .2, rep(.4, 3)), # lower this if noticing subtle batch effects
    sampleCells = 10000,  n.start = 10), 
  varFeatures = 150000, # also can reduce this if noticing subtle batch effects
  dimsToUse = 1:40, force = TRUE)


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
                    resolution = 1, force = TRUE)
proj = saveArchRProject(ArchRProj = proj)

