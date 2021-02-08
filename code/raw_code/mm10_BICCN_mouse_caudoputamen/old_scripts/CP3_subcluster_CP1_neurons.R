suppressMessages(library(ArchR))
library(tidyverse)
ss <- function(x, pattern, slot = 1, ...) { 
 sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
setwd('code/raw_code/mm10_BICCN_mouse_caudoputamen')

##################################
### set Arrow File parameters ####
addArchRGenome("mm10")
PROJDIR='../../../data/raw_data/mm10/BICCN_mouse_caudoputamen'
ARCHDIR=file.path(PROJDIR,'ArchR_BICCN_CP1_Neuron')
proj = loadArchRProject(ARCHDIR)

###################
# add iterative LSI
varFeat = 70; set.seed(1234)
proj <- addIterativeLSI(
 proj, useMatrix = "TileMatrix", 
 name = paste0("IterativeLSI",varFeat),
 LSIMethod = 2, #"tf-logidf", "log(tf-idf)", "logtf-logidf"
 iterations = 4, # increase this if noticing subtle batch effects
 scaleTo = 3000,
 selectionMethod = 'var',
 clusterParams = list( # See Seurat::FindClusters
  resolution = c(.1, .2, .4, .4), # lower this if noticing subtle batch effects
  sampleCells = 10000, n.start = 10), 
 varFeatures = varFeat * 1000, # also can reduce this if noticing subtle batch effects
 dimsToUse = 1:40, force = TRUE)
proj = saveArchRProject(proj)

# add harmony batch correction #
proj <- addHarmony( proj, reducedDims = paste0("IterativeLSI",varFeat), 
                    name = paste0("Harmony",varFeat), groupBy = 'Sample',force = TRUE)

# add imputation
proj <- addImputeWeights(proj, reducedDims = paste0("Harmony",varFeat))
proj = saveArchRProject(proj)


# add umap
proj <- addUMAP( proj, reducedDims = paste0("Harmony",varFeat), 
         name = paste0("UMAP",varFeat), nNeighbors = 30, minDist = 0.5, 
         metric = "cosine", force = TRUE)

# add clusters
proj <- addClusters( proj, reducedDims = paste0("Harmony",varFeat), 
          method = "Seurat", name = paste0("Clusters",varFeat), 
          resolution = 1, force = TRUE)
proj = saveArchRProject(proj)

