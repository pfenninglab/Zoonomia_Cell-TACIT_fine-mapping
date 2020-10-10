suppressMessages(library(ArchR))
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

##################################
### set Arrow File parameters ####
addArchRThreads(threads = 32)
addArchRGenome("hg38")
PROJDIR='../../../data/raw_data/hg38/Corces_2020'
ARCHDIR=file.path(PROJDIR,'ArchR_Corces2020_caudate')

### find the pre-processed arrow files
ArrowFiles = list.files(path = file.path(PROJDIR,'arrow'),
                          pattern = '.arrow', full.names = TRUE)
print(ArrowFiles)

##################################
### make an ArchR Projects ######
proj = ArchRProject( ArrowFiles = ArrowFiles,
                     outputDirectory = ARCHDIR,
                     copyArrows = FALSE, showLogo = FALSE)
proj$Subject = ss(basename(proj$Sample),'\\.', 1)
proj$Region = ss(basename(proj$Sample),'\\.', 2)
proj = filterDoublets( ArchRProj = proj,  cutEnrich = .5)

# add iterative LSI
set.seed(1234)
# add iterative LSI
proj <- addIterativeLSI(
  ArchRProj = proj, useMatrix = "TileMatrix", 
  name = "IterativeLSI",
  LSIMethod = 2, #"tf-logidf","log(tf-idf)", "logtf-logidf"
  iterations = 6, # increase this if noticing subtle batch effects
  scaleTo = 3000,
  selectionMethod = 'var',
  clusterParams = list( # See Seurat::FindClusters
  resolution = c(.1, .2, rep(.4, 3)), # lower this if noticing subtle batch effects
  sampleCells = 10000,  n.start = 10), 
  varFeatures = 150000, # also can reduce this if noticing subtle batch effects
  dimsToUse = 1:30, force = TRUE)


# add harmony batch correction #
proj <- addHarmony( ArchRProj = proj, reducedDims = "IterativeLSI",
                    name = "Harmony", groupBy = "Sample",force = TRUE)

# add imputation
proj <- addImputeWeights(proj, reducedDims = "Harmony")


# add umap
proj <- addUMAP( ArchRProj = proj, reducedDims = "Harmony", 
                 name = "UMAP", nNeighbors = 30, minDist = 0.5, 
                 metric = "cosine", force = TRUE)

# add clusters
proj <- addClusters(input = proj, reducedDims = "Harmony", 
                    method = "Seurat", name = "Clusters", 
                    resolution = .5, force = TRUE)
proj = saveArchRProject(ArchRProj = proj)



