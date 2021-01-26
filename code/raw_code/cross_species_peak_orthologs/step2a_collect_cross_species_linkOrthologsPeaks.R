# to be run in the root github directory
setwd('code/raw_code/cross_species_peak_orthologs')
PROJDIR=file.path('../../../data/raw_data/cross_species_peak_orthologs')
LABEL = 'HumanMergedMultiSpeciesOrthologs'
GENOME = 'hg38'
SOURCE_SPECIES = 'Human'; TARGET_SPECIES = c('Rhesus','Mouse')


#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F, repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(rtracklayer)); suppressMessages(library(ArchR))
suppressMessages(library(Signac)); suppressMessages(library(Seurat))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(EnsDb.Hsapiens.v79))

source('../hal_scripts/narrowPeakFunctions.R')
source('../hal_scripts/seuratFunctions.R')

#################################################################
### collect human peak files that mapped to mouse and macaque ###
ortholog_fn = list.files(file.path(PROJDIR, 'peaks'), full.names = T,
                         pattern = '_orthologPeakList.rds')
names(ortholog_fn) = ss(basename(ortholog_fn), '_orthologPeakList',1)
orthologPeakList = lapply(ortholog_fn, readRDS)
celltypes = Reduce('intersect', lapply(orthologPeakList, names))
names(celltypes) = celltypes
orthologPeakList = lapply(celltypes, function(x){
  lapply(orthologPeakList, `[[`, x)
})


################################################
# interset narrowPeak sets with human-mouse and human-macaque orthologs
linkOrthologPeakList = lapply(orthologPeakList, linkNarrowPeaks, min.gapwidth = 251)
linkOrthologPeaks = linkOrthologPeakList[['Consensus']]
linkPeakNameList = lapply(linkOrthologPeaks, function(x) mcols(x)$name)


##########################################
### read in Corces2020 ArchR project ####
addArchRGenome("hg38")
PROJDIR2='../../../data/raw_data/hg38/Corces_2020'
ARCHDIR2=file.path(PROJDIR2,paste0('ArchR_Corces2020_caudate_labeled'))
humanProj = loadArchRProject(ARCHDIR2, showLogo = F)

# get peak matrix into SummarizedExperiment object
humanPeaks = nameNarrowPeakRanges(getPeakSet(humanProj), genome = 'hg38')
humanOrthPeakMat = getMatrixFromProject(
  ArchRProj = humanProj, useMatrix = "PeakMatrix",  binarize = FALSE)
colData(humanOrthPeakMat)$Species = 'hg38'

# add peak annotations and subset to linked peaks
rowRanges(humanOrthPeakMat) = humanPeaks
humanOrthPeakMat = humanOrthPeakMat[
  match(linkPeakNameList[['Corces_2020_caudate']], 
        rowRanges(humanOrthPeakMat)$name),]

round(table(humanPeaks$peakType)/10)
table(rowRanges(humanOrthPeakMat)$peakType)


#############################################
### read in mouse BICCN CP ArchR project ####
addArchRGenome("mm10")
PROJDIR3='../../../data/raw_data/mm10/BICCN_mouse_caudoputamen'
ARCHDIR3=file.path(PROJDIR3,paste0('ArchR_BICCN_CP_labeled'))
mouseProj = loadArchRProject(ARCHDIR3, showLogo = F)

# get peak matrix into SummarizedExperiment object
mousePeaks = nameNarrowPeakRanges(getPeakSet(mouseProj),genome = 'mm10')
mouseOrthPeakMat = getMatrixFromProject(
  ArchRProj = mouseProj, useMatrix = "PeakMatrix",  binarize = FALSE)
colData(mouseOrthPeakMat)$Species = 'mm10'

# add peak annotations and subset to linked peaks
rowRanges(mouseOrthPeakMat) = mousePeaks
mouseOrthPeakMat = mouseOrthPeakMat[
  match(linkPeakNameList[['BICCN_CP']], rowRanges(mouseOrthPeakMat)$name),]

round(table(mousePeaks$peakType)/10)
table(rowRanges(mouseOrthPeakMat)$peakType)


##############################
### read in ArchR project ####
GENOMEDIR='/home/bnphan/resources/genomes/rheMac10'
load(file.path(GENOMEDIR,'rheMac10_liftoff_GRCh38.p13_ArchR_annotations.rda'))
PROJDIR4='../../../data/raw_data/rheMac10/Stauffer_caudate'
ARCHDIR4=file.path(PROJDIR4,paste0('ArchR_Stauffer_caudate_labeled'))
macaqueProj = loadArchRProject(ARCHDIR4, showLogo = F)

# get peak matrix into SummarizedExperiment object
macaquePeaks = nameNarrowPeakRanges(getPeakSet(macaqueProj),genome = 'rheMac10')
macaqueOrthPeakMat = getMatrixFromProject(
  ArchRProj = macaqueProj, useMatrix = "PeakMatrix",  binarize = FALSE)
colData(macaqueOrthPeakMat)$Species = 'rheMac10'

# add peak annotations and subset to linked peaks
rowRanges(macaqueOrthPeakMat) = macaquePeaks
macaqueOrthPeakMat = macaqueOrthPeakMat[
  match(linkPeakNameList[['Stauffer_caudate']], rowRanges(macaqueOrthPeakMat)$name),]

round(table(macaquePeaks$peakType)/10)
table(rowRanges(macaqueOrthPeakMat)$peakType)

####################################
# combine and save count matrices ##
system(paste('mkdir -p', dirname(peakListRDS_fn)))

matListRDS_fn = file.path(PROJDIR, 'rdas', 'linkOrthologMatrixList.rds')
peakListRDS_fn = file.path(PROJDIR, 'rdas', 'linkOrthologPeakList.rds')

matList = list(humanOrthPeakMat, macaqueOrthPeakMat, mouseOrthPeakMat)
saveRDS(matList, compress = T, file = matListRDS_fn)
saveRDS(linkOrthologPeakList, compress = T, file = peakListRDS_fn)


#####################################################
### combine the data & make Signac/Seurat object ####
matListRDS_fn = file.path(PROJDIR, 'rdas', 'linkOrthologMatrixList.rds')
matList = readRDS(matListRDS_fn)

counts <- do.call('cbind', lapply(matList, function(x) assays(x)[[1]]))
rownames(counts) = GRangesToString(rowRanges(matList[[1]]), sep = c(":", "-"))
keepCols = Reduce('intersect',lapply(matList, function(x) names(colData(x))))
metadata <- do.call('rbind', lapply(matList, function(x) colData(x)[,keepCols]))
metadata = as.data.frame(metadata)

obj_seurat  = makeChomAssaySeurat(
  counts, metadata, EnsDb.Hsapiens.v79, 
  group.by = 'Sample', min.cells = 10, verbose = TRUE)

## save the seurat object
saveRDS_fn = file.path(PROJDIR, 'rdas', 'linkOrthologMatrixSeurat.rds')
saveRDS(obj_seurat, file = saveRDS_fn)


################################
## integration across species ##
obj_seurat.list = SplitObject(obj_seurat, split.by = 'Sample')
obj_seurat.list = lapply(obj_seurat.list, function(x){
  x <- RunTFIDF(x)
  x <- RunSVD(x)
})

# find integration anchors between species
anchors <- FindIntegrationAnchors(
  object.list = obj_seurat.list, reduction = 'cca',
  anchor.features = rownames(obj_seurat),
  assay = rep('peaks', length(obj_seurat.list)),
  k.filter = NA
)

# integrate data and create a new merged object
integrated <- IntegrateData(
  anchorset = anchors,
  dims = 2:30,
  preserve.order = TRUE
)

## save the seurat object
integratedRDS_fn = file.path(PROJDIR,'rdas','linkOrthologIntegratedSeurat.rds')
saveRDS(integrated, file = integratedRDS_fn)



