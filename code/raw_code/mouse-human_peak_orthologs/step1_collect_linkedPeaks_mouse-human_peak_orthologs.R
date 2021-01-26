# to be run in the root github directory
setwd('code/raw_code/mouse-human_peak_orthologs')
PROJDIR=file.path('../../../data/raw_data/mouse-human_peak_orthologs')
LABEL = 'PrimateMergedMultiSpeciesOrthologs'
GENOME = 'hg38'

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
### get the peaks from human experiment and halpered peaks ####
HUMAN_DATADIR=file.path('../../../data/raw_data/hg38/Corces_2020')
humanNarrowPeak_fn = c(
  list.files(file.path(HUMAN_DATADIR,'peak'), full.names = T, 
             pattern ='Consensus.narrowPeak.gz'), 
  list.files(file.path(HUMAN_DATADIR,'halper'), full.names = T,
             pattern ='Consensus_HumanToMouse_HALPER.narrowPeak.gz'))
names(humanNarrowPeak_fn) = ss(basename(humanNarrowPeak_fn), '\\.',2)
humanPeakList = lapply(humanNarrowPeak_fn, import)
humanOrthologPeakRanges = keepOrthologs(humanPeakList, idxReturn = 1)
length(humanOrthologPeakRanges)
# 205236

### get the peaks from mouse experiment and halpered peaks ####
MOUSE_DATADIR=file.path('../../../data/raw_data/mm10/BICCN_mouse_caudoputamen')
mouseNarrowPeak_fn = c(
  # put file w/ hg38 coordinates first
  list.files(file.path(MOUSE_DATADIR,'halper'), full.names = T,
             pattern ='Consensus_MouseToHuman_HALPER.narrowPeak.gz'),
  list.files(file.path(MOUSE_DATADIR,'peak'), full.names = T, 
             pattern ='Consensus.narrowPeak.gz'))
names(mouseNarrowPeak_fn) = ss(basename(mouseNarrowPeak_fn), '\\.',2)
mousePeakList = lapply(mouseNarrowPeak_fn, import)
mouseOrthologPeakRanges = keepOrthologs(mousePeakList, idxReturn = 1)
length(mouseOrthologPeakRanges)
# 144722

### link primate narrowpeaks using peak name from original dataset
orthologPeakList = list(Corces_2020_caudate = humanOrthologPeakRanges,
                        BICCN_CP = mouseOrthologPeakRanges)
linkOrthologPeaks =  linkNarrowPeaks(orthologPeakList, min.gapwidth = 251)
linkPeakNameList = lapply(linkOrthologPeaks, function(x) mcols(x)$name)
lengths(linkOrthologPeaks)
# Corces_2020_caudate            BICCN_CP 
#               54796               54796 


##########################################
### read in Corces2020 ArchR project ####
addArchRGenome("hg38")
ARCHDIR2=file.path(HUMAN_DATADIR, paste0('ArchR_Corces2020_caudate_labeled'))
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

round(table(humanPeaks$peakType)/4)
table(rowRanges(humanOrthPeakMat)$peakType)

##############################
### read in ArchR project ####
addArchRGenome("mm10")
ARCHDIR4 = file.path(MOUSE_DATADIR,paste0('ArchR_BICCN_CP_labeled'))
mouseProj = loadArchRProject(ARCHDIR4, showLogo = F)

# get peak matrix into SummarizedExperiment object
mousePeaks = nameNarrowPeakRanges(getPeakSet(mouseProj),genome = 'mm10')
mouseOrthPeakMat = getMatrixFromProject(
  ArchRProj = mouseProj, useMatrix = "PeakMatrix",  binarize = FALSE)
colData(mouseOrthPeakMat)$Species = 'mm10'

# add peak annotations and subset to linked peaks
rowRanges(mouseOrthPeakMat) = mousePeaks
mouseOrthPeakMat = mouseOrthPeakMat[
  match(linkPeakNameList[['BICCN_CP']], rowRanges(mouseOrthPeakMat)$name),]

round(table(mousePeaks$peakType)/3)
table(rowRanges(mouseOrthPeakMat)$peakType)

####################################
# combine and save count matrices ##
matListRDS_fn = file.path(PROJDIR, 'rdas', 'linkOrthologMatrixList.rds')
peakListRDS_fn = file.path(PROJDIR, 'rdas', 'linkOrthologPeakList.rds')
system(paste('mkdir -p', dirname(peakListRDS_fn)))

matList = list(humanOrthPeakMat, mouseOrthPeakMat)
saveRDS(matList, compress = T, file = matListRDS_fn)
saveRDS(linkOrthologPeaks, compress = T, file = peakListRDS_fn)


#####################################################
### combine the data & make Signac/Seurat object ####
matListRDS_fn = file.path(PROJDIR, 'rdas', 'linkOrthologMatrixList.rds')
matList = readRDS(file = matListRDS_fn)
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
