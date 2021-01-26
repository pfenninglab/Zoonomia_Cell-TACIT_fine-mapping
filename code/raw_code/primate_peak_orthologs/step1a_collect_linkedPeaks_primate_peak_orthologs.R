# to be run in the root github directory
setwd('code/raw_code/primate_peak_orthologs')
PROJDIR=file.path('../../../data/raw_data/primate_peak_orthologs')
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
### collect human peak files that mapped to mouse and macaque ###
DATADIR=file.path('../../../data/raw_data/cross_species_peak_orthologs')


### get the peaks from human experiment and halpered peaks ####
HUMAN_DATADIR=file.path('../../../data/raw_data/hg38/Corces_2020')
humanNarrowPeak_fn = c(
  list.files(file.path(HUMAN_DATADIR,'peak'), full.names = T, 
             pattern ='Consensus.narrowPeak.gz'), 
  list.files(file.path(HUMAN_DATADIR,'halper'), full.names = T,
             pattern ='Consensus_HumanToRhesus_HALPER.narrowPeak.gz'))
names(humanNarrowPeak_fn) = ss(basename(humanNarrowPeak_fn), '\\.',2)
humanPeakList = lapply(humanNarrowPeak_fn, import)
humanOrthologPeakRanges = keepOrthologs(humanPeakList, idxReturn = 1)
length(humanOrthologPeakRanges)

### get the peaks from macaque experiment and halpered peaks ####
MACAQUE_DATADIR=file.path('../../../data/raw_data/rheMac10/Stauffer_caudate')
macaqueNarrowPeak_fn = c(
  # put file w/ hg38 coordinates first
  list.files(file.path(MACAQUE_DATADIR,'halper'), full.names = T,
             pattern ='Consensus_RhesusToHuman_HALPER.narrowPeak.gz'),
  list.files(file.path(MACAQUE_DATADIR,'peak'), full.names = T, 
             pattern ='Consensus.narrowPeak.gz'))
names(macaqueNarrowPeak_fn) = ss(basename(macaqueNarrowPeak_fn), '\\.',2)
macaquePeakList = lapply(macaqueNarrowPeak_fn, import)
macaqueOrthologPeakRanges = keepOrthologs(macaquePeakList, idxReturn = 1)
length(macaqueOrthologPeakRanges)

### link primate narrowpeaks using peak name from original dataset
orthologPeakList = list(Corces_2020_caudate = humanOrthologPeakRanges,
                        Stauffer_caudate = macaqueOrthologPeakRanges)
linkOrthologPeaks =  linkNarrowPeaks(orthologPeakList, min.gapwidth = 251)
linkPeakNameList = lapply(linkOrthologPeaks, function(x) mcols(x)$name)
lengths(linkOrthologPeaks)
# Corces_2020_caudate   \ Stauffer_caudate 
#             136780                136780


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

round(table(humanPeaks$peakType)/4)
table(rowRanges(humanOrthPeakMat)$peakType)

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

round(table(macaquePeaks$peakType)/3)
table(rowRanges(macaqueOrthPeakMat)$peakType)

####################################
# combine and save count matrices ##
matListRDS_fn = file.path(PROJDIR, 'rdas', 'linkOrthologMatrixList.rds')
peakListRDS_fn = file.path(PROJDIR, 'rdas', 'linkOrthologPeakList.rds')
system(paste('mkdir -p', dirname(peakListRDS_fn)))

matList = list(humanOrthPeakMat, macaqueOrthPeakMat)
saveRDS(matList, compress = T, file = matListRDS_fn)
saveRDS(linkOrthologPeakList, compress = T, file = peakListRDS_fn)


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






