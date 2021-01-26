# to be run in the root github directory
setwd('code/raw_code/cross_species_peak_orthologs')
PROJDIR=file.path('../../../data/raw_data/cross_species_peak_orthologs')

#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F, repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(rtracklayer)); suppressMessages(library(ArchR))

source('../../..code/raw_code/hal_scripts/narrowPeakFunctions.R')

#########################################################
### get the peaks from experiment and halpered peaks ####
narrowPeak_fn = list.files(file.path(PROJDIR,'peaks'), full.names = T,
                           pattern ='.narrowPeak.gz')
names(narrowPeak_fn) = ss(basename(narrowPeak_fn), '\\.',2)
halper_narrowPeak_fn = list.files(file.path(PROJDIR,'halper'), full.names = T,
                                  pattern ='HALPER.narrowPeak.gz')

peakList_fn = lapply(narrowPeak_fn, function(n){
  tmp = c(n, grep(ss(basename(n), '\\.',2), halper_narrowPeak_fn, value = T))
  val = grepl('ToHuman', tmp)
  if(any(val))
    # put the human coordinates first
    tmp = c(tmp[val], tmp[!val])
  names(tmp) = ss(basename(tmp), '\\.',2)
  return(tmp)
})

########################################################################
# get the consensus merged human ortholog that mapped to mouse & macaque
peakList = lapply(peakList_fn[['Consensus']], import)

# drop mappings to non-standard chromosomes
peakList = lapply(peakList, keepStandardChromosomes,pruning.mode="coarse")

# lift rheMac8 coordinates from the cactus to rheMac10
chainFile =file.path("/home/bnphan/resources/liftOver_chainz",
                     'rheMac8ToRheMac10.over.chain')
peakList[['Consensus_HumanToRhesus_HALPER']] = liftOver_narrowPeak(
  peakList[['Consensus_HumanToRhesus_HALPER']], chainFile = chainFile)

# reorder orthologs in hg38, mm10, and rheMac10 coordinates
peakList = keepOrthologs(peakList, idxReturn = NA)
lengths(peakList)

##########################################
### read in Corces2020 ArchR project ####
addArchRGenome("hg38")
PROJDIR2='../../../data/raw_data/hg38/Corces_2020'
ARCHDIR2=file.path(PROJDIR2,paste0('ArchR_Corces2020_caudate_labeled'))
humanProj = loadArchRProject(ARCHDIR2, showLogo = F)

# compute accessibility peak matrix in hg38 coordinates
humanProj = addFeatureMatrix( humanProj, matrixName = "OrthologPeakMatrix",
  features = peakList[["Consensus"]], ceiling = 10^9, binarize = F, force = T)
humanProj = saveArchRProject(humanProj)

# get peak matrix into SummarizedExperiment object
humanOrthPeakMat = getMatrixFromProject(
  ArchRProj = humanProj, useMatrix = "OrthologPeakMatrix",  binarize = FALSE)
colData(humanOrthPeakMat)$Species = 'hg38'

#############################################
### read in mouse BICCN CP ArchR project ####
addArchRGenome("mm10")
PROJDIR3='../../../data/raw_data/mm10/BICCN_mouse_caudoputamen'
ARCHDIR3=file.path(PROJDIR3,paste0('ArchR_BICCN_CP_labeled'))
mouseProj = loadArchRProject(ARCHDIR3, showLogo = F)

# compute accessibility peak matrix in mm10 coordinates
mouseProj = addFeatureMatrix( 
  mouseProj, matrixName ="OrthologPeakMatrix", ceiling = 10^9, binarize = F, 
  force = T, features = peakList[["Consensus_HumanToMouse_HALPER"]])
mouseProj = saveArchRProject(mouseProj)

# get peak matrix into SummarizedExperiment object
mouseOrthPeakMat = getMatrixFromProject(
  ArchRProj = mouseProj, useMatrix = "OrthologPeakMatrix",  binarize = FALSE)
colData(mouseOrthPeakMat)$Species = 'mm10'


##############################
### read in ArchR project ####
GENOMEDIR='/home/bnphan/resources/genomes/rheMac10'
load(file.path(GENOMEDIR,'rheMac10_liftoff_GRCh38.p13_ArchR_annotations.rda'))
PROJDIR4='../../../data/raw_data/rheMac10/Stauffer_caudate'
ARCHDIR4=file.path(PROJDIR4,paste0('ArchR_Stauffer_caudate_labeled'))
macaqueProj = loadArchRProject(ARCHDIR4, showLogo = F)

# compute accessibility peak matrix in rheMac10 coordinates
macaqueProj = addFeatureMatrix( 
  macaqueProj, matrixName ="OrthologPeakMatrix", ceiling = 10^9, binarize = F,
  force = T,features = peakList[["Consensus_HumanToRhesus_HALPER"]])
macaqueProj = saveArchRProject(macaqueProj)

# get peak matrix into SummarizedExperiment object
macaqueOrthPeakMat = getMatrixFromProject(
  ArchRProj = macaqueProj, useMatrix = "OrthologPeakMatrix",  binarize = FALSE)
colData(macaqueOrthPeakMat)$Species = 'rheMac10'

matList = list(humanOrthPeakMat, macaqueOrthPeakMat, mouseOrthPeakMat)

matListRDS_fn = file.path(PROJDIR, 'rdas', 'multispeciesMergedMatrixList.rds')
peakListRDS_fn = file.path(PROJDIR, 'rdas', 'multispeciesMergedPeakList.rds')
system(paste('mkdir -p', dirname(peakListRDS_fn)))
saveRDS(matList, compress = T, file = matListRDS_fn)
saveRDS(peakList, compress = T, file = peakListRDS_fn)


