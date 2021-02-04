# to be run in the root github directory
setwd('code/raw_code/cross_species_peak_orthologs')
PROJDIR=file.path('../../../data/raw_data/cross_species_peak_orthologs')

#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F, repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(rtracklayer)); suppressMessages(library(ArchR))

source('../hal_scripts/narrowPeakFunctions.R')

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
chainFile =file.path("/home/bnphan/resources/liftOver_chainz", 'rheMac8ToRheMac10.over.chain')
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

# get peak matrix into SummarizedExperiment object
humanProj$pseudoGroup = paste(humanProj$Sample, humanProj$Clusters2, sep= "#")
humanOrthPeakSE = getGroupSE( ArchRProj = humanProj, divideN = F, scaleTo = NULL,
                              useMatrix = 'OrthologPeakMatrix', groupBy = "pseudoGroup")
colData(humanOrthPeakSE)$Species = 'hg38'
colData(humanOrthPeakSE)$Sample = ss(colnames(humanOrthPeakSE), '#', 1)
colData(humanOrthPeakSE)$Clusters2 = ss(colnames(humanOrthPeakSE), '#', 2)


#############################################
### read in mouse BICCN CP ArchR project ####
addArchRGenome("mm10")
PROJDIR3='../../../data/raw_data/mm10/BICCN_mouse_caudoputamen'
ARCHDIR3=file.path(PROJDIR3,paste0('ArchR_BICCN_CP_labeled'))
mouseProj = loadArchRProject(ARCHDIR3, showLogo = F)

# get peak matrix into SummarizedExperiment object
mouseProj$pseudoGroup = paste(mouseProj$Sample, mouseProj$Clusters2, sep= "#")
mouseOrthPeakSE = getGroupSE( ArchRProj = mouseProj, divideN = F, scaleTo = NULL,
                              useMatrix = 'OrthologPeakMatrix', groupBy = "pseudoGroup")
colData(mouseOrthPeakSE)$Species = 'mm10'
colData(mouseOrthPeakSE)$Sample = ss(colnames(mouseOrthPeakSE), '#', 1)
colData(mouseOrthPeakSE)$Clusters2 = ss(colnames(mouseOrthPeakSE), '#', 2)
colData(mouseOrthPeakSE) = colData(mouseOrthPeakSE)[,keepCols]
mcols(mouseOrthPeakSE) = mcols(humanOrthPeakSE)

##############################
### read in ArchR project ####
GENOMEDIR='/home/bnphan/resources/genomes/rheMac10'
load(file.path(GENOMEDIR,'rheMac10_liftoff_GRCh38.p13_ArchR_annotations.rda'))
PROJDIR4='../../../data/raw_data/rheMac10/Stauffer_caudate'
ARCHDIR4=file.path(PROJDIR4,paste0('ArchR_Stauffer_caudate_labeled'))
macaqueProj = loadArchRProject(ARCHDIR4, showLogo = F)

# get peak matrix into SummarizedExperiment object
macaqueProj$pseudoGroup = paste(macaqueProj$Sample, macaqueProj$Clusters2, sep= "#")
macaqueOrthPeakSE = getGroupSE( ArchRProj = macaqueProj, divideN = F, scaleTo = NULL,
                              useMatrix = 'OrthologPeakMatrix', groupBy = "pseudoGroup")
colData(macaqueOrthPeakSE)$Species = 'rheMac10'
colData(macaqueOrthPeakSE)$Sample = ss(colnames(macaqueOrthPeakSE), '#', 1)
colData(macaqueOrthPeakSE)$Clusters2 = ss(colnames(macaqueOrthPeakSE), '#', 2)
mcols(macaqueOrthPeakSE) = mcols(humanOrthPeakSE)

# merge colData names
matList = list(humanOrthPeakSE, macaqueOrthPeakSE, mouseOrthPeakSE)
counts <- do.call('cbind', lapply(matList, function(x) assays(x)[[1]]))
rowData = rowData(matList[[1]])
keepCols = Reduce('intersect',lapply(matList, function(x) names(colData(x))))
colData <- do.call('rbind', lapply(matList, function(x) colData(x)[,keepCols]))
colData = as.data.frame(colData)
colData$log10nFrags = log10(colData$nFrags)

orthPeakRSE_fn = file.path(PROJDIR, 'rdas', 'orthologPeak_pseudoBulkRSE.rds')
orthPeakRSE = SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = counts), rowRanges=GRanges(rowData), colData=colData)
saveRDS(orthPeakRSE, file = orthPeakRSE_fn)

################################################################
### 