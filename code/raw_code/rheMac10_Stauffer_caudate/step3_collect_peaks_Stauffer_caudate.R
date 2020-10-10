# to be run in the root github directory
setwd('code/raw_code/rheMac10_Stauffer_caudate')

ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(ArchR))
suppressMessages(library(rtracklayer))

source('../hal_scripts/narrowPeakFunctions.R')

##############################
### read in ArchR project ####
PROJDIR='../../../data/raw_data/rheMac10/Stauffer_caudate'
LABEL='Stauffer_caudate'; GENOME = 'rheMac10'; 
SOURCE_SPECIES = 'Rhesus'
TARGET_SPECIES = c('Human','Mouse')
ARCHDIR=file.path(PROJDIR,paste0('ArchR_',LABEL,'_labeled'))

##################################
### set Arrow File parameters ####
# addArchRThreads(threads = 10)
# addArchRGenome(GENOME)
proj = loadArchRProject(ARCHDIR, showLogo = F)

#######################################################
# get the reproducible peaks across Cluster2 cell types
peak_rds_fn = list.files(path = file.path(ARCHDIR,'PeakCalls'), full.names = T,
                         pattern = '.rds')
names(peak_rds_fn) = ss(basename(peak_rds_fn),'-reproduciblePeaks.gr.rds')
peakList = lapply(peak_rds_fn, readRDS)

# get the consensus peak matrix
peakList = c('Consensus' = getPeakSet(proj), peakList)

# label the summit and peak name using rheMac10 coordinates
peakList = lapply(peakList, addSummitCenter)
peakList = lapply(peakList, nameNarrowPeakRanges, genome = GENOME)
peakList = lapply(peakList, sort)

## liftover peaks to rheMac8 (the version in Cactus hal file) ##
chainFile =file.path("/home/bnphan/resources/liftOver_chainz",
                     'rheMac10ToRheMac8.over.chain')
peakList_rheMac8 = lapply(peakList, liftOver_narrowPeak, chainFile = chainFile)

###############################################
# create directory and narrowPeak file names ##
PEAKDIR=file.path(PROJDIR,'peak')
system(paste('mkdir -p',PEAKDIR))
narrowPeak_fn = file.path(PEAKDIR, 
                    paste0(LABEL, '.', names(peakList_rheMac8), '.narrowPeak.gz'))

# write peaks to narrowPeak file
outRanges = mapply(write_GRangesToNarrowPeak,
                   gr = peakList_rheMac8, file = narrowPeak_fn, genome = GENOME)


#####################################
# halLiftOver and HALPER the peaks ##
halmapper_script = '../hal_scripts/halper_map_narrowPeaks.sh'
system('mkdir -p logs')
sbatch = 'sbatch -p pool1'
target_species = paste('-t', TARGET_SPECIES)
source_species = paste('-s', SOURCE_SPECIES)
outdir = paste('-o', file.path(PROJDIR, 'halper'))
peak_files = paste('-p',narrowPeak_fn)

# paste the parameter calls together
thecall = paste(sbatch, halmapper_script, 
                source_species, 
                rep(target_species, each= length(peak_files)), 
                outdir, 
                rep(peak_files, times = length(target_species)))
cat(thecall, file= 'step3b_run_halmaper.sh', sep = '\n')
system('chmod u+x step3b_run_halmaper.sh')




