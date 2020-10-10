# to be run in the root github directory
setwd('code/raw_code/mm10_BICCN_mouse_caudoputamen')

options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(ArchR))
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

source('../hal_scripts/narrowPeakFunctions.R')

### read in ArchR project ####
PROJDIR='../../../data/raw_data/mm10/BICCN_mouse_caudoputamen'
LABEL='BICCN_CP'; GENOME = 'mm10'; 
SOURCE_SPECIES = 'Mouse'
TARGET_SPECIES = c('Human','Rhesus')
ARCHDIR=file.path(PROJDIR,paste0('ArchR_',LABEL,'_labeled'))

### set Arrow File parameters ####
# addArchRThreads(threads = 10)
addArchRGenome(GENOME)
proj = loadArchRProject(ARCHDIR, showLogo = F)

# get the reproducible peaks across Cluster2 cell types
peak_rds_fn = list.files(path = file.path(ARCHDIR,'PeakCalls'), full.names = T,
                         pattern = '.rds')
names(peak_rds_fn) = ss(basename(peak_rds_fn),'-reproduciblePeaks.gr.rds')
peakList = lapply(peak_rds_fn, readRDS)

# get the consensus peak matrix
peakList = c('Consensus' = getPeakSet(proj), peakList)
peakList = lapply(peakList, sort)

# create directory and 
PEAKDIR=file.path(PROJDIR,'peak')
system(paste('mkdir -p',PEAKDIR))
narrowPeak_fn = file.path(PEAKDIR, 
                          paste0(LABEL, '.',names(peakList),'.narrowPeak.gz'))

# write peaks to narrowPeak file
out = mapply(write_GRangesToNarrowPeak,gr = peakList, file = narrowPeak_fn, 
             genome = GENOME)

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


