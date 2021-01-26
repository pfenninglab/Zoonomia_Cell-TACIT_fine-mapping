# to be run in the root github directory
setwd('code/raw_code/hg38_Corces_2020_caudate')

suppressMessages(library(ArchR))
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

source('../hal_scripts/narrowPeakFunctions.R')

### set Arrow File parameters ####
addArchRThreads(threads = 10)
addArchRGenome("hg38")

### read in ArchR project ####
PROJDIR='../../../data/raw_data/hg38/Corces_2020'
LABEL='Corces2020_caudate'; GENOME = 'hg38'; 
SOURCE_SPECIES = 'Homo_sapiens'
TARGET_SPECIES = c('Macaca_mulatta','Mus_musculus')
ARCHDIR=file.path(PROJDIR,paste0('ArchR_',LABEL,'_labeled'))
proj = loadArchRProject(ARCHDIR, showLogo = F)

# get the reproducible peaks across Cluster2 cell types
peak_rds_fn = list.files(path = file.path(ARCHDIR,'PeakCalls'), full.names = T,
                         pattern = '.rds')
names(peak_rds_fn) = ss(basename(peak_rds_fn),'-reproduciblePeaks.gr.rds')
peakList = lapply(peak_rds_fn, readRDS)

# get the consensus peak matrix
peakList = c('Consensus' = getPeakSet(proj), peakList)

# add peak summits, unique peak names, and sort
peakList = lapply(peakList, addSummitCenter)
peakList = lapply(peakList, nameNarrowPeakRanges, genome = GENOME)
peakList = lapply(peakList, sort)

# create directory and 
PEAKDIR=file.path(PROJDIR,'peak')
system(paste('mkdir -p',PEAKDIR))
narrowPeak_fn = file.path(PEAKDIR, 
                          paste0(LABEL, '.',names(peakList),'.narrowPeak.gz'))

# write peaks to narrowPeak file
outList = mapply(write_GRangesToNarrowPeak,gr = peakList, file = narrowPeak_fn, 
             genome = GENOME)

#####################################
# halLiftOver and HALPER the peaks ##
halmapper_script = '../hal_scripts/halper_map_peak_orthologs.sh'
system('mkdir -p logs')
sbatch = 'sbatch -p pfen1 -w compute-1-40'
target_species = paste('-t', TARGET_SPECIES)
source_species = paste('-s', SOURCE_SPECIES)
outdir = paste('-o', file.path(PROJDIR, 'halper'))
peak_files = paste('-b',narrowPeak_fn)

# paste the parameter calls together
thecall = paste(sbatch, halmapper_script, 
                source_species, 
                rep(target_species, each= length(peak_files)), 
                outdir, 
                rep(peak_files, times = length(target_species)))
cat(thecall, file= 'step3b_run_halmaper.sh', sep = '\n')
system('chmod u+x step3b_run_halmaper.sh')



