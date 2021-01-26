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

source('../hal_scripts/narrowPeakFunctions.R')

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
# merge narrowPeak sets with human-mouse and human-macaque orthologs
mergedOrthologPeakList = lapply(orthologPeakList, mergeNarrowPeaks, 
                                width = 501, min.gapwidth = 251)
lengths(mergedOrthologPeakList)

### write orthologous peaks to narrowPeak.gz ###
mergedOrthologNarrowPeak_fn = file.path(PROJDIR,'peaks',
        paste0(LABEL,'.',names(mergedOrthologPeakList), '.narrowPeak.gz'))

outList = mapply(write_GRangesToNarrowPeak, gr = mergedOrthologPeakList, 
                 file = mergedOrthologNarrowPeak_fn, genome = GENOME)


#####################################
# halLiftOver and HALPER the peaks ##
halmapper_script = '../hal_scripts/halper_map_narrowPeaks.sh'
system('mkdir -p logs')
sbatch = 'sbatch -p pool1'
target_species = paste('-t', TARGET_SPECIES)
source_species = paste('-s', SOURCE_SPECIES)
outdir = paste('-o', file.path(PROJDIR, 'halper'))
peak_files = paste('-p',mergedOrthologNarrowPeak_fn)

# paste the parameter calls together
thecall = paste(sbatch, halmapper_script, 
                source_species, 
                rep(target_species, each= length(peak_files)), 
                outdir, 
                rep(peak_files, times = length(target_species)))
cat(thecall, file= 'step1b_run_halmaper.sh', sep = '\n')
system('chmod u+x step1b_run_halmaper.sh')

