# to be run in the root github directory
setwd('code/raw_code/rheMac10_Stauffer_caudate')

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(ArchR))
suppressMessages(library(rtracklayer))

source('../hal_scripts/narrowPeakFunctions.R')

##############################
### read in ArchR project ####
PROJDIR='../../../data/raw_data/rheMac10/Stauffer_caudate'
LABEL='Stauffer_caudate'; GENOME = 'rheMac10'; 
SOURCE_SPECIES = 'Macaca_mulatta'
TARGET_SPECIES = c('Homo_sapiens','Mus_musculus')
ARCHDIR=file.path(PROJDIR,paste0('ArchR_',LABEL,'_labeled'))

##################################
### set Arrow File parameters ####
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

##################################################################
# map rheMac8 UCSC chr to GenBank chromosome names (for halper) ##
library(tidyverse)
CHRDIR = '../../../data/tidy_data/Zoonomia_data/tables/'
chrmap_fn = file.path(CHRDIR, 'GCF_000772875.2_Mmul_8.0.1_assembly_report.txt.gz')
map = read_tsv(chrmap_fn, skip = 30) %>%
  rename_with(gsub, pattern ='# ',replacement =  '') %>%
  rename_with(make.names) %>%
  filter(Assigned.Molecule != 'na') %>%
  filter(GenBank.Accn != 'na') %>% 
  select(GenBank.Accn, UCSC.style.name) %>%
  column_to_rownames(var = "UCSC.style.name")

peakList_rheMac8_renamed = lapply(peakList_rheMac8, function(gr){
  ret = gr %>% as.data.frame() %>%
    filter(seqnames %in% rownames(map)) %>%
    mutate(seqnames = map[as.character(seqnames),'GenBank.Accn']) %>% 
    GRanges()
  return(ret)
}) %>% GRangesList()
seqnames(peakList_rheMac8_renamed[[2]])

###############################################
# create directory and narrowPeak file names ##
PEAKDIR=file.path(PROJDIR,'peak')
system(paste('mkdir -p',PEAKDIR))
narrowPeak_fn = file.path(PEAKDIR, paste0(LABEL, '.', names(peakList_rheMac8), 
                                          '.GenBankRheMac8.narrowPeak.gz'))
# write peaks to narrowPeak file
outRanges = mapply(write_GRangesToNarrowPeak,
                   gr = peakList_rheMac8_renamed, 
                   file = narrowPeak_fn, genome = GENOME)

sapply(outRanges, function(gr){
  tmp = gr %>% as.data.frame() %>%
    mutate(col = paste0(seqnames, start, end, peak, name)) %>% 
    filter(duplicated(col))
  return(nrow(tmp))
})

#####################################
# halLiftOver and HALPER the peaks ##
halmapper_script = '../hal_scripts/halper_map_peak_orthologs.sh'
system('mkdir -p logs')
sbatch = 'sbatch -p pfen1 -w compute-1-40 --mem 10G'
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
# system('./step3b_run_halmaper.sh')




