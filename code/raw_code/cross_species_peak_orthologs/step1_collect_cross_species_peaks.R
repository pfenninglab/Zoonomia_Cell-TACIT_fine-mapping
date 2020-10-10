suppressMessages(library(ArchR))
suppressMessages(library(rtracklayer))
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
source('../hal_scripts/narrowPeakFunctions.R')

PROJDIR='../../../data/raw_data'

LABEL='HumanMergedOrthologs.'; GENOME = 'hg38'; TARGET_SPECIES = c('Rhesus','Mouse')
human_dir = file.path(PROJDIR,'hg38/Corces_2020')
mouse_dir = file.path(PROJDIR,'mm10/BICCN_mouse_caudoputamen')
macaque_dir = file.path(PROJDIR,'rheMac10/Stauffer_caudate')

#################################################################
### collect human peak files that mapped to mouse and macaque ###
hg38_peak_fn = list.files(file.path(human_dir,'peak'),
                          pattern ='narrowPeak.gz', full.names = T)
hg38ToMm10_peak_fn = list.files(file.path(human_dir,'halper'), full.names = T,
                                pattern ='Mouse_HALPER.narrowPeak.gz')
hg38ToMac10_peak_fn = list.files(file.path(human_dir,'halper'), full.names = T,
                                pattern ='Rhesus_HALPER.narrowPeak.gz')

## read in the peakfiles to be mapped
hg38PeakRanges = lapply(hg38_peak_fn,import)
names(hg38PeakRanges) = ss(basename(hg38_peak_fn),'\\.',2)

hg38ToMm10PeakRanges = lapply(hg38ToMm10_peak_fn,import)
names(hg38ToMm10PeakRanges) = ss(basename(hg38ToMm10_peak_fn),'\\.',2)

hg38ToMac10PeakRanges = lapply(hg38ToMac10_peak_fn,import)
names(hg38ToMac10PeakRanges) = ss(basename(hg38ToMac10_peak_fn),'\\.',2)

## get human peaks w/ mouse & macaque orthologs
hg38OrthologPeakRanges = mapply(function(human, mouse, macaque){
  orthologs = keepOrthologs(list(human, mouse, macaque), idxReturn = 1)
}, human = hg38PeakRanges, 
mouse = hg38ToMm10PeakRanges, 
macaque = hg38ToMac10PeakRanges)

names(hg38OrthologPeakRanges) = names(hg38PeakRanges)

## get percent of peaks survived mapping
lengths(hg38OrthologPeakRanges) / lengths(hg38PeakRanges) * 100


######################################################
### collect mouse peak files that mapped to human ###
mm10_peak_fn = list.files(file.path(mouse_dir,'peak'),
                          pattern ='narrowPeak.gz', full.names = T)
mm10ToHg38_peak_fn = list.files(file.path(mouse_dir,'halper'), full.names = T,
                                pattern ='HALPER.narrowPeak.gz')

## read in the peakfiles to be mapped
mm10PeakRanges = lapply(mm10_peak_fn,import)
names(mm10PeakRanges) = ss(basename(mm10_peak_fn),'\\.',2)

mm10ToHg38PeakRanges = lapply(mm10ToHg38_peak_fn,import)
names(mm10ToHg38PeakRanges) = ss(basename(mm10ToHg38_peak_fn),'\\.',2)

## get mouse peaks w/ human orthologs
mm10OrthologPeakRanges = mapply(function(human, mouse){
  orthologs = keepOrthologs(list(human, mouse), idxReturn = 1)
}, human = mm10ToHg38PeakRanges, 
mouse = mm10PeakRanges)
names(mm10OrthologPeakRanges) = names(mm10PeakRanges)

## get percent of peaks survived mapping
lengths(mm10OrthologPeakRanges) / lengths(mm10PeakRanges) * 100


######################################################
### collect macaque peak files that mapped to human ###
rheMac10_peak_fn = list.files(file.path(macaque_dir,'peak'),
                          pattern ='narrowPeak.gz', full.names = T)
rheMac10ToHg38_peak_fn = list.files(file.path(macaque_dir,'halper'), full.names = T,
                                pattern ='HALPER.narrowPeak.gz')

## read in the peakfiles to be mapped
rheMac10PeakRanges = lapply(rheMac10_peak_fn,import)
names(rheMac10PeakRanges) = ss(basename(rheMac10_peak_fn),'\\.',2)

rheMac10ToHg38PeakRanges = lapply(rheMac10ToHg38_peak_fn,import)
names(rheMac10ToHg38PeakRanges) = ss(basename(rheMac10ToHg38_peak_fn),'\\.',2)

## get macaque peaks w/ human orthologs
rheMac10OrthologPeakRanges = mapply(function(human, macaque){
  orthologs = keepOrthologs(list(human, macaque), idxReturn = 1)
}, human = rheMac10ToHg38PeakRanges, 
macaque = rheMac10PeakRanges)
names(rheMac10OrthologPeakRanges) = names(rheMac10PeakRanges)

## get percent of peaks survived mapping
lengths(rheMac10OrthologPeakRanges) / lengths(rheMac10PeakRanges) * 100


######################################################
### merge the peaks w/ orthologs from mouse, macaque, human ###
celltypes = Reduce('intersect',lapply(list(
  hg38OrthologPeakRanges, mm10OrthologPeakRanges, rheMac10OrthologPeakRanges), 
  names))

# merge narrowPeak sets with human-mouse and human-macaque orthologs
mergedOrthologPeaks = mapply(function(human, mouse, macaque){
  mergedOrthologs = mergeNarrowPeaks(
    list(human, mouse, macaque), width = 501, min.gapwidth = 100)}, 
human = hg38OrthologPeakRanges[celltypes], 
mouse = mm10OrthologPeakRanges[celltypes], 
macaque = rheMac10OrthologPeakRanges[celltypes])


################################################
### write orthologous peaks to narrowPeak.gz ###
OUTDIR=file.path(PROJDIR, 'multiSpeciesOrthologs','peaks')
mergedOrthologNarrowPeak_fn = file.path(
  OUTDIR,paste0(LABEL,names(mergedOrthologPeaks), '.narrowPeak.gz'))
outList = mapply(write_GRangesToNarrowPeak, x = mergedOrthologPeaks, 
                 file = mergedOrthologNarrowPeak_fn, genome = GENOME)


#####################################
# halLiftOver and HALPER the peaks ##
halmapper_script = '../hal_scripts/halper_map_narrowPeaks.sh'
system('mkdir -p logs')
sbatch = 'sbatch -p pfen_bigmem'
target_species = paste('-t', TARGET_SPECIES)
source_species = '-s Human'
outdir = paste('-o', file.path(PROJDIR, 'multiSpeciesOrthologs','halper'))
peak_files = paste('-p',mergedOrthologNarrowPeak_fn)

# paste the parameter calls together
thecall = paste(sbatch, halmapper_script, source_species, 
                rep(target_species, each= length(peak_files)), outdir, 
                rep(peak_files, times = length(target_species)))
cat(thecall, file= 'step1b_run_halmaper.sh', sep = '\n')
system('chmod u+x step1b_run_halmaper.sh')

