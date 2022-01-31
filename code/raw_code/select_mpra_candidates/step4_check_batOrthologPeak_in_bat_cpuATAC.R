#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
suppressMessages(library(tidyverse))
suppressMessages(library(GenomicRanges))
suppressMessages(library(here))
library(rtracklayer)

source(here('code/raw_code/hal_scripts/narrowPeakFunctions.R'))
source(here('code/raw_code/hal_scripts/gen_enh_ortholog_sets.R'))
DATADIR='data/raw_data/select_mpra_candidates'

## Map Zoonomia RouAeg4 coordinates to VGP RouAeg1
peak_zooBat = import(here(DATADIR, 'peaks', 'addiction_snps_in_bat_OCRs.bed'))

# map rheMac8 UCSC chr to GenBank chromosome names (for halper) ##
CHRDIR = '/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/data/tidy_data/Zoonomia_data/tables'
chrmap_fn = here(CHRDIR, 'GCA_004024865.1_RouAeg_v1_BIUU_assembly_report.txt.gz')
map = read_tsv(chrmap_fn, skip = 29) %>%
  rename_with(gsub, pattern ='# ',replacement =  '') %>%
  rename_with(make.names) %>% filter(GenBank.Accn != 'na') %>% 
  select(`Sequence.Name`, GenBank.Accn) %>% 
  mutate(GenBank.Accn = ss(GenBank.Accn, '\\.')) %>%
  deframe() 

peak_zooBat2 = peak_zooBat %>% as.data.frame() %>%
  mutate(seqnames = map[as.character(seqnames)]) %>% GRanges()

#########################################################
## liftover peaks to  (the version in Cactus hal file) ##
chainFile =file.path("/home/bnphan/resources/liftOver_chainz", 
                     paste0('HLrouAeg4.HLrouAeg1_NCBIchrs.over.chain'))
ch = import.chain(chainFile, exclude = 'MT')

bat_str_peaks_fn = c(StrC = '/projects/MPRA/Simone/Bats/StrC/atac_out/atac/d73dadff-f69e-4fb3-8f70-652d6fce693d/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak.gz', 
                     StrP = '/projects/MPRA/Simone/Bats/StrP/atac_out/atac/878a0bdd-f8e2-47c5-ac19-a8d89973ae7e/call-reproducibility_idr/execution/idr.optimal_peak.narrowPeak.gz')

## import Bat ATAC StrP and StrC peaks
batAtac_peaks_vgp_gr = lapply(bat_str_peaks_fn, import) %>% GRangesList() %>% unlist()
batAtac_peaks_vgp_gr$tissue = names(batAtac_peaks_vgp_gr)
mcols(batAtac_peaks_vgp_gr)$name = paste0(genome,':', seqnames(batAtac_peaks_vgp_gr),':', 
         start(batAtac_peaks_vgp_gr),'-',end(batAtac_peaks_vgp_gr), 
         ':',mcols(batAtac_peaks_vgp_gr)$peak)

## lift over VGP to Zoonomia RouAeg
batAtac_peaks_vgp_grList = split(batAtac_peaks_vgp_gr, batAtac_peaks_vgp_gr$tissue)
bat_str_peaks_fn = here(DATADIR, 'peaks', 
                        paste0('MEW_bat1_bat3.', names(batAtac_peaks_vgp_grList),'.idr.optimal_peak.narrowPeak.gz' ))
mapply( write_GRangesToNarrowPeak,batAtac_peaks_vgp_grList, bat_str_peaks_fn )

batAtac_peaks_zoo_gr = liftOver_narrowPeak(batAtac_peaks_vgp_gr, chainFile = chainFile)


#############################################################
## find overlap of human -> bat ortholog w/ Bat ATAC peaks ##
oo = findOverlaps(subject = batAtac_peaks_zoo_gr, query = peak_zooBat2)
overlappedSet = batAtac_peaks_zoo_gr[subjectHits(oo)] # each ortholog has 6/7 Bat ATAC hits
overlappedSetList = split(overlappedSet, overlappedSet$tissue)

bat_str_peaks_fn2 = here(DATADIR, 'peaks', 
                        paste0('BNP_D1D2PPRAenhancers_overlapWithMEWBatATAC.', 
                               names(overlappedSetList),'.narrowPeak.gz' ))
mapply( write_GRangesToNarrowPeak,overlappedSetList, bat_str_peaks_fn2 )


