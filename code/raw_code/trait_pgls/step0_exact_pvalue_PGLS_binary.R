ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(here)
library(rcartocolor)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)
library(tidyverse)
library(tidymodels)
library(statmod)
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))

CODEDIR='code/raw_code/trait_pgls'
DATADIR='data/raw_data/trait_pgls'
PLOTDIR='figures/exploratory/trait_pgls'
i_am(file.path(CODEDIR, 'step0_exact_pvalue_PGLS_binary.R'))

###########################################################
## 1a) read in the number of species peak is mappable
mammalTraits_df = here('data/raw_data/ldsc_caudate_zoonomia','traits_pgls',
                       'traitAnnotations_caudate.csv') %>% fread() %>%
  rename_with(make.names) %>% dplyr::rename('Species' =  'Species.Name') %>%
  rowwise(Species) %>% dplyr::select(-contains('_eg')) %>%
  mutate(tmp = any(!is.na(c_across(ActivityPattern:Sleep.Total_daily_sleep_time.adult)))) %>%
  filter(tmp) %>% dplyr::select(c(Species, ActivityPattern)) %>% 
  filter(!is.na(ActivityPattern))


predictions_df = here('data/raw_data/ldsc_caudate_zoonomia','rdas', 
                      paste0('Corces2020.',c('MSN_D1', 'MSN_D2'),
                             '.allPeaksCalibWithSNPs.avgCNN.predictions.rds')) %>%
  lapply(readRDS) %>% data.table::rbindlist(idcol = 'celltype') %>%
  dplyr::rename('peakNames' = 'name' ) %>% 
  pivot_longer(-c(peakNames, celltype), values_to = 'score', names_to = 'Species') %>%
  filter(!is.na(score)) %>% inner_join(mammalTraits_df) %>%
  group_by(peakNames) %>% 
  summarise(
    nDiurnal = sum(ActivityPattern == 1),
    nNocturnal = sum(ActivityPattern == 0)) %>% ungroup()

###########################
## 1) read in pgls results
pgls_fn =  here('data/raw_data/ldsc_caudate_zoonomia/traits_pgls') %>%
  list.files(pattern = 'roundFinal.csv', full.names = T) %>% 
  grep(pattern = 'actPattern', value = T)
names(pgls_fn) = basename(pgls_fn)
pgls_df = pgls_fn %>% lapply(fread) %>% rbindlist(idcol = 'file')   %>%
  rename('peakNames' = 'Enhancer', 'PGLS_Pvalue' = 'Exp_Pvalue') %>%
  inner_join(predictions_df) %>%
  mutate(celltype = case_when(grepl('oligo',file) ~ 'Oligo', 
                              grepl('D1',file) ~ 'MSN_D1',
                              grepl('D2',file) ~ 'MSN_D2',
                              TRUE ~ ''),
         speciesSet = case_when(grepl('eg',file) ~ 'Euarchontoglires', 
                                grepl('boreo|all',file) ~ 'Boreoeutheria'),
         zooTrait = ss(file, '_', 2)) %>%
  group_by(peakNames) %>% 
  mutate( Perm_Pvalue = permp(x = round(PGLS_Pvalue * Trials), 
                            nperm = Trials, n1 = nDiurnal, n2 = nNocturnal)) %>%
  relocate(Perm_Pvalue, .after = PGLS_Pvalue) %>% ungroup()

## annotate peaks w/ nearby gene
pgls_gr = pgls_df %>% mutate(value = gsub('^hg38:|:250$','',peakNames), name = peakNames) %>%
  dplyr::select(c(name, value)) %>% distinct(name, .keep_all = TRUE) %>% deframe() %>% GRanges()
annot_df <- annotatePeak(pgls_gr, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, 
                         annoDb='org.Hs.eg.db') %>% 
  as.GRanges() %>% as.data.frame(row.names = seq(length(.))) %>%
  mutate(annotation = ss(annotation, ' ')) %>%
  rename('Gene.Symbol' = 'SYMBOL') %>% 
  mutate(peakNames = paste0('hg38:', seqnames,':', start, '-', end, ':250')) %>%
  dplyr::select(peakNames, annotation, Gene.Symbol, distanceToTSS)

## filter out peaks not in distal intergenic/introns
pgls_df = right_join(pgls_df, annot_df, by = 'peakNames') %>%
  relocate(c(annotation, Gene.Symbol, distanceToTSS), .after = 'zooTrait') %>%
  filter(annotation %in% c('Distal', 'Intron'))

## annotate peaks w/ nearby gene
pgls_df = pgls_df %>%
  group_by(zooTrait, speciesSet, celltype)  %>% dplyr::select(-file) %>%
  mutate(FDR = p.adjust(Pvalue, 'fdr'), 
         PGLS_FDR = p.adjust(PGLS_Pvalue, 'fdr'),
         Perm_FDR = p.adjust(Perm_Pvalue, 'fdr')) %>%
  ungroup() %>% arrange(PGLS_FDR) %>% 
  relocate(c( PGLS_Pvalue, PGLS_FDR, Perm_Pvalue, Perm_FDR, Pvalue, FDR), .after = Coeff)

sum(pgls_df$PGLS_FDR < 0.05)
sum(pgls_df$Perm_FDR < 0.10)