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
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))

CODEDIR='code/raw_code/trait_pgls'
DATADIR='data/raw_data/trait_pgls'
PLOTDIR='figures/exploratory/trait_pgls'
i_am(file.path(CODEDIR, 'step1_aggregate_trait_pgls_Corces2020_orthologs.R'))

###########################################################
## 1a) read in the number of species peak is mappable
outPeaks_numMappable = here('data/raw_data/ldsc_caudate_zoonomia','rdas',
                            'Corces2020.allPeaks.numMappable.rds')
if(!file.exists(outPeaks_numMappable)){
  outPeaks_rds = here('data/raw_data/ldsc_caudate_zoonomia','rdas') %>%
    list.files(pattern ='.allPeaks.avgCNN.predictions.rds', full.names = T)
  peaksMappable = data.frame()
  for(file in outPeaks_rds){
    print(paste('Working with:', basename(file)))
    input = file %>% readRDS() %>% 
      filter(grepl('hg38', name)) %>%
      mutate_if(is.numeric, ~ !is.na(.x)) %>%
      rowwise(name) %>%
      summarise(numMappable = sum(c_across(Acinonyx_jubatus:Ziphius_cavirostris))) 
    peaksMappable = bind_rows(peaksMappable, input)
  }
  saveRDS(peaksMappable, file = outPeaks_numMappable)
  rm(input); gc()
} 
peaksMappable = readRDS(outPeaks_numMappable) %>% rename('peakNames' = 'name')

## 1b) read in table of fine-mapped SNPs overlapping OCRs ##
allPeaksSNPs_fn = here('data/raw_data/polyfun_caudate/rdas',
                       'polyfun_caudate_finemapped_snpsInAllPeaks_20210518.rds')
snpInPeaks_df = allPeaksSNPs_fn %>% readRDS()
snpInPeaksDedup_df = snpInPeaks_df %>%
  # inner_join(peaksMappable, by = 'peakNames') %>% 
  rowwise(name,celltype) %>%
  mutate(tmp = sum(c_across(MAFbin_lowfreq_1:MAFbin_frequent_10))) %>%
  ungroup() %>%
  pivot_longer(cols = contains('MAFbin'), names_to = 'MAFbin', values_to = 'MAFbin_value') %>%
  filter(MAFbin_value > 0 | tmp ==0) %>% dplyr::select(-tmp) %>%
  distinct(name, celltype, .keep_all = T) %>%
  arrange(label,CHR, POS_hg38, name) %>%
  group_by(peakNames) %>%
  mutate(name = paste(name, collapse = ','),
         label = paste(unique(label), collapse = ','),
         trait = paste(unique(trait), collapse = ','),
         PIP = max(PIP, na.rm = T),
         rank = rev(seq(n()))) %>%
  top_n(1, rank) %>% ungroup() %>% dplyr::select(-rank)


###########################
## 1) read in pgls results
pgls_fn =  here('data/raw_data/ldsc_caudate_zoonomia/traits_pgls') %>%
  list.files(pattern = 'roundFinal.csv', full.names = T)
names(pgls_fn) = basename(pgls_fn)
pgls_df = pgls_fn %>% lapply(fread) %>% rbindlist(idcol = 'file') %>%
  mutate(celltype = case_when(grepl('oligo',file) ~ 'Oligo', 
                              grepl('D1',file) ~ 'MSN_D1',
                              grepl('D2',file) ~ 'MSN_D2',
                              TRUE ~ ''),
         speciesSet = case_when(grepl('eg',file) ~ 'Euarchontoglires', 
                              grepl('boreo|all',file) ~ 'Boreoeutheria'),
         zooTrait = ss(file, '_', 2)) %>%
  rename('peakNames' = 'Enhancer', 'PGLS_Pvalue' = 'Exp_Pvalue') %>%
  inner_join(snpInPeaksDedup_df, by = c('celltype','peakNames')) %>%
  mutate(MAFbin2 = case_when(grepl('low', MAFbin) ~ 'lowFreq',
                             TRUE ~ 'commonFreq'), 
         MAFbin2 = factor(MAFbin2))
  
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
  relocate(c(annotation, Gene.Symbol, distanceToTSS), .after = 'A2') %>%
  dplyr::select(-c(Coding_UCSC_lowfreq:base)) %>%
  filter(annotation %in% c('Distal', 'Intron'))

## annotate peaks w/ nearby gene
pgls_df = pgls_df %>%
  group_by(zooTrait, speciesSet, celltype)  %>% dplyr::select(-file) %>%
  mutate(FDR = p.adjust(Pvalue, 'fdr'), 
         PGLS_FDR = p.adjust(PGLS_Pvalue, 'fdr')) %>%
  ungroup() %>% arrange(PGLS_FDR) %>% 
  relocate(c( PGLS_Pvalue, PGLS_FDR, Pvalue, FDR), .after = Coeff)

sum(pgls_df$PGLS_FDR < 0.05) # 13
sum(pgls_df$PGLS_FDR < 0.10) # 21

sum(pgls_df$FDR < 0.05) # 28
sum(pgls_df$FDR < 0.10) # 36


## write pgls tables to xlsx and rds
pgls_xlsx = here(DATADIR, 'tables', 'trait_pgls_Corces2020_finemapped_snps_20210914.xlsx')
writexl::write_xlsx(pgls_df, pgls_xlsx)
pgls_rds = here(DATADIR, 'rdas', 'trait_pgls_Corces2020_finemapped_snps_20210914.rds')
saveRDS(pgls_df, pgls_rds)

## look at signif enhancers
pgls_df %>% top_n(sum(pgls_df$PGLS_FDR < 0.1), -PGLS_FDR) %>% data.frame() %>%
  dplyr::select(c(label, zooTrait, Gene.Symbol, annotation, PGLS_FDR,
                  top_phyloP, top_phastCons, cCRE_group))


