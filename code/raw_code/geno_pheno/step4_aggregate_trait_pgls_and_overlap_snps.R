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
library(swfdr)

suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))

CODEDIR='code/raw_code/geno_pheno'
DATADIR='data/raw_data/geno_pheno'
PLOTDIR='figures/exploratory/geno_pheno'
i_am(file.path(CODEDIR, 'step4_aggregate_trait_pgls_and_overlap_snps.R'))


###########################################################
## 1) read in table of fine-mapped SNPs overlapping OCRs ##
allPeaksSNPs_fn = here('data/raw_data/polyfun_caudate/rdas',
                       'polyfun_caudate_finemapped_snpsInAllPeaks_20210518.rds')
snpInPeaks_df = allPeaksSNPs_fn %>% readRDS()
snpInPeaksDedup_df = snpInPeaks_df %>%
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
pgls_fn =  here(DATADIR, 'rdas/permulated') %>% list.files(pattern = 'phylolm.rds', full.names = T)
names(pgls_fn) = basename(pgls_fn) %>% ss('.phylolm.rds')
pgls_df = pgls_fn %>% lapply(readRDS) %>% rbindlist(idcol = 'file', fill=TRUE) %>%
  mutate(speciesSet = case_when(grepl('eg',trait) ~ 'Euarchontoglires', 
                                TRUE ~ 'All'), 
         zooTrait = ss(trait, '_', 1)) %>%
  dplyr::rename('peakNames' = 'peak') %>%
  left_join(snpInPeaksDedup_df, by = c('celltype','peakNames'))

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
  dplyr::select(-c(Coding_UCSC_lowfreq:base), -trait.x, -trait.y, ) %>%
  filter(annotation %in% c('Distal', 'Intron'))

with(pgls_df, table(zooTrait, speciesSet, celltype))

#################################################################################
## 2) re-calculate FDR based per trait, species set, and cell type b/c of overlapping
pgls_df = pgls_df %>%
  group_by(zooTrait, speciesSet, celltype)  %>% dplyr::select(-file) %>%
  mutate( FDR = lm_qvalue(p.value, X=numMappable)$qvalues, 
          perm_FDR = lm_qvalue(perm.pvalue, X= numMappable)$qvalues) %>%
  ungroup() %>% arrange(perm_FDR) %>% 
  relocate(c( perm.pvalue, perm_FDR, p.value, FDR), .after = Estimate)


## p-value w/ permulations
sum(pgls_df$perm_FDR < 0.05) # 893
sum(pgls_df$perm_FDR < 0.10) # 1669
with(pgls_df %>% filter(speciesSet == 'All'), table(zooTrait, perm_FDR < 0.05, celltype))
# with(pgls_df %>% filter(speciesSet == 'Euarchontoglires'), table(zooTrait, perm_FDR < 0.05, celltype))

## look if pgls not w/ permulations
sum(pgls_df$FDR < 0.05) # 1038
sum(pgls_df$FDR < 0.10) # 1863
with(pgls_df %>% filter(speciesSet == 'All'), table(zooTrait, FDR < 0.05, celltype))
# with(pgls_df %>% filter(speciesSet == 'Euarchontoglires'), table(zooTrait, FDR < 0.05, celltype))


## write pgls tables to xlsx and rds
pgls_xlsx = here(DATADIR, 'tables', 'geno_pheno_Corces2020_finemapped_snps_20220512.xlsx')
writexl::write_xlsx(pgls_df, pgls_xlsx)
pgls_rds = here(DATADIR, 'rdas', 'geno_pheno_Corces2020_finemapped_snps_20220512.rds')
saveRDS(pgls_df, pgls_rds)

## look at signif enhancers
pgls_df %>% top_n(sum(pgls_df$perm_FDR < 0.05), -perm_FDR) %>% data.frame() %>%
  dplyr::select(c(label, zooTrait, Gene.Symbol, annotation, perm_FDR,
                  top_phyloP, top_phastCons, cCRE_group))


