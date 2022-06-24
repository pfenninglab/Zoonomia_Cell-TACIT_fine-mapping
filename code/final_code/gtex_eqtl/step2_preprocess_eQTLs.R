#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
library(rtracklayer)
library(arrow)
library(GenomicRanges)
library(here)

CODEDIR='code/final_code/gtex_eqtl'
DATADIR='data/tidy_data/gtex_eqtl'
PLOTDIR='figures/exploratory/gtex_eqtl'

dir.create(here(DATADIR, 'rdas'), showWarnings = F)

alpha = 0.05
caudate_eQTL = here(DATADIR, 'GTEx_v8_finemapping_CAVIAR', 
                    'CAVIAR_Results_v8_GTEx_LD_ALL_NOCUTOFF_with_Allele.txt.gz') %>% fread()
caudate_eQTL_gr = caudate_eQTL %>% filter(!duplicated(eQTL)) %>% 
  mutate(tmp = paste0('chr', CHROM, ':', POS)) %>% 
  dplyr::select(eQTL, tmp) %>% deframe() %>% GRanges()
  
caudate_eQTL_gr %>% saveRDS(here(DATADIR, 'rdas', 'GTEx_v8_finemapping_hg38.gr.rds'))
caudate_eQTL %>% saveRDS(here(DATADIR, 'rdas', 'GTEx_v8_finemapping.df.rds'))

##
caudate_eQTL = here(DATADIR, 'GTEx_Analysis_v8_eQTL', 
                    'Brain_Caudate_basal_ganglia.v8.allpairs.parquet') %>% 
  read_parquet(col_select = c('gene_id','variant_id', 'pval_nominal')) %>% 
  mutate(FDR = p.adjust(pval_nominal, 'fdr')) %>% 
  arrange(pval_nominal) %>% filter(!duplicated(variant_id))

table(caudate_eQTL$FDR < alpha)

##
caudate_eQTL_gr = caudate_eQTL %>% 
  mutate(chr = ss(variant_id, '_', 1), pos_19 = ss(variant_id, '_', 2), 
       tmp = paste0(chr, ':', pos_19)) %>% 
  dplyr::select(variant_id, tmp) %>% deframe() %>% GRanges()
seqlevelsStyle(caudate_eQTL_gr) <- "UCSC"

caudate_eQTL_gr %>% saveRDS(here(DATADIR, 'rdas', 'Brain_Caudate_basal_ganglia_hg38.v8.gr.rds'))
caudate_eQTL %>% saveRDS(here(DATADIR, 'rdas', 'Brain_Caudate_basal_ganglia.v8.df.rds'))


