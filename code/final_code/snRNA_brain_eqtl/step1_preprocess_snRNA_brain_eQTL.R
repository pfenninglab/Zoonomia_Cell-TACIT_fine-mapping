#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
library(GenomicRanges)
library(here)

DATADIR='data/tidy_data/snRNA_brain_eqtl'
dir.create(here(DATADIR, 'rdas'))
############################
## 1) read in the SNPs tested 
snp_df = here(DATADIR, 'tables', 'snp_pos.txt') %>% fread() %>%
  dplyr::rename('rsID' = 'SNP') %>% filter(!is.na(chr))
snp_gr = snp_df %>% mutate(tmp = paste0(chr, ':', pos_hg38)) %>% 
  dplyr::select(rsID, tmp) %>% deframe() %>% GRanges()

snp_gr %>% saveRDS(here(DATADIR, 'rdas', 'snRNA_brain_eqtl_hg38.gr.rds'))
snp_df %>% saveRDS(here(DATADIR, 'rdas', 'snRNA_brain_eqtl.df.rds'))


###############################################################################
## 2) read cell type eQTLs for the subset of cell types matching CellTACIT models
celltypes = c('Astrocytes' = 'Astro', 'Microglia' = 'Microglia', 'OPCs' = 'OPC',
              'Inhibitory' = 'INT_Pvalb', 'Oligodendrocytes' = 'Oligo', 
              'Excitatory' = 'EXC_Neur')

cols = c('Gene_id', 'rsID', 'distToTSS', 'p.value', 'beta')

eqtl_fn = here(DATADIR, 'tables') %>% 
  list.files(pattern = '.gz', full.names = T)
names(eqtl_fn) = basename(eqtl_fn) %>% ss('\\.', 1)
eqtl_fn = eqtl_fn[names(eqtl_fn) %in% names(celltypes)]
cell = 'OPCs'

for(cell in names(celltypes)){
  save_fn = here(DATADIR, 'rdas', paste0('snRNA_brain_eqtl.',cell,'.df.rds'))
  if(!file.exists(save_fn)){
  celltype_df = eqtl_fn %>%
    str_subset(cell) %>%
    lapply(fread, header = F, col.names = cols) %>% 
    rbindlist(idcol = 'celltype') %>% 
    mutate(gene_name = ss(Gene_id, '_', 1), gene_id = ss(Gene_id, '_', 2)) %>%
    dplyr::select(-Gene_id)
  
  celltype_df %>% saveRDS(save_fn)
  }
}


