ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
setwd('figures/exploratory/cnn_enhancer_ortholog')
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(arrow))

###################################
### read in the CNN predictions ###
cells = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb',  'Astro', 'Microglia', 'OPC', 'Oligo')
cells = c('MSN_D1', 'MSN_D2', "MSN_SN",'INT_Pvalb')
PROJDIR='../../../data/raw_data/cnn_enhancer_ortholog/predictions'
pred_fn = list.files(path = PROJDIR, pattern = '.feather', full.names = T, recursive = T)
pred = pred_fn %>% lapply(read_feather) %>% bind_rows()
pred = pred %>% mutate(
  group = case_when(grepl('nonEnhBiasAway10x', prefix) ~ 'nonEnhLargeGC',
                    grepl('biasAway10x', prefix) ~ 'largeGC',
                    grepl('nonEnh', prefix) ~ 'nonEnhNeg'),
  fold = ss(prefix, '_fold', 2) %>% ss('_',1), 
  celltype = ss(prefix, '_fold', 1) %>% factor(cells)) %>%
  filter(!is.na(celltype))

table(pred$celltype, pred$group)

pred2 = pred %>% group_by(celltype, fold, group) %>% 
  top_n(1, auPRC) %>% ungroup()

