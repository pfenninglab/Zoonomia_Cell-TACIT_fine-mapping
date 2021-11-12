library(tidyverse)
library(here)
library(ComplexHeatmap)
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

PROJDIR='figures/exploratory/cnn_enhancer_ortholog'
DATADIR='data/raw_data/cnn_enhancer_ortholog'
cells = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb',  'Astro', 'Microglia', 'OPC', 'Oligo')

## read in meme-chip results tables
summary_fn = list.files(here(DATADIR,'meme-chip'), recursive = T, full.names = T, pattern = 'summary.tsv')
names(summary_fn) = ss(summary_fn, 'meme-chip|summary', 2) %>% gsub(pattern = '/', replacement = '')
meme_df = summary_fn %>% lapply(data.table::fread) %>% data.table::rbindlist(idcol = 'file') %>%
  rename_with(make.names) %>%
  ## convert E-values to number and log10
  mutate(genome = ss(file, '\\.', 1), celltype = ss(file, '\\.', 2) %>% factor(cells),
         E.VALUE = parse_number (E.VALUE), logEval = -log10(E.VALUE), 
         logEval = ifelse(is.infinite(logEval), max(logEval[!is.infinite(logEval)]), logEval)) %>%
  ## clean up meme_df motifs common names
  mutate(ALT_ID = ifelse(grepl('^DREME-|^MEME-', ALT_ID), MOST_SIMILAR_MOTIF, ALT_ID),
         ALT_ID = ifelse(grepl('\\(|\\)', ALT_ID), ss(ALT_ID, '\\(|\\)',2), ALT_ID),
         ALT_ID = gsub('_MACMU', '', ALT_ID),
         ALT_ID = toupper(ALT_ID))%>% 
  filter(!grepl('DREME-', ALT_ID), ALT_ID!='') 

  meme_df %>% writexl::write_xlsx(here(PROJDIR, 'tables', 'Data_S13_cross-species_meme-chip_motifs_B.xlsx'))

  

## make a matrix for heatmap
meme_mat = meme_df %>% 
  pivot_wider(id_cols = ALT_ID, names_from = c('celltype', 'genome'), 
              values_from = 'logEval', names_sep = ".",values_fill = NA, values_fn = max) %>%
  column_to_rownames('ALT_ID') %>%
  mutate(tmp = apply(., 1, sd, na.rm = T)) %>%
  filter(tmp > 20) %>% dplyr::select(-tmp)
meme_mat = meme_mat %>% dplyr::select(order(colnames(meme_mat)))


## make into lists of TFs
meme_df2 = meme_df %>% arrange(celltype, ALT_ID) %>%
  ## pick out 1 MOTIF per cell type and species
  group_by(celltype, genome) %>% filter(!duplicated(ALT_ID)) %>% 
  ## pick out MOTIFs shared between species
  group_by(celltype, ALT_ID) %>% filter(n() > 1) %>% 
  ## exclude MOTIFs in many cell types
  group_by(ALT_ID) %>% filter(n() < 7) %>% 
  group_by(celltype) %>% summarise(ALT_ID = paste(unique(ALT_ID), collapse = ', ')) %>%
  ungroup() %>% filter(!duplicated(celltype)) %>% 
  writexl::write_xlsx(here(PROJDIR, 'tables', 'Data_S13_cross-species_meme-chip_motifs_A.xlsx'))




