#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F, bitmapType='cairo')
library(tidyverse)
library(tidymodels)
library(corrr)
library(broom)
library(data.table)
library(RColorBrewer)
library(rtracklayer)
library(ggh4x)
library(here)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# to be run in the root github directory
LABEL='ldsc_atac_and_conservation'
PLOTDIR=file.path('figures/exploratory',LABEL)
sapply(c('plots', 'tables', 'rdas'), function(x)
  dir.create(here(PLOTDIR,x), showWarnings = F, recursive = F))
DATADIR=here('data/raw_data',LABEL)

##############################################################################
## 1) import the mean metric values across metric & Cell TACIT Age
bed_fn1 = list.files(here('data/raw_data/ldsc_zoonomia_meta/peaks'), 
                     full.names = T, pattern = '.predActive.bed.gz') %>% 
  str_subset('Corces') %>% str_subset('all')
peak_df = bed_fn1 %>% set_names() %>% lapply(import) %>% 
  lapply(as.data.frame) %>% rbindlist(idcol = 'file') %>% 
  mutate(file = basename(file), celltype = ss(file, '\\.', 2),
         meta = ss(file, '\\.', 4)) 

peak_wide_df = peak_df %>% dplyr::select(-file) %>% 
  pivot_wider(names_from = 'meta', values_from = 'score')

## annotate peaks
annot_df <- annotatePeak(GRanges(peak_wide_df), annoDb='org.Hs.eg.db', 
                         TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, 
                         tssRegion = c(-5000, 5000)) %>% 
  as.GRanges() %>% as.data.frame(row.names = seq(length(.))) %>%
  mutate(annotation = ss(annotation, ' '))

peak_wide_df$annotation = annot_df$annotation

change_type = c("5'" = 'Other', "3'"  = 'Other', 'Downstream' = 'Other', 
                'Exon' = 'Other', 'UTR' = 'Other', 
                'Promoter' = 'Promoter', 
                'Intron' = 'Enhancer', 'Distal' = 'Enhancer')

peak_wide_df$annotation2 = change_type[peak_wide_df$annotation]
with(peak_wide_df, table(celltype, annotation2))
remap = with(peak_wide_df, setNames(annotation2, name))

# split up the individual peaks & export the split up peaks
peak_df$split = with(peak_df, paste0(celltype,'.', meta, '.', remap[name]))
peak_list = split(peak_df, f = peak_df$split) %>% lapply(GRanges)
lengths(peak_list)
peak_fn = here('data/raw_data/ldsc_atac_and_conservation/peak', 
               paste0('Corces2020.', names(peak_list), '.predActive.bed.gz'))
mapply(export, object =  peak_list, con = peak_fn)

## split up the main peaks
peak_wide_df$split = with(peak_wide_df, paste0(celltype,'.main.', remap[name]))
peak_wide_list = split(peak_wide_df, f = peak_wide_df$split) %>% lapply(GRanges)
lengths(peak_wide_list)
peak_fn2 = here('data/raw_data/ldsc_atac_and_conservation/peak', 
               paste0('Corces2020.', names(peak_wide_list), '.bed.gz'))
mapply(export, object =  peak_wide_list, con = peak_fn2)

