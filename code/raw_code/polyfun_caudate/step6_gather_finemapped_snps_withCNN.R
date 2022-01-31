#######################################
### set up libraries and functions ####
# conda activate r4
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F, bitmapType='cairo')
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
library(rcartocolor)
library(here)

DATADIR = 'data/raw_data/polyfun_caudate'

##############################################
# read in the Zoonomia tree and species list
rda_fn = here('data/tidy_data/Zoonomia_data', 
              'rdas','200_Mammals_Genome_Information.rda')
load(file = rda_fn)
col_clade = df %>% select(Clade, col_clade)%>% filter(!duplicated(Clade)) %>% deframe()
col_order = df %>% select(Order, col_order) %>% filter(!duplicated(Order)) %>% deframe()
celltypes = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb', 'Astro', 'Microglia', 'OPC', 'Oligo')

#############################################
# calibrate the model scores to %validation positives
MODEL_TYPE = 'hgRmMm_nonCelltypeNonEnhBiasAway10x'
DATADIR2=here('data/raw_data/cnn_enhancer_ortholog')
calib_out_fn = here(DATADIR2,paste('rdas/Caudate',MODEL_TYPE,'pos_calibration_ecdf.rds', sep = '_'))
model_calibration = readRDS(file = calib_out_fn)

##############################################
# read in the Zoonomia tree and species list
snp_pred_fn = here(DATADIR,'predictions') %>%
  list.files(pattern = 'SNPs.predictions.txt.gz', recursive = T, full.names = T)
snp_pred_fnList = split(snp_pred_fn, basename(snp_pred_fn) %>% ss('_fold'))

snp_predictions_df = lapply(names(snp_pred_fnList), function(celltype){
  fn = snp_pred_fnList[[celltype]]
  df = lapply(fn, fread) %>% rbindlist()
  df = df %>% group_by(Name) %>% 
    summarize(score = mean(y_pred_logit), score = 1/ (1 + exp(-score)),
              score = model_calibration[[celltype]](score)) %>%
    ungroup() %>% 
    mutate(allele = ifelse(grepl('NonEffect', Name), 'nonEffect', 'Effect'),
           name = gsub(':NonEffect|:Effect','', Name)) %>%
    dplyr::select(-Name) 
  df = df %>% pivot_wider(id_cols = name, values_from = score, names_from = allele) %>%
    mutate(alleleDiff = Effect - nonEffect)
  return(df)
}) %>% rbindlist(idcol = 'celltype')
snp_predictions_df$celltype = names(snp_pred_fnList)[snp_predictions_df$celltype]

summary(snp_predictions_df$alleleDiff)

############################################
## read in the list of fine-mapped SNPs ##
poly_fn = here('data/raw_data/polyfun_caudate/rdas',
               'polyfun_caudate_finemapped_snps_20210518.rds')
snps_df = readRDS(file = poly_fn)

snps_df2 = full_join(snps_df, snp_predictions_df) %>% 
  dplyr::select(-c(population:signif_group, Zoonomia_phastCons.43prim.fdr.05:base)) 

poly_fn2 = here('data/raw_data/polyfun_caudate/rdas',
               'polyfun_caudate_finemapped_snps_with_CNN_predictions_20220119.rds')
saveRDS(snps_df2, file = poly_fn2)

poly_fn3 = here('data/raw_data/polyfun_caudate/tables',
                'polyfun_caudate_finemapped_snps_with_CNN_predictions_20220119.xlsx')
split(snps_df2, snps_df2$celltype) %>% writexl::write_xlsx(poly_fn3)
  










