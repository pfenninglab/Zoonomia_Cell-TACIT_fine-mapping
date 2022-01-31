#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(here))

# to be run in the root github directory
DATADIR='data/raw_data/select_mpra_candidates'

######################################################
# calibrate the model scores to %validation positives
MODEL_TYPE = 'hgRmMm_nonCelltypeNonEnhBiasAway10x'
DATADIR2='data/raw_data/cnn_enhancer_ortholog'
calib_out_fn = here(DATADIR2,paste('rdas/Caudate',MODEL_TYPE,'pos_calibration_ecdf.rds', sep = '_'))
model_calibration = readRDS(file = calib_out_fn)
celltypes = names(model_calibration)

pred_files_prefix = c('BNP_PPRA_enhancers_20220126', 
                      'IMK_PPRA_SequencesAll', 'MEW_VL_PPRA2Sequences')

cell = 'MSN_D1'
for(pref in pred_files_prefix){
  for( cell in celltypes){
    in_pred_tsv = here(DATADIR,'predictions',paste0(pref,'.',cell, '.predictions.txt.gz'))
    out_pred_tsv = here(DATADIR,'predictions',paste0(pref,'.',cell, '.normPred.txt.gz'))
    # if(any(!file.exists(c(outCalib_tsv, outCalib_rds)))){
    df = fread(in_pred_tsv, header = F) %>% 
      mutate_if(is.numeric, ~ model_calibration[[cell]](.)) %>% 
      mutate_if(is.numeric, ~ ifelse(. == Inf, 1, .))
    # write out the predictions over all mapped peaks
    write_tsv(df, file = out_pred_tsv)
    # }
  }
}


