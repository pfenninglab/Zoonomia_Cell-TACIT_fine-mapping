library(tidyverse)
library(data.table)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(reshape2)
library(here)
options(stringsAsFactors = F)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

# to be run in the root github directory
DATADIR=here('data/raw_data/cnn_enhancer_ortholog')
celltypes = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb', 'Astro', 'Microglia', 'OPC', 'Oligo')
MODEL_TYPE = 'hgRmMm_nonCelltypeNonEnhBiasAway10x'
calibPos_fn = list.files(path = here(DATADIR,'predictions'), full.names = T, recursive = T,
                         pattern = '.validPositiveForCalibration.predictions.txt.gz')
calibPos_fn = calibPos_fn %>% 
  grep(pattern = MODEL_TYPE, value = T) %>% grep(pattern = 'DO0.25.valid', value = T)
names(calibPos_fn) = calibPos_fn %>% basename() %>% ss('_fold')

input = calibPos_fn %>% map(fread) %>% rbindlist(idcol = 'celltype') 

predScores = input %>% select(-y_pred_class) %>% group_by(Name, celltype) %>%
  summarise(y_pred_score = mean(y_pred_score)) %>% ungroup() %>%
  mutate(celltype = factor(celltype, celltypes))

model_calibration = lapply(split(predScores$y_pred_score, predScores$celltype), ecdf)
calib_out_fn = here(DATADIR,paste('rdas/Caudate',MODEL_TYPE,'pos_calibration_ecdf.rds', sep = '_'))
saveRDS(model_calibration, file = calib_out_fn)

