#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F, bitmapType='cairo')
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(ggtree))
suppressMessages(library(ggtreeExtra))
suppressMessages(library(GenomicRanges))
suppressMessages(library(rcartocolor))
suppressMessages(library(ggrepel))
suppressMessages(library(ggforce))
suppressMessages(library(ggpubr))
suppressMessages(library(ggsci))
suppressMessages(library(here))
suppressMessages(library(rmeta))
library(Matrix)

# to be run in the root github directory
LABEL='Zoonomia_data'
PROJDIR=here('data/raw_data/ldsc_caudate_zoonomia')
SETDIR=here('figures/exploratory/ldsc_caudate_zoonomia')

##############################################
# read in the Zoonomia tree and species list
rda_fn = here('data/tidy_data/Zoonomia_data', 
              'rdas','200_Mammals_Genome_Information.rda')
load(file = rda_fn)
col_clade = df %>% select(Clade, col_clade)%>% filter(!duplicated(Clade)) %>% deframe()
col_order = df %>% select(Order, col_order) %>% filter(!duplicated(Order)) %>% deframe()
celltypes = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb', 'Astro', 'Microglia', 'OPC', 'Oligo')

## read in human_peakList and human_enhList
load(here('data/raw_data/cnn_enhancer_ortholog', 'rdas', 
            'caudate_positive_sequences.hg38..rda'))
human_peakList = human_peakList[celltypes] %>% as.list() %>% lapply(as.data.frame)
human_enhList = human_enhList[celltypes] %>% as.list() %>% lapply(as.data.frame)

######################################################
# read in the LDSC partitioned heritability estimation
enrich_fn =here(PROJDIR,'predictions') %>% 
  list.files(path = ., pattern = '.avgCNN.predictions.txt.gz', full.names = T)
names(enrich_fn) = enrich_fn %>% basename() %>% ss('.avgCNN.predictions.txt.gz')

cell = 'MSN_D1'
dir.create(here(PROJDIR,'tables'), showWarnings = F)
dir.create(here(PROJDIR,'rdas'), showWarnings = F)
for( cell in celltypes){
  outPeaks_tsv = here(PROJDIR,'tables',paste0('Corces2020.', cell, '.allPeaks.avgCNN.predictions.txt.gz'))
  outPeaks_rds = here(PROJDIR,'rdas',paste0('Corces2020.',cell, '.allPeaks.avgCNN.predictions.rds'))
  outEnh_tsv = here(PROJDIR,'tables',paste0('Corces2020.',cell, '.enhPeaks.avgCNN.predictions.txt.gz'))
  outEnh_rds = here(PROJDIR,'rdas',paste0('Corces2020.',cell, '.enhPeaks.avgCNN.predictions.rds'))
  
  if(any(!file.exists(c(outPeaks_tsv, outEnh_tsv, outPeaks_rds, outEnh_rds)))){
  enrich_fn2 = enrich_fn %>% grep(pattern = cell, value = T)
  input = enrich_fn2 %>% lapply(fread, header = F, col.names = c('name', 'score'))  %>%
    data.table::rbindlist(idcol = 'file') %>% 
    mutate(species = file %>% ss('\\.',3)) 
  
  predScores_peaks = input %>% select(-file) %>% 
    pivot_wider(values_from = score, names_from = species, values_fn = mean)

  # write out the predictions over all mapped peaks
  write_tsv(predScores_peaks, file = outPeaks_tsv)
  saveRDS(predScores_peaks, outPeaks_rds)
  
  # write out the predictions over mapped enhancers peaks
  predScores_enh = predScores_peaks %>% filter(name %in% human_enhList[[cell]]$name)
  write_tsv(predScores_enh, file = outEnh_tsv)
  saveRDS(predScores_enh, outEnh_rds)
}
}



######################################################
# calibrate the model scores to %validation positives
MODEL_TYPE = 'hgRmMm_nonCelltypeNonEnhBiasAway10x'
DATADIR=here('data/raw_data/cnn_enhancer_ortholog')
calib_out_fn = here(DATADIR,paste('rdas/Caudate',MODEL_TYPE,'pos_calibration_ecdf.rds', sep = '_'))
model_calibration = readRDS(file = calib_out_fn)

cell = 'MSN_D1'
for( cell in celltypes){
  outEnh_rds = here(PROJDIR,'rdas',paste0('Corces2020.',cell, '.enhPeaks.avgCNN.predictions.rds'))
  outCalib_tsv = here(PROJDIR,'tables',paste0('Corces2020.',cell, '.calibPeaks.avgCNN.predictions.txt.gz'))
  outCalib_rds = here(PROJDIR,'rdas',paste0('Corces2020.',cell, '.calibPeaks.avgCNN.predictions.rds'))
  # if(any(!file.exists(c(outCalib_tsv, outCalib_rds)))){
    df = readRDS(outEnh_rds) %>% 
      mutate_if(is.numeric, ~ model_calibration[[cell]](.)) %>% 
      mutate_if(is.numeric, ~ ifelse(. == Inf, 1, .))
  
    # write out the predictions over all mapped peaks
    write_tsv(df, file = outCalib_tsv)
    saveRDS(df, outCalib_rds)
  # }
  }





