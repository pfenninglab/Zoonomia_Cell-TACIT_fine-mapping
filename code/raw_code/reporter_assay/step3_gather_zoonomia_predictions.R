#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F, bitmapType='cairo')
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(here))
suppressMessages(library(rmeta))
library(Matrix)

# to be run in the root github directory
DATADIR='data/raw_data/reporter_assay'
PLOTDIR='figures/exploratory/reporter_assay'

##############################################
# read in the Zoonomia tree and species list
rda_fn = here('data/tidy_data/Zoonomia_data', 
              'rdas','200_Mammals_Genome_Information.rda')
load(file = rda_fn)
col_clade = df %>% select(Clade, col_clade)%>% filter(!duplicated(Clade)) %>% deframe()
col_order = df %>% select(Order, col_order) %>% filter(!duplicated(Order)) %>% deframe()
celltypes = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb', 'Astro', 'Microglia', 'OPC', 'Oligo')


######################################################
# read in the LDSC partitioned heritability estimation
enrich_fn =here(DATADIR,'predictions') %>% 
  list.files(path = ., pattern = '.avgCNN.predictions.txt.gz', full.names = T)
names(enrich_fn) = enrich_fn %>% basename() %>% ss('.avgCNN.predictions.txt.gz')

cell = 'Oligo'
dir.create(here(DATADIR,'tables'), showWarnings = F)
dir.create(here(DATADIR,'rdas'), showWarnings = F)

for( cell in celltypes){
  outPeaks_tsv = here(DATADIR,'tables',paste0('Zoonomia_reporter_assay.', cell, '.allPeaks.avgCNN.predictions.txt.gz'))
  outPeaks_rds = here(DATADIR,'rdas',paste0('Zoonomia_reporter_assay.',cell, '.allPeaks.avgCNN.predictions.rds'))
  enrich_fn2 = enrich_fn %>% grep(pattern = cell, value = T)
  input = enrich_fn2 %>% lapply(fread, header = F, col.names = c('name', 'score'))  %>%
    data.table::rbindlist(idcol = 'file') %>% 
    mutate(species = file %>% ss('\\.',3)) 
  
  predScores_peaks = input %>% select(-file) %>% 
    pivot_wider(values_from = score, names_from = species, values_fn = mean) %>% 
    filter(grepl('hg38', name))
  
  # write out the predictions over all mapped peaks
  write_tsv(predScores_peaks, file = outPeaks_tsv)
  saveRDS(predScores_peaks, outPeaks_rds)
}



######################################################
# calibrate the model scores to %validation positives
MODEL_TYPE = 'hgRmMm_nonCelltypeNonEnhBiasAway10x'
calib_out_fn = here('data/raw_data/cnn_enhancer_ortholog',
                    paste('rdas/Caudate',MODEL_TYPE,'pos_calibration_ecdf.rds', sep = '_'))
model_calibration = readRDS(file = calib_out_fn)

cell = 'MSN_D1'
for( cell in celltypes){
  outEnh_rds = here(DATADIR,'rdas',paste0('Zoonomia_reporter_assay.',cell, '.allPeaks.avgCNN.predictions.rds'))
  outCalib_tsv = here(DATADIR,'tables',paste0('Zoonomia_reporter_assay.',cell, '.calibPeaks.avgCNN.predictions.txt.gz'))
  outCalib_rds = here(DATADIR,'rdas',paste0('Zoonomia_reporter_assay.',cell, '.calibPeaks.avgCNN.predictions.rds'))
  # if(any(!file.exists(c(outCalib_tsv, outCalib_rds)))){
  df2 = readRDS(outEnh_rds) %>% 
    mutate_if(is.numeric, ~ model_calibration[[cell]](.)) %>% 
    mutate_if(is.numeric, ~ ifelse(. == Inf, 1, .))
  
  # write out the predictions over all mapped peaks
  write_tsv(df2, file = outCalib_tsv)
  saveRDS(df2, outCalib_rds)
  # }
}





########################################################################
# CellTACIT score: compute calibrated CNN score then compute cross-species 
dir.create(here(DATADIR,'CellTACIT'), showWarnings = F)
cell = 'MSN_D1'
for( cell in celltypes){
  outGroupMeta_rds = here(DATADIR,'rdas',paste0('Zoonomia_reporter_assay.',cell, '.ZooMeta.rds'))
  outCellTACIT_rds = here(DATADIR,'rdas',paste0('Zoonomia_reporter_assay.',cell, '.CellTACIT.mean.rds'))
  outCellTACIT_bed = here(DATADIR,'CellTACIT',paste0('Zoonomia_reporter_assay.',cell, '.CellTACIT.mean.bed.gz'))

  ## read in the peak x species mtx -> peak x group long
  fn = here(DATADIR,'rdas',paste('Zoonomia_reporter_assay',cell, 'allPeaks.avgCNN.predictions.rds', sep = '.'))
  df_allMeta = readRDS(fn) %>% 
    filter(grepl('hg38', name)) %>% 
    # for peaks w/ predicted open, get the calibrated positive score
    mutate_if(is.numeric, ~ model_calibration[[cell]](.)) %>% 
    pivot_longer(cols = !name, names_to = 'Species', values_to = 'score') %>%
    right_join(df %>% select(c(Species, group_meta)), by = 'Species') %>% 
    group_by(group_meta, name) %>% summarise(score = mean(score, na.rm = T)) %>%
    mutate(score = ifelse(is.na(score),0, score)) %>% 
    ungroup()
  df_allMeta %>% saveRDS(outGroupMeta_rds)
  
  df_allMeta2 = df_allMeta %>%
    mutate(MYA = group_meta %>% as.character() %>% ss('#', 2), 
           MYA = as.numeric(MYA) + 1) %>%
    filter(!is.na(name)) %>% group_by(name) %>%    
    summarise(score = sum(score * MYA) / n()) %>%
    mutate( score = ifelse(is.na(score), 1, score))
  
  gr_allMeta2 = df_allMeta2 %>%
    mutate(seqnames = ss(name, ":", 2), 
           start = ss(name, ":", 3) %>% ss('-', 1), 
           end = ss(name, ":", 3) %>% ss('-', 2)) %>%
    column_to_rownames(var ='name') %>% GRanges()
  
  gr_allMeta2 %>% saveRDS(outCellTACIT_rds)
  export(gr_allMeta2, outCellTACIT_bed)
}

