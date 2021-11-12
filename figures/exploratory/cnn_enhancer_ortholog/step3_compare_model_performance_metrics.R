ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(arrow))
library(here)

###################################
### 1) overall model performances ###
PROJDIR='figures/exploratory/cnn_enhancer_ortholog'
cells = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb',  'Astro', 'Microglia', 'OPC', 'Oligo')
groups = c('nonEnhNeg','largeGC','nonEnhLargeGC','nonCellEnhLargeGC')
pred_fn = list.files(path = here('data/raw_data','cnn_enhancer_ortholog','predictions'), 
                     pattern = '.feather', full.names = T, recursive = T) %>%
  grep(pattern = 'DO0.25.performance|DO0.25.None.performance', value = TRUE)
input = pred_fn %>% lapply(read_feather) %>% rbindlist(fill=TRUE)
pred = input %>% mutate(
  trainingSet = case_when(grepl('hg_', model)~ 'HgOnly',  TRUE ~ 'HgRmMm'),
  group = case_when(grepl('nonCelltypeNonEnhBiasAway10x', prefix) ~ 'nonCellEnhLargeGC',
                    grepl('nonEnhBiasAway10x', prefix) ~ 'nonEnhLargeGC',
                    grepl('biasAway10x', prefix) ~ 'largeGC',
                    grepl('nonEnh', prefix) ~ 'nonEnhNeg') %>% factor(groups), 
  fold = ss(prefix, '_fold', 2) %>% ss('_',1), 
  celltype = ss(prefix, '_fold', 1) %>% factor(cells), 
  Pos_to_Neg_ratio = (tp + fn) / (tn + fp),
  NPos = tp + fn, NNeg = tn + fp,
  tpr = tp /(tp + fn), fpr = fp /(fp + tn), tnr = tn /(tn + fp),
  npv = tn / (tn + fn),
  auPRC.adj = auPRC - min(NPos/(NPos + NNeg),NNeg/(NPos + NNeg))) %>%
  filter(!is.na(celltype)) %>% select(-predict_fasta)

table(pred$celltype, pred$group)

####################################
### cross-species group averages ###
df_long = pred %>% pivot_longer(cols = c('auROC','auPRC.adj','f1_score','fpr',
                                         'tpr', 'npv', 'NPos', 'NNeg', 'Pos_to_Neg_ratio'), 
                                values_to ='value',  names_to = 'variable')
table_out = here(PROJDIR, 'tables', 'Data_S11_Cell-TACIT_model_performance_values_A.xlsx')
performance_out = df_long %>% group_by(celltype, variable, group, trainingSet) %>% 
  summarize(mean = mean(value)) %>%
  pivot_wider(names_from = 'variable', values_from = 'mean') %>%
  writexl::write_xlsx(table_out)

## negative to positive ratios across various negative sets
pred %>% group_by(group, trainingSet) %>% 
  mutate(nneg = tn + fp, npos = tp + fn, 
         P2N_ratio = (tn + fp)/ (tp + fn)) %>%
  summarize(
    mean_nneg = mean(nneg), mean_npos = mean(npos),
    mean_N2P_ratio = mean(P2N_ratio))

## negative to positive ratios across various negative sets
df_long %>% filter(group == 'nonCellEnhLargeGC', trainingSet == 'HgRmMm') %>%
  group_by(variable, celltype) %>% 
  summarize(mean1 = mean(value)) %>%
  group_by(variable) %>%
  summarize(mean = mean(mean1), se = sd(mean1))


## negative to positive ratios across various negative sets
df_long %>% filter(group == 'nonCellEnhLargeGC', trainingSet == 'HgOnly') %>%
  group_by(variable, celltype) %>% 
  summarize(mean1 = mean(value)) %>%
  group_by(variable) %>%
  summarize(mean = mean(mean1), se = sd(mean1))




###############################################
## 2) cell type specific model performances ###
pred2_fn = list.files(path = here('data/raw_data','cnn_enhancer_ortholog','predictions'), 
                     pattern = 'CelltypeOnly_NonCelltype_valid.performance.feather', 
                     full.names = T, recursive = T)
names(pred2_fn) = pred2_fn
input2 = pred2_fn %>% lapply(read_feather) %>% data.table::rbindlist(idcol = 'file')
pred2 = input2 %>% mutate(
  group = case_when(grepl('nonCelltypeNonEnhBiasAway10x', prefix) ~ 'nonCellEnhLargeGC',
                    grepl('nonEnhBiasAway10x', prefix) ~ 'nonEnhLargeGC',
                    grepl('biasAway10x', prefix) ~ 'largeGC',
                    grepl('nonEnh', prefix) ~ 'nonEnhNeg') %>% factor(groups), 
  genome = basename(file) %>% ss('\\.', 7) %>% ss('CelltypeOnly', 1),
  trainingSet = case_when(grepl('hgRmMm_',file) ~ 'HgRmMm', TRUE ~ 'HgOnly'),
  fold = ss(prefix, '_fold', 2) %>% ss('_',1), 
  celltype = ss(prefix, '_fold', 1) %>% factor(cells), 
  Pos_to_Neg_ratio = (tp + fn) / (tn + fp),
  NPos = tp + fn, NNeg = tn + fp,
  tpr = tp /(tp + fn), fpr = fp /(fp + tn), tnr = tn /(tn + fp),
  npv = tn / (tn + fn),
  auPRC.adj = auPRC - min(NPos/(NPos + NNeg),NNeg/(NPos + NNeg))) %>%
  filter(!is.na(celltype)) %>% select(-predict_fasta)


## negative to positive ratios across various negative sets
pred2 %>% group_by(group, trainingSet) %>% 
  mutate(nneg = tn + fp, npos = tp + fn, 
         P2N_ratio = (tn + fp)/ (tp + fn)) %>%
  summarize(
    mean_nneg = mean(nneg),
    mean_npos = mean(npos),
    mean_N2P_ratio = mean(P2N_ratio))

## negative to positive ratios across various negative sets
df_long2 = pred2 %>% pivot_longer(cols = c('auPRC.adj','auROC','npv'), 
                                  values_to ='value',  names_to = 'variable')
## write performance metrics out to table
table_out2 = here(PROJDIR, 'tables', 'Data_S11_Cell-TACIT_model_performance_values_B.xlsx')
pred2 %>% pivot_longer(cols = c('auROC','auPRC.adj','f1_score','fpr',
                        'tpr', 'npv', 'NPos', 'NNeg', 'Pos_to_Neg_ratio'), 
               values_to ='value',  names_to = 'variable') %>% 
  group_by(celltype, variable, group, trainingSet) %>% 
  summarize(mean = mean(value)) %>%
  pivot_wider(names_from = 'variable', values_from = 'mean') %>%
  writexl::write_xlsx(table_out2)


df_long2 %>% filter(group == 'nonCellEnhLargeGC', trainingSet == 'HgRmMm') %>%
  group_by(variable, celltype) %>% 
  summarize(mean1 = mean(value)) %>%
  group_by(variable) %>%
  summarize(mean = mean(mean1), se = sd(mean1))

df_long2 %>% filter(group == 'nonCellEnhLargeGC', trainingSet == 'HgOnly') %>%
  group_by(variable, celltype) %>% 
  summarize(mean1 = mean(value)) %>%
  group_by(variable) %>%
  summarize(mean = mean(mean1), se = sd(mean1))



##############################################
### 3) species specific model performances ###
pred3_fn = list.files(path = here('data/raw_data','cnn_enhancer_ortholog','predictions'), 
                     pattern = 'Only_nonEnh_valid.performance.feather', full.names = T, recursive = T)
names(pred3_fn) = pred3_fn
input = pred3_fn %>% lapply(read_feather) %>% data.table::rbindlist(idcol = 'file')
pred3 = input %>% mutate(
  group = case_when(grepl('nonCelltypeNonEnhBiasAway10x', prefix) ~ 'nonCellEnhLargeGC',
                         grepl('nonEnhBiasAway10x', prefix) ~ 'nonEnhLargeGC',
                         grepl('biasAway10x', prefix) ~ 'largeGC',
                         grepl('nonEnh', prefix) ~ 'nonEnhNeg') %>% factor(groups), 
  genome = basename(file) %>% ss('\\.', 7) %>% ss('Only', 1),
  genome = factor(genome, c('hg38', 'rheMac10', 'mm10')),
  trainingSet = case_when(grepl('hgRmMm_',file) ~ 'HgRmMm', TRUE ~ 'HgOnly'),
  fold = ss(prefix, '_fold', 2) %>% ss('_',1), 
  celltype = ss(prefix, '_fold', 1) %>% factor(cells), 
  Pos_to_Neg_ratio = (tp + fn) / (tn + fp),
  NPos = tp + fn, NNeg = tn + fp,
  tpr = tp /(tp + fn), fpr = fp /(fp + tn), tnr = tn /(tn + fp),
  npv = tn / (tn + fn),
  auPRC.adj = auPRC - min(NPos/(NPos + NNeg),NNeg/(NPos + NNeg))) %>%
  filter(!is.na(celltype)) %>% select(-predict_fasta)

## negative to positive ratios across various negative sets
pred3 %>% group_by(group, trainingSet) %>% 
  mutate(nneg = tn + fp, npos = tp + fn, 
         P2N_ratio = (tn + fp)/ (tp + fn)) %>%
  summarize(
    mean_nneg = mean(nneg),
    mean_npos = mean(npos),
    mean_N2P_ratio = mean(P2N_ratio))

## negative to positive ratios across various negative sets
df_long3 = pred3 %>% pivot_longer(cols = c('auPRC','auROC'),
                                values_to ='value',  names_to = 'variable')
## write performance metrics out to table
table_out3 = here(PROJDIR, 'tables', 'Data_S11_Cell-TACIT_model_performance_values_C.xlsx')
pred3 %>% pivot_longer(cols = c('auROC','auPRC.adj','f1_score','fpr',
                                'tpr', 'npv', 'NPos', 'NNeg', 'Pos_to_Neg_ratio'), 
                       values_to ='value',  names_to = 'variable') %>% 
  group_by(celltype, variable, group, trainingSet) %>% 
  summarize(mean = mean(value)) %>%
  pivot_wider(names_from = 'variable', values_from = 'mean') %>%
  writexl::write_xlsx(table_out3)


df_long3 %>% filter(group == 'nonCellEnhLargeGC', trainingSet == 'HgRmMm') %>%
  group_by(variable, celltype, genome) %>% 
  summarize(mean1 = mean(value)) %>%
  group_by(variable, genome) %>%
  summarize(mean = mean(mean1), se = sd(mean1))


## negative to positive ratios across various negative sets
df_long3 %>% filter(group == 'nonCellEnhLargeGC', trainingSet == 'HgOnly') %>%
  group_by(variable, celltype, genome) %>% 
  summarize(mean1 = mean(value)) %>%
  group_by(variable, genome) %>%
  summarize(mean = mean(mean1), se = sd(mean1))



############################################################
### write cell type=specific model performances to table ###
out_fn2 = here(PROJDIR, 'tables', 'Cell-TACIT_celltype-specific_performance_by_species.xlsx')
celltype_specific = 
  df_long2 %>% filter(group == 'nonCellEnhLargeGC', trainingSet == 'HgRmMm') %>%
  group_by(trainingSet, variable, celltype, genome) %>% 
  summarize(mean = mean(value), se = sd(value)) %>%
  writexl::write_xlsx(out_fn2)
  
  
out_fn3 = here(PROJDIR, 'tables', 'Cell-TACIT_species-specific_performance_by_species.xlsx')
species_specific = 
 df_long3 %>% filter(group == 'nonCellEnhLargeGC', trainingSet == 'HgRmMm') %>%
  group_by(trainingSet, variable, celltype, genome) %>% 
  summarize(mean = mean(value), se = sd(value)) %>%
  writexl::write_xlsx(out_fn3)

  