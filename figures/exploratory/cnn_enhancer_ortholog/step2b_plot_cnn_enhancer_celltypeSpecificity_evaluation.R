ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
setwd('figures/exploratory/cnn_enhancer_ortholog')
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(arrow))
library(rcartocolor)
library(here)

###################################
### read in the CNN predictions ###
PROJDIR='cnn_enhancer_ortholog'
cells = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb',  'Astro', 'Microglia', 'OPC', 'Oligo')
groups = c('nonEnhNeg','largeGC','nonEnhLargeGC','nonCellEnhLargeGC')
pred_fn = list.files(path = here('data/raw_data',PROJDIR,'predictions'), 
                     pattern = 'CelltypeOnly_NonCelltype_valid.performance.feather', 
                     full.names = T, recursive = T)
names(pred_fn) = pred_fn
input = pred_fn %>% lapply(read_feather) %>% data.table::rbindlist(idcol = 'file')
pred = input %>% mutate(
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
  filter(!is.na(celltype)) %>% select(-predict_fasta) %>%
  pivot_longer(cols = all_of(c('tpr', 'tnr')),
               names_to = 'eval_acc', values_to = 'value') %>%
  mutate(eval_acc = factor(eval_acc, c('tpr', 'tnr')))

table(pred$celltype, pred$group)
table(pred$celltype, pred$genome)
table(pred$celltype, pred$genome, pred$trainingSet)

###################################
### plot validation performance ###
library(rcartocolor)
height_ppt = 5; width_ppt = 8;
height_fig = 1.75; width_fig = 4.75; font_fig = 7

pdf('plots/cross_celltype_specific_performance_20210510.pdf', 
    width = width_ppt, height = height_ppt)
ggplot(data = pred %>% filter(trainingSet == 'HgOnly'), aes(x = eval_acc, y = value)) + 
  geom_boxplot(aes(fill = group)) + 
  scale_fill_carto_d(name = "Cell type:", palette = "Vivid") +
  facet_grid(genome~celltype, scales = 'free',space = 'free_x') + 
  theme_bw(base_size = 11) + ylim(c(0, 1)) + 
  guides(fill = guide_legend( nrow = 1)) + 
  ggtitle('HgOnly Cell Type Specific Validation Performance') + 
  xlab('Cell Type Specific Performance') + ylab('Accuracy') + 
  theme(legend.position = "bottom", legend.text=element_text(size=10),
        legend.title=element_text(size=10), legend.key.height=unit(.5,"line"), 
        legend.key.width=unit(.5,"line"), legend.margin=margin(-15, 0, 0, 0))

ggplot(data = pred %>% filter(trainingSet == 'HgRmMm'), aes(x = eval_acc, y = value)) + 
  geom_boxplot(aes(fill = group)) + 
  scale_fill_carto_d(name = "Cell type:", palette = "Vivid") +
  facet_grid(genome~celltype, scales = 'free',space = 'free_x') + 
  theme_bw(base_size = 11) + ylim(c(0, 1)) + 
  guides(fill = guide_legend( nrow = 1)) + 
  ggtitle('HgRmMm Cell Type Specific Validation Performance') + 
  xlab('') + ylab('Accuracy') + 
  theme(legend.position = "bottom", legend.text=element_text(size=10),
        legend.title=element_text(size=10), legend.key.height=unit(.5,"line"), 
        legend.key.width=unit(.5,"line"), legend.margin=margin(-15, 0, 0, 0))
dev.off()







pred2 = input %>% mutate(
  group = case_when(grepl('nonCelltypeNonEnhBiasAway10x', prefix) ~ 'nonCellEnhLargeGC',
                    grepl('nonEnhBiasAway10x', prefix) ~ 'nonEnhLargeGC',
                    grepl('biasAway10x', prefix) ~ 'largeGC',
                    grepl('nonEnh', prefix) ~ 'nonEnhNeg') %>% factor(groups), 
  genome = basename(file) %>% ss('\\.', 7) %>% ss('CelltypeOnly', 1),
  genome = factor(genome, c('hg38', 'rheMac10', 'mm10')),
  trainingSet = case_when(grepl('hgRmMm_',file) ~ 'HgRmMm', TRUE ~ 'HgOnly'),
  fold = ss(prefix, '_fold', 2) %>% ss('_',1), 
  celltype = ss(prefix, '_fold', 1) %>% factor(cells), 
  tpr = tp /(tp + fn), fpr = fp /(fp + tn), tnr = tn /(tn + fp),
  auPRC.adj = auPRC - (tp + fn) / (tn + fp)) %>%
  filter(!is.na(celltype)) %>% select(-predict_fasta) %>%
  pivot_longer(cols = all_of(c('auROC', 'auPRC.adj')),names_to = 'eval_acc', values_to = 'value')%>%
  mutate(eval_acc = factor(eval_acc, c('auROC', 'auPRC.adj'))) %>% 
  filter(group == 'nonCellEnhLargeGC') 

pdf('plots/cross_celltype_specific_performance_auPRC_auROC_20211109.pdf', 
    width = width_ppt, height = height_ppt)
ggplot(data = pred2, aes(x = eval_acc, y = value)) + 
  geom_boxplot(aes(fill = genome, color = genome)) + 
  scale_fill_carto_d(name = "Cell type:", palette = "Bold") +
  scale_colour_carto_d(name = "Cell type:", palette = "Bold") +
  facet_grid(trainingSet~celltype, scales = 'free',space = 'free_x') + 
  theme_bw(base_size = 11) + ylim(c(0, 1)) + 
  ggtitle('Cell Type-Specific Validation Performance') + 
  xlab('') + ylab('Accuracy') + 
  guides(fill = guide_legend( nrow = 1)) + 
  theme(legend.position = "bottom", legend.text=element_text(size=10),
        legend.title=element_text(size=10), legend.key.height=unit(.5,"line"), 
        legend.key.width=unit(.5,"line"), legend.margin=margin(-15, 0, 0, 0),
        axis.text.x=element_text(angle = -30, hjust = 0))
dev.off()






