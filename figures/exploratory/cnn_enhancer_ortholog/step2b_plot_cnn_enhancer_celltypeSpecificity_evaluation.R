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
model_types = c('nonEnhNeg','largeGC','nonEnhLargeGC','nonCellEnhLargeGC')
pred_fn = list.files(path = here('data/raw_data',PROJDIR,'predictions'), 
                     pattern = 'CelltypeOnly_NonCelltype_valid.performance.feather', full.names = T, recursive = T)
names(pred_fn) = pred_fn
input = pred_fn %>% lapply(read_feather) %>% data.table::rbindlist(idcol = 'file')
pred = input %>% mutate(
  model_type = case_when(grepl('nonCelltypeNonEnhBiasAway10x', prefix) ~ 'nonCellEnhLargeGC',
                    grepl('nonEnhBiasAway10x', prefix) ~ 'nonEnhLargeGC',
                    grepl('biasAway10x', prefix) ~ 'largeGC',
                    grepl('nonEnh', prefix) ~ 'nonEnhNeg') %>% factor(model_types), 
  genome = basename(file) %>% ss('\\.', 7) %>% ss('CelltypeOnly', 1),
  fold = ss(prefix, '_fold', 2) %>% ss('_',1), 
  celltype = ss(prefix, '_fold', 1) %>% factor(cells), 
  tpr = tp /(tp + fn), fpr = fp /(fp + tn), tnr = tn /(tn + fp),
  auPRC.adj = auPRC - (tp + fn) / (tn + fp)) %>%
  filter(!is.na(celltype)) %>% select(-predict_fasta) %>%
  pivot_longer(cols = all_of(c('tpr', 'tnr')),names_to = 'eval_acc', values_to = 'value') %>%
  mutate(eval_acc = factor(eval_acc, c('tpr', 'tnr')))

table(pred$celltype, pred$model_type)
table(pred$celltype, pred$genome)

###################################
### plot validation performance ###
library(rcartocolor)
height_ppt = 5; width_ppt = 8;
height_fig = 1.75; width_fig = 4.75; font_fig = 7

pdf('plots/cross_celltype_specific_performance_20210427.pdf', 
    width = width_ppt, height = height_ppt)
ggplot(data = pred, aes(x = eval_acc, y = value)) + 
  geom_boxplot(aes(fill = model_type)) + 
  # geom_jitter(width = .25, pch = 21, aes(fill = model_type)) + 
  scale_fill_carto_d(name = "Cell type:", palette = "Vivid") +
  facet_grid(genome~celltype, scales = 'free',space = 'free_x') + 
  theme_bw(base_size = 11) + ylim(c(0, 1)) + 
  guides(fill = guide_legend( nrow = 1)) + 
  xlab('Cell Type Specific Performance') + ylab('Accuracy') + 
  theme(legend.position = "bottom", legend.text=element_text(size=10),
        legend.title=element_text(size=10),
        legend.key.height=unit(.5,"line"), 
        legend.key.width=unit(.5,"line"), 
        ) 
dev.off()





