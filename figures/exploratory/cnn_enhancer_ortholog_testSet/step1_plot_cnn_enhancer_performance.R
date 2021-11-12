ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
setwd('figures/exploratory/cnn_enhancer_ortholog_testSet')
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(arrow))
library(here)

###################################
### read in the CNN predictions ###
PROJDIR='cnn_enhancer_ortholog'
cells = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb',  'Astro', 'Microglia', 'OPC', 'Oligo')
groups = c('nonEnhNeg','largeGC','nonEnhLargeGC','nonCellEnhLargeGC')
pred_fn = list.files(path = here('data/raw_data',PROJDIR,'predictions'), 
                     pattern = '.feather', full.names = T, recursive = T) %>%
  grep(pattern = 'nonCelltypeNonEnhBiasAway10x_test', value = TRUE) 
input = pred_fn %>% lapply(read_feather) %>% rbindlist(fill=TRUE)
pred = input %>% mutate(
  trainingSet = case_when(grepl('hg_', model)~ 'HgOnly',  TRUE ~ 'HgRmMm'),
  trainingSet = factor(trainingSet, c('HgOnly', 'HgRmMm')),
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

df_long = pred %>% pivot_longer(cols = c('auROC','auPRC.adj','f1_score','fpr','tpr'), 
                              values_to ='value',  names_to = 'variable')


###################################
### plot validation performance ###
library(rcartocolor)
height_ppt = 5; width_ppt = 8;
height_fig = 1.75; width_fig = 4.75; font_fig = 7

pdf('plots/cross_species_test_performance_20211109.pdf', 
    width = width_ppt, height = height_ppt/3*2)
ggplot(data = df_long, aes(x = trainingSet, y = value)) + 
  geom_bar(stat = 'summary', fun = 'mean', aes(fill = trainingSet)) + 
  geom_jitter(width = .25, pch = 21, aes(fill = trainingSet)) + 
  stat_summary(fun.data = mean_se,  geom = "errorbar",  width = .5) +
  scale_fill_carto_d(name = "Cell type:", palette = "Vivid") +
  facet_grid(variable~celltype, scales = 'free',space = 'free_x') + 
  theme_bw(base_size = font_fig) + ylim(c(0,NA)) + 
  ggtitle('nonCellEnhLargeGC Test Performance') + 
  guides(fill = guide_legend( nrow = 1)) + 
  theme(legend.position = "bottom", legend.text=element_text(size=10),
        legend.title=element_text(size=font_fig), legend.key.height=unit(.5,"line"), 
        legend.key.width=unit(.5,"line"),  axis.title.x=element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        legend.margin=margin(-10, 0, 0, 0)) 
dev.off()






