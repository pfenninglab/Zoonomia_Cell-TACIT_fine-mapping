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

# to be run in the root github directory
LABEL='Zoonomia_data'
DATADIR=here('data/raw_data/ldsc_caudate_zoonomia')

##############################################
# read in the Zoonomia tree and species list
rda_fn = here('data/tidy_data/Zoonomia_data', 
              'rdas','200_Mammals_Genome_Information.rda')
load(file = rda_fn)
col_clade = df %>% select(Clade, col_clade)%>% filter(!duplicated(Clade)) %>% deframe()
col_order = df %>% select(Order, col_order) %>% filter(!duplicated(Order)) %>% deframe()
celltypes = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb', 'Astro', 'Microglia', 'OPC', 'Oligo')


###################################
### plot validation performance ###
cell = 'MSN_D1'
outEnh_fn = here(DATADIR,'rdas',paste0('Corces2020.',celltypes, '.enhPeaks.avgCNN.predictions.rds'))
outCalib_fn = here(DATADIR,'rdas',paste0('Corces2020.',celltypes, '.calibPeaks.avgCNN.predictions.rds'))
names(outEnh_fn) = names(outCalib_fn) = celltypes

df_enh = outEnh_fn %>% map(readRDS) %>% rbindlist(idcol = 'celltype') %>% 
  pivot_longer(-c(name, celltype), names_to = "Species", values_to = 'scores') %>%
  replace(is.na(.), 0) %>% group_by(Species, celltype) %>% 
  summarise(avg_score = mean(scores, na.rm = T), 
            sd_score = sd(scores, na.rm = T)) %>%
  mutate(celltype = factor(celltype, celltypes)) %>% inner_join(df, by = 'Species')

df_calib = outCalib_fn %>% map(readRDS) %>% rbindlist(idcol = 'celltype') %>% 
  pivot_longer(-c(name, celltype), names_to = "Species", values_to = 'scores') %>%
  replace(is.na(.), 0) %>% group_by(Species, celltype) %>% 
  summarise(avg_score = mean(scores, na.rm = T), 
            sd_score = sd(scores, na.rm = T)) %>%
  mutate(celltype = factor(celltype, celltypes)) %>% inner_join(df, by = 'Species')

###################################
### plot validation performance ###
library(rcartocolor)
height_ppt = 5; width_ppt = 8;
height_fig = 1.75; width_fig = 4.75; font_fig = 7

pdf(here('figures/exploratory/cnn_enhancer_ortholog', 
         'plots/cross_species_predictions_withNA_performance_20210508.pdf'), 
    width = width_ppt, height = height_ppt)
ggplot(data = df_enh, aes(x = Time.Since.Split.from.Human.TimeTree.median, y = avg_score)) + 
  geom_point(pch = 21, aes(fill = Order)) + 
  scale_fill_manual(values = col_order) + 
  facet_wrap(~celltype, scales = 'fixed', nrow = 2) + 
  theme_bw(base_size = 9) + 
  guides(fill = guide_legend( nrow = 4)) + 
  xlab(paste0('MY from Human (TimeTree median est.)')) + 
  ylab(paste0('Average Cell type CNN Score')) + 
  theme(legend.position = "bottom", legend.text=element_text(size=10),
        legend.title=element_text(size=10),
        legend.key.height=unit(.5,"line"), 
        legend.key.width=unit(.5,"line"), 
  ) 
ggplot(data = df_calib, aes(x = Time.Since.Split.from.Human.TimeTree.median, y = avg_score)) + 
  geom_point(pch = 21, aes(fill = Order)) + 
  scale_fill_manual(values = col_order) + 
  facet_wrap(~celltype, scales = 'fixed', nrow = 2) + 
  theme_bw(base_size = 9) + 
  guides(fill = guide_legend( nrow = 4)) + 
  xlab(paste0('MY from Human (TimeTree median est.)')) + 
  ylab(paste0('Average Cell type CNN Calibrated Score')) + 
  theme(legend.position = "bottom", legend.text=element_text(size=10),
        legend.title=element_text(size=10),
        legend.key.height=unit(.5,"line"), 
        legend.key.width=unit(.5,"line"), 
  ) 
dev.off()





