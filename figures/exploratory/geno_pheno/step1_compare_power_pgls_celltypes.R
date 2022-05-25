ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(here)
library(rcartocolor)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(tidymodels)
library(viridis)

CODEDIR='code/raw_code/geno_pheno'
DATADIR='data/raw_data/geno_pheno'
PLOTDIR='figures/exploratory/geno_pheno'

#######################################
## 0) read in the complete PGLS analyses
celltypes = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb',  'Astro',  'Microglia', 'OPC', 'Oligo')
cell_cols = setNames( carto_pal(8, "Safe"), celltypes)
celltypes = celltypes[-c(4:6)]; cell_cols = cell_cols[-c(4:6)]

pgls_df = here(DATADIR, 'rdas', 'geno_pheno_Corces2020_finemapped_snps_20220512.rds') %>% readRDS()



## 1. plot the p-value histograms
pdf(here(PLOTDIR,'plots', 'plot_phistogram_pgls_all_mammals_permulations.pdf'),
    height = 3, width = 6)
ggplot(pgls_df, aes(x = perm.pvalue, fill = celltype)) + 
  geom_density(adjust = .5) +  facet_grid(zooTrait~celltype, scales = 'free_y') + 
  scale_fill_manual(values = cell_cols) + 
  theme_bw(base_size = 10) + xlab('PGLS P-value') + 
  theme(legend.position = 'none')
dev.off()



## 2. plot the p-values from euarchontoglires or all placental mammals
signif_df = pgls_df %>% group_by(peakNames) %>% filter(any(perm_FDR < 0.05)) %>% 
  ungroup() %>%  mutate(perm.pvalue = -log10(perm.pvalue)) %>%
  pivot_wider(id_cols = c('zooTrait', 'celltype', 'peakNames'), 
              names_from = 'speciesSet', values_from = 'perm.pvalue', values_fill = 0)

pdf(here(PLOTDIR,'plots', 'plot_pvalue_changes_with_more_species.pdf'),
    height = 3, width = 6)
ggplot(signif_df, aes(x = Euarchontoglires, y =All)) + 
  geom_abline( slope= 1, intercept = 0, color = 'red', linetype = 'dashed') + 
  geom_point(pch = 20, alpha= .3) +  facet_grid(zooTrait~celltype) + 
  theme_bw(base_size = 10) +
  xlab('-log10(permulated P-value) Euarchontoglires') + 
  ylab('-log10(permulated P) Placental mammals') + 
  ggtitle('Change in power to detect significant associations with more species')
dev.off()

## 3. plot the p-values from pgls or permulations
signifPerm_df = pgls_df %>% filter(speciesSet == 'All') %>%
  group_by(peakNames) %>% filter(any(perm_FDR < 0.05 | FDR < 0.05)) %>% 
  ungroup() %>%  mutate(perm.pvalue = -log10(perm.pvalue), p.value = -log10(p.value))

pdf(here(PLOTDIR,'plots', 'plot_pvalue_changes_with_permulations.pdf'),
    height = 3, width = 6)
ggplot(signifPerm_df, aes(x = p.value, y =perm.pvalue)) + 
  geom_abline( slope= 1, intercept = 0, color = 'red', linetype = 'dashed') + 
  geom_point(pch = 20, alpha= .3) +  facet_grid(zooTrait~celltype) + 
  theme_bw(base_size = 10) +
  xlab('-log10(P-value) from PGLS') + ylab('-log10(P-value) from permulations')
dev.off()

### 4. plot the overall # of significant peaks using all species ###
signifSummary_df = pgls_df %>% filter(perm_FDR < 0.05) %>%
  group_by(speciesSet, celltype, zooTrait) %>% summarise(num = n()) %>% ungroup() %>%
  pivot_wider(id_cols = c('zooTrait', 'celltype'), 
              names_from = 'speciesSet', values_from = 'num', values_fill = 0)

signifSummary_df %>% writexl::write_xlsx(here(PLOTDIR,'tables', 'numPeak_changes_with_more_species.xlsx'))

pdf(here(PLOTDIR,'plots', 'plot_numPeak_changes_with_more_species.pdf'),
    height = 3, width = 6)
ggplot(signifSummary_df, aes(x = Euarchontoglires, y =All)) + 
  geom_abline( slope= 1, intercept = 0, color = 'red', linetype = 'dashed') + 
  geom_point(pch = 21, size = 2, aes(fill = celltype)) + 
  facet_wrap(~zooTrait, scales = 'free') + 
  scale_fill_manual(values = cell_cols) + 
  theme_bw(base_size = 12) + 
  theme(legend.position = 'bottom',
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-5,-5)) + 
  xlab('Number of peaks FDR < 0.05 using only Euarchontoglires') + 
  ylab('# signif peaks with \nPlacental mammals')
dev.off()


### 5. plot the # signif peaks using permulations vs.  ###
signifWithPerm_df = pgls_df %>% filter(speciesSet == 'All') %>%
  group_by(zooTrait, celltype) %>% 
  summarise(num_pgls = sum(FDR < 0.05), num_perm = sum(perm_FDR < 0.05)) %>% ungroup()


pdf(here(PLOTDIR,'plots', 'plot_numPeak_changes_with_permulations.pdf'),
    height = 3, width = 6)
ggplot(signifWithPerm_df, aes(x = num_pgls, y =num_perm)) + 
  geom_abline( slope= 1, intercept = 0, color = 'red', linetype = 'dashed') + 
  geom_point(pch = 21, size = 2, aes(fill = celltype)) + 
  facet_wrap(~zooTrait, scales = 'free') + 
  scale_fill_manual(values = cell_cols) + 
  theme_bw(base_size = 12) + 
  theme(legend.position = 'bottom',
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-5,-5)) + 
  xlab('Number of peaks FDR < 0.05 using PGLS') + 
  ylab('# signif peaks w/ permulations')
dev.off()
