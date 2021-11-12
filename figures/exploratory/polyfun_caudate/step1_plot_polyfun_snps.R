ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(here)
library(tidyverse)
library(rcartocolor)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(GenomicRanges)
library(rtracklayer)

CODEDIR='code/raw_code/polyfun_caudate'
DATADIR='data/raw_data/polyfun_caudate'
PLOTDIR='figures/exploratory/polyfun_caudate'

####################################
## 1) read in table of fine-mapped SNPs ##
dir.create(here('data/raw_data/polyfun_caudate/rdas'), showWarnings = F)
poly_fn = here('data/raw_data/polyfun_caudate/rdas',
               'polyfun_caudate_finemapped_snps_20210518.rds')
snps_df = readRDS(file = poly_fn) %>%
  mutate(SNP = paste(CHR,POS_hg38,SNP,A1,A2, sep = ':'))


nrow(snps_df) # 98468
snps_df %>% filter(PIP > 0.1) %>% tally() #7423
snps_df %>% filter(PIP > 0.5) %>% tally() #1796
snps_df %>% filter(PIP > 0.95) %>% tally() #748

snps_df %>% group_by(SNP) %>% mutate(tmp = n()) %>% filter(tmp > 1) %>%
  ungroup() %>% filter(!duplicated(SNP)) %>% nrow() # 9744

snps_df %>% group_by(SNP) %>%  mutate(tmp = n()) %>% ungroup() %>% 
  filter(!duplicated(SNP)) %>% pull(tmp) %>% summary()




#################################
## make plots for presentation ##
dir.create(here(file.path(PLOTDIR, 'plots')), showWarnings = F, recursive = T)
height_ppt = 4; width_ppt = 8;
height_fig = 1.75; width_fig = 2.25; font_fig = 7

plot_fn = here(PLOTDIR,'plots',
               paste0('polyfun_caudate_finemapped_snps_overlap_HAR_CHAR_20210518.ppt.pdf'))
# pdf(plot_fn, height = height_ppt, width = width_ppt)
for (cutoff in c(.95, .9, .5, .25, .10)){
  plot_fn = here(PLOTDIR,'plots',
                 paste0('polyfun_zoonomia_finemapping_PIP',cutoff,'_20210513.ppt.pdf'))
  pdf(plot_fn, height = height_ppt, width = width_ppt)
  
  pp = ggplot(data = snps_df %>% filter(PIP >= cutoff) , aes(x = label)) +
    geom_bar(aes(fill = top_phyloP)) + 
    geom_text(stat='count', aes(label = ..count.. ), vjust=-1, size= 2.5)+
    scale_fill_manual(values = top_phyloP_cols) + 
    facet_grid(~group, scales = 'free', space = 'free') + 
    scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
    xlab('Fine-mapping group') + ylab(paste0('Number of SNPs w/ PIP >',cutoff)) +
    theme_classic(base_size = 6 ) + guides(fill = guide_legend(nrow = 1)) + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    theme(legend.key.size = unit(.5, 'cm'), legend.position = 'bottom') 
  print(pp)
  
  pp = ggplot(data = snps_df %>% filter(PIP >= cutoff) , aes(x = label)) +
    geom_bar(aes(fill = top_phastCons)) + 
    geom_text(stat='count', aes(label = ..count.. ), vjust=-1, size= 2.5)+
    scale_fill_manual(values = top_phastCons_cols) + 
    facet_grid(~group, scales = 'free', space = 'free') + 
    scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
    xlab('Fine-mapping group') + ylab(paste0('Number of SNPs w/ PIP >',cutoff)) +
    theme_classic(base_size = 6 ) + guides(fill = guide_legend(nrow = 1)) + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    theme(legend.key.size = unit(.5, 'cm'), legend.position = 'bottom') 
  print(pp)
  
  pp = ggplot(data = snps_df %>% filter(PIP >= cutoff) , aes(x = label)) +
    geom_bar(aes(fill = cCRE_group)) + 
    geom_text(stat='count', aes(label = ..count.. ), vjust=-1, size= 2.5)+
    scale_fill_manual(values = ENCODE3_cCRE_cols) + 
    facet_grid(~group, scales = 'free', space = 'free') + 
    scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
    xlab('Fine-mapping group') + ylab(paste0('Number of SNPs w/ PIP >',cutoff)) +
    theme_classic(base_size = 6 ) + guides(fill = guide_legend(nrow = 1)) + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    theme(legend.key.size = unit(.5, 'cm'), legend.position = 'bottom') 
  print(pp)
  dev.off()
}



#######################################
## plot numFine-mapped SNPs by power ##
snps_df2 = snps_df %>% filter(PIP >= .1) %>%
  group_by(trait, h2, h2_Z, group) %>%
  summarize(numSNP = n())

lm(numSNP ~ h2 +h2_Z  , data = snps_df2) %>% summary()
snps_df2 %>% group_by(group) %>% summarize(mean(numSNP))


  
