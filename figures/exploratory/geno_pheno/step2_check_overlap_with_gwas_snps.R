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
library(ChIPseeker)
library(tidyverse)
library(tidymodels)
library(viridis)

CODEDIR='code/raw_code/geno_pheno'
DATADIR='data/raw_data/geno_pheno'
PLOTDIR='figures/exploratory/geno_pheno'

celltypes = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb',  'Astro',  'Microglia', 'OPC', 'Oligo')
cell_cols = setNames( carto_pal(8, "Safe"), celltypes)

# didn't use interneurons, astrocyte, and microglia
celltypes = celltypes[-c(4:6)]; cell_cols = cell_cols[-c(4:6)]

###########################################################
## 1) read in table of fine-mapped SNPs overlapping OCRs ##
allPeaksSNPs_fn = here('data/raw_data/polyfun_caudate/rdas',
                       'polyfun_caudate_finemapped_snps_20210518.rds')
snpInPeaks_df = allPeaksSNPs_fn %>% readRDS() %>% 
  mutate(trait = droplevels(trait)) 

## traits that are coded as sleep traits
snpInPeaks_sleep_df = snpInPeaks_df %>% filter(group == 'Sleep')
snpInPeaks_sleep_gr = snpInPeaks_sleep_df %>% mutate(tmp = paste0('chr', CHR, ':', POS_hg38)) %>%
  dplyr::select(SNP, tmp) %>% deframe() %>% GRanges()

## traits signif correlated w/ BrainVol reported by Jansen et al.
brainVol_traits = c('BrainVol_E', 'Intelligence_E', 'Neuroticsim_E', 'EduAttain_E',
                    'BMI_E', 'ADHD_E', 'Insomnia_E')
snpInPeaks_brain_df = snpInPeaks_df %>% filter(trait %in% brainVol_traits)
snpInPeaks_brain_gr = snpInPeaks_brain_df %>% mutate(tmp = paste0('chr', CHR, ':', POS_hg38)) %>%
  dplyr::select(SNP, tmp) %>% deframe() %>% GRanges()

#######################################
## 2) read in the complete PGLS analyses
pgls_df = here(DATADIR, 'rdas', 'geno_pheno_Corces2020_finemapped_snps_20220512.rds') %>% readRDS() %>%
  filter(speciesSet == 'All')
alpha = 0.1

## 3. test whether fine-mapped sleep trait SNPs near PGLS signif total daily sleep enhancers
pgls_sleep_df = pgls_df %>% filter(zooTrait == 'Sleep.Total')
pgls_sleep_gr = pgls_sleep_df %>% mutate(tmp = gsub('^hg38:|:250$','', peakNames)) %>%
  pull(tmp) %>% GRanges()
dd = distanceToNearest(pgls_sleep_gr, snpInPeaks_sleep_gr)
distToSNP = sapply(split(mcols(dd)$distance,queryHits(dd)), min)
distToSNP = distToSNP[seq(length(pgls_sleep_gr))]
pgls_sleep_df = pgls_sleep_df %>% mutate(
  disToSNP = )

pdf()
ggplot(pgls_sleep_df, aes(x = disToSNP)) + 
  geom_density(aes(fill = perm_FDR < alpha), alpha = 0.4)+ 
  theme_bw()
dev.off()

     
pgls_sleep_df %>% nest(data = -c(celltype))   %>%
  mutate(fit = map(data, ~ fisher.test(x = .x$hasFineMappedSNP>0, y = .x$perm_FDR < alpha,
                                       alternative = 'greater')),
         tidied = map(fit, tidy)
  ) %>% unnest(tidied) %>% 
  select(-data, -fit) %>% arrange(p.value)


## 3. test whether fine-mapped sleep trait SNPs near PGLS signif total daily sleep enhancers
pgls_brain_df = pgls_df %>% filter(zooTrait == 'Brain.resid')
pgls_brain_gr = pgls_brain_df %>% mutate(tmp = gsub('^hg38:|:250$','', peakNames)) %>%
  pull(tmp) %>% GRanges()




pgls_snp_count_df = pgls_df %>%
  nest(data = -c(celltype, zooTrait))   %>%
  mutate(fit = map(data, ~ fisher.test(x = !is.na(.x$name), y = .x$perm_FDR < alpha,
                                       alternative = 'greater')),
         tidied = map(fit, tidy)
  ) %>% unnest(tidied) %>% 
  select(-data, -fit) %>% arrange(p.value)




