ss <- function(x, pattern, slot = 1, ...) { 
  sapply(base::strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(here)
library(rcartocolor)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(tidyverse)
library(tidymodels)
# library(conflicted)
# rename <- dplyr::rename; select <- dplyr::select
tidymodels_prefer()
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))

CODEDIR='code/raw_code/trait_pgls'
DATADIR='data/raw_data/trait_pgls'
PLOTDIR='figures/e# install.packages("devtools")
devtools::install_github("r-lib/conflicted")xploratory/trait_pgls'
i_am(file.path(CODEDIR, 'step2_chisqTest_sleepTrait_zooTrait_overlap.R'))

#### read in the PGLS runs ####
pgls_rds = here(DATADIR, 'rdas', 'trait_pgls_Corces2020_finemapped_snps_20210824.rds')
pgls_df = readRDS(pgls_rds)

##########################
# read in the GWAS traits
load(here('data/tidy_data/ldsc_gwas','rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% select( -file) %>% 
  mutate(label = ss(as.character(trait), '_'))

# grab the labels that match sleep GWASs
matchSleep = pheno %>% filter(group == 'Sleep') %>% 
  filter(grepl('Morningness|Insomnia|Gettingup|Sleep', label)) %>%
  pull(label) %>% paste(collapse = '|')

## check hypergeometric overlap
alpha = .05
pgls_df %>% filter(celltype != 'Oligo') %>%
  group_by(peakNames, celltype, zooTrait) %>%
  summarise(inSleep = grepl(matchSleep, label),
         isSignif = any(PGLS_Pvalue < alpha)) %>%
  ungroup() %>% 
  distinct(peakNames, celltype, zooTrait, inSleep, isSignif) %>%
  nest(data = -c(celltype, zooTrait)) %>%
  mutate(
    # overlap of OCR w/ SNP from sleep GWAS w/ marginal signif PGLS pvalue
    test = map(data, ~ fisher.test(.x$inSleep, .x$isSignif)), # S3 list-col
    tidied = map(test, tidy)
  ) %>% 
  unnest(cols = tidied) %>% select(-data, -test) %>% arrange(p.value)


# grab the labels that match neuro GWASs
matchNeuro = pheno %>% filter(group == 'Neuro') %>% 
  filter(grepl('Morningness|Insomnia|Gettingup|Sleep', label)) %>%
  pull(label) %>% paste(collapse = '|')
matchNeuro = "EduAttain|Income|BrainVol|Neuroticsim|Intelligence|Headache"
matchNeuro = "EduAttain|BrainVol|Intelligence"
matchNeuro = "BrainVol"


## check hypergeometric overlap of brain residual
alpha = .1
pgls_df %>% filter(celltype == 'Oligo') %>%
  group_by(peakNames, celltype, zooTrait) %>%
  summarise(inSleep = grepl(matchNeuro, label),
            isSignif = any(PGLS_Pvalue < alpha)) %>%
  ungroup() %>% 
  distinct(peakNames, celltype, zooTrait, inSleep, isSignif) %>%
  nest(data = -c(celltype, zooTrait)) %>%
  mutate(
    # overlap of OCR w/ SNP from sleep GWAS w/ marginal signif PGLS pvalue
    test = map(data, ~ fisher.test(.x$inSleep, .x$isSignif)), # S3 list-col
    tidied = map(test, tidy)
  ) %>% 
  unnest(cols = tidied) %>% select(-data, -test) %>% arrange(p.value)

