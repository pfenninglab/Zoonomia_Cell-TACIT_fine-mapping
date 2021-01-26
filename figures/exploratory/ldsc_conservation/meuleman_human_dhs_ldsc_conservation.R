#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F)
suppressMessages(library(tidyverse))

# to be run in the root github directory
LABEL='meuleman_dhs_index'
setwd('figures/exploratory/ldsc_conservation')
PROJDIR=file.path('../../../data/raw_data/',LABEL)

#########################################
# read in the LDSC heritability estimation
enrich_fn =file.path(PROJDIR,'enrichments') %>% 
  list.files(path = ., pattern = 'results', full.names = T)
names(enrich_fn) = basename(enrich_fn)
enrichments = lapply(enrich_fn, read_tsv) %>% bind_rows(.id = 'file')

## format groupings and calculate conditional cell type enrichment p-value
enrichments = enrichments %>% 
  mutate(
    celltype = ss(file, '\\.', 1), 
    gwas = ss(file, '\\.', 2), 
    trait = ss(file, '-', 1), 
    ref = ss(file, '-', 2), 
    tmp_col1 = ss(Categories,'\\.', 1),
    group = case_when(
      tmp_col1 == 'DHS_Index_and_Vocabulary' ~ 'DHS_hg38Peaks',
      tmp_col1 == 'DHS_Index_and_Vocabulary' ~ 'DHS_mappedToMm10',
      TRUE ~ 'baselineLD_v2.1'
    ), 
    peaktype = ss(Categories,'\\.', 2),
    Coefficient_p = 2*pnorm(q=abs(Coefficient_z), lower.tail=FALSE), 
    Coefficient_fdr = p.adjust(Coefficient_p, 'fdr')) %>%
  select(-tmp_col1)

## truncate SE range of `Observed_scale_h2` to be non-negative
## this effects Prop_of_h2g and Enrichment estimates
enrichments = enrichments %>% 
  # filter(grepl('DHS_Index_and_Vocabulary', Categories)) %>% 
  group_by(file) %>%
  mutate(Total_h2g = mean(Observed_scale_h2 / Proportion_of_h2g)) %>% 
  ungroup() %>% mutate(
    Observed_scale_max = pmax(Observed_scale_h2 + Observed_scale_h2_SE, 0),
    Observed_scale_min = pmax(Observed_scale_h2 + Observed_scale_h2_SE, 0),
    Observed_scale_trun = pmax(Observed_scale_h2, 0), 
    Proportion_of_h2g_trun = Observed_scale_trun / Total_h2g, 
    Enrichments_trun = Proportion_of_h2g_trun / Proportion_of_SNPs 
  )

enrichments %>% filter(Coefficient_fdr < .05) %>% 
  filter(group != 'baselineLD_v2.1') %>% 
  arrange(Coefficient_fdr) %>% data.frame() %>% head(.,20)




