#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F)
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
library(ggrepel)
library(ggpubr)
library(here)
library(tidymodels)

# to be run in the root github directory
LABEL='caudate_conservation_ldsc'
PROJDIR= 'figures/exploratory/ldsc_conservation'
DATADIR=here('data/raw_data/',LABEL)

##########################
# read in the GWAS traits
load(here('data/tidy_data/ldsc_gwas','rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% select( -file) %>% mutate(label = ss(as.character(trait), '_'))

#################################
## save enrichments table to RDS
save_fn = here('figures/exploratory/ldsc_conservation','rdas',
               'caudate_conservation_ldsc_coeff_hg_rm_mm.rds')
enrichments = readRDS(file = save_fn)
save_fn2 = here('figures/exploratory/ldsc_conservation','rdas',
                'caudate_conservation_ldsc_coeff-wide_hg_rm_mm.rds')
enrich_wide = readRDS(file = save_fn2)

enrichments2 = enrichments %>% 
  dplyr::select(-c(reference:group_col, h2:h2_perSNP, Name, file, annot_group, subgroup:label)) %>%
  relocate(celltype:Coef_norm_se, .before = everything()) %>%
  relocate(c(peaktype:model_species, Coefficient_P_value), .after = 'cell_group') %>%
  arrange(celltype,model_species,peaktype) %>% split(., f = .$trait)
save_fn3 = here(PROJDIR,'tables','Table_S2_ldsc_coefficient_caudate_celltype_cross-species.xlsx')
# writexl::write_xlsx(enrichments2, save_fn3)


#################################
## make plots for presentation ##
save_fn4 = here(PROJDIR,'tables','Table_S3_cross-species_cross-trait_LDSC_effect_size_increase_regression.xlsx')
tidy_mod = enrich_wide %>% filter(p.signif != 'NS') %>%
  mutate(model_species = factor(model_species, c('rheMac10', 'mm10'))) %>%
  nest(data = -c(cell_group, celltype)) %>%
  mutate( test = map(data, ~ lm(norm_coeff_diff ~ model_species/norm_coeff_mean + peaktype, data = .x)),
          tidied = map(test, tidy)
  )  %>% 
  unnest(cols = tidied) %>% 
  filter(term != '(Intercept)') %>%
  # select(-data, -test, -term) %>%
  select(-data, -test) %>%
  arrange(term, desc(cell_group), celltype) %>%
  mutate(FDR = p.adjust(p.value, 'fdr'))

writexl::write_xlsx(tidy_mod, save_fn4)




