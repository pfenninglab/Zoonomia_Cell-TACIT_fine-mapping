#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F, bitmapType='cairo')
suppressMessages(library(tidyverse))
suppressMessages(library(tidymodels))
suppressMessages(library(data.table))
suppressMessages(library(here))

#########################################
# to be run in the root github directory
LABEL='ldsc_zoonomia_meta'
PLOTDIR=here('figures/exploratory/', LABEL)
DATADIR=here('data/raw_data/',LABEL)

#################################
## load enrichments table from RDS
enrichments_fn = here(PLOTDIR,'rdas','zoonomia_meta_prop_heritability_Corces2020.rds')
enrichments = readRDS(file = enrichments_fn)

zoo_meta_lm = enrichments %>%
  group_by(trait, celltype) %>% nest() %>%
  mutate(
    fit = map(data, ~ lm(Enrichment ~ Time.Since.Split.from.Human.TimeTree.median:peaktype , data = .x)),
    tidied = map(fit, tidy)
  ) %>% 
  unnest(tidied) %>% ungroup() %>%
  select(-data, -fit) %>% filter(term != '(Intercept)') %>%
  mutate(term = gsub('Time.Since.Split.from.Human.TimeTree.median', 'hMYA', term),
         term = gsub('peaktype', '', term)) %>%
  pivot_wider(values_from = c(estimate, std.error, statistic, p.value), 
              names_from = term)
  
zoo_meta_lm = zoo_meta_lm %>%
  arrange(desc(`estimate_hMYA:predActive` > `estimate_hMYA:mappable`),
          desc(`estimate_hMYA:predActive` > 0), 
          desc(`estimate_hMYA:mappable` > 0),
          sign(`estimate_hMYA:predActive`) * log10(`p.value_hMYA:predActive`), 
          sign(`estimate_hMYA:mappable`) * log10(`p.value_hMYA:mappable`))
          
dir.create(here(PLOTDIR,'tables'))
lm_tsv_fn = here(PLOTDIR,'tables','zoonomia_meta_prop_heritability_Corces2020_linRegCoefficients.xlsx')
writexl::write_xlsx(zoo_meta_lm, lm_tsv_fn)
