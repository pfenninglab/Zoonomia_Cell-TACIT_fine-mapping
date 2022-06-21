#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
library(here)

DATADIR=

alpha = 0.05
caudate_eQTL = here('data/tidy_data/gtex_eqtl', 'GTEx_Analysis_v8_eQTL', 
                    'Brain_Caudate_basal_ganglia.v8.egenes.txt.gz') %>%
  fread() %>% filter(qval <= alpha )

summary(caudate_eQTL$qval)
unique(caudate_eQTL$gene_name)
