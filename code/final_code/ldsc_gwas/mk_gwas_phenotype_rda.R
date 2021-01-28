#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
library(RColorBrewer)
library(rcartocolor)
library(wesanderson)

cols_darjeeling = c(wesanderson::wes_palette('Darjeeling1', 5), 
                    wesanderson::wes_palette('Darjeeling2', 5))

# to be run in the root github directory
LABEL='ldsc_gwas'
setwd('code/final_code/ldsc_gwas')
PROJDIR=file.path('../../../data/tidy_data',LABEL)

# read in the SNP heritability estimates 
pheno = file.path(PROJDIR,'gwas_list_sumstats.tsv') %>% read_tsv() %>%
  filter(!grepl('XX', fig_group)) %>%
  mutate(trait = factor(trait, unique(trait)), 
         group = factor(group, unique(group)))
group_col = cols_darjeeling[seq_along(levels(pheno$group))]
names(group_col) = levels(pheno$group)
pheno = pheno %>% mutate(group_col = group_col[group])

# read in the SNP heritability estimates 
h2_fn = file.path(PROJDIR,'heritability') %>% 
  list.files(path = ., pattern = 'log',full.names = T)
names(h2_fn) = ss(basename(h2_fn),'\\.')
h2 = h2_fn %>% lapply(fread, skip = 'two-step', col.names ='log', nrows = 5, sep = '\n') %>% 
  bind_rows(.id = 'match') %>%
  separate(log, sep = ':', into = c('var', 'val')) %>% 
  mutate(var = make.names(var), val = str_trim(val)) %>%
  pivot_wider(names_from = var, values_from = val) %>%
  separate(Total.Observed.scale.h2, sep = ' ', into = c('h2', 'h2_se')) %>%
  separate(Intercept, sep = ' ', into = c('Inter', 'Inter_se')) %>% 
  mutate_if(is.character, str_replace_all, pattern = "\\(", replacement = "") %>% 
  mutate_if(is.character, str_replace_all, pattern = "\\)", replacement = "") %>%
  select(match:Inter_se) %>% mutate(across(!match, as.numeric)) %>% arrange(desc(h2))

# read in how many SNPs used to estimate heritability 
M = h2_fn %>% lapply(fread, skip = 'After merging with reference panel LD', 
                     col.names ='log', nrows = 1, sep = '\n') %>% 
  bind_rows(.id = 'match') %>%
  separate(log, sep = ',', into = c('var', 'num_SNP')) %>% 
  mutate(var = make.names(var), num_SNP = ss(str_trim(num_SNP),' ')) %>%
  select(match,num_SNP) %>% mutate(across(!match, as.numeric)) %>%
  arrange(desc(num_SNP))

# merge gwas phenotypes and SNP heritabilities
pheno = left_join(h2, M, by = 'match') %>% left_join(pheno, ., by = 'match') %>%
  filter(h2 > 0) %>%  mutate(h2_perSNP = h2 / num_SNP) %>% data.frame()

# export gwas phenotypes w/ SNP heritabilites
system(paste('mkdir -p',file.path(PROJDIR,'rdas')))
save(pheno, group_col, file = file.path(PROJDIR,'rdas','gwas_list_sumstats.rda'))
system(paste('mkdir -p',file.path(PROJDIR,'tables')))
write.table(pheno, file =file.path(PROJDIR,'tables','gwas_list_sumstats_with_h2.tsv'), 
            quote = F, row.names = F, sep = '\t')


