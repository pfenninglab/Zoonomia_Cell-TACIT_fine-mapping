## conda activate r4
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(here)
library(ape)     #Phylogenetic tree processing
library(data.table)
library(ggplot2)
library(tidyverse)
library(swfdr)
library(optparse)

CODEDIR='code/raw_code/geno_pheno'
DATADIR='data/raw_data/geno_pheno'
PLOTDIR='figures/exploratory/geno_pheno'
dir.create(here(DATADIR, 'rdas', 'phylolm_out'), showWarnings = F, recursive = T) 
dir.create(here(PLOTDIR, 'plots'), showWarnings = F, recursive = T) 

# Read in tree and traits
traits = read.csv(file= here('data/raw_data/ldsc_caudate_zoonomia/traits_pgls',
                             'traitAnnotations_caudate_phylolm.csv'))

## cell types to do associations over
celltypes = c('MSN_D1', 'MSN_D2', 'MSN_SN', 'Oligo', 'OPC')
perm_types = c( rep("NULL", 2), 'character', rep("NULL", 3), 'numeric')
perm_cols = c(  'peak', 'p.value')

## command line options
option_list = list(
  make_option(c("-c", "--cell"), type="numeric", default=NULL, 
              help="index of cell type to run", metavar="character"),
  make_option(c("-t", "--trait"), type="numeric", default=NULL, 
              help="index of trait to run", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# celltype = celltypes[1]
# idx_trait = names(traits)[-c(1:2)][2]

celltype = celltypes[opt$c]
idx_trait = names(traits)[-c(1:2)][opt$t]

dir.create(here(DATADIR, 'rdas', 'permulated'))
pgls_fn = here(DATADIR, 'rdas', 'phylolm_out', paste0(celltype, '.', idx_trait, '.phylolm.rds'))
print(paste0('Adding permulations for celltype: ',celltype, ', Trait: ', idx_trait))
pgls_df = readRDS(pgls_fn)

perm_fn = list.files(here(DATADIR, 'tables', 'permulations_out'), full.names = T,
                     pattern =  paste0(celltype, '.', idx_trait,'.permulations.'))

## gather permulation p values
pvalue_list = perm_fn %>% 
  lapply(fread, header= F, colClasses = perm_types, col.names = perm_cols) %>%
  rbindlist() %>% split(x = .$p.value, f = .$peak )
ecdf_list =lapply(pvalue_list, ecdf)

## add the empirical p values to permuations
pgls_df = pgls_df %>% filter(peak %in% names(ecdf_list)) %>%
  mutate(numPerm =lengths(pvalue_list)[peak], 
    perm.pvalue = map2_dbl(.x = peak, .y =p.value, .f = ~ecdf_list[[.x]](.y)), 
    perm.pvalue = pmax(1/numPerm, perm.pvalue)) %>%
  arrange(perm.pvalue) %>% 
  # calculate a smarter FDR correction w/ swfdr
  mutate( FDR = lm_qvalue(p.value, X=numMappable)$qvalues, 
          perm_FDR = lm_qvalue(perm.pvalue, X=numMappable)$qvalues)

pgls_fn2 = here(DATADIR, 'rdas', 'permulated', paste0(celltype, '.', idx_trait, '.phylolm.rds'))
pgls_df %>% saveRDS(pgls_fn2)

print(paste('num signif FDRs',  sum(pgls_df$FDR < 0.05)))
print(paste('num signif permulation FDRs', sum(pgls_df$perm_FDR < 0.05)))

qqplot_fn = here(PLOTDIR, 'plots', 
                 paste0(celltype, '.', idx_trait, '.qqplot.pdf'))
tmp_df = data.frame(method = rep(c('pgls', 'permulations'), each = nrow(pgls_df)), 
                    pvalue = c(pgls_df$p.value, pgls_df$perm.pvalue)) %>%
  group_by(method) %>%
  mutate(exp.pvalues= (rank(pvalue, ties.method="first")+.5)/(length(pvalue)+1)) %>% 
  ungroup() %>%
  mutate(log.exp = -log10(exp.pvalues), log.obs = -log10(pvalue))

pdf(qqplot_fn)
# p1 = ggplot(tmp_df, aes(pvalue, fill = method)) +
#   # geom_density(alpha = 0.75) + theme_bw()+
#   geom_histogram(alpha = 0.75, bins = 40) + theme_bw()+
#   theme(legend.position = 'bottom')
# print(p1)
p2 = ggplot(tmp_df, aes(log.exp, log.obs, shape = method, color = method)) +
  geom_point( size = 3) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.5)  + 
  theme_bw()+
  theme(legend.position = 'bottom')
print(p2)
dev.off()
