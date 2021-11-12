#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F, bitmapType='cairo')

library(rcartocolor)
library(GenomicRanges)
library(ggtree)
library(ggtreeExtra)
library(tidyverse)
library(tidytree)
library(treeio)
library(here)
library(ape)
library(aplot)
library(ggfun)

DATADIR='data/raw_data/trait_pgls'
PLOTDIR='figures/explanatory/figure4_trait_pgls_loci'
i_am(file.path(PLOTDIR, 'figure4_plot_loci_pgls_tracks_MSN_D1.R'))

##########################
# read in the GWAS traits
load(here('data/tidy_data/ldsc_gwas','rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% select( -file) %>% 
  mutate(label = ss(as.character(trait), '_'))

#  read in the Zoonomia tree and phenotypes annotations
load(here('data/tidy_data/Zoonomia_data', 'rdas','200_Mammals_Genome_Information.rda'))
trait_df$Sleep.total_daily_sleep_time.adult[trait_df$Species== 'Homo_sapiens'] = 8

## read in PGLS results
pgls_rds = here(DATADIR, 'rdas', 'trait_pgls_Corces2020_finemapped_snps_20210824.rds')
pgls_df = readRDS(pgls_rds) %>% filter(zooTrait != "Brain Residual") %>%
  group_by(peakNames, zooTrait) %>% mutate(tmp = sum(PGLS_FDR < .2)) %>% ungroup() %>%
  filter(tmp == 2, grepl(pattern = 'Sleep',group), speciesSet =='Boreoeutheria', celltype == 'MSN_D1') %>% 
  dplyr::select(-tmp)

pgls_gr = pgls_df %>% pull(peakNames) %>% gsub(pattern = '^hg38:|:250$', replacement = '') %>%
  GRanges()
loci = pgls_gr %>% resize( 4*5e3, fix="start") %>% resize( 4e4, fix="end") %>% 
  GenomicRanges::reduce()

##########################
# read in the GWAS traits
celltype = 'MSN_D1'
peaksPred_df = here('data/raw_data/ldsc_caudate_zoonomia','rdas',
                    paste0('Corces2020.',celltype,'.allPeaks.avgCNN.predictions.rds')) %>%
  readRDS() %>% filter(grepl('hg38', name), name %in% pgls_df$peakNames) 

peaksPred_df2 = peaksPred_df %>%
  gather(key = Species, value = value, 2:ncol(peaksPred_df))

###############################################
## trim the tree for species w/o annotations
to_drop1 = tree %>% as_tibble() %>% filter(!is.na(label)) %>%
  filter(!label %in%trait_df$Species ) %>% pull(label)

to_drop2 = trait_df %>% 
  filter(is.na(Sleep.total_daily_sleep_time.adult)) %>%
  pull(Species)

tree2 = drop.tip(tree, c(to_drop1, to_drop2))
tree_dt = full_join(tree2 %>% as_tibble(), 
                    df %>% rename('label' = 'Species')) %>% 
  filter(!is.na(node)) %>% as.treedata()

trait_df$name = 'Sleep'

pdf('tmp.pdf')
g = ggtree(tree_dt, aes(color = Clade, fill = Clade)) + 
  geom_tiplab(aes(label = label),  color = 'black', align=TRUE, size = 2) +
  xlim_tree(.55) + theme_tree2(legend.position = 'bottom')
p2 <- ggplot(trait_df, aes(x = name, y=Species)) + 
  geom_tile(aes(fill=Sleep.total_daily_sleep_time.adult)) + 
  scale_fill_viridis_c(option = "plasma") + 
  theme_tree2() 
p3 <- ggplot(peaksPred_df2, aes(x=name, y=Species)) + 
  geom_tile(aes(fill=value)) + scale_fill_viridis_c() + 
  theme_tree2() 
p2 %>% insert_left(g) %>% insert_right(p3)
dev.off()










  