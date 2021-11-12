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
library(rcartocolor)

DATADIR='data/raw_data/trait_pgls'
PLOTDIR='figures/explanatory/figure4_trait_pgls_loci'
i_am(file.path(PLOTDIR, 'figure4_plot_pgls_matrix_tracks_MSN.R'))

##########################
# read in the GWAS traits
load(here('data/tidy_data/ldsc_gwas','rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% dplyr::select( -file) %>% 
  mutate(label = ss(as.character(trait), '_'))

#  read in the Zoonomia tree and phenotypes annotations
load(here('data/tidy_data/Zoonomia_data', 'rdas','200_Mammals_Genome_Information.rda'))
trait_df$Sleep.total_daily_sleep_time.adult[trait_df$Species== 'Homo_sapiens'] = 8
col_clade = df %>% dplyr::select(c(Clade, col_clade))%>% filter(!duplicated(Clade)) %>% deframe()
clade_list = with(df, split(Species, Clade))

## read in PGLS results
alpha = 0.10
pgls_rds = here(DATADIR, 'rdas', 'trait_pgls_Corces2020_finemapped_snps_20210914.rds')
pgls_df = readRDS(pgls_rds) %>% filter(zooTrait != "Brain Residual") %>%
  filter(PGLS_FDR < alpha, grepl('MSN', celltype)) %>%
  arrange(PGLS_FDR) %>%
  mutate(name = gsub('^hg38:|:250$', '', peakNames) %>% ss('-', 1),
         name = factor(name, unique(name)))
with(pgls_df, table(celltype, zooTrait))

pgls_gr = pgls_df %>% pull(peakNames) %>% gsub(pattern = '^hg38:|:250$', replacement = '') %>%
  GRanges()
loci = pgls_gr %>% resize( 4*5e3, fix="start") %>% resize( 4e4, fix="end") %>% 
  GenomicRanges::reduce()

##########################
# read in the GWAS traits
celltypes = c('MSN_D1' = 'MSN_D1', 'MSN_D2' = 'MSN_D2')
peaksPred_df =celltypes %>% 
  lapply(function(x) here('data/raw_data/ldsc_caudate_zoonomia','rdas',
                    paste0('Corces2020.',x,'.allPeaks.avgCNN.predictions.rds')) %>% readRDS()) %>%
  data.table::rbindlist(idcol = 'celltype') %>% 
  filter(grepl('hg38', name), name %in% pgls_df$peakNames) 

peaksPred_df2 = peaksPred_df %>%
  pivot_longer(!c('name', 'celltype'), names_to = 'Species', values_to = 'value') %>%
  right_join(df[,c('Species', 'Clade')]) %>% rename('peakNames' = 'name') %>%
  right_join(pgls_df[,c('peakNames', 'speciesSet')]) %>%
  group_by(peakNames) %>%
  mutate(value = ifelse(speciesSet == 'Euarchontoglires' & 
                        !(Clade %in% c('Euarchonta', 'Glires')), 
                      NA, value)) %>% ungroup() %>%
  mutate(name = gsub('^hg38:|:250$', '', peakNames) %>% ss('-', 1),
          name = factor(name, levels(pgls_df$name)))

###############################################
## trim the tree for species w/o annotations
tree = root(tree, 'Homo_sapiens') #reroot to have humans a s
tree = root(tree, df%>% filter(Order =='Primates') %>% pull(Species) %>% as.character()) #reroot to have humans a s
tree = root(tree, df%>% filter(Clade %in% c('Euarchonta', 'Glires')) %>% pull(Species) %>% as.character()) #reroot to have humans a s

to_drop1 = tree %>% as_tibble() %>% filter(!is.na(label)) %>%
  filter(!label %in%trait_df$Species ) %>% pull(label)

to_drop2 = trait_df %>% 
  filter(is.na(Sleep.total_daily_sleep_time.adult)) %>%
  pull(Species)

tree2 = drop.tip(tree, c(to_drop1, to_drop2))
tree_dt = full_join(tree2 %>% as_tibble(), 
                    df %>% rename('label' = 'Species')) %>% 
  filter(!is.na(node)) %>% as.treedata()
tree_dt2 =  groupOTU(tree_dt, clade_list, 'Clade2')

Clade = tree_dt %>% as_tibble() %>% filter(!is.na(Clade)) %>% pull(Clade) %>% unique()
col_clade = col_clade[Clade]
trait_df$name1 = 'Act. Pat'
trait_df$name2 = 'Sleep'

height_ppt = 4; width_ppt = 8;
height_fig = 2; width_fig = 2.25; font_fig = 7

pdf(here(PLOTDIR, 'plots','fig4_tree_pgls_all_peaksMSN_tilePlot.pdf'), 
    height = height_fig * 3, width = width_fig * 3)
g = ggtree(tree_dt2, aes(color = Clade2, fill = Clade2), size = 1) +
  geom_tiplab(aes(label = label),  color = 'black', align=TRUE, size = 2) +
  scale_color_manual(values = col_clade, name = 'Clade') + 
  xlim_tree(1.2) + theme_tree2() 

act.pat.lvls = c('NOCTURNAL','CATHEMERAL','CREPUSCULAR','DIURNAL')
p1 <- ggplot(trait_df, aes(x = name1, y=Species)) + 
  geom_tile(aes(fill=factor(Activity.pattern, act.pat.lvls))) + 
  scale_fill_viridis_d(option = 'plasma', name = 'Act. Pattern') + 
  theme_bw(base_size = font_fig- 3) +
  theme_tree2()

p2 <- ggplot(trait_df, aes(x = name2, y=Species)) + 
  geom_tile(aes(fill=Sleep.total_daily_sleep_time.adult)) + 
  scale_fill_viridis_c( option = 'inferno', name = 'Daily Sleep (hr)') + 
  guides(fill = guide_legend(title.position = "right")) +
  theme_tree2() + theme(legend.title = element_text(angle = -90))

p3 <- ggplot(peaksPred_df2, aes(x=name, y=Species)) + 
  geom_tile(aes(fill=value)) + scale_fill_viridis_c(option = 'inferno', name ='Pred OCR Act.') + 
  facet_grid(.~ celltype, scales = 'free', space = 'free') +
  guides(fill = guide_legend(title.position = "right")) +
  scale_x_discrete(position = "bottom")+ theme_tree2() +  
  theme(legend.title = element_text(angle = -90)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

pg = p2 %>% insert_left(g, width = 3.3) %>%
  insert_right(p3, width = 5) %>%
  insert_right(p1, width = 1) 
print(pg) & theme(legend.position='bottom')
dev.off()




height_fig = 1; width_fig = 2.25; font_fig = 7
pdf(here(PLOTDIR, 'plots','fig4_tree_pgls_signif_MSN_OCR_FDR.pdf'), 
    height = height_fig/2 * 3, width = width_fig/2 * 3)
ggplot(pgls_df, aes(x = name, y=-log10(PGLS_FDR))) + 
  geom_bar(stat = 'identity', aes(fill = celltype)) + 
  geom_hline(yintercept = -log10(alpha), color = 'red') +
  facet_grid(.~ celltype, scales = 'free', space = 'free') +
  scale_fill_carto_d(palette = "Safe", guide = 'none') + ylab('-log10(FDR)') + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))+ 
    theme_bw()
dev.off()


pdf(here(PLOTDIR, 'plots','fig4_tree_pgls_signif_MSN_OCR_Coeff.pdf'), 
    height = height_fig/2 * 3.2, width = width_fig/2 * 3)
ggplot(pgls_df, aes(x = name, y= Coeff)) + 
  geom_bar(stat = 'identity', aes(fill = celltype)) + 
  geom_hline(yintercept = 0, color = 'black') +
  facet_grid(.~ celltype, scales = 'free', space = 'free') +
  scale_fill_carto_d(palette = "Safe", guide = 'none') + ylab('PGLS Coef.') + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))+ 
  theme_bw()
dev.off()











  