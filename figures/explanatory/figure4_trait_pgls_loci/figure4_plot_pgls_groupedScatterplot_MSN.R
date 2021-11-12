#######################################
### set up libraries and functions ####
# conda activate r4
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F, bitmapType='cairo')
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
library(rcartocolor)
library(here)
library(rtracklayer)
library(trackViewer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

PLOTDIR='figures/explanatory/figure4_trait_pgls_loci'

##############################################
# 1a) read in the Zoonomia tree and group_meta list
rda_fn = here('data/tidy_data/Zoonomia_data', 
              'rdas','200_Mammals_Genome_Information.rda')
load(file = rda_fn)
df_meta = df_meta %>% mutate(Order2 = ifelse(Time.Since.Split.from.Human.TimeTree.median <=94, as.character(Order), as.character(Clade)))

col_clade = df %>% dplyr::select(Clade, col_clade)%>% filter(!duplicated(Clade)) %>% deframe()
col_meta = df_meta %>% mutate(value = Order2, name = col_meta) %>% 
  filter(!duplicated(value)) %>% dplyr::select(value,name) %>% deframe()

col_celltypes = rcartocolor::carto_pal(n = 8, 'Safe')
names(col_celltypes) = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb',  'Astro', 
                         'Microglia', 'OPC', 'Oligo')  

# 1b) read in the GWAS traits
load(here('data/tidy_data/ldsc_gwas','rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% dplyr::select( -file) %>% 
  mutate(label = ss(as.character(trait), '_'))

##############################################
## read in PGLS results
DATADIR='data/raw_data/trait_pgls'
alpha = 0.10
pgls_rds = here(DATADIR, 'rdas', 'trait_pgls_Corces2020_finemapped_snps_20210914.rds')
pgls_df = readRDS(pgls_rds) %>% filter(zooTrait != "Brain Residual") %>%
  filter(PGLS_FDR < alpha, grepl('MSN', celltype)) %>%
  top_n(3, -PGLS_FDR) %>%
  mutate(name = gsub('^hg38:|:250$', '', peakNames) %>% ss('-', 1),
         name = factor(name, unique(name)))
with(pgls_df, table(celltype, zooTrait))

pgls_gr = pgls_df %>% pull(peakNames) %>% gsub(pattern = '^hg38:|:250$', replacement = '') %>%
  GRanges()
loci = pgls_gr %>% resize( 4e4, fix="start") %>% resize( 2*4e4, fix="end") %>% 
  GenomicRanges::reduce()

pgls_df2 = readRDS(pgls_rds) %>% filter(zooTrait != "Brain Residual") %>%
  filter(grepl('MSN', celltype), grepl('Cross|Morn', label)) %>%
  mutate(name = gsub('^hg38:|:250$', '', peakNames) %>% ss('-', 1),
         name = factor(name, unique(name)))

##########################
# read in the GWAS traits
celltypes = c('MSN_D1' = 'MSN_D1', 'MSN_D2' = 'MSN_D2')
peaksPred_df =celltypes %>% 
  lapply(function(x) here('data/raw_data/ldsc_caudate_zoonomia','rdas',
                          paste0('Corces2020.',x,'.allPeaks.avgCNN.predictions.rds')) %>% readRDS()) %>%
  data.table::rbindlist(idcol = 'celltype') %>%
  filter(grepl('hg38', name), name %in% pgls_df$peakNames)

peaksPred_gr = peaksPred_df %>% mutate(name = gsub('^hg38:|:250$', '', name)) %>%
  pull(name) %>% GRanges()

height_fig = .75; width_fig = 2.2; font_fig = 7

idx = 1

oo = findOverlaps(subject = peaksPred_gr, query = loci)
peaksPred_df2 = peaksPred_df[subjectHits(oo),] %>%
  pivot_longer(!c('name', 'celltype'), names_to = 'Species', values_to = 'value') %>%
  right_join(df[,c('Species', 'Clade')]) %>% rename( 'name' = 'peakNames') %>%
  right_join(trait_df[,c('Species', 'Sleep.total_daily_sleep_time.adult')]) %>% 
  mutate(name = gsub('^hg38:|:250$', '', peakNames) ) %>%
  filter(!is.na(value), !is.na(Sleep.total_daily_sleep_time.adult)) %>% 
  filter(Clade %in% c('Euarchonta', 'Glires'))

with(peaksPred_df2, table(peakNames))
numRow = ceiling(length(unique(peaksPred_df2$name))/6)

plot_fn = with(as.data.frame(loci),here(PLOTDIR, 'plots',paste0('fig4_locus_totalDailySleep.MSN_D2.', unique(seqnames), ':', min(start), '-', max(start), '.pdf')))
pdf(plot_fn, height = height_fig, width = width_fig)
p1 = ggplot(peaksPred_df2, aes(y = value, x = Sleep.total_daily_sleep_time.adult)) + 
  geom_smooth(method = 'lm', color = 'black', size = 1) + 
  geom_smooth(method = 'lm', aes(color = Clade), size = 1) + 
  geom_point(pch = 21, aes(fill = Clade), size = 1) + 
  scale_fill_manual(values = col_clade[c('Euarchonta', 'Glires')], guide = 'none') + 
  scale_color_manual(values = col_clade[c('Euarchonta', 'Glires')], guide = 'none') + 
  ylab('CellTACIT OCR Prediction') + xlab('Daily Sleep (hr)') + 
  facet_wrap(.~ name, scales = 'free_x', nrow = numRow) +  
  ylim(range(peaksPred_df2$value)) +
  theme_bw(base_size = font_fig -3 ) + 
  theme(axis.title=element_text(size=font_fig-3))
print(p1)
dev.off()






