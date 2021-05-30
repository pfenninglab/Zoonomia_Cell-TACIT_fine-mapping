#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F, bitmapType='cairo')

library(rcartocolor)
library(ggrepel)
library(ggforce)
library(here)
library(tidyverse)
library(rcartocolor)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(ggtree)
library(ggtreeExtra)

LABEL='polyfun_caudate'
CODEDIR='code/raw_code/polyfun_caudate'
DATADIR='data/raw_data/polyfun_caudate'
PLOTDIR='figures/exploratory/polyfun_caudate'
i_am(file.path(PLOTDIR, 'step2_plot_trait_enhancer_pgls_associations.R'))

##########################
# read in the GWAS traits
load(here('data/tidy_data/ldsc_gwas','rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% select( -file) %>% 
  mutate(label = ss(as.character(trait), '_'))

##########################
# read in the Zoonomia tree and species list
rda_fn = here('data/tidy_data/Zoonomia_data', 
              'rdas','200_Mammals_Genome_Information.rda')
load(file = rda_fn)
df = df %>% arrange(Clade) %>%
  mutate(Clade = ifelse(grepl('Xen|Afro', as.character(Clade)), 
                        'Xenarthra & Afrotheria', as.character(Clade)),
         Clade = factor(Clade, unique(Clade)))
col_clade = df %>% dplyr::select(c(Clade, col_clade))%>% filter(!duplicated(Clade)) %>% deframe()
col_order = df %>% dplyr::select(c(Order, col_order)) %>% filter(!duplicated(Order)) %>% deframe()

##########################
# read in the PGLS results
alpha = 0.05
pgls_rds = here(DATADIR, 'rdas', 'trait_pgls_Corces2020_finemapped_snps_brainResid.rds')
pgls_df = readRDS(pgls_rds) %>% filter(PGLS_FDR < alpha)

mammalTraits_df = here('data/raw_data/ldsc_caudate_zoonomia','traits_pgls',
                       'traitAnnotations_caudate.csv') %>% fread() %>%
  rename_with(make.names) %>% dplyr::rename( 'Species' = 'species.binomial') %>%
  dplyr::select(c(Species,contains('eg'))) %>% rename_with(~ gsub('_eg','',.x)) %>%
  rowwise(Species) %>% mutate(tmp = any(!is.na(c_across(ActivityPattern:Sleep.Total_daily_sleep_time.adult)))) %>%
  filter(tmp) %>% dplyr::select(-tmp)


pgls_df2 =  pgls_df %>% filter(!duplicated(celltype)) %>% 
  arrange(PGLS_FDR) %>%
  mutate(peakNames = factor(peakNames, rev(unique(peakNames)))) %>%
  mutate(tmp = here('data/raw_data/ldsc_caudate_zoonomia','rdas',
                    paste0('Corces2020.',celltype,'.allPeaksCalibWithSNPs.avgCNN.predictions.rds'))) %>%
           pull(tmp) %>% lapply(readRDS) %>% rbindlist() %>%
  rename('peakNames' = 'name' ) %>%
  inner_join(pgls_df, by = 'peakNames') %>%
  pivot_longer(cols = Acinonyx_jubatus:Ziphius_cavirostris, names_to = 'Species',
               values_to = 'score') %>%
  left_join(df, by = 'Species') %>%
  inner_join(mammalTraits_df, by = 'Species') %>%
  filter(!is.na(Brain.resid))

#################################
## make plots for presentation ##
height_ppt = 4; width_ppt = 8;
height_fig = 1.75; width_fig = 2.25; font_fig = 7

gene_text <- pgls_df2 %>% group_by(peakNames, celltype) %>% 
  summarise(label = paste(paste0('FDR=',signif(unique(PGLS_FDR), 4)),
                          unique(Gene.Symbol), sep = '\n'))

pgls_pdf = here(PLOTDIR, 'plots', 'trait_pgls_Corces2020_finemapped_snps_brainResid.ppt.pdf')
pdf(pgls_pdf, height = height_ppt, width = width_ppt)
ggplot(pgls_df2, aes(Brain.resid, score)) + 
  geom_smooth(method = 'lm') + geom_point(pch = 21, aes(fill = Order)) + 
  scale_fill_manual(values = col_order) + 
  facet_grid(celltype ~ peakNames) + 
  geom_text(mapping = aes(x = +Inf, y = -Inf, label = label),
            data = gene_text, hjust   = 1, vjust   = -.5, size = 5 ) +
  guides(fill = guide_legend(nrow =1, title.position="left", override.aes = list(size = 4))) + 
  theme_bw(base_size = 14) + 
  xlab('Brain Volume Adj. Body Size') + ylab('Predicted OCR Activity') + 
  theme(legend.position = "bottom", legend.text=element_text(size=font_fig),
        legend.title=element_text(size=font_fig), legend.key.height=unit(.5,"line")) 

dev.off()







