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
library(rcartocolor)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(ggtree)
library(ggtreeExtra)
library(tidyverse)

LABEL='trait_pgls'
CODEDIR='code/raw_code/trait_pgls'
DATADIR='data/raw_data/trait_pgls'
PLOTDIR='figures/exploratory/trait_pgls'
i_am(file.path(PLOTDIR, 'step1_plot_trait_enhancer_pgls_associations_sleepTraits.R'))

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
alpha = .10
pgls_rds = here(DATADIR, 'rdas', 'trait_pgls_Corces2020_finemapped_snps_20210914.rds')
pgls_df = readRDS(pgls_rds) %>% filter(zooTrait != "Brain Residual") %>%
  group_by(peakNames) %>% mutate(tmp = sum(PGLS_FDR < .2)) %>% ungroup() %>%
  filter(tmp == 2) %>% dplyr::select(-tmp)

mammalTraits_df = here('data/raw_data/ldsc_caudate_zoonomia','traits_pgls',
                       'traitAnnotations_caudate.csv') %>% fread() %>%
  rename_with(make.names) %>% dplyr::rename('Species' =  'Species.Name') %>%
  rowwise(Species) %>% dplyr::select(-contains('_eg')) %>%
  mutate(tmp = any(!is.na(c_across(ActivityPattern:Sleep.Total_daily_sleep_time.adult)))) %>%
  filter(tmp) %>% dplyr::select(-tmp)

predictions_df = pgls_df %>% filter(!duplicated(celltype)) %>%
  dplyr::select(celltype) %>%
  mutate(tmp = here('data/raw_data/ldsc_caudate_zoonomia','rdas',
                    paste0('Corces2020.',celltype,'.allPeaksCalibWithSNPs.avgCNN.predictions.rds')),
         predList = map(tmp, readRDS)) %>% 
  unnest(cols = predList) %>% dplyr::select(-tmp) %>%
  rename('peakNames' = 'name' )

pgls_df2 = predictions_df  %>%
  inner_join(pgls_df, by = c('celltype', 'peakNames')) %>%
  pivot_longer(cols = Acinonyx_jubatus:Ziphius_cavirostris, names_to = 'Species',
               values_to = 'score') %>%
  left_join(df, by = 'Species') %>%
  inner_join(mammalTraits_df, by = 'Species') %>%
  mutate(Clade = droplevels(Clade)) %>%
  filter(!is.na(Sleep.Total_daily_sleep_time.adult), !is.na(score)) %>%
  group_by(peakNames, speciesSet) %>%
  mutate(
    tmpInd = ifelse(speciesSet == 'Boreoeutheria', Clade != 'Xenarthra & Afrotheria', 
                    Clade %in% c('Euarchonta','Glires'))
  ) %>% filter(tmpInd) %>% dplyr::select(-tmpInd) %>%
  ungroup() %>%
  mutate(
    peakNames = peakNames %>% as.character(), 
    peakNames = gsub('hg38:', '', peakNames), 
    peakNames = gsub(':250', '', peakNames)
  )

#################################
## make plots for presentation ##
height_ppt = 4; width_ppt = 8;
height_fig = 1.75; width_fig = 2.25; font_fig = 7

gene_text <- pgls_df2 %>% filter(zooTrait == 'dailySleepAdult') %>%
  group_by(peakNames, celltype, speciesSet) %>% 
  summarise(label = paste(paste0('FDR=',signif(unique(PGLS_FDR), 4)),
                          paste('near', unique(Gene.Symbol)), sep = ', '))

pgls_pdf = here(PLOTDIR, 'plots', 'trait_pgls_Corces2020_finemapped_snps_totalSleepAdult.20210910.ppt.pdf')
pdf(pgls_pdf,  width = width_ppt, height = width_ppt)
ggplot(pgls_df2 %>% filter(zooTrait == 'dailySleepAdult'),
       aes(y = Sleep.Total_daily_sleep_time.adult, x = score)) + 
  geom_smooth(method = 'lm', color = 'black') + 
  geom_smooth(method = 'lm', aes(color = Clade)) + 
  geom_point(pch = 21, aes(fill = Clade)) + 
  scale_fill_manual(values = col_clade, name = '' ) + 
  facet_grid(speciesSet~peakNames, scales = 'free_x') + 
  geom_text(mapping = aes(x = -Inf, y = +Inf, label = label),
            data = gene_text, hjust   = -.1, vjust   = 1.1, size = 3 ) +
  guides(fill = guide_legend(nrow =1, title.position="left", override.aes = list(size = 4))) + 
  theme_classic(base_size = 14) + 
  ylab('Total daily sleep in adults (hrs)') + xlab('Predicted D1 MSN OCR Activity') + 
  theme(legend.position = "bottom", legend.text=element_text(size=font_fig-2),
        legend.title=element_text(size=font_fig), legend.key.height=unit(.5,"line"),
        strip.text.y = element_text(size = font_fig)) 
dev.off()







