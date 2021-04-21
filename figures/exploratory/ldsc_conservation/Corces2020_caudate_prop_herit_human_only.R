#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F, bitmapType='cairo')
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
library(rcartocolor)
library(ggrepel)
library(ggforce)
library(ggpubr)
library(ggsci)
library(here)

# to be run in the root github directory
LABEL='caudate_conservation_ldsc'
PROJDIR=here('data/raw_data/',LABEL)

##########################
# read in the GWAS traits
load(here('data/tidy_data/ldsc_gwas','rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% select( -file) %>% 
  mutate(label = ss(as.character(trait), '_'))

## load enrichments table from RDS
enrichments_fn = here('figures/exploratory/ldsc_conservation','rdas',
               'caudate_conservation_prop_herit_hg_rm_mm.rds')
enrichments = readRDS(file = enrichments_fn)
enrichments = enrichments %>% filter(model_species =='hg38') %>%
  arrange(celltype, peaktype) %>%
  mutate(celltype = ifelse(as.character(celltype)=='Zoonomia', 
                           as.character(peaktype), as.character(celltype)),
         celltype = factor(celltype, unique(celltype)))
enrichments$celltype %>% table()
enrichments$model_species %>% table()


#################################
## make plots for presentation ##
system(paste('mkdir -p', file.path( 'plots')))
height_ppt = 4; width_ppt = 8;
height_fig = 1.75; width_fig = 2.25; font_fig = 7
# plot_celltypes = enrichments %>% filter(Padj < alpha) %>% pull(celltype) %>% unique() 
plot_traits = sort(unique(enrichments$trait))
plot_group = sort(unique(enrichments$group))

# make plots
dir.create(file.path('plots','prop_herit_plots'))
plot_fn = here('figures/exploratory/ldsc_conservation','plots', 
               paste('HgOnly_Corces2020_caudate_prop_herit_enrichments.ppt.pdf'))
pdf(width = width_ppt, height = height_ppt, file = plot_fn, onefile = T)
for (plot in plot_group){
  ## first row - proportion
  pp1 = ggplot(data = enrichments %>% filter(group %in% plot ), 
               aes(y = Enrichment, x = celltype, fill = celltype, color = p.signif)) +
    geom_bar(stat = 'identity') + 
    geom_errorbar(aes(ymin=Enrichment_min, ymax=Enrichment_max),width=.5) +
    scale_fill_carto_d(name = "Cell type:", palette = "Safe") +
    scale_color_manual(values = c(alpha('black', .4), alpha('black',1))) + 
    facet_wrap( ~ trait , scales = 'free_y') +  
    xlab('') +  ylab('Heritability Enrichment') + 
    theme_bw(base_size = 9) + ggtitle(paste("Trait Group:", plot)) + 
    guides(fill = guide_legend(nrow = 2, title.position = 'left'),
           colour = guide_legend(override.aes = list(fill = 'white'))) + 
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 40, hjust = 1),
          legend.text=element_text(size=7), legend.title=element_text(size=8),
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank())

  ## make plots
  plot_fn = here('figures/exploratory/ldsc_conservation','plots','prop_herit_plots',
                      paste('HgOnly_Corces2020_caudate_prop_herit_enrichments',
                            plot,'ppt.pdf', sep = '.'))
  pdf(width = width_ppt, height = height_ppt, file = plot_fn, onefile = F)
  print(pp1)
  dev.off()
  print(pp1)
}
dev.off()



## plot GWAS that are enriched for accelerated or conserved
tmp = enrichments %>% group_by(trait) %>%
  filter(any(grepl("PhyloP|Phast", celltype) & Padj < 0.05)) %>%
  ungroup()
pp2 = ggplot(data = tmp, aes(y = Enrichment, x = celltype, fill = celltype, 
                             color = p.signif)) + geom_bar(stat = 'identity') + 
  geom_errorbar(aes(ymin=Enrichment_min, ymax=Enrichment_max),width=.5) +
  scale_fill_carto_d(name = "Cell type:", palette = "Safe") +
  scale_color_manual(values = c(alpha('black', .4), alpha('black',1))) + 
  facet_wrap( ~ trait , scales = 'free_y') +  
  xlab('') +  ylab('Heritability Enrichment') + 
  theme_bw(base_size = 9) + 
  guides(fill = guide_legend(nrow = 2, title.position = 'left'),
         colour = guide_legend(override.aes = list(fill = 'white'))) + 
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 40, hjust = 1),
        legend.text=element_text(size=7), legend.title=element_text(size=8),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank())

## make plots
plot_fn = here('figures/exploratory/ldsc_conservation','plots','prop_herit_plots',
               'HgOnly_Corces2020_caudate_prop_herit_enrichments_Zoonomia_ppt.pdf')
pdf(width = width_ppt, height = height_ppt, file = plot_fn, onefile = F)
print(pp2)
dev.off()

