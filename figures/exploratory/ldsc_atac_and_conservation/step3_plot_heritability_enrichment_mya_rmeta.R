#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F, bitmapType='cairo')
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
library(rcartocolor)
library(ggrepel)
library(ggforce)
library(ggpubr)
library(ggsci)
library(here)
library(metafor)
library(tidymeta)
library(ggplot2)
library(dplyr)
library(broom)

#########################################
# 0) to be run in the root github directory
LABEL='ldsc_atac_and_conservation'
PLOTDIR=here('figures/exploratory/', 'ldsc_atac_and_conservation')
DATADIR=here('data/raw_data/',LABEL)

# read in the GWAS traits
load(here('data/tidy_data/ldsc_gwas','rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% select( -file) %>% 
  mutate(label = ss(as.character(trait), '_'))

##############################################
# 1) read in the Zoonomia tree and group_meta list
rda_fn = here('data/tidy_data/Zoonomia_data', 
              'rdas','200_Mammals_Genome_Information.rda')
load(file = rda_fn)
df_meta = df_meta %>% mutate(Order2 = ifelse(Time.Since.Split.from.Human.TimeTree.median <=94, as.character(Order), as.character(Clade)))
col_meta = df_meta %>% mutate(value = Order2, name = col_meta) %>% 
  filter(!duplicated(value)) %>% dplyr::select(value,name) %>% deframe()

## load in the saved enrichments table to RDS
dir.create(here(PLOTDIR,'rdas'), showWarnings = F)
save_fn = here(PLOTDIR,'rdas','zoo_meta_atac_and_conservation_prop_heritability_Corces2020.rds')
enrichments = readRDS(file = save_fn) %>% 
  dplyr::rename('TimeSplit' = 'Time.Since.Split.from.Human.TimeTree.median')

##########################################################
## 2) compute meta-analyses of the heritability enrichment
enrichments_meta = enrichments %>% 
  nest(data = -c(group, celltype, peaktype, constype, TimeSplit)) %>% 
  mutate(meta = map(data, ~ rma(method = 'HS', slab = file, data = .x, 
                                yi = Enrichment, sei = Enrichment_max - Enrichment)), 
         tidied = map(meta, tidy)) %>% 
  unnest(tidied) %>% dplyr::select(-c(data))


#################################
## 3) make plots for presentation
dir.create(here(PLOTDIR,'plots'))
dir.create(here(PLOTDIR,'plots','prop_herit'))
height_ppt = 4; width_ppt = 8;
height_fig = 1.75; width_fig = 2.25; font_fig = 7
plot_group = sort(unique(enrichments$group))

enrichments_meta = enrichments_meta %>% mutate(`Cons + ATAC` = constype)

#### make plots
plot_fn = here(PLOTDIR,'plots',paste('meta_group_All_prop_herit_enrichments.ppt.pdf', sep = '_'))
pdf(width = width_ppt, height = 6, file = plot_fn, onefile = T)
for(cell in unique(enrichments_meta$celltype)){
  tmp = enrichments_meta %>% filter(celltype %in% cell)
  plot_fn = here(PLOTDIR,'plots', 'prop_herit', paste('meta_group', 
                                       cell, 'prop_herit_enrichments.ppt.pdf', sep = '_'))
  pdf(width = width_ppt, height = 6, file = plot_fn)
  ## the SNP enrichments ##
  pp1 = ggplot(data = tmp, aes(x = TimeSplit + 1, y = estimate)) +
    # add line fit
    geom_smooth(aes(color = `Cons + ATAC`), method = "lm", 
                formula = y ~ x + I(x^2), size = 1, se = TRUE) +  
    geom_errorbar(aes(fill = `Cons + ATAC`, shape = `Cons + ATAC`, 
                      ymin = pmax(0, estimate-std.error), 
                      ymax = estimate+std.error), width = 0.1) + 
    geom_point(aes(fill = `Cons + ATAC`, shape = `Cons + ATAC`)) +
    scale_fill_brewer(palette = "Set1") + scale_shape_manual(values = c(21:24)) + 
    scale_color_brewer(palette = "Set1") + scale_x_continuous(trans='log2') + 
    facet_grid(peaktype ~ group, scales = 'free_y') + ylim(c(0,NA)) +
    theme_bw(base_size = 10) + ggtitle(cell) +
    xlab(paste0('Million years from Human')) + ylab(paste0('Heritability Enrichment')) + 
    guides(shape = guide_legend(nrow =1, title.position="top", override.aes = list(size = 3))) + 
    theme(legend.position = "bottom", legend.text=element_text(size=font_fig),
          legend.title=element_text(size=font_fig), legend.key.height=unit(.5,"line"))
  print(pp1)
  dev.off()
  print(pp1)
}
dev.off()



#### make plots
for(peak in unique(enrichments_meta$peaktype)){
  tmp = enrichments_meta %>% filter(peaktype %in% peak)
  plot_fn = here(PLOTDIR,'plots', paste('meta_group', peak, 'prop_herit_enrichments.ppt.pdf', sep = '_'))
  pdf(width = width_ppt, height = 12, file = plot_fn)
  ## the SNP enrichments ##
  pp1 = ggplot(data = tmp, aes(x = TimeSplit + 1, y = estimate)) +
    # add line fit
    geom_smooth(aes(color = `Cons + ATAC`), method = "lm", 
                formula = y ~ x + I(x^2), size = 1, se = TRUE) +  
    geom_errorbar(aes(fill = `Cons + ATAC`, shape = `Cons + ATAC`, 
                      ymin = pmax(0, estimate-std.error), 
                      ymax = estimate+std.error), width = 0.1) + 
    geom_point(aes(fill = `Cons + ATAC`, shape = `Cons + ATAC`)) +
    scale_fill_brewer(palette = "Set1") + scale_shape_manual(values = c(21:24)) + 
    scale_color_brewer(palette = "Set1") + scale_x_continuous(trans='log2') + 
    facet_grid(group ~ celltype, scales = 'free_y') + ylim(c(0,NA)) +
    theme_bw(base_size = 10) + ggtitle(peak) +
    xlab(paste0('Million years from Human')) + ylab(paste0('Heritability Enrichment')) + 
    guides(shape = guide_legend(nrow =1, title.position="top", override.aes = list(size = 3))) + 
    theme(legend.position = "bottom", legend.text=element_text(size=font_fig),
          legend.title=element_text(size=font_fig), legend.key.height=unit(.5,"line"))
  print(pp1)
  dev.off()
}





