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
library(ggh4x)

# to be run in the root github directory
PLOTDIR=here('figures/explanatory/','figure2_ldsc_cross-species_grid')

# read in the GWAS traits
load(here('data/tidy_data/ldsc_gwas','rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% select( -file) %>% 
  mutate(label = ss(as.character(trait), '_'))

################################################
## 1) Gather the OCR orthologs + phyloP alone
keepTraits = pheno %>% pull(trait)

## OCR orthologs
keepPeakTypes =  c('Human' = 'Hg peak', 'Human + Macaque' = 'hgRmOrth', 'Human + Mouse' = 'hgMmOrth')
keepPeakTypes2 = names(keepPeakTypes); names(keepPeakTypes2) = keepPeakTypes
enrichments_orthologs = here('figures/exploratory/ldsc_conservation','rdas',
                             'caudate_conservation_prop_herit_hg_rm_mm.rds') %>%
  readRDS() %>% mutate(datatype = 'ATAC-seq') %>%
  filter(trait %in% keepTraits, peaktype %in%keepPeakTypes) %>%
  mutate(peaktype = keepPeakTypes2[as.character(peaktype)],
         celltype = droplevels(celltype), celltype2 = celltype)
# enrichments_orthologs %>% count(celltype, celltype2)

## 
keepPeakTypesPhyloP =  c('PhyloP.accl' = 'PhyloP.accl', 'PhyloP.cons' = 'PhyloP.cons')
keepPeakTypesPhyloP2 = names(keepPeakTypesPhyloP); names(keepPeakTypesPhyloP2) = keepPeakTypesPhyloP
celltypes = enrichments_orthologs %>% pull(celltype) %>% levels()

cell_cols =c(carto_pal(8, "Safe"), 'darkgray', 'black','darkseagreen','darkslategray')
names(cell_cols) = c(celltypes,'PhyloP', 'Mappable', 'M1ctx', 'Liver')

save_fn = here(PLOTDIR,'rdas','figure2_heritability_enrichments_Corces2020_phyloP_CellTACIT-ML.rds')
enrichments = readRDS(save_fn)

################################################
## 3) make multipage cell type and trait enrichment plots
height_ppt = 4; width_ppt = 8;
height_fig = 1.25; width_fig = 1.75; font_fig = 4

## plot the traits
plot_fn = here(PLOTDIR,'plots', paste('Corces2020_caudate.sfig_neuronal_cell_types.pdf', sep = '.'))
pdf(width = 4.75, height = 6, file = plot_fn, onefile = T)

## plot the traits
for (plotTrait in split(keepTraits, rep(1:8, each = 8))) {
  tmp = enrichments %>% filter(trait %in% plotTrait) %>% 
    filter(celltype %in% c('PhyloP', 'MSN_D1','MSN_D2','MSN_SN','INT_Pvalb','M1ctx'))
  if(nrow(tmp) > 0){
    pp1 = ggplot(data = tmp, aes(y = Enrichment, x = peaktype, fill = celltype2,
                                 color = celltype, ymin=Enrichment_min, ymax=Enrichment_max)) +
      geom_bar(stat = 'identity', position = position_dodge2(), size = .2) +
      geom_errorbar(position = position_dodge2(width = 0.5, padding = 0.25), 
                    size = .25, color = 'black') +
      scale_fill_manual(values = cell_cols, guide = 'none') +
      scale_color_manual(values = cell_cols, guide = 'none') +
      facet_nested( label ~ celltype + datatype , scales = 'free',space = 'free_x') +
      xlab('Species/Order') +  ylab('Heritability Enrichment') +
      theme_bw(base_size = font_fig) +
      theme(legend.position = "bottom", panel.spacing=unit(0,"lines"),
            legend.text=element_text(size=font_fig), legend.title=element_text(size=font_fig+1),
            axis.text.x = element_text(angle = 40, hjust = 1, face = 'bold'),
            axis.title.x = element_blank(),
            strip.text.x = element_text(face = 'bold'),
            strip.text.y = element_text(face = 'bold'))
    print(pp1)
  }
}
dev.off()

## plot the glia and liver
plot_fn = here(PLOTDIR,'plots', paste('Corces2020_caudate.sfig_nonneuronal_cell_types.pdf', sep = '.'))
pdf(width = 4.75, height = 6, file = plot_fn, onefile = T)

## plot the traits
for (plotTrait in split(keepTraits, rep(1:8, each = 8))) {
  tmp = enrichments %>% filter(trait %in% plotTrait) %>% 
    filter(celltype %in% c('PhyloP', 'Astro','Microglia','OPC','Oligo','Liver'))
  if(nrow(tmp) > 0){
    pp1 = ggplot(data = tmp, aes(y = Enrichment, x = peaktype, fill = celltype2,
                                 color = celltype, ymin=Enrichment_min, ymax=Enrichment_max)) +
      geom_bar(stat = 'identity', position = position_dodge2(), size = .2) +
      geom_errorbar(position = position_dodge2(width = 0.5, padding = 0.25), 
                    size = .25, color = 'black') +
      scale_fill_manual(values = cell_cols, guide = 'none') +
      scale_color_manual(values = cell_cols, guide = 'none') +
      facet_nested( label ~ celltype + datatype , scales = 'free',space = 'free_x') +
      xlab('Species/Order') +  ylab('Heritability Enrichment') +
      theme_bw(base_size = font_fig) +
      theme(legend.position = "bottom", panel.spacing=unit(0,"lines"),
            legend.text=element_text(size=font_fig), legend.title=element_text(size=font_fig+1),
            axis.text.x = element_text(angle = 40, hjust = 1, face = 'bold'),
            axis.title.x = element_blank(),
            strip.text.x = element_text(face = 'bold'),
            strip.text.y = element_text(face = 'bold'))
    print(pp1)
  }
}
dev.off()

