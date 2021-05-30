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
library(ggpattern)
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
keepTraits = c('AD_E', 'SmokInitiation_E')
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

## phylop only
keepPeakTypesPhyloP =  c('PhyloP.accl' = 'PhyloP.accl', 'PhyloP.cons' = 'PhyloP.cons')
keepPeakTypesPhyloP2 = names(keepPeakTypesPhyloP); names(keepPeakTypesPhyloP2) = keepPeakTypesPhyloP
celltypes = enrichments_orthologs %>% pull(celltype) %>% levels()
enrichments_phyloP = here('figures/exploratory/ldsc_conservation','rdas',
                             'caudate_conservation_prop_herit_hg_rm_mm.rds') %>%
  readRDS() %>% mutate(datatype = 'PhyloP') %>%
  filter(trait %in% keepTraits, peaktype %in%keepPeakTypesPhyloP) %>%
  mutate(peaktype = keepPeakTypesPhyloP2[as.character(peaktype)],
         celltype = 'PhyloP', celltype2 = 'PhyloP')
# enrichments_phyloP %>% count(celltype, celltype2)

## mappable to species
keepSpecies = c('Macaque' = 'Macaca_mulatta',
                'Mouse' = 'Mus_musculus',
                'Fruit bat' = 'Rousettus_aegyptiacus')
keepSpecies2 = names(keepSpecies)
names(keepSpecies2) = keepSpecies
enrichments_mappable = here('figures/exploratory/ldsc_caudate_zoonomia',
               'rdas','caudate_conservation_prop_herit_Zoonomia.rds') %>%
  readRDS() %>% mutate(datatype = 'Mappable to') %>%
  filter(trait %in% keepTraits, Species %in% keepSpecies) %>%
  mutate(peaktype = keepSpecies2[as.character(Species)],
         celltype2 = 'Mappable')
# enrichments_mappable %>% count(celltype, celltype2)

## mappable and overlaps + phyloP
enrichments_mapPhylop = here('figures/exploratory/ldsc_caudate_zoonomia',
                             'rdas','caudate_conservation_prop_herit_phyloP_Zoonomia.rds') %>%
  readRDS() %>% mutate(datatype = 'Mappable to') %>%
  filter(trait %in% keepTraits, Species %in% keepSpecies) %>%
  mutate(peaktype = keepSpecies2[as.character(Species)],
         celltype2 = 'PhyloP')
# enrichments_mapPhylop %>% count(celltype, celltype2)

## CellTACIT-ML
keepGroupMeta = c('Human' = 'Primates#0', 'OW Monkey, 20MY' = 'Primates#19.8', 
                  'NW Monkey, 29MY' = 'Primates#28.81', 'Primates, 74MY' = 'Primates#74.1', 
                  'Rodentia, 89MY' = 'Rodentia#89', 'Chiroptera, 94MY' ='Chiroptera#94')
keepGroupMeta2 = names(keepGroupMeta)
names(keepGroupMeta2) = keepGroupMeta
enrichments_meta = here('figures/exploratory/', 'ldsc_zoonomia_meta',
                        'rdas','zoonomia_meta_prop_heritability_Corces2020.rds') %>%
  readRDS() %>% mutate(datatype = 'CellTACIT-ML') %>%
  filter(trait %in% keepTraits, group_meta %in% keepGroupMeta) %>%
  mutate(celltype2 = ifelse(peaktype != 'predActive', 'Mappable', as.character(celltype)),
         peaktype = keepGroupMeta2[as.character(group_meta)])
# enrichments_meta %>% count(celltype, celltype2)


################################################
## 2) Combine dataframes and reset x labels
peaktype2 = c(keepPeakTypesPhyloP2, keepPeakTypes2, keepSpecies2, keepGroupMeta2) %>% unique()
datatype2 = c('PhyloP', 'ATAC-seq', 'Mappable to', 'CellTACIT-ML')
celltypes2 =  c('Mappable', 'PhyloP', celltypes)
enrichments = list(enrichments_phyloP, enrichments_orthologs,  enrichments_mappable, 
                   enrichments_mapPhylop, enrichments_meta) %>% 
  rbindlist(fill=TRUE) %>% 
  mutate(peaktype = factor(peaktype, peaktype2),
         datatype = factor(datatype, datatype2), 
         celltype2 = factor(celltype2, celltypes2),
         label = trait %>% as.character() %>% ss('_'))

cell_cols =c(carto_pal(9, "Safe"), 'black')
names(cell_cols) = c(celltypes,'PhyloP', 'Mappable')
enrichments %>% count(datatype, peaktype)
enrichments %>% count(datatype, celltype2)

dir.create(here(PLOTDIR,'rdas'), showWarnings = F)
save_fn = here(PLOTDIR,'rdas','figure2_heritability_enrichments_Corces2020_phyloP_CellTACIT-ML.rds')
saveRDS(enrichments, file = save_fn)


################################################
## 3) make plots for each trait cell type grid
height_ppt = 4; width_ppt = 8;
height_fig = 1.25; width_fig = 1.75; font_fig = 4

dir.create(here(PLOTDIR,'plots','prop_herit_bar_plots'), showWarnings = F, recursive = T)
dir.create(here(PLOTDIR,'plots','prop_herit_bar_plots_grids'), showWarnings = F, recursive = T)
plotTrait='DrinksPerWeek_E'; cell = 'MSN_D1'

## plot grid selection of traits and cell types
plotTraits = c('AD_E', 'BrainVol_E', 'Gettingup_E', 'Schizophrenia_E', 
               'SmokInitiation_E')
plotCells = c('PhyloP', 'MSN_D1','MSN_D2','INT_Pvalb','Astro','Oligo','Microglia')
tmp = enrichments %>% filter(trait %in% plotTraits, datatype == 'CellTACIT-ML') %>%
  filter(celltype %in% plotCells) %>%
  mutate(trait = factor(trait, plotTraits), celltype = factor(celltype, plotCells),
         label = trait %>% as.character() %>% ss('_')) %>%
  arrange(trait, peaktype) %>% 
  mutate(label = factor(label, unique(label)), 
         peaktype =peaktype  %>% as.character() %>% ss(','), 
         peaktype = factor(peaktype, unique(peaktype)))

plot_fn = here(PLOTDIR,'plots','Corces2020_caudate.CellTACIT-ML.fig2.pdf')
pdf(width = 4.75, height = 4, file = plot_fn, onefile = T)
pp1 = ggplot(data = tmp, aes(y = Enrichment, x = peaktype, fill = celltype2,
                             color = celltype, ymin=Enrichment_min, ymax=Enrichment_max)) +
  geom_bar(stat = 'identity', position = position_dodge2(), size = .2) +
  geom_errorbar(position = position_dodge2(width = 0.5, padding = 0.25), 
                size = .25, color = 'black') +
  scale_fill_manual(values = cell_cols, guide = 'none') +
  scale_color_manual(values = cell_cols, guide = 'none') +
  facet_grid( label ~ celltype , scales = 'free',space = 'free_x') +
  xlab('Species/Order') +  ylab('Heritability Enrichment') +
  theme_bw(base_size = 8) +
  # guides(fill = guide_legend(nrow = 2, title.position = 'left'),
  #        colour = guide_legend(override.aes = list(fill = 'white'))) +
  theme(legend.position = "bottom",
        legend.text=element_text(size=8), legend.title=element_text(size=8+1),
        axis.text.x = element_text(angle = 40, hjust = 1, size = 8 - 2),
        axis.title.x = element_blank(), panel.spacing=unit(0.2,"lines"),
        strip.text.x = element_text(face = 'bold'),
        strip.text.y = element_text(face = 'bold'))
print(pp1)
dev.off()


## plot the traits
for (plotTrait in keepTraits) {
  tmp = enrichments %>% filter(trait %in% plotTrait)
  if(nrow(tmp) > 0){
    plot_fn = here(PLOTDIR,'plots','prop_herit_bar_plots',
                   paste('Corces2020_caudate',plotTrait,'fig2.pdf', sep = '.'))
    pdf(width = 8, height = height_fig, file = plot_fn, onefile = T)
    pp1 = ggplot(data = tmp, aes(y = Enrichment, x = peaktype, fill = celltype2,
                                 color = celltype, ymin=Enrichment_min, ymax=Enrichment_max)) +
      geom_bar(stat = 'identity', position = position_dodge2(), size = .2) +
      geom_errorbar(position = position_dodge2(width = 0.5, padding = 0.25), 
                    size = .25, color = 'black') +
      scale_fill_manual(values = cell_cols, guide = 'none') +
      scale_color_manual(values = cell_cols, guide = 'none') +
      facet_nested( label ~ celltype + datatype , scales = 'free_x',space = 'free_x') +
      xlab('Species/Order') +  ylab('Heritability Enrichment') +
      theme_bw(base_size = font_fig) +
      theme(legend.position = "bottom", panel.spacing=unit(0,"lines"),
            legend.text=element_text(size=font_fig), legend.title=element_text(size=font_fig+1),
            axis.text.x = element_text(angle = 40, hjust = 1, face = 'bold'),
            axis.title.x = element_blank(),
            strip.text.x = element_text(face = 'bold'),
            strip.text.y = element_text(face = 'bold'))
    print(pp1)
    dev.off()
  }
}


## plot the traits and cell types
for (plotTrait in keepTraits) {
  for (cell in celltypes){
    tmp = enrichments %>% filter(trait %in% plotTrait, celltype %in% c(cell,'PhyloP'))
    if(nrow(tmp) > 0){
      pp1 = ggplot(data = enrichments %>% 
                     filter(trait %in% plotTrait, celltype %in% c(cell,'PhyloP')), 
                   aes(y = Enrichment, x = peaktype, fill = celltype2,
                       color = celltype, ymin=Enrichment_min, ymax=Enrichment_max)) +
        geom_bar(stat = 'identity', position = position_dodge2(), size = .2) +
        geom_errorbar(position = position_dodge2(width = 0.5, padding = 0.25), 
                      size = .25, color = 'black') +
        scale_fill_manual(values = cell_cols, guide = 'none') +
        scale_color_manual(values = cell_cols, guide = 'none') +
        facet_grid( label ~  datatype , scales = 'free_x',space = 'free_x') +
        xlab('Species/Order') +  ylab('Heritability Enrichment') +
        theme_bw(base_size = font_fig) +
        # guides(fill = guide_legend(nrow = 2, title.position = 'left'),
        #        colour = guide_legend(override.aes = list(fill = 'white'))) +
        theme(legend.position = "bottom", 
              legend.text=element_text(size=font_fig), legend.title=element_text(size=font_fig+1),
              axis.text.x = element_text(angle = 40, hjust = 1, face = 'bold'),
              axis.title.x = element_blank(),
              strip.text.x = element_text(face = 'bold'),
              strip.text.y = element_text(face = 'bold'))

      ## make plots
      plot_fn2 = here(PLOTDIR,'plots','prop_herit_bar_plots_grids',
                      paste('Corces2020_caudate',plotTrait,cell,'fig.pdf', sep = '.'))
      pdf(width = width_fig, height = height_fig, file = plot_fn2, onefile = F)
      print(pp1); dev.off()
    }
  }
}