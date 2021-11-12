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

## 
keepPeakTypesPhyloP =  c('PhyloP.accl' = 'PhyloP.accl', 'PhyloP.cons' = 'PhyloP.cons')
keepPeakTypesPhyloP2 = names(keepPeakTypesPhyloP); names(keepPeakTypesPhyloP2) = keepPeakTypesPhyloP
celltypes = enrichments_orthologs %>% pull(celltype) %>% levels()

cell_cols =c(carto_pal(8, "Safe"), 'darkgray', 'black','darkseagreen','darkslategray')
names(cell_cols) = c(celltypes,'PhyloP', 'Mappable', 'M1ctx', 'Liver')

save_fn = here(PLOTDIR,'rdas','figure2_heritability_enrichments_Corces2020_phyloP_CellTACIT-ML.rds')
enrichments = readRDS(save_fn)

## export heritability enrichments to table
table_fn = here(PLOTDIR,'tables','Table_S4_heritability_enrichments_Corces2020_phyloP_CellTACIT-ML.xlsx')
enrichments2 = enrichments %>% 
  relocate(c(datatype,  celltype, celltype2, peaktype, model_species, group_meta), .after = label) %>%
  relocate(c(Padj, logPadj, p.signif), .after = group_meta) %>%
  relocate(c(starts_with('Coefficient')), .after = group_meta) %>%
  relocate(c(starts_with('Enrichment')), .after = group_meta) %>%
  dplyr::select(-c(group:signif_group, file,annot_group, tmpcol, Categories,
                   Zoonomia.Index:col_meta, ends_with(c('.y','.x')),
                   Observed_scale_h2_z:numSpecies)) %>%
  arrange(datatype,  celltype, celltype2, peaktype, group_meta) %>%
  split(., .$label)

enrichments2 %>% writexl::write_xlsx(table_fn)


################################################
## 3) make plots for each trait cell type grid
height_ppt = 4; width_ppt = 8;
height_fig = 1.25; width_fig = 1.75; font_fig = 4

dir.create(here(PLOTDIR,'plots','prop_herit_bar_plots'), showWarnings = F, recursive = T)
dir.create(here(PLOTDIR,'plots','prop_herit_bar_plots_grids'), showWarnings = F, recursive = T)
plotTrait='DrinksPerWeek_E'; cell = 'MSN_D1'

## plot grid selection of traits and cell types
plotTraits = c('TotCholesterol_E', 'AD_E', 'BrainVol_E', 'Schizophrenia_E', 
               'SmokInitiation_E')
plotCells = c('PhyloP', 'MSN_D1','MSN_D2','MSN_SN','Oligo','Microglia','Liver')
tmp = enrichments %>% filter(trait %in% plotTraits, datatype %in% c('CellTACIT-ML','PhyloP')) %>%
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

plot_fn2 = here(PLOTDIR,'plots','Corces2020_caudate.CellTACIT-ML.ppt2.pdf')
pdf(width = 7, height = 4, file = plot_fn2, onefile = T)
print(pp1)
dev.off()

enrichments %>% filter(trait %in% c('AD_E', 'Schizophrenia_E', 'SmokInitiation_E')) %>%
  filter(datatype == 'PhyloP', peaktype == 'PhyloP.cons') %>%
  group_by(trait, celltype, peaktype) %>%
  summarise(maxEnrich = max(Enrichment), 
            minFDR = min(Padj))

enrichments %>% filter(trait %in% c('AD_E', 'Schizophrenia_E', 'SmokInitiation_E')) %>%
  filter(celltype %in% c('Microglia','MSN_D1', 'MSN_D2')) %>%
  filter(datatype == 'CellTACIT-ML') %>%
  group_by(trait, celltype, celltype2) %>%
  summarise(minEnrich = min(Enrichment), 
            maxEnrich = max(Enrichment), 
            minFDR = min(Padj), 
            maxFDR = max(Padj))



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