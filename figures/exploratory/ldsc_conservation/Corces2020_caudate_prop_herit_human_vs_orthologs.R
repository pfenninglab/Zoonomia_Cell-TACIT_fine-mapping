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

# to be run in the root github directory
LABEL='caudate_conservation_ldsc'
DATADIR=here('data/raw_data/',LABEL)

##########################
# read in the GWAS traits
load(here('data/tidy_data/ldsc_gwas','rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% select( -file) %>% 
  mutate(label = ss(as.character(trait), '_'))

#########################################
# read in the LDSC partitioned heritability estimation
enrich_fn =here(DATADIR,'prop_herit') %>% 
  list.files(path = ., pattern = '.agg.gz', full.names = T)
names(enrich_fn) = ss(basename(enrich_fn), '.agg.gz')
input = lapply(enrich_fn, read_tsv, col_type = cols()) %>% 
  rbindlist(fill = T, idcol='file') %>% select(file:Coefficient_SE)
input %>% data.frame() %>% head()
input %>% pull(Categories) %>% ss('\\.', 1) %>% table()
input %>% pull(Categories) %>% ss('\\.', 2) %>% table()

#########################################
## format groupings and calculate conditional cell type enrichment p-value
enrichments = input %>% 
  mutate(
    annot_group = ss(file, '\\.', 1),
    # get the 2nd to last slot between '.'
    match = ss(file, '\\.', 2), 
    tmpcol = ss(Categories, '\\.', 1), 
    peaktype = case_when(
      grepl('PhastCons_241mammals', Categories) ~ 'Phast.mam', 
      grepl('PhastCons_43primates', Categories) ~ 'Phast.pri', 
      grepl('phyloPam_accl', Categories) ~ 'PhyloP.accl', 
      grepl('phyloPam_cons', Categories) ~ 'PhyloP.cons', 
      tmpcol == 'Corces2020_caudate' ~ 'Hg peak',
      grepl('Corces2020_caudate_mappedToMm10', tmpcol) ~ 'Hg -> Mm', 
      grepl('Corces2020_caudate_mappedToRheMac10', tmpcol) ~ 'Hg -> Rm', 
      grepl('hgRmOrth',tmpcol)  ~ 'hgRmOrth', 
      grepl('hgMmOrth',tmpcol)  ~ 'hgMmOrth', 
      grepl('Rm|Rhe|Stauffer',tmpcol) ~ 'Rm -> Hg',
      TRUE ~ 'Mm -> Hg'),
    peaktype = factor(peaktype, c('Hg peak', 'Hg -> Rm', 'Hg -> Mm', 
                                  'Rm -> Hg','Mm -> Hg', 'hgRmOrth', 'hgMmOrth',
                                  'PhyloP.accl','PhyloP.cons', 'Phast.mam','Phast.pri')),
    model_species = case_when(
      grepl('PhastCons|phyloP', Categories) ~ 'hg38', 
      tmpcol == 'Corces2020_caudate' ~ 'hg38',
      grepl('Rm|Rhe|Stauffer',tmpcol) ~ 'rheMac10',
      TRUE ~ 'mm10'),
    model_species = factor(model_species, c('hg38', 'rheMac10', 'mm10')),
    annot_group = case_when(
      grepl('binary', annot_group) ~ 'binary annotation', 
      grepl('phyloP', annot_group) ~ 'phyloP score + peak', 
    ),
    celltype = ifelse(grepl('Phast|phylo',Categories), 'Zoonomia', 
                      ss(Categories,'\\.', 2)),
    celltype = celltype %>% 
      factor(c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb',  'Astro', 
               'Microglia', 'OPC', 'Oligo','Zoonomia')), 
    cell_group = case_when(
      grepl('MSN|INT|Zoonomia', celltype) ~ 'Neuron', 
      TRUE ~ 'Glia'
    )) %>% inner_join(x = pheno, by = 'match') %>% 
  group_by(file) %>% type_convert() %>%
  filter(!grepl('CP.MSN_D|CP.INT', Categories), !is.na(celltype))

enrichments$peaktype %>% table()
enrichments$celltype %>% table()
enrichments$model_species %>% table()

enrichments = enrichments %>% mutate(
  ## Truncate heritability predictions to be nonnative
  Observed_scale_h2 = pmax(Observed_scale_h2, 0),
  Observed_scale_h2_min = pmax(Observed_scale_h2 - Observed_scale_h2_SE, 0), 
  Observed_scale_h2_max = pmax(Observed_scale_h2 + Observed_scale_h2_SE, 0), 

  # percent of total h2g heritability
  Proportion_of_h2g = Observed_scale_h2 / h2,
  Proportion_of_h2g_min = Observed_scale_h2_min / h2,
  Proportion_of_h2g_max = Observed_scale_h2_max / h2,
  
  # SNP enrichment
  Enrichment = Proportion_of_h2g / Proportion_of_SNPs,
  Enrichment_min = Proportion_of_h2g_min / Proportion_of_SNPs,
  Enrichment_max = Proportion_of_h2g_max / Proportion_of_SNPs,
)

enrichments %>% data.frame() %>% head()
summary(enrichments$Enrichment)

## normalize the coefficients by per SNP heritability
# compute Padj w/ bonferroni family mutliple hypothesis correction
alpha = 0.05;
enrichments = enrichments %>% 
  mutate(
  # normalize conditional coefficients
  Coef_norm = Coefficients / h2_perSNP,
  Coef_norm_se = Coefficient_SE / h2_perSNP, 

  # calculate 2-sided p-value for non-zero h2
  Observed_scale_h2_z = Observed_scale_h2 / Observed_scale_h2_SE,
  Observed_scale_h2_P_value = pnorm(- Observed_scale_h2_z ),

  # calculate 2-sided p-value for non-significant contribution
  Coef_norm_z = Coef_norm / Coef_norm_se,
  Coef_norm_z_P_value = pnorm(- Coef_norm_z ),
  Padj = p.adjust(Coef_norm_z_P_value, 'fdr'), 
  logPadj = -log10(Padj),
  p.signif = ifelse(Padj < alpha, paste('FDR <',alpha), 'NS'),
  p.signif = factor(p.signif, levels = c('NS', paste('FDR <',alpha)))
  ) %>% ungroup() %>% filter(complete.cases(Coef_norm_z_P_value))

# look at top signif celltypes
enrichments %>% filter(Padj < alpha) %>% 
  select(match, celltype, annot_group, Padj) %>%
  group_by(match) %>% top_n(1, -Padj) %>% ungroup() %>% 
  data.frame() %>% arrange(Padj)

table(enrichments$p.signif)


#################################
## save enrichments table to RDS
save_fn = here('figures/exploratory/ldsc_conservation','rdas',
               'caudate_conservation_prop_herit_hg_rm_mm.rds')
saveRDS(enrichments, file = save_fn)


#################################
## make plots for presentation ##
system(paste('mkdir -p', here( 'plots')))
height_ppt = 4; width_ppt = 8;
height_fig = 1.75; width_fig = 2.25; font_fig = 7
# plot_celltypes = enrichments %>% filter(Padj < alpha) %>% pull(celltype) %>% unique() 
dir.create(here('figures/exploratory/ldsc_conservation', 'plots','prop_herit_plots'))
plot_traits = sort(unique(enrichments$trait))

# make plots
plot_fn = here('figures/exploratory/ldsc_conservation','plots', 
               paste('Corces2020_caudate_prop_herit_enrichments.ppt.pdf'))
pdf(width = width_ppt, height = height_ppt, file = plot_fn, onefile = T)
for (lab in plot_traits){
  plot_cells = enrichments %>% filter(trait %in% lab) %>% 
    filter(Padj < alpha) %>% pull(celltype) %>% unique()
  plot_cells = unique(c(plot_cells, c('MSN_D1', 'MSN_D2', 'MSN_SN', 'INT_Pvalb','Zoonomia')))
  ## first row - enrichment

  pp1 = ggplot(data = enrichments %>% filter(trait %in% lab) %>% 
                 filter(celltype %in% plot_cells), 
               aes(y = Enrichment, x = peaktype, fill = celltype, color = p.signif)) +
    geom_bar_pattern( aes(pattern = model_species, pattern_colour = p.signif), 
                      stat = 'identity', pattern_angle = 60, pattern_spacing = 0.05, 
                      pattern_density = .1, pattern_fill = 'black') + 
    scale_pattern_manual(values = c(hg38 = "none", mm10 = 'stripe', rheMac10 = 'circle'), name = '') +
    scale_pattern_colour_manual(values = c('grey20', 'black'), guide = 'none') + 
    geom_errorbar(aes(ymin=Enrichment_min, ymax=Enrichment_max), 
                  width=.2, position=position_dodge(.9)) +
    scale_fill_carto_d(name = "Cell type:", palette = "Safe") +
    scale_color_manual(values = c(alpha('black', .4), alpha('black',1))) + 
    facet_grid( ~ celltype , scales = 'free', space = 'free_x' ) +  
    xlab('') +  ylab('Enrichment') + 
    theme_bw(base_size = 10) + ggtitle(paste("Trait:", lab)) + 
    guides(fill = FALSE, colour = guide_legend(override.aes = list(fill = 'white')),
           pattern = guide_legend(
             override.aes = list(fill = 'white', color = 'black', pattern_size = .2,
                                 pattern_spacing = 0.01, pattern_density = 0.2))) + 
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1),
          legend.text=element_text(size=7), legend.title=element_text(size=8),
          axis.title.x = element_blank())
  
  
  ## 2nd row - proportions 
  pp2 = ggplot(data = enrichments %>% filter(trait %in% lab ) %>%
                 filter(celltype %in% plot_cells), 
               aes(y = Proportion_of_h2g, x = peaktype, fill = celltype, color = p.signif)) +
    geom_bar_pattern( aes(pattern = model_species, pattern_colour = p.signif), 
                      stat = 'identity', pattern_angle = 60, pattern_spacing = 0.05, 
                      pattern_density = .1, pattern_fill = 'black') + 
    scale_pattern_manual(values = c(hg38 = "none", mm10 = 'stripe', rheMac10 = 'circle'), name = '') +
    scale_pattern_colour_manual(values = c('grey20', 'black'), guide = 'none') + 
    geom_errorbar(aes(ymin=Proportion_of_h2g_min, ymax=Proportion_of_h2g_max), 
                  width=.2, position=position_dodge(.9)) +
    scale_fill_carto_d(name = "Cell type:", palette = "Safe") +
    scale_color_manual(values = c(alpha('black', .4), alpha('black',1))) + 
    facet_grid( ~ celltype , scales = 'free_x', space = 'free_x' ) +  
    xlab('') +  ylab('Prop. Heritability') + 
    theme_bw(base_size = 10) + 
    guides(fill = FALSE, colour = guide_legend(override.aes = list(fill = 'white')),
           pattern = guide_legend(
             override.aes = list(fill = 'white', color = 'black', pattern_size = .2,
                                 pattern_spacing = 0.01, pattern_density = 0.2))) + 
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1),
          legend.text=element_text(size=7), legend.title=element_text(size=8),
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank())
  
  ## make plots
  plot_fn = here('figures/exploratory/ldsc_conservation','plots','prop_herit_plots',
                      paste('Corces2020_caudate_prop_herit_enrichments',
                            lab,'ppt.pdf', sep = '.'))
  pdf(width = width_ppt, height = height_ppt, file = plot_fn, onefile = F)
  pp = ggarrange(pp1, pp2, nrow = 2, common.legend = TRUE, 
                 align = 'h', legend="bottom")
  print(pp)
  dev.off()
  print(pp)
}
dev.off()


# make plots
plot_fn2 = here('figures/exploratory/ldsc_conservation', 'plots',
                'Corces2020_caudate_prop_herit_enrichments.fig.pdf')
pdf(width = width_fig, height = height_fig, file = plot_fn2)
pp = ggplot(data = enrichments %>% filter(trait %in% plot_traits), 
            aes(y = logPadj, x = celltype, fill = celltype, color = p.signif)) +
  geom_hline(yintercept = -log10(alpha), color = 'red',linetype = 'dashed') + 
  geom_point(aes(alpha = p.signif, shape = cell_group)) + 
  scale_color_manual(values = c('grey', 'black'), guide = 'none') + 
  scale_alpha_manual(values = c(.5, 1), guide = 'none') + 
  scale_fill_carto_d(name = "", palette = "Safe") +
  scale_shape_manual(values = c(24, 21), breaks = c('Neuron', 'Glia'), guide = 'none') + 
  facet_row(. ~ group, scales = 'free_x', space = 'free') +  
  xlab('GWAS Trait Conditional Cell Type Enrichment') + 
  ylab('-log10(Bonf. adj. P)') + 
  theme_bw(base_size = font_fig) + 
  guides(fill = guide_legend(
    override.aes = list(shape = c(24, 24, 24, 24, 21, 21, 21, 21)), nrow = 1)) + 
  theme(legend.position = "bottom", legend.text=element_text(size=font_fig-2),
        legend.title=element_text(size=font_fig-2),
        legend.key.height=unit(.5,"line"), 
        legend.key.width=unit(.5,"line")) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(pp)
dev.off()



