#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F)
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
library(rcartocolor)
library(ggrepel)

# to be run in the root github directory
LABEL='caudate_conservation_ldsc'
setwd('figures/exploratory/ldsc_conservation')
PROJDIR=file.path('../../../data/raw_data/',LABEL)

##########################
# read in the GWAS traits
load(file.path('../../../data/tidy_data/ldsc_gwas','rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% select( -file) %>% mutate(label = ss(as.character(trait), '_'))

#########################################
# read in the LDSC partitioned heritability estimation
enrich_fn =file.path(PROJDIR,'enrichments') %>% 
  list.files(path = ., pattern = '.results.gz', full.names = T)
names(enrich_fn) = ss(basename(enrich_fn), '.results.gz')
input = lapply(enrich_fn, read_tsv) %>% rbindlist(fill = T, idcol='file') %>%
  select(file:Coefficient_SE)
input %>% data.frame() %>% head()

#########################################
## format groupings and calculate conditional cell type enrichment p-value
enrichments = input %>% 
  mutate(
    annot_group = ss(file, '\\.', 2),
    match = ss(file, '\\.', 3), 
    tmpcol = ss(file, '\\.', 1), 
    peaktype = case_when(
      tmpcol == 'Corces2020_caudate' ~ 'hgPeak',
      grepl('peakInMouse', tmpcol) ~ 'hgToMm', 
      grepl('hgMmOrth', tmpcol) ~ 'hgMmOrth', 
      TRUE ~ 'mmToHg',
    ),
    peak_group = case_when(
      peaktype == 'hgPeak' ~ 'hgPeak' ,
      TRUE ~ 'inMm'
    ),
    annot_group = case_when(
      grepl('binary', annot_group) ~ 'binary annotation', 
      grepl('phyloP', annot_group) ~ 'phyloP score + peak', 
    ),
    celltype = ss(Categories,'\\.', 2),
    celltype = celltype %>% 
      factor(c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb', 'Astro', 
               'Microglia', 'OPC', 'Oligo')), 
    cell_group = case_when(
      grepl('MSN|INT', celltype) ~ 'Neuron', 
      TRUE ~ 'Glia'
    )) %>% inner_join(x = pheno, by = 'match') %>% 
  filter(! group %in% c('Other', 'Metabolism') ) %>%
  group_by(file) %>% type_convert() %>%
  filter(!grepl('CP.MSN_D|CP.INT', Categories)) %>%
  filter(grepl('Corces2020', Categories)) %>%
  filter(Categories != 'Corces2020_caudate.All_Roadmap_DHS_cCREs.BG')

enrichments$file %>% ss('\\.') %>% table()

## Truncate heritability preditions to be positive ##
enrichments = enrichments %>% mutate(
  Observed_scale_h2 = pmax(Observed_scale_h2, 0),
  Proportion_of_h2g = Observed_scale_h2 / h2,
  Enrichment = Proportion_of_h2g / Proportion_of_SNPs,
  Observed_scale_h2_min = pmax(Observed_scale_h2 - Observed_scale_h2_SE, 0), 
  Observed_scale_h2_max = pmax(Observed_scale_h2 + Observed_scale_h2_SE, 0), 
)

enrichments %>% data.frame() %>% head()
summary(enrichments$Enrichment)

## normalize the coefficients by per SNP heritability
# compute Padj w/ bonferroni family mutliple hypothesis correction
enrichments = enrichments %>% group_by(annot_group, peak_group) %>% 
  mutate(
  # normalize conditional coefficients
  Coef_norm = Coefficients / h2_perSNP,
  Coef_norm_se = Coefficient_SE / h2_perSNP, 
  # calculate 2-sided p-value for non-zero effect
  Coef_z = Coef_norm / Coef_norm_se,
  Coefficient_P_value = pnorm(- Coef_z ),
  Padj = p.adjust(Coefficient_P_value, 'fdr'), 
  logPadj = -log10(Padj)
) %>% ungroup() %>% filter(complete.cases(Coefficient_P_value))
# look at top signif celltypes
enrichments %>% filter(Padj < 0.05) %>% 
  select(match, celltype, annot_group, Padj) %>%
  group_by(match) %>% top_n(1, -Padj) %>% ungroup() %>% 
  data.frame() %>% arrange(Padj)

###################################################################
## pivot table to pair up peak_group with each celltype and trait  ##
alpha = 0.05; 
enrich_wide = enrichments %>% 
  # these id columns make up all combinations aside form peakgroup
  pivot_wider(id_col = c(eval(names(pheno)), celltype, cell_group,  annot_group), 
              names_from = c(peak_group), values_from = where(is.numeric)) %>% 
  mutate(
    p.signif = case_when(
      Padj_hg38 < alpha & Padj_mappableToMm10 < alpha ~ 'both',
      Padj_hg38 < alpha ~ 'hg38', 
      Padj_mappableToMm10 < alpha ~ 'mappableToMm10',
      TRUE ~ 'NS'
    ),
    pointsize = max((Padj_hg38 + Padj_mappableToMm10)^2, .1),
    norm_coeff_diff = Coef_norm_mappableToMm10 - Coef_norm_hg38, 
    norm_coeff_mean = (Coef_norm_mappableToMm10 + Coef_norm_hg38)/2, 
    label = ifelse(p.signif != 'NS',as.character(celltype),NA),
    ) %>% arrange(trait) %>%
mutate(trait = as.character(trait) %>% ss('_') %>% factor(., unique(.)))

# look at conditionally enriched cell types
# both 8        hg38 6      mappableToMm10 3
enrich_wide%>% pull(p.signif) %>% table()

#################################
## make plots for presentation ##
system(paste('mkdir -p', file.path( 'plots')))
height_ppt = 5; width_ppt = 8
height_fig = 6; width_fig = 2.25; font_fig = 7
# plot_celltypes = enrichments %>% filter(Padj < alpha) %>% pull(celltype) %>% unique() 
plot_traits = enrich_wide %>% filter(p.signif !='NS') %>% pull(trait) %>% unique() 

# make plots
plot_fn = file.path('plots','Corces2020_caudate_conditional_herit_enrichments.ppt.pdf')
pdf(width = width_ppt, height = height_ppt, file = plot_fn)
  pp = ggplot(data = enrich_wide %>% filter(trait %in% plot_traits), 
              aes(y = norm_coeff_diff, x = norm_coeff_mean, 
                  fill = celltype, color = p.signif)) +
    geom_hline(yintercept = 0, color = 'black') + 
    geom_vline(xintercept = 0, color = 'black') + 
    geom_point(aes(alpha = p.signif != 'NS', shape = annot_group), size =2.5) + 
    scale_color_manual(values = c('black', 'darkred', 'blue','grey50'), 
                       name = paste('P_bonf <',alpha)) + 
    scale_shape_manual(values = c(21,22), name = '') + 
    scale_fill_carto_d(name = "Cell type:", palette = "Safe") +
    scale_alpha_manual(values = c(.2, 1), guide = 'none') + 
    coord_cartesian(ylim=c(-2, 4), xlim=c(0, 4)) +
    facet_wrap( ~ trait, scales = 'fixed', nrow =2 ) +  
    geom_label_repel(aes(label = label), box.padding = .1, label.size = .01,
                     max.overlaps = 40, size = 3, show.legend = F,na.rm = T,
                     point.padding = .1, segment.color = 'grey50', max.time = 2,
                     min.segment.length = .1, alpha = .7, 
                     label.padding = .1, force_pull = 1, force = 40) +
    xlab('Avg. Normalized Herit. Contribution (mappedToMm10 + hg38)/2') + 
    ylab('Norm. Contrib. Diff. (mappedToMm10 - hg38)') + 
    theme_bw(base_size = 12) + 
    guides(colour = guide_legend(nrow = 2), 
           fill = guide_legend(override.aes = list(shape = 21), nrow = 2),
           shape = guide_legend(nrow = 2)) + 
    theme(legend.position = "bottom", 
          legend.text=element_text(size=7),
          legend.title=element_text(size=8))
  print(pp)
dev.off()


# make plots
plot_fn2 = file.path('plots','Corces2020_caudate_conditional_herit_enrichments.fig.pdf')
pdf(width = width_fig, height = height_fig, file = plot_fn2)
for(cell in c("Neuron", "Glia")){
  pp = ggplot(data = enrich_wide %>% filter(cell_group == cell), 
              aes(y = norm_coeff_diff, x = norm_coeff_mean, 
                  fill = group, color = p.signif)) +
    geom_point(pch = 21, aes(alpha = p.signif != 'NS',), size =1) + 
    scale_color_manual(values = c('black', 'darkred', 'blue','grey50'), 
                       name = paste('P_bonf <',alpha)) + 
    scale_fill_manual(values = group_col, guide = 'none') + 
    scale_alpha_manual(values = c(.2, 1), guide = 'none') + 
    geom_abline( slope = 0, intercept = 0, color = 'black', linetype = 'dashed') + 
    coord_cartesian(ylim=c(-2, 4), xlim=c(0, 4)) +
    facet_grid(celltype ~ annot_group, scales = 'fixed') +  
    geom_label_repel(aes(label = label), box.padding = .08, label.size = .01,
                     max.overlaps = 40, size = 1.2, show.legend = F,na.rm = T,
                     point.padding = .2, segment.color = 'grey50', max.time = 5,
                     min.segment.length = .08, alpha = .8,
                     label.padding = .08, force_pull = 1, force = 40) +
    xlab('Avg. Norm. Contrib. (mappedToMm10 + hg38)/2') + 
    ylab('Difference Normalized Heritability Contribution (mappedToMm10 - hg38)') + 
    theme_bw(base_size = font_fig) + guides(color = guide_legend(nrow = 2)) + 
    theme(legend.position = "bottom", legend.text=element_text(size=font_fig-2),
          legend.title=element_text(size=font_fig-2))
  print(pp)
}
dev.off()


