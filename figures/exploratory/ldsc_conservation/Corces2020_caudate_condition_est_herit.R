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
library(ggforce)
# to be run in the root github directory
LABEL='caudate_conservation_ldsc'
setwd('figures/exploratory/ldsc_conservation')
PROJDIR=file.path('../../../data/raw_data/',LABEL)

##########################
# read in the GWAS traits
load(file.path('../../../data/tidy_data/ldsc_gwas','rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% select( -file) %>% 
  mutate(label = ss(as.character(trait), '_'))

#########################################
# read in the LDSC partitioned heritability estimation
enrich_fn =file.path(PROJDIR,'est_herit_cond') %>% 
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
      tmpcol == 'Corces2020_caudate' ~ 'Human peak',
      grepl('Corces2020_caudate_mappedTo', tmpcol) ~ 'Hg -> Model Org', 
      grepl('Orth',tmpcol)  ~ 'Enhancer Ortholog', 
      TRUE ~ 'Model Org -> Hg',
    ),
    model_species = case_when(
      tmpcol == 'Corces2020_caudate' ~ 'hg38',
      grepl('Rm|Rhe|Stauffer',tmpcol) ~ 'rheMac10',
      TRUE ~ 'mm10'
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
  group_by(file) %>% type_convert() %>%
  filter(peaktype == 'Human peak') %>%
  filter(!grepl('CP.MSN_D|CP.INT', Categories)) %>%
  filter(grepl('Corces2020', Categories))

enrichments$file %>% ss('\\.') %>% table()
enrichments$Categories %>% ss('\\.',2) %>% table()
enrichments$celltype %>% table()
enrichments$model_species %>% ss('\\.') %>% table()

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
alpha = 0.05;
enrichments = enrichments %>% group_by(annot_group, match) %>% 
  mutate(
  # normalize conditional coefficients
  Coef_norm = Coefficients / h2_perSNP,
  Coef_norm_se = Coefficient_SE / h2_perSNP, 
  # calculate 2-sided p-value for non-zero effect
  Coef_z = Coef_norm / Coef_norm_se,
  Coefficient_P_value = pnorm(- Coef_z ),
  Padj = p.adjust(Coefficient_P_value, 'bonferroni'), 
  logPadj = -log10(Padj),
  p.signif = Padj < alpha
) %>% ungroup() %>% filter(complete.cases(Coefficient_P_value))

enrichments = enrichments %>%
  filter(Categories != 'Corces2020_caudate.All_Roadmap_DHS_cCREs_BICCN_Stauffer.BG')

# look at top signif celltypes
enrichments %>% filter(Padj < 0.05) %>% 
  select(match, celltype, annot_group, Padj) %>%
  group_by(match) %>% top_n(1, -Padj) %>% ungroup() %>% 
  data.frame() %>% arrange(Padj)

#################################
## make plots for presentation ##
system(paste('mkdir -p', file.path( 'plots')))
height_ppt = 5; width_ppt = 8;
height_fig = 1.75; width_fig = 4.75; font_fig = 7
# plot_celltypes = enrichments %>% filter(Padj < alpha) %>% pull(celltype) %>% unique() 
plot_traits = enrichments %>% filter(p.signif) %>% pull(trait) %>% unique() 

# make plots
plot_fn = file.path('plots','Corces2020_caudate_conditional_herit_enrichments.ppt.pdf')
pdf(width = width_ppt, height = height_ppt, file = plot_fn)
  pp = ggplot(data = enrichments %>% filter(trait %in% plot_traits), 
              aes(y = logPadj, x = celltype, fill = celltype, color = p.signif)) +
    geom_bar(aes(alpha = p.signif != 'NS'), stat = 'identity') + 
    scale_color_manual(values = c(NA, 'black'), guide = 'none') + 
    scale_fill_carto_d(name = "Cell type:", palette = "Safe") +
    scale_alpha_manual(values = c(1,.2), guide = 'none') + 
    facet_wrap( ~ trait, scales = 'fixed', nrow =4 ) +  
    geom_hline(yintercept = -log10(alpha), color = 'red',linetype = 'dashed') + 
    xlab('Enriched cell type | other cell types & annotations') + 
    ylab('-log10(Bonferroni P-value)') + 
    theme_bw(base_size = 10) + coord_flip() + 
    guides(fill = guide_legend(override.aes = list(shape = 21), nrow = 1)) + 
    theme(legend.position = "bottom", 
          legend.text=element_text(size=7),
          legend.title=element_text(size=8))
  print(pp)
dev.off()


# make plots
plot_fn2 = file.path('plots','Corces2020_caudate_conditional_herit_enrichments.fig.pdf')
pdf(width = width_fig, height = height_fig, file = plot_fn2)
pp = ggplot(data = enrichments %>% filter(trait %in% plot_traits), 
            aes(y = logPadj, x = label, fill = celltype, color = p.signif)) +
  geom_hline(yintercept = -log10(alpha), color = 'red',linetype = 'dashed') + 
  geom_point(aes(alpha = p.signif, shape = cell_group)) + 
  scale_color_manual(values = c('grey', 'black'), guide = 'none') + 
  scale_fill_carto_d(name = "", palette = "Safe") +
  scale_alpha_manual(values = c(.5, 1), guide = 'none') + 
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



