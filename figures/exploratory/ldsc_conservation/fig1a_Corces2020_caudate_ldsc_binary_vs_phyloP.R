#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F)
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
library(ggrepel)
library(ggpubr)
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
  list.files(path = ., pattern = '.cell_type_results.txt', full.names = T)
names(enrich_fn) = ss(basename(enrich_fn), '.cell_type_results.txt')
input = lapply(enrich_fn, read_tsv) %>% bind_rows(.id = 'file')

#########################################
## format groupings and calculate conditional cell type enrichment p-value
enrichments = input %>% 
  filter(!grepl('CP.MSN_D|CP.INT', Name)) %>%
  mutate(
    file = gsub('Corces2020_', '', file), 
    annot_group = ss(file, '\\.', 1),
    match = ss(file, '\\.', 2), 
    tmpcol = ss(Name,'\\.', 1),
    peaktype = case_when(
      tmpcol == 'Corces2020_caudate' ~ 'Human peak',
      tmpcol == 'Corces2020_caudate_mappedToMm10' ~ 'Hg peak -> Mm', 
      tmpcol == 'Corces2020_caudate_hgMmOrth' ~ 'HgMm Ortholog', 
      TRUE ~ 'Mm peak -> Hg',
    ),
    peak_group = case_when(
      peaktype == 'hgPeak' ~ 'hgPeak' ,
      TRUE ~ 'inMm'
    ),
    annot_group = case_when(
      grepl('binary', annot_group) ~ 'binary annotation', 
      grepl('phyloP', annot_group) ~ 'phyloP score + peak', 
    ),
    celltype = ss(Name,'\\.', 2),
    celltype = gsub('Consensus', 'Caudate',celltype) %>% 
      factor(c('Caudate', 'MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb', 'Astro', 
               'Microglia', 'OPC', 'Oligo')), 
    cell_group = case_when(
      celltype =='Caudate' ~ 'Neuron', 
      grepl('MSN|INT', celltype) ~ 'Neuron', 
      TRUE ~ 'Glia'
    )) %>%
  select(-tmpcol) %>%
  inner_join(x = pheno, by = 'match') %>%
  filter(celltype != 'stromal_b')
enrichments %>% data.frame() %>% head()
enrichments %>% pull(peaktype) %>% table()

## normalize the coefficients by per SNP heritability
# compute Padj within a particular peak_group (hg38 or inMm) 
# and peaktype (full, core), these features are overlapping, nested
enrichments = enrichments %>% #group_by(peaktype) %>%
  mutate(Padj = p.adjust(Coefficient_P_value, 'bonferroni'), 
         logPadj = -log10(Padj), 
         Coef_norm = Coefficient / h2_perSNP, 
         Coef_norm_se = Coefficient_std_error / h2_perSNP) %>% 
  ungroup()

# look at top signif celltypes
enrichments %>% filter(Padj < .1) %>% select(match, celltype, annot_group, Padj) %>%
  group_by(match) %>% top_n(1, -Padj) %>% ungroup() %>% 
  data.frame() %>% arrange(Padj) 
enrichments %>% data.frame() %>% head(2)

###################################################################
## pivot table to pair up peak_group with each celltype and trait  ##
alpha = 0.05; 
enrich_wide = enrichments %>% 
  # these id columns make up all combinations aside form peakgroup
  group_by(celltype, annot_group, match) %>% 
  mutate(
    Padj_inMm = Padj,
    Coef_norm_inMm = Coef_norm,
    ind = which(peaktype == 'Human peak'),
    Padj_hgPeak = Padj[ind],
    Coef_norm_hgPeak = Coef_norm[ind],
    ) %>% ungroup() %>% select(-ind) %>% 
  filter(peaktype != 'Human peak') %>%
  mutate(
    p.signif = case_when(
      Padj_hgPeak < alpha & Padj_inMm < alpha ~ 'both',
      Padj_hgPeak < alpha ~ 'hg38', 
      Padj_inMm < alpha ~ 'inMm',
      TRUE ~ 'NS'
    ),
    peaktype = factor(peaktype, levels = 
                        c('HgMm Ortholog', 'Hg peak -> Mm','Mm peak -> Hg')),
    norm_coeff_diff = Coef_norm_inMm - Coef_norm_hgPeak, 
    norm_coeff_mean = (Coef_norm_inMm + Coef_norm_hgPeak)/2, 
    label = ifelse(p.signif !='NS' & annot_group=='phyloP score + peak' |
                     p.signif =='both',label,NA))

# look at the top values
enrich_wide %>% data.frame() %>% head(2)
enrich_wide%>% pull(p.signif) %>% table()
# both hg38 inMm 
# 128  163    7 

#################################
## make plots for presentation ##
system(paste('mkdir -p', file.path( 'plots')))
height_ppt = 4.5; width_ppt = 8
height_fig = 5; width_fig = 2.25; font_fig = 5
plot_celltypes = enrichments %>% filter(Padj < alpha) %>% pull(celltype) %>% unique() 

# make plots
plot_fn = file.path('plots','Corces2020_caudate_ldsc_hum_v_hum2mouse.ppt.pdf')
pdf(width = width_ppt, height = height_ppt, file = plot_fn)
for(annot in c('binary annotation', 'phyloP score + peak')){
for(cell in c("Neuron", "Glia")){
  pp = ggplot(data = enrich_wide %>% filter(cell_group == cell) %>%
                filter(annot_group == annot), 
              aes(y = norm_coeff_diff, x = norm_coeff_mean, 
                  fill = group, color = p.signif)) +
    geom_hline(yintercept = 0, color = 'darkgrey') + 
    geom_vline(xintercept = 0, color = 'darkgrey') + 
    geom_point(pch = 21, aes(alpha = p.signif != 'NS'), size =1.5) + 
    scale_color_manual(values = c('black', 'darkred', 'blue','grey50'), 
                       breaks = c("both", 'hg38', 'inMm', 'NS'),
                       name = paste('P_bonf <',alpha)) + 
    scale_fill_manual(values = group_col, name = 'Trait Group') + 
    # scale_shape_manual(values = c(21, 22, 24), name = '', 
    #                    breaks = c('hgToMm', 'mmToHg', 'hgMmOrth')) + 
    scale_alpha_manual(values = c(.2, 1), guide = 'none') + 
    coord_cartesian(ylim=c(-4, 10), xlim=c(0, 8)) +
    facet_grid(peaktype ~ celltype, scales = 'fixed') +  
    geom_label_repel(aes(label = label), box.padding = .1, label.size = .01,
                   max.overlaps = 40, size = 2, show.legend = F,na.rm = T,
                   point.padding = .1, segment.color = 'grey50', max.time = 2,
                   min.segment.length = .1, alpha = .4, 
                   label.padding = .1, force_pull = 1, force = 40) +
    xlab('Normalized Heritability Coefficient Average (Orth + hgPeak)/2') + 
    ylab('Norm. Coef. Diff. (Orth - hgPeak)') + 
    theme_bw(base_size = 11) + 
    guides(fill = guide_legend(override.aes = list(
             shape = 22, color = 'white', size = 3))) + 
    theme(legend.position = "right", 
          legend.text=element_text(size=8),
          legend.title=element_text(size=9))
    print(pp)
  }}
dev.off()


# make plots
plot_fn2 = file.path('plots','Corces2020_caudate_ldsc_hum_v_hum2mouse.fig.pdf')
pdf(width = width_fig, height = height_fig, file = plot_fn2)
for(annot in c('binary annotation', 'phyloP score + peak')){
  for(cell in c("Neuron", "Glia")){
  pp = ggplot(data = enrich_wide %>% filter(cell_group == cell) %>% 
                filter(annot_group == annot), 
              aes(y = norm_coeff_diff, x = norm_coeff_mean, 
                  fill = group, color = p.signif)) +
    geom_hline(yintercept = 0, color = 'darkgrey') + 
    geom_vline(xintercept = 0, color = 'darkgrey') + 
    geom_point(pch = 21, stroke =.3, aes(alpha = p.signif != 'NS'), size =.5) + 
    scale_color_manual(values = c('black', 'darkred', 'blue','grey50'), 
                       breaks = c("both", 'hg38', 'inMm', 'NS'),
                       name = paste('P_bonf <',alpha)) + 
    scale_fill_manual(values = group_col, name = 'Trait Group') + 
    scale_alpha_manual(values = c(.2, 1), guide = 'none') + 
    coord_cartesian(ylim=c(-4, 10), xlim=c(0, 8)) +
    facet_grid(celltype ~ peaktype, scales = 'fixed') +  
    geom_label_repel(aes(label = label), box.padding = .1, label.size = .01,
                     max.overlaps = 40, size = 1, show.legend = F,na.rm = T,
                     point.padding = .1, segment.color = 'grey50', max.time = 2,
                     min.segment.length = .1, alpha = .5, 
                     label.padding = .1, force_pull = 1, force = 40) +
    xlab('Norm Herit. Coef. Avg. (Orth + hgPeak)/2') + 
    ylab('Normalized Contribution Coefficient Diff. (Orth - hgPeak)') + 
    theme_bw(base_size = font_fig) +     
    guides(colour = guide_legend(ncol = 1, title.position="top",override.aes = list(size = 1)), 
           fill = guide_legend(ncol = 2, title.position="top", override.aes = list(size = 1))) + 
    theme(legend.position = "bottom", legend.text=element_text(size=font_fig),
          legend.title=element_text(size=font_fig),
          legend.key.height=unit(.5,"line"))
  print(pp)
}}
dev.off()


