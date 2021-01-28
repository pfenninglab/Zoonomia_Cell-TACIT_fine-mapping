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

# to be run in the root github directory
LABEL='meuleman_dhs_index'
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

## format groupings and calculate conditional cell type enrichment p-value
enrichments = input %>% 
  mutate(
    file = gsub('Meuleman_DHS_', '', file), 
    annot_group = ss(file, '\\.', 1),
    match = ss(file, '\\.', 2), 
    tmpcol = ss(Name,'\\.', 1),
    peak_group = case_when(
      grepl('mappedToMm10', tmpcol) ~ 'mappableToMm10', 
      TRUE ~ 'hg38'
    ),
    annot_group = case_when(
      annot_group == 'binary' ~ 'binary peak annotation', 
      annot_group == 'phyloP' ~ 'phyloP score + peak', 
    ),
    tissue = ss(Name,'\\.', 2),
    peaktype = ss(Name,'\\.', 3)) %>%
  select(-tmpcol) %>%
  inner_join(x = pheno, by = 'match') %>%
  filter(tissue != 'stromal_b')

## normalize the coefficients by per SNP heritability
# compute Padj within a particular peak_group (hg38 or mappableToMm10) 
# and peaktype (full, core), these features are overlapping, nested
enrichments = enrichments %>% group_by(peaktype) %>%
  mutate(Padj = p.adjust(Coefficient_P_value, 'bonferroni'), 
         logPadj = -log10(Padj), 
         Coef_norm = Coefficient / h2_perSNP, 
         Coef_norm_se = Coefficient_std_error / h2_perSNP) %>% 
  ungroup()

# look at top signif tissues
enrichments %>% filter(Padj < .1) %>% select(match, tissue, Padj) %>%
  group_by(match) %>% top_n(1, -Padj) %>% ungroup() %>% 
  data.frame() %>% arrange(Padj) 

enrichments %>% data.frame() %>% head(15)

###################################################################
## pivot table to pair up peak_group with each tissue and trait  ##
alpha = .05; 
enrich_wide = enrichments %>% 
  # these id columns make up all combinations aside form peakgroup
  pivot_wider(id_col = c(eval(names(pheno)), tissue, peaktype, annot_group), 
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
     label = ifelse(p.signif != 'NS',label,NA))

enrich_wide %>% data.frame() %>% head(2)


#################################
## make plots for presentation ##
system(paste('mkdir -p', file.path( 'plots')))
height_ppt = 5; width_ppt = 8
plot_tissues = enrichments %>% filter(Padj < alpha) %>% pull(tissue) %>% unique() 

# make plots
plot_fn = file.path('plots','meuleman_dhs_index_ldsc_hum_v_hum2mouse.ppt.pdf')
pdf(width = width_ppt, height = height_ppt, file = plot_fn)
for(peak in c('full', 'core')){
  pp = ggplot(data = enrich_wide %>% filter(tissue %in% plot_tissues & peaktype ==peak), 
         aes(y = norm_coeff_diff, x = norm_coeff_mean, fill = group, color = p.signif)) +
    geom_hline(yintercept = 0, color = 'black') + 
    geom_vline(xintercept = 0, color = 'black') + 
    geom_point(pch = 21, aes(alpha = p.signif != 'NS'), size =2) + 
    scale_color_manual(values = c('black', 'darkred', 'blue','grey50'), 
                       breaks = c('both', 'hg38','mappableToMm10','NS'),
                       name = paste('P_bonf <',alpha)) + 
    scale_fill_manual(values = group_col, name = 'GWAS Trait') + 
    scale_alpha_manual(values = c(.2, 1), guide = 'none') + 
    coord_cartesian(ylim=c(-2, 5), xlim=c(0, 7)) +
    facet_grid(annot_group~tissue, scales = 'fixed') +  
    geom_label_repel(aes(label = label), box.padding = .25, label.size = .01,
                     max.overlaps = 30, size = 1.5, show.legend = F,na.rm = T,
                     point.padding = .1, segment.color = 'grey50', max.time = 2,
                     min.segment.length = .1, alpha = .7, 
                     label.padding = .1, force_pull = 1, force = 2) +
    xlab('Avg. Herit. Normalized Contribution (mappedToMm10 + hg38)/2') + 
    ylab('Norm. Contrib. Diff. (mappedToMm10 - hg38)') + 
    theme_bw(base_size = 12) +
    guides(colour = guide_legend(nrow = 2),  fill = guide_legend(nrow = 2)) + 
    theme(legend.position = "bottom", 
          legend.text=element_text(size=8), legend.title=element_text(size=9)) + 
    ggtitle(paste('Meuleman DHS Index',peak,'peak'))
  print(pp)
}
dev.off()



# make plots for figure, tall
height_fig = 4; width_fig = 4.75; font_fig = 7
plot_fn2 = file.path('plots','meuleman_dhs_index_ldsc_hum_v_hum2mouse.fig.pdf')
pdf(width = width_fig, height = height_fig, file = plot_fn2)
  pp = ggplot(data = enrich_wide %>% filter(tissue %in% plot_tissues), 
              aes(y = norm_coeff_diff, x = norm_coeff_mean,
                  fill = group, color = p.signif)) +
    geom_hline(yintercept = 0, color = 'black') + 
    geom_vline(xintercept = 0, color = 'black') + 
    geom_point(pch = 21, aes(alpha = p.signif != 'NS'), size =2) + 
    scale_color_manual(values = c('black', 'darkred', 'blue','grey50'), 
                       breaks = c('both', 'hg38','mappableToMm10','NS'),
                       name = paste('P_bonf <',alpha)) + 
    scale_fill_manual(values = group_col, name = 'GWAS Trait') + 
    scale_alpha_manual(values = c(.2, 1), guide = 'none') + 
    coord_cartesian(ylim=c(-2, 7), xlim=c(0, 7)) +
    facet_grid(tissue~peaktype + annot_group, scales = 'fixed') +  
    geom_label_repel(aes(label = label), box.padding = .25, label.size = .01,
                     max.overlaps = 40, size = 1, show.legend = F,na.rm = T,
                     point.padding = .1, segment.color = 'grey50', max.time = 2,
                     min.segment.length = .1,alpha = .7, 
                     label.padding = .1, force_pull = 1, force = 2) +
    xlab('Avg. Normalized Herit. Contribution (mappedToMm10 + hg38)/2') + 
    ylab('Norm. Contrib. Diff. (mappedToMm10 - hg38)') + 
    theme_bw(base_size = font_fig) + guides(color = guide_legend(nrow = 2)) + 
    theme(legend.position = "bottom", legend.text=element_text(size=font_fig-2),
          legend.title=element_text(size=font_fig-2)) + 
    guides(colour = guide_legend(nrow = 2),  fill = guide_legend(nrow = 2)) 
  print(pp)
dev.off()


