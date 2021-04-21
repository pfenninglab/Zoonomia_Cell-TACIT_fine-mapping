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
library(here)

# to be run in the root github directory
LABEL='caudate_conservation_ldsc'
PROJDIR= 'figures/exploratory/ldsc_conservation'
DATADIR=here('data/raw_data/',LABEL)

##########################
# read in the GWAS traits
load(here('data/tidy_data/ldsc_gwas','rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% select( -file) %>% mutate(label = ss(as.character(trait), '_'))

#########################################
# read in the LDSC partitioned heritability estimation
enrich_fn =here(DATADIR,'enrichments') %>% 
  list.files(path = ., pattern = '.cell_type_results.txt', full.names = T)
names(enrich_fn) = ss(basename(enrich_fn), '.cell_type_results.txt')
input = lapply(enrich_fn, read_tsv, col_types = cols()) %>% bind_rows(.id = 'file')

#########################################
## format groupings and calculate conditional cell type enrichment p-value
enrichments = input %>% 
  filter(!grepl('CP.MSN_D|CP.INT', Name)) %>%
  mutate(
    file = gsub('caudate_conservation_', '', file), 
    annot_group = ss(file, '\\.', 1),
    match = ss(file, '\\.', 2), 
    tmpcol = ss(Name,'\\.', 1),
    peaktype = case_when(
      tmpcol == 'Corces2020_caudate' ~ 'Human peak',
      grepl('Corces2020_caudate_mappedTo', tmpcol) ~ 'Hg -> Model Org', 
      grepl('Orth',tmpcol)  ~ 'OCR Ortholog', 
      TRUE ~ 'Model Org -> Hg',
    ),
    model_species = case_when(
      tmpcol == 'Corces2020_caudate' ~ 'hg38',
      grepl('Rm|Rhe|Stauffer',tmpcol) ~ 'rheMac10',
      TRUE ~ 'mm10'
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
  filter(celltype != 'Caudate') %>%
  select(-tmpcol) %>% inner_join(x = pheno, by = 'match')
enrichments %>% data.frame() %>% head()
enrichments %>% pull(peaktype) %>% table()
enrichments %>% pull(model_species) %>% table()

## normalize the coefficients by per SNP heritability
# compute Padj within a particular peak_group (hg38 or inModelOrg) 
# and peaktype (full, core), these features are overlapping, nested
enrichments = enrichments %>% 
  mutate(Padj = p.adjust(Coefficient_P_value, 'fdr'), 
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
  group_by(celltype, match, annot_group) %>% 
  mutate(
    Padj_inModel = Padj,
    Coef_norm_inModel = Coef_norm,
    ind = which(peaktype == 'Human peak'),
    Padj_hgPeak = Padj[ind],
    Coef_norm_hgPeak = Coef_norm[ind],
    ) %>% ungroup() %>% select(-ind) %>% 
  filter(peaktype != 'Human peak') %>%
  mutate(
    p.signif = case_when(
      Padj_hgPeak < alpha & Padj_inModel < alpha ~ 'both',
      Padj_hgPeak < alpha ~ 'inHg38', 
      Padj_inModel < alpha ~ 'inModelOrg',
      TRUE ~ 'NS'
    ),
    peaktype = factor(peaktype, levels = 
                       rev(c('OCR Ortholog', 'Hg -> Model Org','Model Org -> Hg'))),
    norm_coeff_diff = Coef_norm_inModel - Coef_norm_hgPeak, 
    norm_coeff_mean = (Coef_norm_inModel + Coef_norm_hgPeak)/2, 
    label = ifelse(p.signif !='NS' & annot_group=='phyloP score + peak' |
                     p.signif =='both',label,NA))

# look at the top values
enrich_wide %>% data.frame() %>% head(2)
enrich_wide%>% pull(p.signif) %>% table()
# both     inHg38 inModelOrg         NS 
#  577        257         96       2142 

#################################
## save enrichments table to RDS
save_fn = here('figures/exploratory/ldsc_conservation','rdas',
               'caudate_conservation_ldsc_coeff_hg_rm_mm.rds')
saveRDS(enrichments, file = save_fn)
save_fn2 = here('figures/exploratory/ldsc_conservation','rdas',
               'caudate_conservation_ldsc_coeff-wide_hg_rm_mm.rds')
saveRDS(enrich_wide, file = save_fn2)

to_map = enrichments %>% arrange(celltype) %>% group_by(match) %>% 
  filter(Padj < alpha) %>% summarise(celltypes = paste(unique(celltype), collapse = ','))
save_fn3 = here(DATADIR,'rdas','caudate_celltype_gwas_enriched.tsv')
write_tsv(to_map, file = save_fn3)


#################################
## make plots for presentation ##
system(paste('mkdir -p', here( PROJDIR, 'plots')))
height_ppt = 4.5; width_ppt = 8
height_fig = 3.5; width_fig = 4.75; font_fig = 5
plot_celltypes = enrichments %>% filter(Padj < alpha) %>% pull(celltype) %>% unique() 

# make plots
for(annot in unique(enrich_wide$annot_group)){
  for(cell in c("Neuron", "Glia")){
    for(species in c('rheMac10','mm10')){
      plot_fn = here(PROJDIR, 'plots',paste('Corces2020_caudate_ldsc',
                                annot,species,cell,'ppt.pdf', sep = '_'))
      pdf(width = width_ppt, height = height_ppt, file = plot_fn)
      pp = ggplot(data = enrich_wide %>% filter(cell_group == cell)  %>%
                    arrange(peaktype) %>%
                    mutate(peaktype = gsub('Model Org',species,as.character(peaktype)),
                           peaktype = factor(peaktype, rev(unique(peaktype)))) %>%
                filter(annot_group == annot) %>% filter(model_species == species), 
              aes(y = norm_coeff_diff, x = norm_coeff_mean, 
                  fill = group, color = p.signif)) +
      geom_hline(yintercept = 0, color = 'darkgrey') + 
      geom_vline(xintercept = 0, color = 'darkgrey') + 
      geom_point(pch = 21, aes(alpha = p.signif != 'NS'), size =1.5) + 
      scale_color_manual(values = c('black', 'darkred', 'blue','grey50'), 
                       breaks = c("both", 'inHg38', 'inModelOrg', 'NS'),
                       name = paste('FDR <',alpha)) + 
      scale_fill_manual(values = group_col, name = 'Trait Group') + 
      # scale_shape_manual(values = c(21, 22, 24), name = '', 
      #                    breaks = c('hgToMm', 'mmToHg', 'hgMmOrth')) + 
      scale_alpha_manual(values = c(.2, 1), guide = 'none') + 
      coord_cartesian(ylim=c(-6, 10), xlim=c(0, 12)) +
      facet_grid( peaktype ~  celltype, scales = 'fixed') +  
      geom_label_repel(aes(label = label), box.padding = .1, label.size = .01,
                   max.overlaps = 40, size = 1.75, show.legend = F,na.rm = T,
                   point.padding = .1, segment.color = 'grey50', max.time = 2,
                   min.segment.length = .1, alpha = .4, 
                   label.padding = .1, force_pull = 1, force = 40) +
      xlab(paste0('Normalized Heritability Coefficient Average (',species,' Peak + hg38 Peak)/2')) + 
      ylab(paste0('Norm. Coef. Diff. (',species,' Peak - hg38 Peak)')) + 
      theme_bw(base_size = 11) + 
      guides(fill = guide_legend(override.aes = list(
             shape = 22, color = 'white', size = 3))) + 
      theme(legend.position = "right", 
          legend.text=element_text(size=8),
          legend.title=element_text(size=9))
      print(pp)
      dev.off()
    }}}


# make plots
for(annot in unique(enrich_wide$annot_group)){
  for(cell in c("Neuron", "Glia")){
    plot_fn2 = here(PROJDIR, 'plots',paste('Corces2020_caudate_ldsc',
                                       annot,cell,'fig.pdf', sep = '_'))
    pdf(width = width_fig, height = height_fig, file = plot_fn2, onefile = F)
    ppList = lapply(c('rheMac10','mm10'), function(species){
      pp = ggplot(data = enrich_wide %>% filter(cell_group == cell) %>%
                    arrange(peaktype) %>%
                    mutate(peaktype = gsub('Model Org',species,as.character(peaktype)),
                           peaktype = factor(peaktype, unique(peaktype))) %>%
                    filter(annot_group == annot) %>% filter(model_species == species),
                  aes(y = norm_coeff_diff, x = norm_coeff_mean, 
                      fill = group, color = p.signif)) +
        geom_hline(yintercept = 0, color = 'darkgrey') + 
        geom_vline(xintercept = 0, color = 'darkgrey') + 
        geom_point(pch = 21, stroke =.3, aes(alpha = p.signif != 'NS'), size =.5) + 
        scale_color_manual(values = c('black', 'darkred', 'blue','grey50'), 
                           breaks = c("both", 'inHg38', 'inModelOrg', 'NS'),
                           name = paste('FDR <',alpha)) + 
        scale_fill_manual(values = group_col, name = 'Trait Group') + 
        scale_alpha_manual(values = c(.2, 1), guide = 'none') + 
        coord_cartesian(ylim=c(-6, 10), xlim=c(0, 12)) +
        facet_grid(celltype ~ peaktype, scales = 'fixed') +  
        geom_label_repel(aes(label = label), box.padding = .1, label.size = .01,
                         max.overlaps = 40, size = .8, show.legend = F,na.rm = T,
                         point.padding = .1, segment.color = 'grey50', max.time = 2,
                         min.segment.length = .1, alpha = .5, 
                         label.padding = .1, force_pull = 1, force = 40) +
        xlab(paste0('Norm Herit. Coef. Avg. (',species,' Peak + hg38 Peak)/2')) + 
        ylab(paste0('Normalized Contribution Coefficient Diff. (',species,' Peak - hg38 Peak)')) + 
        theme_bw(base_size = font_fig) +     
        guides(colour = guide_legend(nrow = 2, title.position="left",override.aes = list(size = 1)), 
               fill = guide_legend(nrow =2, title.position="left", override.aes = list(size = 1))) + 
        theme(legend.position = "bottom", legend.text=element_text(size=font_fig),
              legend.title=element_text(size=font_fig),
              legend.key.height=unit(.5,"line"))
      return(pp)
    })
    pp = ggarrange(ppList[[1]], ppList[[2]], ncol = 2, 
                   common.legend = TRUE, align = 'h', legend="bottom")
    print(pp)
    dev.off()
  }}


