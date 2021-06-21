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

#########################################
# to be run in the root github directory
LABEL='ldsc_zoonomia_meta'
PLOTDIR=here('figures/exploratory/', LABEL)
DATADIR=here('data/raw_data/',LABEL)

##########################
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

####################################################
# read in the enrichments from phyloP scores cons/accl
phyloP_fn = here('figures/exploratory/ldsc_conservation','rdas',
                 'caudate_conservation_prop_herit_hg_rm_mm.rds')
hgEnrich = readRDS(file = phyloP_fn) %>% filter(model_species =='hg38') %>%
  arrange(celltype, peaktype) %>%
  group_by(trait) %>%
  mutate(
    ind1 = which(peaktype =='PhyloP.accl'), 
    ind2 = which(peaktype =='PhyloP.cons'),
    phyloPacc.Enr     = Enrichment[ind1],
    phyloPacc.Enr_min = Enrichment_min[ind1],
    phyloPacc.Enr_max = Enrichment_max[ind1],
    phyloPcon.Enr     = Enrichment[ind2],
    phyloPcon.Enr_min = Enrichment_min[ind2],
    phyloPcon.Enr_max = Enrichment_max[ind2],
    phyloPacc.Prop     = Proportion_of_h2g[ind1],
    phyloPacc.Prop_min = Proportion_of_h2g_min[ind1],
    phyloPacc.Prop_max = Proportion_of_h2g_max[ind1],
    phyloPcon.Prop     = Proportion_of_h2g[ind2],
    phyloPcon.Prop_min = Proportion_of_h2g_min[ind2],
    phyloPcon.Prop_max = Proportion_of_h2g_max[ind2]) %>% 
  filter(peaktype %in% c('PhyloP.cons', 'PhyloP.accl')) %>%
  ungroup() %>% select(-(file:model_species)) %>% 
  select(-(Observed_scale_h2_min:ind2))

hgEnrich$peaktype %>% table()
hgEnrich$trait %>% droplevels() %>% table()


######################################################
# read in the LDSC partitioned heritability estimation
enrich_fn =here(DATADIR,'prop_herit_M1ctx_liver') %>% 
  list.files(path = ., pattern = '.agg.gz', full.names = T)
names(enrich_fn) = ss(basename(enrich_fn), '.agg.gz')
input = lapply(enrich_fn, fread) %>% 
  rbindlist(fill = T, idcol='file') %>% select(file:Coefficient_SE)
input %>% data.frame() %>% head()

#########################################
## format groupings and calculate conditional cell type enrichment p-value
enrichments = input %>% 
  mutate(
    peaktype = case_when(grepl('mappable',Categories) ~ 'mappable', TRUE ~ 'predActive'),
    group_meta = Categories %>% ss('enhPeaks.',2) %>% 
      gsub(pattern = '\\.mappable|\\.predActive', replacement = ''),
    match = ss(file, '\\.', 3), 
    celltype = ss(file, '\\.', 5) %>% factor(c('M1ctx', 'Liver'))) %>% 
  # GWAS  signif in humans
  inner_join(x = hgEnrich %>% select(-celltype), by = c('match')) %>%
  # zoonomia group_meta data
  distinct(peaktype, celltype, group_meta, match, .keep_all = TRUE) %>%
  inner_join(x = df_meta, by = 'group_meta') %>%
  group_by(file) %>% type_convert()

enrichments %>% pull(celltype) %>% table()
enrichments %>% pull(group_meta) %>% table()

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
enrichments %>% count(peaktype)

## normalize the coefficients by per SNP heritability
# compute Padj w/ bonferroni family mutliple hypothesis correction
alpha = 0.05;
enrichments = enrichments %>% group_by(match, celltype) %>% 
  mutate(
    # normalize conditional coefficients
    Coef_norm = Coefficients / h2_perSNP,
    Coef_norm_se = Coefficient_SE / h2_perSNP, 
    # calculate 2-sided p-value for non-zero effect
    Coef_z = Coef_norm / Coef_norm_se,
    Coefficient_P_value = pnorm(- Coef_z ),
    Padj = p.adjust(Coefficient_P_value, 'fdr'), 
    logPadj = -log10(Padj),
    p.signif = ifelse(Padj < alpha, paste('FDR <',alpha), 'NS'),
    p.signif = factor(p.signif, levels = c('NS', paste('FDR <',alpha)))
  ) %>% ungroup() %>% filter(complete.cases(Coefficient_P_value))

to_label = c('Primates#0', 'Primates#28.81', 'Rodentia#89')


#################################
## save enrichments table to RDS
dir.create(here(PLOTDIR,'rdas'), showWarnings = F)
save_fn = here(PLOTDIR,'rdas','zoonomia_meta_prop_heritability_M1ctx_liver.rds')
saveRDS(enrichments, file = save_fn)

save_excel = here(PLOTDIR,'tables','zoonomia_meta_prop_heritability_M1ctx_liver.xlsx')
enrichments %>% writexl::write_xlsx(save_excel)

#################################
## make plots for presentation ##
dir.create(here(PLOTDIR,'plots'), showWarnings = F)
dir.create(here(PLOTDIR,'plots','prop_herit_M1ctx_liver'), showWarnings = F)
height_ppt = 4; width_ppt = 8;
height_fig = 1.75; width_fig = 2.25; font_fig = 7
plot_traits = sort(unique(enrichments$trait))
cell = 'M1ctx'; lab = 'LDL_E'


###############
#### make plots
for(cell in unique(enrichments$celltype)) {
  plot_fn = here(PLOTDIR,'plots',paste('zoonomia_bulk_ATAseq', cell, 'prop_herit_enrichments.ppt.pdf', sep = '_'))
  pdf(width = width_ppt, height = height_ppt, file = plot_fn, onefile = T)
  for(lab in plot_traits) {
    tmp = enrichments %>% filter(celltype %in% cell, trait %in% lab)
    if(nrow(tmp) > 0) {
      pp1 = ggplot(data = tmp, aes(x = Time.Since.Split.from.Human.TimeTree.median, y = Enrichment)) +
        geom_hline(aes(yintercept = mean(tmp$phyloPcon.Enr), 
                       linetype = "phyloP.Con"),color = "#882255") + 
        geom_hline(aes(yintercept = mean(tmp$phyloPacc.Enr), 
                       linetype = "phyloP.Acc"),color = "#888888")  +
        geom_smooth(aes(color = peaktype), method = "lm", formula = y ~ x + I(x^2), size = 1, se = TRUE) +  
        geom_errorbar(aes(fill = Order2, shape = peaktype, ymin = Enrichment_min, ymax = Enrichment_max),
                      width = 0.1, position=position_dodge(width=2)) + 
        geom_jitter(aes(fill = Order2, shape = peaktype), position=position_dodge(width=2)) +
        scale_fill_manual(values = col_meta) + 
        scale_shape_manual(values = c(21:23)) + 
        scale_linetype_manual(name = "Zoonomia", values = c(1, 1), 
                              guide = guide_legend(title.position="top", nrow =2,override.aes = 
                                                     list(size = 1,color = c("#882255", "#888888")))) + 
        theme_bw(base_size = 10) + ggtitle(paste(cell, lab, sep = ', ')) +
        xlab(paste0('MY from Human (TimeTree median est.)')) + ylab(paste0('Heritability Enrichment')) + 
        guides(fill = guide_legend(nrow =2, title.position="top", override.aes = list(pch = 21 ,size = 4)),
               shape = guide_legend(nrow =2, title.position="top", override.aes = list(size = 4))) + 
        theme(legend.position = "bottom", legend.text=element_text(size=font_fig),
              legend.title=element_text(size=font_fig), legend.key.height=unit(.5,"line"))
      
      pp2 = ggplot(data = tmp, aes(x = Time.Since.Split.from.Human.TimeTree.median, 
                                   y = Proportion_of_h2g)) +
        geom_hline(aes(yintercept = mean(tmp$phyloPcon.Prop), 
                       linetype = "phyloP.Con"),color = "#882255") + 
        geom_hline(aes(yintercept= mean(tmp$phyloPacc.Prop), 
                       linetype = "phyloP.Acc"),color = "#888888")  +
        geom_smooth(aes(color = peaktype), method = "lm", formula = y ~ x + I(x^2), size = 1, se = TRUE) +  
        geom_errorbar(aes(fill = Order2, shape = peaktype, ymin = Proportion_of_h2g_min, ymax = Proportion_of_h2g_max),
                      width = 0.1, position=position_dodge(width=2)) + 
        geom_jitter(aes(fill = Order2, shape = peaktype), position=position_dodge(width=2)) +
        scale_fill_manual(values = col_meta) + 
        scale_shape_manual(values = c(21:23)) + 
        scale_linetype_manual(name = "Zoonomia", values = c(1, 1), 
                              guide = guide_legend(title.position="top", nrow =2,override.aes = 
                                                     list(size = 1,color = c("#882255", "#888888")))) + 
        theme_bw(base_size = 10) + ggtitle(paste(cell, lab, sep = ', ')) +
        xlab(paste0('MY from Human (TimeTree median est.)')) + ylab(paste0('Proportion SNP Heritability')) + 
        guides(fill = guide_legend(nrow =2, title.position="top", override.aes = list(pch = 21 ,size = 4)),
               shape = guide_legend(nrow =2, title.position="top", override.aes = list(size = 4))) + 
        theme(legend.position = "bottom", legend.text=element_text(size=font_fig),
              legend.title=element_text(size=font_fig), legend.key.height=unit(.5,"line"))
      
      #### 
      plot_fn1 = here(PLOTDIR,'plots', 'prop_herit_M1ctx_liver', 
                      paste('zoonomia_bulk_ATAseq',lab, cell, 'prop_herit_enrichments.ppt.pdf', sep = '.'))
      pdf(width = width_ppt, height = height_ppt, file = plot_fn1, onefile = F)
      pp = ggarrange(pp1, pp2, nrow = 1, common.legend = TRUE, 
                     align = 'h', legend="bottom")
      print(pp)
      dev.off()
      print(pp)
    }
  }
  dev.off()
}


