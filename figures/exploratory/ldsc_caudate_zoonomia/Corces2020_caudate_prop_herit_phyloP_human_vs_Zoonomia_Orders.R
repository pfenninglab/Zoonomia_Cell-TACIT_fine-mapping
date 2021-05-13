#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F, bitmapType='cairo')
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(ggtree))
suppressMessages(library(ggtreeExtra))
library(rcartocolor)
library(ggrepel)
library(ggforce)
library(ggpubr)
library(ggsci)
library(here)
library(rmeta)

# to be run in the root github directory
LABEL='Zoonomia_data'
PROJDIR=here('data/raw_data/ldsc_caudate_zoonomia')
SETDIR=here('figures/exploratory/ldsc_caudate_zoonomia')

##########################
# read in the GWAS traits
load(here('data/tidy_data/ldsc_gwas','rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% select( -file) %>% 
  mutate(label = ss(as.character(trait), '_'))

##########################
# read in the Zoonomia tree and species list
rda_fn = here('data/tidy_data/Zoonomia_data', 
              'rdas','200_Mammals_Genome_Information.rda')
load(file = rda_fn)
col_clade = df %>% select(Clade, col_clade)%>% filter(!duplicated(Clade)) %>% deframe()
col_order = df %>% select(Order, col_order) %>% filter(!duplicated(Order)) %>% deframe()

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
  ungroup() %>% select(-(file:model_species)) %>%
  select(-(Observed_scale_h2_min:ind2))

hgEnrich$celltype %>% table()
hgEnrich$trait %>% droplevels() %>% table()

######################################################
# read in the LDSC partitioned heritability estimation
enrich_fn =here(PROJDIR) %>% 
  list.files(path = ., pattern = '.agg.gz', full.names = T, recursive = T)
names(enrich_fn) = ss(basename(enrich_fn), '.agg.gz')
input = lapply(enrich_fn, fread) %>% 
  rbindlist(fill = T, idcol='file') %>% select(file:Coefficient_SE)
input %>% data.frame() %>% head()

#########################################
## format groupings and calculate conditional cell type enrichment p-value
enrichments = input %>% 
  mutate(
    match = ss(file, '\\.', 2), 
    Species = ss(Categories, '\\.', 4),
    Species = ifelse(is.na(Species), ss(Categories, 'PredActiveIn', 2), Species), 
    Species = ifelse(is.na(Species), 'Homo_sapiens', Species), 
    plotGroup = case_when(
      grepl('phyloP',file) ~  'mappable + phyloP',
      grepl('predActive', file) ~ 'mappable + predActive',
      TRUE ~ 'mappable'
    ), 
    plotGroup = factor(plotGroup, c('mappable', 'mappable + phyloP', 'mappable + predActive')),
    celltype = ss(file, '\\.', 4) %>% 
      factor(c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb',  'Astro', 
               'Microglia', 'OPC', 'Oligo')), 
    cell_group = case_when(
      grepl('MSN|INT', celltype) ~ 'Neuron', 
      TRUE ~ 'Glia'
    )) %>% 
  filter(! grepl('PredActiveInHomo_sapiens', Categories)) %>%
  # GWAS  signif in humans
  inner_join(x = hgEnrich, by = c('match', 'celltype')) %>%
  # zoonomia species data
  inner_join(x = df, by = 'Species') %>%
  group_by(file) %>% type_convert()

enrichments$celltype %>% table()
table(enrichments$Species, enrichments$plotGroup)

### meta-analyze traits across orders
enrichments_meta = enrichments %>%
  mutate(Enrichment_SE = Observed_scale_h2_SE / h2 / Proportion_of_SNPs) %>%
  group_by(plotGroup, celltype, trait, group, Order) %>%
  group_split() %>%
  lapply(function(df){
    # meta-analyze the enrichments
    tmp = meta.summaries(df$Enrichment, df$Enrichment_SE, method="random")
    tmp2 = meta.summaries(df$Coefficients, df$Coefficient_SE, method="random")
    tmp3 = meta.summaries(df$Observed_scale_h2, df$Observed_scale_h2_SE, method="random")
    
    tab = data.frame( 
      Enrichment = tmp$summary, Enrichment_SE = tmp$se.summary,
      Coefficients = tmp2$summary, Coefficient_SE = tmp2$se.summary,
      Observed_scale_h2 = tmp3$summary, Observed_scale_h2_SE = tmp3$se.summary, 
      Time.Since.Split.from.Human.TimeTree.median = mean(df$Time.Since.Split.from.Human.TimeTree.median))
    
    ret = df[1,] %>% select(c(match:label, plotGroup,celltype, 
                              trait, group, Order, Proportion_of_SNPs)) %>%
      bind_cols(tab)
    return(ret)
  }) %>% rbindlist()


df_meta = enrichments_meta %>% 
  mutate(
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
  Enrichment_max = Proportion_of_h2g_max / Proportion_of_SNPs)


## normalize the coefficients by per SNP heritability
# compute Padj w/ bonferroni family mutliple hypothesis correction
alpha = 0.05;
df_meta = df_meta %>% group_by(match, celltype) %>% 
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

to_label = c('Homo_sapiens', 'Mus_musculus', 'Macaca_mulatta')

table(df_meta$Order, df_meta$plotGroup )


#################################
## save enrichments table to RDS
dir.create(here(SETDIR,'rdas'), showWarnings = FALSE)
save_fn = here(SETDIR,'rdas','caudate_conservation_prop_herit_phyloP_Zoonomia_OrdersMeta.rds')
saveRDS(df_meta, file = save_fn)

#################################
## make plots for presentation ##
dir.create(here(SETDIR,'plots'))
dir.create(here(SETDIR,'plots','zoonomia_caudate_Orders'))
height_ppt = 4; width_ppt = 8;
height_fig = 1.75; width_fig = 2.25; font_fig = 7
plot_traits = sort(unique(enrichments$trait))

# make plots
for(cell in unique(enrichments$celltype)){
plot_fn = here(SETDIR,'plots',paste('zoonomia_caudate_Orders', 
                    cell, 'prop_herit_enrichments_Orders.ppt.pdf', sep = '.'))
pdf(width = width_ppt, height = height_ppt, file = plot_fn, onefile = T)
for(lab in plot_traits){
  ## subset data to cell type and GWAS
  tmp = df_meta %>% filter(celltype %in% cell, trait %in% lab)
  if(nrow(tmp) > 0){
  ## the SNP enrichments ##
  pp1 = ggplot(data = tmp, aes(x = Time.Since.Split.from.Human.TimeTree.median, y = Enrichment)) +
    geom_smooth(aes(color = plotGroup), method = "lm", formula = y ~ x, size = 1, se = TRUE) +  
    geom_errorbar(aes(ymin = Enrichment_min, ymax = Enrichment_max, 
                      fill = Order, shape = plotGroup),
                  width = 0.1, position=position_dodge(width=2)) + 
    geom_jitter(aes(fill = Order, shape = plotGroup), position=position_dodge(width=2)) +
    scale_fill_manual(values = col_order) + 
    scale_shape_manual(values = 21:23) + 
    scale_linetype_manual(name = "Zoonomia", values = c(1, 1), 
    guide = guide_legend(title.position="top", nrow =2,override.aes = 
                          list(size = 1,color = c("#882255", "#888888")))) + 
    theme_bw(base_size = 10) + ggtitle(paste(cell, lab, sep = ', ')) +
    xlab(paste0('MY from Human (TimeTree median est.)')) + ylab(paste0('Heritability Enrichment')) + 
    guides(fill = guide_legend(nrow =3, title.position="top", override.aes = list(size = 3, pch = 21)),
           color = guide_legend(nrow =3, title.position="top", override.aes = list(size = 3)),
           shape = guide_legend(nrow =3, title.position="top", override.aes = list(size = 3))) + 
    theme(legend.position = "bottom", legend.text=element_text(size=font_fig),
          legend.title=element_text(size=font_fig), legend.key.height=unit(.5,"line"))
  
  ## the proportion of SNP heritability ##
  pp2 = ggplot(data = tmp, aes(x = Time.Since.Split.from.Human.TimeTree.median, 
                   y = Proportion_of_h2g)) +
    geom_smooth(aes(color = plotGroup), method = "lm", formula = y ~ x + I(x^2), size = 1, se = TRUE) +  
    geom_errorbar(aes(ymin = Proportion_of_h2g_min, ymax = Proportion_of_h2g_max, 
                      fill = Order, shape = plotGroup),
                  width = 0.1, position=position_dodge(width=2)) + 
    geom_jitter(aes(fill = Order, shape = plotGroup), position=position_dodge(width=2)) +
    scale_fill_manual(values = col_order) + 
    scale_shape_manual(values = 21:23) + 
    scale_linetype_manual(name = "Zoonomia", values = c(1, 1), 
    guide = guide_legend(title.position="top", nrow =2,override.aes = 
                          list(size = 1,color = c("#882255", "#888888")))) + 
    theme_bw(base_size = 10) + ggtitle(paste(cell, lab, sep = ', ')) +
    xlab(paste0('MY from Human (TimeTree median est.)')) + ylab(paste0('Proportion SNP Heritability')) + 
    guides(fill = guide_legend(nrow =3, title.position="top", override.aes = list(size = 3, pch = 21)),
           color = guide_legend(nrow =3, title.position="top", override.aes = list(size = 3)),
           shape = guide_legend(nrow =3, title.position="top", override.aes = list(size = 3))) + 
    theme(legend.position = "bottom", legend.text=element_text(size=font_fig),
          legend.title=element_text(size=font_fig), legend.key.height=unit(.5,"line"))
  
  ## Plot the proportion of SNP heritability ##
  plot_fn1 = here(SETDIR,'plots', 'zoonomia_caudate_Orders',
                       paste('zoonomia_caudate_Orders',lab, cell,
                             'prop_herit_enrichments.ppt.pdf', sep = '.'))
  pdf(width = width_ppt, height = height_ppt, file = plot_fn1, onefile = F)
  pp = ggarrange(pp1, pp2, nrow = 1, common.legend = TRUE, 
                 align = 'h', legend="bottom")
  print(pp)
  dev.off()
  print(pp)
  }}
dev.off()
}



