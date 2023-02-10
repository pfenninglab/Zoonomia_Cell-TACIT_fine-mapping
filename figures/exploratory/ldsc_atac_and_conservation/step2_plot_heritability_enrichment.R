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
LABEL='ldsc_atac_and_conservation'
PLOTDIR=here('figures/exploratory/', 'ldsc_atac_and_conservation')
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
  ungroup() %>% select(-(file:model_species)) %>%
  select(-(Observed_scale_h2_min:ind2))

hgEnrich$celltype %>% table()
hgEnrich$trait %>% droplevels() %>% table()


######################################################
# read in the LDSC partitioned heritability estimation
enrich_fn =here(DATADIR,'prop_herit') %>% 
  list.files(path = ., pattern = '.agg.gz', full.names = T)
names(enrich_fn) = ss(basename(enrich_fn), '.agg.gz')

# exclude older enrichments
enrich_fn = enrich_fn[!grepl('Corces', ss(names(enrich_fn), '\\.', 2))]
input = lapply(enrich_fn, fread) %>% 
  rbindlist(fill = T, idcol='file') %>% select(file:Coefficient_SE)
input %>% data.frame() %>% head()

celltypes = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb',  
              'Astro', 'Microglia', 'OPC', 'Oligo')

cols = c( '91mam_GERP', '200mam_PhyloP', '43prim_PhastCons','predActive')
rename= c('91 Mammals GERP', '240 Mammals PhyloP', '43 Primate PhastCons', 
          'Cell-TACIT') %>% setNames(cols)

#########################################
## format groupings and calculate conditional cell type enrichment p-value
enrichments = input %>% 
  mutate(
    match = ss(file, '\\.', 4),
    constype = case_when(grepl('91mam_GERP',Categories) ~ '91mam_GERP', 
                         grepl('200mam_PhyloP',Categories) ~ '200mam_PhyloP', 
                         grepl('43prim_PhastCons',Categories) ~ '43prim_PhastCons', 
                         grepl('predActive',Categories) ~ 'predActive'),
    peaktype = ss(file, '\\.', 2) %>% factor(c('Enhancer', 'Promoter', 'Other')),
    celltype = ss(Categories, '\\.', 2) %>% factor(celltypes), 
    cell_group = case_when( grepl('MSN|INT', celltype) ~ 'Neuron', TRUE ~ 'Glia'), 
    Categories = Categories %>% gsub(paste(unique(match), collapse = '|'), '', .),
    Categories = Categories %>% gsub(paste(unique(peaktype), collapse = '|'), '', .),
    Categories = Categories %>% gsub(paste(unique(constype), collapse = '|'), '', .),
    Categories = Categories %>% gsub(paste(celltypes, collapse = '|'), '', .),
    group_meta = Categories %>% ss('\\.', 3)
 ) %>% dplyr::select(-Categories) %>% 
  # GWAS  signif in humans
  inner_join(x = hgEnrich, by = c('match', 'celltype')) %>%
  # zoonomia group_meta data
  inner_join(x = df_meta %>% 
               mutate(group_meta = ss(as.character(group_meta), '\\.', 1) %>% factor), 
             by = 'group_meta') %>%
  group_by(file) %>% type_convert()%>%
  mutate(constype = rename[constype], constype = factor(constype, rename),
         celltype = plyr::revalue(celltype,c("MSN_SN"="MSN_D1.D2H")),
         cell_group = factor(cell_group.x, c('Neuron', 'Glia')))

enrichments$celltype %>% table()
enrichments$group_meta %>% table()
enrichments$constype %>% table()
enrichments$peaktype %>% table()

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
save_fn = here(PLOTDIR,'rdas','zoo_meta_atac_and_conservation_prop_heritability_Corces2020.rds')
saveRDS(enrichments, file = save_fn)

#################################
## make plots for presentation ##
dir.create(here(PLOTDIR,'plots'))
dir.create(here(PLOTDIR,'plots','prop_herit'))
height_ppt = 4; width_ppt = 8;
height_fig = 1.75; width_fig = 2.25; font_fig = 7
plot_traits = sort(unique(enrichments$trait))

enrichments = enrichments %>% mutate(`Cons + ATAC` = constype)

###############
#### make plots
for(cell in unique(enrichments$celltype)){
  plot_fn = here(PLOTDIR,'plots',paste('zoonomia_caudate', 
                                       cell, 'prop_herit_enrichments.ppt.pdf', sep = '_'))
  pdf(width = width_ppt, height = 6, file = plot_fn, onefile = T)
  for(lab in plot_traits){
    ## subset data to cell type and GWAS
    tmp = enrichments %>% filter(celltype %in% cell, trait %in% lab)
    if(nrow(tmp) > 0){
      ## the SNP enrichments ##
      pp1 = ggplot(data = tmp, aes(x = Time.Since.Split.from.Human.TimeTree.median + 1, 
                                   y = Enrichment)) +
        # add phyloP enrichments lines
        geom_hline(aes(yintercept = mean(tmp$phyloPcon.Enr), 
                       linetype = "phyloP.Con"),color = "#882255") + 
        geom_hline(aes(yintercept = mean(tmp$phyloPacc.Enr), 
                       linetype = "phyloP.Acc"),color = "#888888")  +
        # add line fit
        geom_smooth(aes(color = `Cons + ATAC`), method = "lm", 
                    formula = y ~ x, size = 1, se = TRUE) +  
        geom_errorbar(aes(fill = `Cons + ATAC`, shape = `Cons + ATAC`, 
                          ymin = Enrichment_min, ymax = Enrichment_max), width = 0.1) + 
        geom_point(aes(fill = `Cons + ATAC`, shape = `Cons + ATAC`)) +
        scale_fill_brewer(palette = "Set1") + scale_shape_manual(values = c(21:24)) + 
        scale_color_brewer(palette = "Set1") + scale_x_continuous(trans='log2') + 
        facet_wrap( ~ peaktype , scales = 'free_y') + 
        scale_linetype_manual(name = "200mam_PhyloP alone", values = c(1, 1), 
                              guide = guide_legend(nrow =1,override.aes = 
                                                     list(size = 1,color = c("#882255", "#888888")))) + 
        theme_bw(base_size = 9) + ggtitle(paste(cell, lab, sep = ', ')) +
        xlab(paste0('Million years from Human')) + ylab(paste0('Heritability Enrichment')) + 
        guides(shape = guide_legend(nrow =1, title.position="top", override.aes = list(size = 3))) + 
        theme(legend.position = "bottom", legend.text=element_text(size=font_fig),
              legend.title=element_text(size=font_fig), legend.key.height=unit(.5,"line"))
      
      ## the proportion of SNP heritability ##
      pp2 = ggplot(data = tmp, aes(x = Time.Since.Split.from.Human.TimeTree.median + 1, 
                                   y = Proportion_of_h2g)) +
        # add phyloP enrichments lines
        geom_hline(aes(yintercept = mean(tmp$phyloPcon.Prop), 
                       linetype = "phyloP.Con"),color = "#882255") + 
        geom_hline(aes(yintercept= mean(tmp$phyloPacc.Prop), 
                       linetype = "phyloP.Acc"),color = "#888888")  +
        # add line fit
        geom_smooth(aes(color = `Cons + ATAC`), method = "lm", 
                    formula = y ~ x, size = 1, se = TRUE) +  
        geom_errorbar(aes(fill = `Cons + ATAC`, shape = `Cons + ATAC`, 
                          ymin = Proportion_of_h2g_min, ymax = Proportion_of_h2g_max),
                      width = 0.1) + 
        geom_point(aes(fill = `Cons + ATAC`, shape = `Cons + ATAC`)) +
        scale_fill_brewer(palette = "Set1") + scale_shape_manual(values = c(21:24)) + 
        scale_color_brewer(palette = "Set1") + scale_x_continuous(trans='log2') + 
        facet_wrap( ~ peaktype , scales = 'free_y') + 
        scale_linetype_manual(name = "200mam_PhyloP alone", values = c(1, 1), 
                              guide = guide_legend(nrow =1,override.aes = 
                                                     list(size = 1,color = c("#882255", "#888888")))) + 
        theme_bw(base_size = 9) +
        xlab(paste0('Million years from Human')) + ylab(paste0('Proportion SNP Heritability')) + 
        guides(shape = guide_legend(title= 'Conservation + ATAC', nrow =1, 
                                    title.position="top", override.aes = list(size = 3))) + 
        theme(legend.position = "bottom", legend.text=element_text(size=font_fig),
              legend.title=element_text(size=font_fig), legend.key.height=unit(.5,"line"))
      
      ## Plot the proportion of SNP heritability ##
      plot_fn1 = here(PLOTDIR,'plots', 'prop_herit',
                      paste('zoonomia_caudate',lab, cell, 'prop_herit_enrichments.ppt.pdf', sep = '_'))
      pdf(width = width_ppt, height = 6, file = plot_fn1, onefile = F)
      pp = ggarrange(pp1, pp2, nrow = 2, common.legend = TRUE, 
                     align = 'v', legend="bottom")
      print(pp)
      dev.off()
      print(pp)
    }}
  dev.off()
}



