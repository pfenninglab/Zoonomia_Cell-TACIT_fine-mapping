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
library(ggh4x)
library(ggsci)
library(here)
library(metafor)
library(tidymeta)
library(ggplot2)
library(dplyr)
library(broom)
library(swfdr)
library(qvalue)

#########################################
# 0) to be run in the root github directory
LABEL='ldsc_atac_and_conservation'
PLOTDIR=file.path('figures/exploratory',LABEL)
DATADIR=here('data/raw_data/',LABEL)

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


## load in the saved enrichments table to RDS
dir.create(here(PLOTDIR,'rdas'), showWarnings = F)
save_fn = here(PLOTDIR,'rdas','zoo_meta_atac_and_conservation_prop_heritability_Corces2020.rds')
enrichments = readRDS(file = save_fn) %>% 
  dplyr::rename('TimeSplit' = 'Time.Since.Split.from.Human.TimeTree.median')

##########################################################################################
## 2) compute meta-analyses of the heritability enrichment identifying effect of Time split
alpha = 0.05
cols = c( '91mam_GERP', '200mam_PhyloP', '43prim_PhastCons','predActive')
rename= c('91 Mammals GERP', '240 Mammals PhyloP', '43 Primate PhastCons', 
          'Cell-TACIT') %>% setNames(cols)

## compute meta analyses on the time-split coefficient, across traits
enrichments_group = enrichments %>% 
  nest(data = -c(group, celltype, cell_group.x, peaktype, constype)) %>% 
  mutate(meta = map(data, ~ rma(method = 'REML', slab = file, data = .x, 
                                weights = h2_Z, mods = ~ TimeSplit + match,
                                yi = Enrichment, 
                                sei = Enrichment_max - Enrichment)),
         tidied = map(meta, tidy)) %>% 
  unnest(c(tidied)) %>% dplyr::select(-c(data)) %>% 
  filter(term == 'TimeSplit') %>% 
  mutate(fdr = qvalue(p.value, fdr.level = alpha)$qvalues) %>% 
  dplyr::rename('herit_enrich_per_MYA' = 'estimate', 
                'herit_enrich_per_MYA_SE' = 'std.error')

## look at the meta analyses regression
enrichments_group %>% filter(celltype =='MSN_D2') %>% filter(constype =='predActive') %>% pull(meta)

save_tbl_fn = here(PLOTDIR,'tables','ldsc_atac_and_conservation_gain_enrichment_group_TimeSplit.xlsx')
enrichments_group %>% dplyr::select(-meta) %>% writexl::write_xlsx(save_tbl_fn)
save_rds_fn = here(PLOTDIR,'rdas','ldsc_atac_and_conservation_gain_enrichment_group_TimeSplit.rds')
enrichments_group %>% dplyr::select(-meta) %>% saveRDS(save_rds_fn)


#################################
## 3) make plots for presentation
dir.create(here(PLOTDIR,'plots'))
height_ppt = 4; width_ppt = 8;


enrichments_group2 = enrichments_group %>%
  mutate(constype = rename[constype], constype = factor(constype, rename)) %>% 
  mutate(`Cons + ATAC` = constype, 
         cell_group = factor(cell_group.x, c('Neuron', 'Glia')))

theme_ppt =  theme_bw(base_size = 10) + 
  theme(legend.position = "top", axis.title.x = element_blank(), 
        legend.key.height=unit(.5,"line"), legend.key.width=unit(.5,"line"), 
        axis.text.x=element_text(angle = 30, hjust = 1))

plot_fn = here(PLOTDIR,'plots', 'meta_group_herit_enrichment_per_mya.ppt.pdf')
pdf(width = width_ppt, height = height_ppt, file = plot_fn)
ggplot(enrichments_group2, aes(x = celltype, y = herit_enrich_per_MYA)) + 
  geom_errorbar(aes(fill = `Cons + ATAC`, shape = `Cons + ATAC`, alpha = fdr < 0.05, 
                    ymin = pmax(0, herit_enrich_per_MYA - herit_enrich_per_MYA_SE), 
                    ymax = herit_enrich_per_MYA + herit_enrich_per_MYA_SE), 
                width = .3) + 
  geom_point(aes(fill = `Cons + ATAC`, shape = `Cons + ATAC`, alpha = fdr < 0.05)) +
  geom_hline(yintercept = 0, color = 'black') +
  scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1") + 
  scale_shape_manual(values = c(21:24)) +
  facet_nested(peaktype ~ cell_group + group, scales = 'free') + 
  ylab(paste0('Heritability Enrichment per MYA')) + theme_ppt
dev.off()


plot_fn = here(PLOTDIR,'plots', 'meta_group_herit_enrichment_per_mya.Neuron.ppt.pdf')
pdf(width = width_ppt, height = height_ppt, file = plot_fn)
ggplot(enrichments_group2 %>% filter(cell_group == 'Neuron'), 
       aes(x = celltype, y = herit_enrich_per_MYA)) + 
  geom_errorbar(aes(fill = `Cons + ATAC`, shape = `Cons + ATAC`, alpha = fdr < 0.05, 
                    ymin = pmax(0, herit_enrich_per_MYA - herit_enrich_per_MYA_SE), 
                    ymax = herit_enrich_per_MYA + herit_enrich_per_MYA_SE), 
                width = .3) + 
  geom_point(size = 2, aes(fill = `Cons + ATAC`, shape = `Cons + ATAC`, alpha = fdr < 0.05)) +
  geom_hline(yintercept = 0, color = 'black') +
  scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1") + 
  scale_shape_manual(values = c(21:24)) +
  facet_nested(peaktype ~ group, scales = 'free') + 
  ylab(paste0('Heritability Enrichment per MYA')) + theme_ppt
dev.off()


plot_fn = here(PLOTDIR,'plots', 'meta_group_herit_enrichment_per_mya.Glia.ppt.pdf')
pdf(width = width_ppt, height = height_ppt, file = plot_fn)
ggplot(enrichments_group2 %>% filter(cell_group == 'Glia'), 
       aes(x = celltype, y = herit_enrich_per_MYA)) + 
  geom_errorbar(aes(fill = `Cons + ATAC`, shape = `Cons + ATAC`, alpha = fdr < 0.05, 
                    ymin = pmax(0, herit_enrich_per_MYA - herit_enrich_per_MYA_SE), 
                    ymax = herit_enrich_per_MYA + herit_enrich_per_MYA_SE), 
                width = .3) + 
  geom_point(size = 2, aes(fill = `Cons + ATAC`, shape = `Cons + ATAC`, alpha = fdr < 0.05)) +
  geom_hline(yintercept = 0, color = 'black') +
  scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1") + 
  scale_shape_manual(values = c(21:24)) +
  facet_nested(peaktype ~ group, scales = 'free') + 
  ylab(paste0('Heritability Enrichment per MYA')) + theme_ppt
dev.off()


#################################
## 4) make plots for figures
in2mm = 25.4
height_fig = 2; width_fig = 3.464567; font_fig = 5

theme_fig =  theme_bw(base_size = 5) + 
  theme(legend.position = "top", legend.title=element_text(size=3),
        legend.text=element_text(size=3),
        axis.title.x = element_blank(), legend.box.spacing = unit(.05,"line"),
        legend.key.height=unit(.05,"line"), legend.key.width=unit(.05,"line"), 
        axis.text.x=element_text(angle = 45, hjust = 1))


plot_fn = here(PLOTDIR,'plots', 'meta_group_herit_enrichment_per_mya.Neuron.fig.pdf')
pdf(width = width_fig, height = height_fig, file = plot_fn)
ggplot(enrichments_group2 %>% filter(cell_group == 'Neuron'), 
       aes(x = celltype, y = herit_enrich_per_MYA)) + 
  geom_hline(yintercept = 0, color = 'black', size = .15) +
  geom_errorbar(size = .4, aes(fill = `Cons + ATAC`, shape = `Cons + ATAC`, alpha = fdr < 0.05, 
                              ymin = pmax(0, herit_enrich_per_MYA - herit_enrich_per_MYA_SE), 
                              ymax = herit_enrich_per_MYA + herit_enrich_per_MYA_SE), 
                width = .3) + 
  geom_point(size = .75, stroke = .4, 
             aes(fill = `Cons + ATAC`, shape = `Cons + ATAC`, alpha = fdr < 0.05)) +
  scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1") + 
  scale_shape_manual(values = c(21:24)) + 
  facet_nested(peaktype ~ group, scales = 'free') + 
  ylab(paste0('Heritability Enrichment per MYA')) + theme_fig +
  guides(shape = guide_legend(nrow =1, title.position="top", override.aes = list(size = 1)))
  
dev.off()


plot_fn = here(PLOTDIR,'plots', 'meta_group_herit_enrichment_per_mya.Glia.fig.pdf')
pdf(width = width_fig, height = height_fig, file = plot_fn)
ggplot(enrichments_group2 %>% filter(cell_group == 'Glia'), 
       aes(x = celltype, y = herit_enrich_per_MYA)) + 
  geom_hline(yintercept = 0, color = 'black', size = .15) +
  geom_errorbar(size = .4, aes(fill = `Cons + ATAC`, shape = `Cons + ATAC`, alpha = fdr < 0.05, 
                    ymin = pmax(0, herit_enrich_per_MYA - herit_enrich_per_MYA_SE), 
                    ymax = herit_enrich_per_MYA + herit_enrich_per_MYA_SE), 
                width = .3) + 
  geom_point(size = .75, stroke = .4, 
             aes(fill = `Cons + ATAC`, shape = `Cons + ATAC`, alpha = fdr < 0.05)) +
  scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1") + 
  scale_shape_manual(values = c(21:24)) + 
  facet_nested(peaktype ~ group, scales = 'free') + 
  ylab(paste0('Heritability Enrichment per MYA')) + theme_fig
  
dev.off()



###############################################################
## 5) compute meta-analyses per trait using fixed effect model

## compute meta analyses on the time-split coefficient per trait
enrichments_match = enrichments %>% 
  nest(data = -c(match, celltype, cell_group.x, peaktype, constype)) %>%
  mutate(meta = map(data, ~ rma(method = 'EE', slab = file, data = .x, 
                                weights = h2_Z, mods = ~ TimeSplit,
                                yi = Enrichment, 
                                sei = Enrichment_max - Enrichment)),
         h2_Z = map_dbl(data, ~ mean(.x$h2_Z)),
         tidied = map(meta, tidy)) %>% dplyr::select(-c(data)) %>% 
  unnest(c(tidied)) %>% filter(term == 'TimeSplit') %>% 
  mutate(fdr = lm_qvalue(p.value, X=h2_Z)$q) %>% 
  dplyr::rename('herit_enrich_per_MYA' = 'estimate', 
                'herit_enrich_per_MYA_SE' = 'std.error')

save_tbl_fn = here(PLOTDIR,'tables','ldsc_atac_and_conservation_gain_enrichment_trait_TimeSplit.xlsx')
enrichments_match %>% dplyr::select(-meta) %>% writexl::write_xlsx(save_tbl_fn)
save_rds_fn = here(PLOTDIR,'rdas','ldsc_atac_and_conservation_gain_enrichment_trait_TimeSplit.rds')
enrichments_match %>% dplyr::select(-meta) %>% saveRDS(save_rds_fn)



enrichments_match2 = enrichments_match %>%
  mutate(constype = rename[constype], constype = factor(constype, rename)) %>% 
  mutate(`Cons + ATAC` = constype, 
         celltype = plyr::revalue(celltype,c("MSN_SN"="MSN_D1D2H")),
         cell_group = factor(cell_group.x, c('Neuron', 'Glia')), 
         peaktype = factor(peaktype, c('Enhancer', 'Promoter', 'Other'))) %>% 
  arrange(peaktype, constype)



############################################################
## 6) plot the trait-cell type concordance matrix per trait
enrichments_match3 = enrichments_match2 %>% 
  arrange(desc(p.value)) %>% inner_join(pheno) %>% 
  group_by(celltype, match, peaktype) %>% 
  mutate(ind = which(constype == 'Cell-TACIT'), 
         keep = fdr[ind] < 0.05 & sum(fdr <0.05) >=2) %>% 
  ungroup() %>% filter(keep) %>% dplyr::select(-c(ind, keep)) %>% 
  mutate(match = factor(match, unique(as.character(match))))


plot_fn = here(PLOTDIR,'plots', 'meta_trait_signif_heatmap.Enhancer.fig.pdf')
pdf(width = width_fig, height = height_fig, file = plot_fn)
ggplot(enrichments_match3 %>% filter(peaktype == 'Enhancer'), 
       aes(y = ss(as.character(trait), '_'), x =constype, fill = fdr < 0.05)) + 
  geom_tile(aes(alpha = fdr < 0.05)) +
  scale_fill_manual(values = c('black', '#00000000')) +
  scale_y_discrete( limits = rev) +
  scale_x_discrete( position = "top") +
  facet_nested( group~ celltype , scales = 'free', space = 'free') + 
  theme_fig + theme(axis.text.x=element_text(angle = -45, hjust = 1), 
                    axis.title = element_blank())
dev.off()



plot_fn = here(PLOTDIR,'plots', 'meta_trait_signif_heatmap.Promoter.fig.pdf')
pdf(width = width_fig, height = height_fig/3*2, file = plot_fn)
ggplot(enrichments_match3 %>% filter(peaktype == 'Promoter'), 
       aes(y = ss(as.character(trait), '_'), x =constype, fill = fdr < 0.05)) + 
  geom_tile(aes(alpha = fdr < 0.05)) +
  scale_fill_manual(values = c('black', '#00000000')) +
  scale_y_discrete( limits = rev) +
  scale_x_discrete( position = "top") +
  facet_nested( group~ celltype , scales = 'free', space = 'free') + 
  theme_fig + theme(axis.text.x=element_text(angle = -45, hjust = 1), 
                    axis.title = element_blank())
dev.off()


## look at some interesting cell types and traits
enrichments_match2 %>% filter(grepl('^Intell',match) & celltype == 'MSN_D1') %>% 
  dplyr::select(peaktype, constype, herit_enrich_per_MYA) %>% 
  pivot_wider(names_from = 'peaktype', values_from = 'herit_enrich_per_MYA')


enrichments_match2 %>% filter(grepl('^Intell',match) & celltype == 'MSN_D1') %>% 
  filter(peaktype =='Enhancer') %>% dplyr::select(-meta) %>% data.frame()


