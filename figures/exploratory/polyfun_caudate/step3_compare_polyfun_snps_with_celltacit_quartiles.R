ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(here)
library(qvalue)
library(tidyverse)
library(tidymodels)
library(broom)
library(MASS)
library(rcartocolor)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(GenomicRanges)
library(ChIPseeker)
library(rtracklayer)

CODEDIR='code/raw_code/polyfun_caudate'
DATADIR='data/raw_data/polyfun_caudate'
PLOTDIR='figures/exploratory/polyfun_caudate'

suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))

##########################
# read in the GWAS traits
load(here('data/tidy_data/ldsc_gwas','rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% dplyr::select( -file) %>% 
  mutate(label = ss(as.character(trait), '_'))

celltypes = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb',  'Astro', 
  'Microglia', 'OPC', 'Oligo')

#####################
## 1) read in table of fine-mapped SNPs ##
dir.create(here('data/raw_data/polyfun_caudate/rdas'), showWarnings = F)
poly_fn = here('data/raw_data/polyfun_caudate/rdas',
               'polyfun_caudate_finemapped_snps_20210518.rds')
snps_brain_df = readRDS(file = poly_fn) %>%
  mutate(SNP = paste(CHR,POS_hg38,SNP,A1,A2, sep = ':')) %>%
  filter(PIP > 0.1)

## read in the fine-mapped SNPs from Flagship I
exclude_traits = c('BMI', 'Intelligence score', 'Neuroticism','Morning Person')
snps_fsii_df = readRDS(file = here('../zoonomia_finemapping/data/tidy_data/polyfun/polyfun_finemapped_snps_zoonomia_20220517.rds')) %>%
  mutate(SNP = paste(CHR,POS_hg38,SNP,A1,A2, sep = ':')) %>%
  ## only keep the final fine-mapping, exclude brain-related traits
  filter(group == 'base + ZooAnnot + cCRE', !TRAIT %in% exclude_traits ) %>%
  mutate(group = 'Other', TRAIT = droplevels(TRAIT)) %>%
  dplyr::rename('label' = 'TRAIT') %>%
  filter(PIP > 0.1)

## combine the two set of fine-mapped SNPs
common_names = intersect(names(snps_brain_df),names(snps_fsii_df))
snps_df = bind_rows(snps_brain_df %>% dplyr::select(all_of(common_names)), 
                    snps_fsii_df %>% dplyr::select(all_of(common_names))) %>%
  dplyr::select(label:PIP, group) %>% 
  filter(!duplicated(name))
table(snps_df$group)

## get the gnomic ranges position of the SNPs
snps_gr = snps_df %>% mutate(tmp = paste0('chr', CHR, ':', POS_hg38)) %>%
  dplyr::select(name, tmp) %>% deframe() %>% GRanges()


##########################################
## 2) read in the Cell TACIT Age quartiles
peak_fn = list.files('data/raw_data/ldsc_celltacit_age_decile/peaks', 
                     pattern = '.bed.gz', full.names = T) %>% str_subset('quartile')
peak_fn = peak_fn[!grepl('rand', peak_fn)]
names(peak_fn) = basename(peak_fn) %>% gsub('Corces2020.|CellTACIT.|.bed.gz', '', .)
peak_list = peak_fn %>% lapply(import)

## annotate each SNP with the overlapping peak
overlapMat = lapply(peak_list, countOverlaps, query = snps_gr) %>% 
  lapply(enframe) %>% rbindlist(idcol = 'peak') %>% 
  mutate(celltype = ss(peak, '\\.', 1), quartile = ss(peak, '\\.', 2)) %>%
  filter(value ==1) %>%
  pivot_wider(id_cols = c('name'), names_from = celltype, 
              values_from = quartile, values_fill = 'non-overlap')

## add the annotations to the SNP data frame
snps_annotated_df = inner_join(snps_df, overlapMat) %>%
  pivot_longer(cols = names(overlapMat)[-1], names_to = "celltype", values_to = "quantile")

snps_summary_df = snps_annotated_df %>% 
  group_by(group, celltype, quantile) %>%
  summarize(count = n()) %>%
  mutate(group = factor(group), group = relevel(group, ref = 'Other'), 
         quantile = factor(quantile), quantile = relevel(quantile, ref = 'non-overlap'))


################################
## 4) plot the log odds ratio
quartiles = c('0-25%', '25-50%', '50-75%', '75-100%')
names(quartiles) = paste0('quartile', 1:4); 
groups = names(group_col)[1:6]; names(groups) = names(group_col)[1:6]


snp_enrichment_df = lapply(groups, function(y){
  snps_annotated_df %>% 
    filter(group %in% c(y, 'Other')) %>%
    mutate(group = factor(group), group = relevel(group, ref = 'Other'), 
           quantile = factor(quantile), quantile = relevel(quantile, ref = 'non-overlap')) %>%
    nest(data = c(-celltype)) %>%
    mutate(    test = map(data, ~glm(group ~ quantile, data = .x,family=binomial())),
               tidied = map(test, tidy)
    ) %>% unnest(cols = tidied) %>% 
    dplyr::select(-data, -test)
}) %>% rbindlist(idcol = 'group') %>%
  filter(term != '(Intercept)') %>%
  mutate(OR = exp(estimate), FDR = p.adjust(p.value, 'fdr'), 
         OR_max = exp(estimate + std.error), 
         OR_min = exp(estimate - std.error), 
         term = gsub('quantile', '', term),
         celltype = factor(celltype, celltypes)) %>%
  dplyr::rename('quantile' = 'term') %>%
  arrange(log10(FDR) * OR ) %>% 
  mutate(group = factor(group), 
         quantile = factor(quantile), quantile = recode(quantile, !!!quartiles))

dir.create( here(PLOTDIR, 'tables'), showWarnings = F)
save_fn = here(PLOTDIR, 'tables', 'fine-mapped_SNP_enrichmentByGroup_quartiles_logistic_regression.xlsx')
snp_enrichment_df %>% writexl::write_xlsx(save_fn)


################################
## 4) plot the log odds ratio
plot_fn = here(PLOTDIR, 'plots', 'fine-mapped_SNP_logEnrichmentByGroup_quartiles_main.pdf')
pdf(plot_fn, height = 2.5, width = 5)
ggplot(snp_enrichment_df %>% filter(celltype %in% celltypes[c(1:3,6)]), 
       aes(x = quantile, y = estimate, fill = group, alpha = FDR < 0.05)) + 
  geom_hline(yintercept = 0, color = 'red', linetype = 'dashed') + 
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), 
                width=1, position=position_dodge(width=0.75)) + 
  geom_point(pch =21, position=position_dodge(width=0.75), size = 2) +
  scale_fill_manual(values = group_col[c(1:3,5:6)]) + 
  scale_alpha_manual(values = c(0.3, 1)) + 
  facet_wrap(~celltype, nrow = 1, scales = 'fixed') + 
  theme_bw(base_size = 7) + ylim(c(-3, 3.25)) + 
  xlab('Cell-TACIT Age Quantile') + 
  ylab('Log Enrichment Ratio\n(Neuropsych SNPs/ Other SNPs)') +
  guides(fill = guide_legend(override.aes = list(size = 1),nrow = 1, title=""),
         alpha = guide_legend(nrow = 1)) + 
  theme(legend.position = 'top',
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))+ 
  theme(legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-7,-5), 
        legend.spacing.x = unit(0.05, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'))
dev.off()


plot_fn2 = here(PLOTDIR, 'plots', 'fine-mapped_SNP_logEnrichmentByGroup_quartiles_other.pdf')
pdf(plot_fn2, height = 2.5, width = 5)
ggplot(snp_enrichment_df %>% filter(celltype %in% celltypes[c(4,5,7,8)]), 
       aes(x = quantile, y = estimate, fill = group, alpha = FDR < 0.05)) + 
  geom_hline(yintercept = 0, color = 'red', linetype = 'dashed') + 
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), 
                width=1, position=position_dodge(width=0.75)) + 
  geom_point(pch =21, position=position_dodge(width=0.75), size = 2) +
  scale_fill_manual(values = group_col[c(1:3,5:6)]) + 
  scale_alpha_manual(values = c(0.3, 1)) + 
  facet_wrap(~celltype, nrow = 1, scales = 'fixed') + 
  theme_bw(base_size = 7) + ylim(c(-2.3, 2.5)) + 
  xlab('Cell-TACIT Age Quantile') + 
  ylab('Log Enrichment Ratio\n(Neuropsych SNPs/ Other SNPs)') +
  guides(fill = guide_legend(override.aes = list(size = 1),nrow = 1, title=""),
         alpha = guide_legend(nrow = 1)) + 
  theme(legend.position = 'top',
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))+ 
  theme(legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-7,-5), 
        legend.spacing.x = unit(0.05, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'))
dev.off()


#################################
## 5) plot the regular odds ratio
plot_fn = here(PLOTDIR, 'plots', 'fine-mapped_SNP_enrichmentByGroup_quartiles_main.pdf')
pdf(plot_fn, height = 2.5, width = 4.25)
ggplot(snp_enrichment_df %>% filter(celltype %in% celltypes[c(1:3,6)], !is.infinite(OR_max)), 
       aes(x = quantile, y = OR, fill = group, alpha = FDR < 0.05)) + 
  geom_hline(yintercept = 1, color = 'red', linetype = 'dashed') + 
  geom_errorbar(aes(ymin=OR_min, ymax=OR_max), 
                width=1, position=position_dodge(width=0.75)) + 
  geom_point(pch =21, position=position_dodge(width=0.75), size = 2) +
  scale_fill_manual(values = group_col[c(1:3,5:6)]) + 
  scale_alpha_manual(values = c(0.3, 1)) + 
  facet_wrap(~celltype, nrow = 1, scales = 'fixed') + 
  theme_bw(base_size = 7) + ylim(c(0, 26)) + 
  xlab('Cell-TACIT Age Quantile') + 
  ylab('Odds Ratio (Brain SNPs/Other SNPs)') +
  guides(fill = guide_legend(override.aes = list(size = 1),nrow = 1, title=""),
         alpha = guide_legend(nrow = 1)) + 
  theme(legend.position = 'top',
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))+ 
  theme(legend.margin=margin(1,1,2,1),
        legend.box.margin=margin(-5,-5,-7,-5), 
        legend.key.size = unit(0.25, "cm"))
dev.off()


plot_fn2 = here(PLOTDIR, 'plots', 'fine-mapped_SNP_enrichmentByGroup_quartiles_other.pdf')
pdf(plot_fn2, height = 2.5, width = 4.25)
ggplot(snp_enrichment_df %>% filter(celltype %in% celltypes[c(4,5,7,8)], !is.infinite(OR_max)), 
       aes(x = quantile, y = OR, fill = group, alpha = FDR < 0.05)) + 
  geom_hline(yintercept = 1, color = 'red', linetype = 'dashed') + 
  geom_errorbar(aes(ymin=OR_min, ymax=OR_max), 
                width=1, position=position_dodge(width=0.75)) + 
  geom_point(pch =21, position=position_dodge(width=0.75), size = 2) +
  scale_fill_manual(values = group_col[c(1:3,5:6)]) + 
  scale_alpha_manual(values = c(0.3, 1)) + 
  facet_wrap(~celltype, nrow = 1, scales = 'fixed') + 
  theme_bw(base_size = 7) + ylim(c(0, 11)) + 
  xlab('Cell-TACIT Age Quantile') + 
  ylab('Odds Ratio (Brain SNPs/Other SNPs)') +
  guides(fill = guide_legend(override.aes = list(size = 1),nrow = 1, title=""),
         alpha = guide_legend(nrow = 1)) + 
  theme(legend.position = 'top',
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))+ 
  theme(legend.margin=margin(1,1,2,1),
        legend.box.margin=margin(-5,-5,-7,-5), 
        legend.key.size = unit(0.25, "cm"))
dev.off()











###########################################################################
## 6) calculate the enrichment at the trait level
traits = snps_annotated_df %>% filter(!duplicated(label), group != 'Other') %>% pull(label)
names(traits) = traits

snp_enrichment_df2 = parallel::mclapply(traits, function(y){
  snps_annotated_df %>% 
    filter(group %in% c('Other') | label %in% y) %>%
    mutate(group = factor(group), group = relevel(group, ref = 'Other'), 
           quantile = factor(quantile), quantile = relevel(quantile, ref = 'non-overlap')) %>%
    nest(data = c(-celltype)) %>%
    mutate(    test = map(data, ~glm(group ~ quantile, data = .x,family=binomial())),
               tidied = map(test, tidy)
    ) %>% unnest(cols = tidied) %>% 
    dplyr::select(-data, -test)
}, mc.cores = 8) %>% rbindlist(idcol = 'label') %>%
  filter(term != '(Intercept)') %>%
  mutate(OR = exp(estimate), FDR = p.adjust(p.value, 'fdr'), 
         OR_max = exp(estimate + std.error), 
         OR_min = exp(estimate - std.error), 
         term = gsub('quantile', '', term),
         celltype = factor(celltype, celltypes)) %>%
  dplyr::rename('quantile' = 'term') %>%
  arrange(log10(FDR) * OR ) %>% 
  inner_join(pheno %>% dplyr::select(label, group)) %>%
  mutate(group = factor(group), 
         quantile = factor(quantile), quantile = recode(quantile, !!!quartiles))

save_fn = here(PLOTDIR, 'tables', 'fine-mapped_SNP_enrichmentByTrait_quartiles_logistic_regression.xlsx')
snp_enrichment_df2 %>% writexl::write_xlsx(save_fn)


################################
## 7) plot the log odds ratio
plot_fn = here(PLOTDIR, 'plots', 'fine-mapped_SNP_enrichmentByTrait_quartiles_neurons.pdf')
pdf(plot_fn, height = 2.5, width = 5)
ggplot(snp_enrichment_df2 %>% filter(celltype %in% celltypes[1:4], !is.infinite(OR_max)), 
       aes(x = quantile, y = OR, fill = group, alpha = FDR < 0.05)) + 
  geom_hline(yintercept = 1, color = 'red', linetype = 'dashed') + 
  geom_errorbar(aes(ymin=OR_min, ymax=OR_max), 
                width=1, position=position_dodge(width=0.75)) + 
  geom_point(pch =21, position=position_dodge(width=0.75), size = 2) +
  scale_fill_manual(values = group_col[c(1:3,5:6)]) + 
  scale_alpha_manual(values = c(0.1, 1)) + 
  facet_wrap(~celltype, nrow = 1, scales = 'fixed') + 
  theme_bw(base_size = 7) + ylim(c(-0, 100)) + 
  xlab('Cell-TACIT Age Quantile') + 
  ylab('Enrichment Odds Ratio\n(Neuropsych SNPs/ Other SNPs)') +
  guides(fill = guide_legend(override.aes = list(size = 1),nrow = 1, title=""),
         alpha = guide_legend(nrow = 1)) + 
  theme(legend.position = 'top',
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))+ 
  theme(legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-7,-5), 
        legend.spacing.x = unit(0.05, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'))
dev.off()


plot_fn2 = here(PLOTDIR, 'plots', 'fine-mapped_SNP_enrichmentByTrait_quartiles_glia.pdf')
pdf(plot_fn2, height = 2.5, width = 5)
ggplot(snp_enrichment_df2 %>% filter(celltype %in% celltypes[5:8], !is.infinite(OR_max)), 
       aes(x = quantile, y = estimate, fill = group, alpha = FDR < 0.05)) + 
  geom_hline(yintercept = 0, color = 'red', linetype = 'dashed') + 
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), 
                width=1, position=position_dodge(width=0.75)) + 
  geom_point(pch =21, position=position_dodge(width=0.75), size = 2) +
  scale_fill_manual(values = group_col[c(1:3,5:6)]) + 
  scale_alpha_manual(values = c(0.1, 1)) + 
  facet_wrap(~celltype, nrow = 1, scales = 'fixed') + 
  theme_bw(base_size = 7) + ylim(c(-1, 5)) + 
  xlab('Cell-TACIT Age Quantile') + 
  ylab('Log Enrichment Ratio\n(Neuropsych SNPs/ Other SNPs)') +
  guides(fill = guide_legend(override.aes = list(size = 1),nrow = 1, title=""),
         alpha = guide_legend(nrow = 1)) + 
  theme(legend.position = 'top',
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))+ 
  theme(legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-7,-5), 
        legend.spacing.x = unit(0.05, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'))
dev.off()


