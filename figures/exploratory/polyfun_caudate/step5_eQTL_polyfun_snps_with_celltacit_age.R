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
library(RRHO)

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




#######################################################################
# 2) import the Cell Tacit tracks, subset to just intronic/distal peaks
bed_fn = list.files('data/raw_data/ldsc_caudate_zoonomia/CellTACIT', full.names = T, 
                    pattern = '.CellTACIT.mean.bed.g')
names(bed_fn) = basename(bed_fn) %>% ss('Corces2020.|.mean.bed.gz', 2)

cell_tacit_peaks = lapply(bed_fn, import) %>%
  lapply(function(gr){
    ## annotate peaks
    annot_df <- annotatePeak(gr, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, 
                             annoDb='org.Hs.eg.db') %>% 
      as.GRanges() %>% as.data.frame(row.names = seq(length(.))) %>%
      mutate(annotation = ss(annotation, ' '))
    # keep only non-coding distal/intronic peaks
    gr = gr[annot_df$annotation %in% c('Distal', 'Intron')]
    return(gr)
  })


## annotate each SNP with the overlapping peak
snps_celltacit_df = lapply(cell_tacit_peaks, function(gr){
  ret =snps_gr
  oo = findOverlaps(subject = gr, query = ret)
  mcols(ret)$score = 0
  mcols(ret)$score[queryHits(oo)] = mcols(gr)$score[subjectHits(oo)]
  return(data.frame(name = names(ret), CellTACIT = mcols(ret)$score))
}) %>% rbindlist(idcol = 'peak') %>% 
  mutate(celltype = ss(peak, '\\.', 1)) %>% inner_join(snps_df) 


compare_SNPs = snps_celltacit_df %>% 
  filter(PIP > 0.10 ) %>%
  nest(data = -c(group, celltype)) %>%
  mutate(test = map(data, ~cor.test(.x$PIP, .x$CellTACIT, method = "spearman")),
         tidied = map(test, tidy)
  ) %>% unnest(cols = tidied) %>% 
  dplyr::select(-data, -test) %>%
  mutate(FDR = p.adjust(p.value, 'fdr')) %>%
  arrange(FDR)






snps_celltacitMax_df = snps_celltacit_df %>%
  filter(group== 'SU') %>%
  filter(grepl('MSN', celltype)) %>%
  arrange(desc(CellTACIT)) %>% group_by(name, group) %>% 
  mutate(CellTACIT_max = max(CellTACIT)) %>%
  filter(!duplicated(name)) %>% ungroup()

pip.list1 = snps_celltacitMax_df %>% dplyr::select(name, PIP) %>% as.data.frame()
celltacit.list2 = snps_celltacitMax_df %>% dplyr::select(name, CellTACIT_max)%>% as.data.frame()

RRHO.celltacitAge <-  RRHO(pip.list1, celltacit.list2, BY=TRUE, alternative='enrichment')
pdf()
lattice::levelplot(RRHO.celltacitAge$hypermat.by)
dev.off()

pval.testing <- pvalRRHO(RRHO.celltacitAge, 1000)
pval.testing


###########################################################################
## 3) calculate the enrichment of fine-mapped SNPs in cell tacit enhancers
quartiles = c('0-25%', '25-50%', '50-75%', '75-100%')
names(quartiles) = paste0('quartile', 1:4); 
groups = names(group_col)[1:6]; names(groups) = names(group_col)[1:6]

plot_fn = here(PLOTDIR, 'plots', 'gtex_finemapped_eQTL_PIP_quartiles_main.pdf')
pdf(plot_fn, height = 2.5, width = 5)
ggplot(snps_annotated_df%>% filter(celltype %in% celltypes[c(1:3,6)]), 
                                   aes(x = quantile, y = eQTL_PIP, color = quantile)) + 
  geom_boxplot() +
  facet_wrap(~celltype, nrow = 1, scales = 'fixed') + 
  theme_bw(base_size = 7) + 
  ylab('Fraction Noncoding Variants\nin Cell-TACIT Age OCR') + 
  xlab('Posterior Probability Threshold') +
  guides(color = guide_legend(override.aes = list(size = 1),nrow = 1, title=""),
         alpha = guide_legend(nrow = 1)) + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1), 
        legend.position = 'top',
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-7,-5), 
        legend.spacing.x = unit(0.05, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'))
dev.off()


plot_fn2 = here(PLOTDIR, 'plots', 'fine-mapped_SNP_ecdfByGroup_quartiles_other.pdf')
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
  ylab('Log Enrichment Ratio\n(Brain SNPs/ Other SNPs)') +
  guides(fill = guide_legend(override.aes = list(size = 1),nrow = 1, title=""),
         alpha = guide_legend(nrow = 1)) + 
  theme(legend.position = 'top',
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))+ 
  theme(legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-5,-5,-7,-5), 
        legend.spacing.x = unit(0.05, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'))
dev.off()