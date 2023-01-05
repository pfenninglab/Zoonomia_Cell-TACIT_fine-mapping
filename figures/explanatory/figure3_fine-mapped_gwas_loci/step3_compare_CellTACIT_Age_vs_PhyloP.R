ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(here)
library(tidyverse)
library(rcartocolor)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(GenomicRanges)
library(rtracklayer)
library(trackViewer)
library(tidymodels)
library(precrec)
library(lmerTest)
library(lme4)

library(ChIPseeker)
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))

# transformation function for axis sig figs
scaleFUN <- function(x) sprintf("%.1f", x)

PLOTDIR='figures/explanatory/figure3_fine-mapped_gwas_loci'
DATADIR='data/raw_data'

i_am(file.path(PLOTDIR, 'step3_compare_CellTACIT_Age_vs_PhyloP.R'))

###############################
## read in the CellTACIT Scores
save_track_fn= here('figures/explanatory/figure3_fine-mapped_gwas_loci/rdas',
                    paste0('fig3_trackViewer.CellTACIT_track.rds'))
CellTACIT_track = readRDS(file = save_track_fn)
CellTACIT_peakList = lapply(CellTACIT_track, function(x) x@dat)
CellTACIT_peakList = CellTACIT_peakList[rev(names(CellTACIT_peakList))]
CellTACIT_peakList = mapply(function(gr, lab) {
  gr$celltype = lab
  
  ## annotate peaks
  annot_df <- annotatePeak(gr, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, 
                           annoDb='org.Hs.eg.db') %>% 
    as.GRanges() %>% as.data.frame(row.names = seq(length(.))) %>%
    mutate(annotation = ss(annotation, ' '))
  # keep only non-coding distal/intronic peaks
  gr = gr[annot_df$annotation %in% c('Distal', 'Intron')]
  
  return(gr)
  }, CellTACIT_peakList, names(CellTACIT_peakList)) %>% GRangesList()
CellTACIT_peaks = unlist(CellTACIT_peakList)

###############################
## read in the PhyloP in OCRs
save_fn = here('data/raw_data/hg38/Corces_2020/rdas', 
               'Corces2020_caudate.peakList.phyloPconsFrac.rds')
consFDR_peakList = readRDS(file = save_fn)
consFDR_peakList = consFDR_peakList[names(CellTACIT_peakList)]
consFDR_peakList = mapply(function(gr1, gr2){
  oo = findOverlaps(query = gr1, subject = gr2)
  return(gr2[unique(subjectHits(oo))])
}, gr1 =CellTACIT_peakList, gr2 = consFDR_peakList )
all.equal(lengths(consFDR_peakList), lengths(CellTACIT_peakList))

consFDR_peakList = mapply(function(gr, lab) {
  gr$celltype = lab
  return(gr)
}, consFDR_peakList, names(consFDR_peakList)) %>% GRangesList()
consFDR_peaks = unlist(consFDR_peakList)


#################################
## read in the saved polyfun SNPs
finemap_df = here(DATADIR, 'polyfun_caudate/rdas/polyfun_caudate_finemapped_snps_20220718.rds') %>%
  readRDS() %>% mutate(start = POS_hg38, end = POS_hg38) %>%
  dplyr::select(-c(POS_hg19,POS_hg38, SNPVAR:P, BETA_MEAN:MAF, index:runPolyfun, h2:base))
finemap_gr = finemap_df %>% GRanges()
seqlevelsStyle(finemap_gr) = 'UCSC'

## annotate fine-mapped SNPs in distal and intronic regions
annot_df2 <- annotatePeak(finemap_gr, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, 
                         annoDb='org.Hs.eg.db') %>% 
  as.GRanges() %>% as.data.frame(row.names = seq(length(.))) %>%
  mutate(annotation = ss(annotation, ' '))
# keep only non-coding distal/intronic SNPs
finemap_gr = finemap_gr[annot_df2$annotation %in% c('Distal', 'Intron')]
finemap_df = finemap_df[annot_df2$annotation %in% c('Distal', 'Intron'),]

## add the single base phyloP to the SNPs
bws = list.files('/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/phyloP/human-centered-200m-Feb2021', pattern = '.bigWig', full.names = T)
names(bws) = basename(bws) %>% ss('\\.', 2)
phyloP_tracks <- lapply(bws, importScore, format="BigWig", ranges=finemap_gr)
phyloP_peaks = lapply(phyloP_tracks, function(x) x@dat) %>% GRangesList() %>% unlist()

ooPhylop= findOverlaps(subject = finemap_gr, query = phyloP_peaks)
finemap_df$phyloP = 0
finemap_df$phyloP[subjectHits(ooPhylop)] = phyloP_peaks$score[queryHits(ooPhylop)]
finemap_df$phyloPcons = finemap_df$phyloP>=2.270

## gather conservation scores per cell type
finemap_df2 = lapply(names(consFDR_peakList), function(celltype){
  tmp_df = finemap_df; tmp_gr = finemap_gr
  tmp_df$celltype = celltype
  
  ## add the ,max phyloP enhancer fraction conserved
  ooCons = findOverlaps(subject = tmp_gr, query = consFDR_peakList[[celltype]])
  scores = split(consFDR_peakList[[celltype]]$score[queryHits(ooCons)], subjectHits(ooCons))
  tmp_df$consFrac = 0
  tmp_df$consFrac[as.numeric(names(scores))] = sapply(scores, max)
  
  ## add the max CellTACIT scores
  ooCellTACIT = findOverlaps(subject = tmp_gr, query = CellTACIT_peakList[[celltype]])
  scores2 = split(CellTACIT_peakList[[celltype]]$score[queryHits(ooCellTACIT)], subjectHits(ooCellTACIT))
  tmp_df$CellTACIT_Age = 0
  tmp_df$CellTACIT_Age[as.numeric(names(scores2))] = sapply(scores2, max)
  
  return(tmp_df)
}) %>% rbindlist()

finemap_df2$celltype = factor(finemap_df2$celltype, names(CellTACIT_peakList))
finemap_df2 %>% filter(celltype == 'MSN_D2', SNP == 'rs7933981') %>% pull(CellTACIT_Age)


# #################################
## output fine-mapped SNPs to table
dir.create(here(PLOTDIR,'tables'), showWarnings = F)
# table_fn = here(PLOTDIR,'tables',
#                 'Data_S9_polyfun_caudate_finemapped_snps_20220718.xlsx')
# finemap_df %>% dplyr::select(-c(population:end)) %>%
#   split(., .$label) %>% writexl::write_xlsx(path = table_fn)
# finemap_df %>% saveRDS(here(PLOTDIR,'rdas',
#                             'Data_S9_polyfun_caudate_finemapped_snps_20220718.rds'))

table2_fn = here(PLOTDIR,'tables',
                'Data_S10_CellTACIT_Age_caudate_finemapped_snps_20220718.xlsx')
finemap_df2 %>% dplyr::select(-c(population:end)) %>%
  filter(CellTACIT_Age > 0) %>%
  split(., .$celltype) %>% writexl::write_xlsx(path = table2_fn)
finemap_df2 %>% saveRDS(here(PLOTDIR,'rdas',
                            'Data_S10_CellTACIT_Age_caudate_finemapped_snps_20220718.rds'))


# #################################
## aggregate cell tacit Age across cell types
finemap_df3 = finemap_df2 %>%
  group_by(match, name) %>%
  mutate(celltype = celltype[which.max(CellTACIT_Age)], 
         CellTACIT_Age = sum(CellTACIT_Age), 
         consFrac = max(consFrac)) %>%
  distinct(match, name, CellTACIT_Age, consFrac, celltype, .keep_all = T)

finemap_df3 %>% filter(CellTACIT_Age > 0, PIP > .1) %>%
  filter(!duplicated(name)) %>% nrow() # 946 PIP>.1

finemap_df3 %>% filter(phyloPcons, PIP > .1) %>%
  filter(!duplicated(name)) %>% nrow() # 693 PIP>.1

with(finemap_df3, table(CellTACIT_Age > 10, PIP > .1)) 
with(finemap_df3, table(phyloPcons, PIP > .1)) 

## relationship b/t PIP and CellTACIT Age across all traits
lm1 = finemap_df3 %>% filter(PIP > 0.1) %>% 
  mutate(CellTACIT_Age = CellTACIT_Age/100) %>%
  lm(PIP ~ CellTACIT_Age + phyloPcons + consFrac, data = .)
summary(lm1)

## relationship b/t PIP and CellTACIT Age regress group
lm2 = finemap_df3 %>% filter(PIP > 0.1) %>% 
  mutate(CellTACIT_Age = CellTACIT_Age/100) %>%
  lm(PIP ~  group + CellTACIT_Age + phyloPcons  + consFrac, data = .) 
summary(lm2)

## relationship b/t PIP and CellTACIT Age regress trait
lm3 = finemap_df3 %>% filter(PIP > 0.1) %>% 
  mutate(CellTACIT_Age = CellTACIT_Age/100) %>%
  lm(PIP ~  trait + CellTACIT_Age + phyloPcons  + consFrac, data = .)
summary(lm3)

## nested effects w/ group and trait
lme1 =  finemap_df3 %>% filter(PIP > 0.1, CellTACIT_Age >0) %>% 
  mutate(CellTACIT_Age = CellTACIT_Age/100) %>%
  lme4::lmer(PIP ~  CellTACIT_Age +phyloPcons+ (1|group) + (1|trait) , data = .)
summary(lme1)

lme2 =  finemap_df3 %>% filter(PIP > 0.1, CellTACIT_Age >0) %>% 
  mutate(CellTACIT_Age = CellTACIT_Age/100) %>%
  lme4::lmer(PIP ~  phyloPcons + (1|group) + (1|trait) , data = .)

lme3 =  finemap_df3 %>% filter(PIP > 0.1, CellTACIT_Age >0) %>% 
  mutate(CellTACIT_Age = CellTACIT_Age/100) %>%
  lme4::lmer(PIP ~  CellTACIT_Age + (1|group) + (1|trait) , data = .)

anova(lme1, lme2) # effect of Cell TACIT Age
# Models:
#   lme2: PIP ~ phyloPcons + (1 | group) + (1 | trait)
# lme1: PIP ~ CellTACIT_Age + phyloPcons + (1 | group) + (1 | trait)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
# lme2    5 293.24 317.50 -141.62   283.24                       
# lme1    6 289.50 318.61 -138.75   277.50 5.7415  1    0.01657 *
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



anova(lme1, lme3) # effect of in PhyloPcons bas
# Models:
#   lme3: PIP ~ CellTACIT_Age + (1 | group) + (1 | trait)
# lme1: PIP ~ CellTACIT_Age + phyloPcons + (1 | group) + (1 | trait)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# lme3    5 307.01 331.27 -148.50   297.01                         
# lme1    6 289.50 318.61 -138.75   277.50 19.514  1  9.986e-06 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# #################################
## make some plots
pdf('tmp.pdf')
ggplot(finemap_df2 %>% filter(consFrac > 0 & CellTACIT_Age >0, PIP > .1), 
       aes(x = consFrac, y = CellTACIT_Age)) + 
  geom_abline( slope = 100, intercept = 0, color = 'red') + 
  geom_point(aes(fill = PIP, alpha = PIP), pch = 21)+ 
  scale_fill_viridis_c(option = 'inferno') + 
  facet_grid(celltype~group, scales = 'free') + 
  xlim(c(0,1)) +  theme_bw()

ggplot(finemap_df2 %>% filter(consFrac > 0 & CellTACIT_Age >0, PIP > .1, phyloP > 0), 
       aes(x = phyloP, y = CellTACIT_Age)) + 
  geom_abline( slope = 10, intercept = 0, color = 'red') + 
  geom_point(aes(fill = PIP, alpha = PIP), pch = 21)+ 
  scale_fill_viridis_c(option = 'inferno') + 
  facet_grid(celltype~group, scales = 'free') + 
  xlim(c(0,10)) +theme_bw()

dev.off()

plot_fn1= here(PLOTDIR, 'plots/sfig_CellTACIT_vs_consFrac_score_20220718.pdf')
pdf(plot_fn1,width = 2.25*4,  height = 1.5*4)
ggplot(finemap_df3 %>% filter(consFrac > 0 & CellTACIT_Age >0, PIP > .1), 
       aes(x = consFrac, y = CellTACIT_Age)) + 
  geom_point(aes(fill = PIP, alpha = PIP, size = PIP,), pch = 21)+ 
  geom_abline( slope = 100, intercept = 0, color = 'red') + 
  scale_fill_viridis_c() + 
  xlab('Max SNP OCR Fraction of Conserved Bases') + ylab('Max SNP CellTACIT Age') +
  facet_grid(celltype~group, scales = 'free_x') + 
  xlim(c(0,1)) + theme_bw() + scale_x_continuous(labels=scaleFUN)
dev.off()

plot_fn2= here(PLOTDIR, 'plots/sfig_CellTACIT_vs_phyloP_score_20220718.pdf')
pdf(plot_fn2,width = 2.25*4,  height = 1.5*4)
ggplot(finemap_df3 %>% filter(phyloP > 0 & CellTACIT_Age >0, PIP > .1), 
       aes(x = phyloP, y = CellTACIT_Age)) + 
  geom_point(aes(fill = PIP, alpha = PIP, size = PIP, shape = celltype), pch = 21)+ 
  geom_abline( slope = 10, intercept = 0, color = 'red') + 
  scale_fill_viridis_c() + 
  xlab('Mammals PhyloP @ SNP') + ylab('Max SNP CellTACIT Age') +
  facet_grid(celltype~group, scales = 'free_x') + 
  xlim(c(0,10)) + theme_bw() + scale_x_continuous(labels=scaleFUN)
dev.off()


# #################################
## read in the saved polyfun SNPs
# finemap_df3 = finemap_df3 %>% filter(CellTACIT_Age > 0)
finemap_df3$label = with(finemap_df3, ifelse(PIP > .90, 'High PIP', 'Low PIP'))

scores = with(finemap_df3, join_scores(phyloP, consFrac, CellTACIT_Age))
labels = with(finemap_df3, join_labels(label, label, label) )
mmmdat1 <- mmdata(scores, labels, modnames= c("phyloP", "consFrac", 'CellTACITage'), dsids = c(1, 2, 3))
mmcurves <- evalmod(mmmdat1)
mmcurves
precrec_dat =   as.data.frame(mmcurves)

pdf('tmp2.pdf',width = 2.25*4,  height = 2.25*4)
ggplot(precrec_dat %>% filter(type =='ROC'), aes(x = x, y = y, color = modname)) +
  geom_line() +
  theme_bw()+ scale_x_continuous(labels=scaleFUN)

ggplot(precrec_dat %>% filter(type =='PRC'), aes(x = x, y = y, color = modname))+
  geom_line() +
  theme_bw()+ scale_x_continuous(labels=scaleFUN)
dev.off()


finemap_df2 %>% filter(PIP > 0.1) %>%
  nest(data = -c(celltype, group)) %>% 
  mutate(
    test = map(data, ~ lm(PIP ~ phyloP + CellTACIT_Age + consFrac, data = .x)), # S3 list-col
    tidied = map(test, tidy)
  ) %>% 
  unnest(cols = tidied) %>% 
  filter(term != '(Intercept)') %>%
  dplyr::select(-data, -test) %>%
  writexl::write_xlsx('tmp.xlsx')

