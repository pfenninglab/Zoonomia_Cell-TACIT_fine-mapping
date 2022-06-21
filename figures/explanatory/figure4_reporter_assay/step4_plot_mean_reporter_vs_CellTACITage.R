library(tidyverse)
library(ggplot2)
library(arrow)
library(data.table)
library(here)
library(ArchR)
library(Hmisc)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

PLOTDIR='figures/explanatory/figure4_reporter_assay'
DATADIR='data/raw_data/reporter_assay'
DATADIR2='data/raw_data/select_mpra_candidates'

##################################################################
celltypes = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb','Astro', 'Microglia', 'OPC', 'Oligo') 
## pre-calculate the proportion of each cell type
proj = loadArchRProject(here('data/raw_data/hg38/Corces_2020/ArchR_Corces2020_caudate_labeled'))
proportion_df = getCellColData(proj) %>% as_tibble() %>%
  group_by(Sample, Clusters2) %>% 
  summarise(isNeuronal = case_when(grepl('MSN|INT', Clusters2) ~ 'NeuN+', 
                                   TRUE ~ 'NeuN-'), 
            num = n()) %>% 
  distinct(Sample, Clusters2, isNeuronal, .keep_all = TRUE) %>%
  group_by(Sample, isNeuronal) %>% mutate(prop = num/sum(num)) %>% 
  group_by(Clusters2) %>% summarise(prop = mean(prop)) %>%
  dplyr::rename('Celltype' = 'Clusters2')

## load in the 
save_fn2 = here(PLOTDIR, 'rdas', 'reporter_assay_summary_nuclei.rds')
seg_summary_df = readRDS(save_fn2) %>%
  mutate(Enhancer = ss(as.character(Condition), 'HS|MM', 2) %>% factor(), 
         Species = ss(as.character(Condition), 'A|B', 1) %>% factor()) %>%
  group_by(isNeuronal, Enhancer) %>%
  mutate(
    indH = which(Species =='HS'), 
    indM = which(Species =='MM'), 
    mCherry_diff = mCherry_mean[indH] - mCherry_mean[indM],
    mCherry_mean = mean(mCherry_mean), 
    mCherry_sem = sqrt((mCherry_sem * sqrt(3) )^2 + (mCherry_sem * sqrt(3))^2) / sqrt(3)
  ) %>%
  ungroup() %>% dplyr::select(-c(Species, indH, indM, Condition)) %>% 
  distinct(isNeuronal, Enhancer, .keep_all = T)

enhancers = c('hg38:chr11:113567061-113567561:250' = 'A', 
              'hg38:chr11:113577668-113578168:250' = 'B')

#################################
### get the CellTACITage values ###
CellTACITage_fn = here(DATADIR, 'rdas') %>% 
  list.files(pattern = 'CellTACIT.mean.rds', full.names = T)

names(CellTACITage_fn) = basename(CellTACITage_fn) %>% ss('\\.', 2)

pred_df = CellTACITage_fn %>% lapply(readRDS) %>% 
  lapply(function(gr) { gr$name = names(gr); return(gr)}) %>%
  lapply(as.data.frame) %>%
  rbindlist(idcol = 'Celltype') %>% 
  mutate(Enhancer =enhancers[name] %>% factor(),
         isNeuronal = case_when(grepl('MSN|INT', Celltype) ~ 'NeuN+', 
                                TRUE ~ 'NeuN-'), 
         isNeuronal = factor(isNeuronal, c('NeuN+', 'NeuN-'))) %>% 
  inner_join(seg_summary_df) %>%
  inner_join(proportion_df)

pred_df2 = pred_df %>% 
  group_by(Enhancer, isNeuronal) %>%
  mutate(pred_mean = sum(score * prop), 
         pred_sem = sqrt(wtd.var(score, prop*n()))/sqrt(n())) %>%
  dplyr::select(-c(Celltype)) %>% 
  distinct(Enhancer, isNeuronal, .keep_all = T)


mod = lm(data = pred_df2, formula =  mCherry_mean ~ pred_mean)
summary(mod)

mod2 = lm(data = pred_df2, formula =  mCherry_diff ~ pred_mean)
summary(mod2)
cor_p = with(pred_df2, cor.test(pred_mean, mCherry_mean, method="pearson"))
cor_s = with(pred_df2, cor.test(pred_mean, mCherry_mean, method="spearman"))

the_string = 
  paste0('R = ', signif(cor_p$estimate, 3), ', p= ', signif(cor_p$p.value, 2),'\n', 
         bquote(rho), ' = ', signif(cor_s$estimate, 3), ', p= ', signif(cor_s$p.value, 2))


pdf(here(PLOTDIR, 'plots', 'correlate_CellTACITage_mCherryDiff.pdf'), height = 1.7, width = 2.375)
ggplot(pred_df2,   aes(x=pred_mean, y=abs(mCherry_diff), fill=Enhancer)) + 
  geom_abline(slope = mod2$coefficients[2], intercept = mod$coefficients[1], 
              alpha = .75, linetype = 'dashed') + 
  geom_errorbar(aes(ymin=abs(mCherry_diff)-mCherry_sem, 
                    ymax=abs(mCherry_diff)+mCherry_sem), width = 2) + 
  geom_errorbarh(aes(xmin=pred_mean-pred_sem, 
                    xmax=pred_mean+pred_sem), height = .5) + 
  geom_point(aes(shape= isNeuronal), size = 2) +
  guides(fill = guide_legend(override.aes = list(pch = 21, size = 2), nrow = 1 ), 
         shape = guide_legend(override.aes = list(size = 2), nrow = 1 )) +
  # xlim(c(0, NA))+ ylim(c(0, NA)) +
  # geom_text(x=.33, y=10, label=the_string, hjust = 1, size = 2.5) + 
  ylab(bquote('Human-mouse abs. diff. mCherry+ nuc./'~mm^2)) +
  xlab('Cell-TACIT Age') +
  scale_shape_manual(values = c('NeuN+' = 24, 'NeuN-' = 21), name="") +
  scale_fill_manual(values = c('A' = 'blue', 'B' = 'green')) + 
  theme_bw(base_size = 5) + 
  theme(plot.margin = margin(2, 1, 2, 1),
        legend.position = 'none')
dev.off()










