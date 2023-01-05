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

cols = setNames(RColorBrewer::brewer.pal(4, 'Paired'),
                c('HSA', 'MMA', 'HSB', 'MMB'))


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
seg_summary_df = readRDS(save_fn2)

#################################
### get the prediction values ###
prediction_fn = here(DATADIR2, 'predictions') %>% 
  list.files(pattern = 'normPred.txt.gz', full.names = T) %>% 
  str_subset('reporter_enhancers')
names(prediction_fn) = basename(prediction_fn) %>% ss('\\.', 3)

pred_df = prediction_fn %>% lapply(fread, col.names = c('sequence', 'pred')) %>%
  rbindlist(idcol = 'Celltype') %>% 
  mutate(Condition = ss(sequence, '_', 1), 
         Condition = factor(Condition, c('HSA',  'HSB','MMA', 'MMB')),
         isNeuronal = case_when(grepl('MSN|INT', Celltype) ~ 'NeuN+', 
                                TRUE ~ 'NeuN-'), 
         isNeuronal = factor(isNeuronal, c('NeuN+', 'NeuN-'))) %>% 
  inner_join(seg_summary_df) %>%
  inner_join(proportion_df)

pred_df2 = pred_df %>% 
  group_by(Condition, isNeuronal) %>%
  mutate(pred_mean = sum(pred * prop), 
         pred_sem = sqrt(wtd.var(pred, prop*n()))/sqrt(n())) %>%
  distinct(Condition, isNeuronal, .keep_all = T)

pred_df2 %>% writexl::write_xlsx(here(PLOTDIR, 'tables', 
                                        'table_Sx_reporter_assay_CellTACIT-score_hs_mm.xlsx') )


pdf(here(PLOTDIR, 'plots', 'reporter_assay_CellTACITpred_barplot.pdf'), height = 1.5, width = 1.4)
ggplot(pred_df2,   aes(x=Condition, y=pred_mean, fill=Condition)) + 
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(ymin=pred_mean-pred_sem, 
                    ymax=pred_mean+pred_sem), width=.6) + 
  facet_grid(~isNeuronal)+
  ylab('Cell-TACIT Calibrated Probability Score') + xlab('') +
  scale_fill_manual(values = cols) + theme_bw(base_size = 5) + 
  theme(legend.position = 'none',
        plot.margin = margin(2, 2, -4, 1))
dev.off()


mod = lm(data = pred_df2, formula =  mCherry_mean ~ pred_mean)
summary(mod)
cor_p = with(pred_df2, cor.test(pred_mean, mCherry_mean, method="pearson"))
cor_s = with(pred_df2, cor.test(pred_mean, mCherry_mean, method="spearman"))

the_string = 
  paste0('R = ', signif(cor_p$estimate, 3), ', p= ', signif(cor_p$p.value, 2),'\n', 
         bquote(rho), ' = ', signif(cor_s$estimate, 3), ', p= ', signif(cor_s$p.value, 2))


pdf(here(PLOTDIR, 'plots', 'correlate_CellTACITpred_mCherryNum.pdf'), height = 1.5, width = 1.8)
ggplot(pred_df2,   aes(x=pred_mean, y=mCherry_mean, fill=Condition)) + 
  geom_abline(slope = mod$coefficients[2], intercept = mod$coefficients[1], 
              alpha = .5, linetype = 'dashed') + 
  geom_errorbar(aes(ymin=mCherry_mean-mCherry_sem, 
                    ymax=mCherry_mean+mCherry_sem), width = .015) + 
  geom_errorbarh(aes(xmin=pred_mean-pred_sem, 
                    xmax=pred_mean+pred_sem), height = 3) + 
  geom_point(aes(shape= isNeuronal), size = 1) +
  guides(fill = guide_legend(override.aes = list(pch = 21, size = 1), nrow = 1 ), 
         shape = guide_legend(override.aes = list(size = 1), nrow = 1 )) +
  # xlim(c(0, NA))+ ylim(c(0, NA)) +
  geom_text(x=.34, y=10, label=the_string, hjust = 1, size = 1.7) + 
  ylab(bquote('mCherry+ nuclei per'~mm^2)) +
  xlab('Cell-TACIT Calibrated Probability Score') +
  scale_shape_manual(values = c('NeuN+' = 24, 'NeuN-' = 21), name="") +
  scale_fill_brewer(palette="Paired", name="") + theme_bw(base_size = 5) + 
  theme(plot.margin = margin(2, 2, 2, 2), legend.position = 'bottom', 
        legend.spacing.x = unit(.02, 'cm'), legend.spacing.y = unit(.3, 'cm'), 
        legend.key.size = unit(.2, "cm"), 
        legend.box.margin=margin(-5,-10,-5,-10))
dev.off()










