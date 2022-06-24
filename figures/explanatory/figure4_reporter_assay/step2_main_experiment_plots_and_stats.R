library(tidyverse)
library(ggplot2)
library(arrow)
library(data.table)
library(here)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(RColorBrewer)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

PLOTDIR='figures/explanatory/figure4_reporter_assay'
DATADIR='data/raw_data/reporter_assay'

dir.create(here(PLOTDIR, 'tables'), showWarnings = F, recursive = T)

save_fn1 = here(PLOTDIR, 'rdas', 'reporter_assay_segmented_nuclei_per_nuclei.rds')
seg_filtered_df = readRDS(save_fn1)

save_fn = here(PLOTDIR, 'rdas', 'reporter_assay_segmented_nuclei_per_image.rds')
seg_per_img_df2 = readRDS(save_fn) %>%
  mutate(Condition = factor(Condition, c('HSA',  'HSB','MMA', 'MMB')))

cols = setNames(RColorBrewer::brewer.pal(4, 'Paired'),
                c('HSA', 'MMA', 'HSB', 'MMB'))

###########################################
### make the main plots for the figures ###

## get averages and std. error for barplots
seg_summary_df = seg_per_img_df2 %>% 
  dplyr::select(-contains(c('number', 'File')), -c('NeuN', 'GFP')) %>% 
  group_by(Animal_ID,isNeuronal) %>% 
  mutate(mCherry = mean(mCherry), 
         mCherryPerGFP = mean(mCherryPerGFP)) %>% 
  group_by(Condition,isNeuronal) %>% 
  summarise(mCherry_mean = mean(mCherry),
            mCherry_sem = sd(mCherry) / sqrt(n())) %>% 
  ungroup()

save_fn2 = here(PLOTDIR, 'rdas', 'reporter_assay_summary_nuclei.rds')
seg_summary_df %>% saveRDS(save_fn2)

pdf(here(PLOTDIR, 'plots', 'reporter_assay_mCherryNum_barplot.pdf'), height = 1.5, width = 1.6)
ggplot(seg_summary_df,   aes(x=Condition, y=mCherry_mean, fill=Condition)) + 
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(ymin=mCherry_mean-mCherry_sem, 
                    ymax=mCherry_mean+mCherry_sem), width=.6) + 
  facet_grid(~isNeuronal)+ 
  # ylim(c(0, 68)) +
  ylab(bquote('mCherry+ nuclei per'~mm^2)) + xlab('') +
  scale_fill_manual(values = cols) + theme_bw(base_size = 5) + 
  theme(legend.position = 'none',
        plot.margin = margin(2, 2, -4, 1))
dev.off()




##############################################################
### test significance using number mCherry positive nuclei ###
modCondMCherryNeuronal1 =  lmer(mCherry ~ Enhancer + Species + Cage + GFP + (1 | Animal_ID), 
                        data = seg_per_img_df2 %>% filter(isNeuronal == 'NeuN+'))
summary(modCondMCherryNeuronal1)


modCondMCherryGlia1 =  lmer(mCherry ~ Enhancer + Species + Cage + GFP + (1 | Animal_ID), 
                               data = seg_per_img_df2 %>% filter(isNeuronal == 'NeuN-'))
summary(modCondMCherryGlia1)


modCondMCherryNeuronal2 =  lmer(mCherry ~ Condition + Cage + GFP + (1 | Animal_ID), 
                        data = seg_per_img_df2 %>% filter(isNeuronal == 'NeuN+'))
summary(modCondMCherryNeuronal2)


modCondMCherryGlia2 =  lmer(mCherry ~ Condition + Cage + GFP + (1 | Animal_ID), 
                                data = seg_per_img_df2 %>% filter(isNeuronal == 'NeuN+'))
summary(modCondMCherryGlia2)


modList = list('model1_NeuN+' = modCondMCherryNeuronal1, 
               'model1_NeuN-' = modCondMCherryGlia1, 
               'model2_NeuN+' = modCondMCherryNeuronal2, 
               'model2_NeuN-' = modCondMCherryGlia2) %>% lapply(tidy)

outList = c(modList, list('per_image_data' = seg_per_img_df2 %>% 
                            dplyr::rename('GFP_per_mm2' = 'GFP', 
                                          'NeuN_per_mm2' = 'NeuN',
                                          'mCherry_per_mm2'= 'mCherry'), 
                          'per_condition_data' = seg_summary_df)) %>%
  writexl::write_xlsx(here(PLOTDIR, 'tables', 'reporter_assay_stats_and_data.xlsx'))




