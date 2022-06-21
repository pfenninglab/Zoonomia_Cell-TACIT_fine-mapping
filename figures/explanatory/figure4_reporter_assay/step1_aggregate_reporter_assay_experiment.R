library(tidyverse)
library(ggplot2)
library(arrow)
library(data.table)
library(here)
library(lme4)
library(lmerTest)
library(broom.mixed)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

PLOTDIR='figures/explanatory/figure4_reporter_assay'
DATADIR='data/raw_data/reporter_assay'

dir.create(here(PLOTDIR, 'plots'), showWarnings = F, recursive = T)
dir.create(here(PLOTDIR, 'tables'), showWarnings = F, recursive = T)
dir.create(here(PLOTDIR, 'rdas'), showWarnings = F, recursive = T)

### set some thresholds
mpp = 0.365 # micron per pixel scaling factor
min_nuclei_diam = 5 # microns
max_nuclei_diam = 25 # microns
neun_frac_thres = 0.05

###################################################################
### read in the blinding key table with experimental conditions ###
df = readxl::read_xlsx(here(DATADIR, 'tables', 'Reporter_Assay_Imaging_Blind_Key.xlsx')) %>%
  mutate(Blind_number = as.numeric(str_trim(Blind_number)))

### read in the per-cell segmentation data for each image
## the denoised data seemed the best segmentation performance, so ignore other images
segmentation_fn = here(DATADIR, 'tables') %>%
  list.files(pattern = '_segmented_objects.feather', full.names = T) %>%
  str_subset('Denoised')
names(segmentation_fn) = basename(segmentation_fn) %>% ss(' ', 1)

seg_df = segmentation_fn %>% lapply(read_feather) %>%
  rbindlist(idcol = 'ImageFile')  %>%
  mutate(Blind_number = ImageFile %>% ss('Slide|\\.', 2) %>% as.numeric())

## scale the area from pixel units to micron units, get nuclei object diameter 
seg_df = seg_df %>% 
  mutate(
    x = (x - mean(x)) * mpp, y = (y - mean(y)) * mpp,
    dist_to_center = sqrt(x^2 + y^2),
    area = area* (mpp)^2, diam = sqrt(area/pi) * 2)

table(seg_df$Segmentation)
summary(seg_df$dist_to_center)
summary(seg_df$x)


######################################################################
### make some QC plots to find good filters of noisy segmentations ###
seg_filtered_df = seg_df %>% 
  filter(diam < max_nuclei_diam, diam > min_nuclei_diam, 
         solidity > .925, dist_to_center < 0.8 * max(dist_to_center)) %>%
  inner_join(df)

save_fn1 = here(PLOTDIR, 'rdas', 'reporter_assay_segmented_nuclei_per_nuclei.rds')
seg_filtered_df %>% saveRDS(save_fn1)

with(seg_filtered_df %>% filter(Segmentation =='mCherry'), summary(diam))

pdf(here(PLOTDIR, 'plots', 'segmentation_qc_measures.pdf'), height = 8, width = 11)
ggplot(seg_filtered_df, aes(y = solidity, x = Segmentation)) + 
  geom_violin() + facet_wrap(~ImageFile) + theme_bw()

ggplot(seg_filtered_df, aes(y = diam, x = Segmentation)) + 
  geom_violin() + facet_wrap(~ImageFile) + theme_bw()
dev.off()


seg_df3 = seg_filtered_df  %>% filter(Segmentation != 'NeuN')
pdf(here(PLOTDIR, 'plots', 'fluorescence_intensity.pdf'), height = 8, width = 11)
ggplot(seg_df3, aes(y = mean_mCherry, x = mean_NeuN)) + 
  geom_hex() + facet_wrap(~ImageFile) + 
  scale_fill_binned(type = "viridis") + theme_bw()

ggplot(seg_df3, aes(y = mean_GFP, x = mean_NeuN)) + 
  geom_hex() + facet_wrap(~ImageFile) + 
  scale_fill_binned(type = "viridis") + theme_bw()

ggplot(seg_df3, aes(y = mean_GFP, x = mean_mCherry)) + 
  geom_hex() + facet_wrap(~ImageFile) + 
  scale_fill_binned(type = "viridis") + theme_bw()

ggplot(seg_df3, aes(y = percent_NeuN, x = max_NeuN)) + 
  geom_hex() + facet_wrap(~ImageFile) + 
  scale_fill_binned(type = "viridis") + theme_bw()

ggplot(seg_df3, aes(y = percent_NeuN, x = max_NeuN ))+ 
  geom_hex() + facet_wrap(~ImageFile) + 
  scale_fill_binned(type = "viridis") + theme_bw()
dev.off()


######################################################################
### make some QC plots to find good filters of noisy segmentations ###
seg_per_img_df = seg_filtered_df %>% 
  mutate(isNeuronal = ifelse(percent_NeuN > neun_frac_thres, 'NeuN+', 'NeuN-')) %>% 
  group_by(ImageFile, Blind_number, Segmentation, isNeuronal) %>% 
  summarise(num_nuclei = n(),
            measured_area_mm2 = (2 * pi * max(x)^2) / 1e6,
            nuclei_per_mm2 = num_nuclei / measured_area_mm2 ) %>%
  ungroup() %>% inner_join(df)

pdf(here(PLOTDIR, 'plots', 'reporter_assay_conditionByChannel.pdf'), height = 8, width = 11)
ggplot(seg_per_img_df, aes(y = num_nuclei, x = Condition)) + 
  geom_boxplot() + geom_point(aes(fill = Condition), pch= 21) + 
  facet_wrap(~Segmentation + isNeuronal, scale = 'free') + 
  scale_colour_viridis_d() + theme_bw()
dev.off()


seg_per_img_df2 = seg_per_img_df %>% 
  dplyr::select(-c(measured_area_mm2, num_nuclei)) %>%
  pivot_wider(values_from = nuclei_per_mm2, names_from = 'Segmentation') %>%
  group_by(ImageFile) %>% 
  mutate(mCherryPerGFP = mCherry/GFP, 
         mCherryPerGFPPerNeuN = mCherry/GFP/NeuN, 
         species = Condition %>% ss('A|B', 1),
         enhancer = Condition %>% ss('HS|MM', 2)) %>%
  group_by(Condition) %>% 
  mutate(Rep = sample(LETTERS[1:3])[as.numeric(factor(Animal_ID))]) %>%
  ungroup()

save_fn = here(PLOTDIR, 'rdas', 'reporter_assay_segmented_nuclei_per_image.rds')
seg_per_img_df2 %>% saveRDS(save_fn)


pdf(here(PLOTDIR, 'plots', 'reporter_assay_condition_neuronal.pdf'), height = 3, width = 4.5)
ggplot(seg_per_img_df2, aes(y = mCherryPerGFP, x = Condition )) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(shape = Cage), position=position_dodge(width = .5)) + 
  facet_wrap(~isNeuronal, scale = 'free_y') + 
  theme_bw() + theme(legend.position = 'none')

ggplot(seg_per_img_df2, aes(y = mCherry, x = Condition )) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(shape = Cage), position=position_dodge(width = .5)) + 
  facet_grid(~isNeuronal) + 
  theme_bw() + theme(legend.position = 'none')
dev.off()


pdf(here(PLOTDIR, 'plots', 'reporter_assay_condition_neuronal.pdf'), height = 3, width = 4.5)
ggplot(seg_per_img_df2, aes(y = mCherryPerGFP, x = Condition )) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(shape = Cage), position=position_dodge(width = .5)) + 
  facet_wrap(~isNeuronal, scale = 'free_y') + 
  theme_bw() + theme(legend.position = 'none')

ggplot(seg_per_img_df2, aes(y = mCherry, x = Condition )) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(shape = Cage), position=position_dodge(width = .5)) + 
  facet_grid(~isNeuronal) + 
  theme_bw() + theme(legend.position = 'none')
dev.off()

pdf(here(PLOTDIR, 'plots', 'reporter_assay_proportion_neuronal.pdf'), height = 3, width = 4.5)
ggplot(seg_per_img_df2, aes(y = mCherry, x = Animal_ID)) + 
  geom_bar(aes(fill = isNeuronal), stat = "summary", fun = "mean", position = 'fill') + 
  facet_grid(~Condition, space = 'free', scale = 'free') + 
  theme_bw() + theme(legend.position = 'none') + 
  ylab('Proportion of mCherry+ nuclei') + 
  xlab('Animal Replicate') + 
  theme(axis.text.x=element_blank())

ggplot(seg_per_img_df2, aes(y = GFP, x = Animal_ID)) + 
  geom_bar(aes(fill = isNeuronal), stat = "summary", fun = "mean", position = 'fill') + 
  facet_grid(~Condition, space = 'free', scale = 'free') + 
  theme_bw() + theme(legend.position = 'none') + 
  ylab('Proportion of GFP+ nuclei') + 
  xlab('Animal Replicate') + 
  theme(axis.text.x=element_blank())
dev.off()