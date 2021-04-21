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
library(here)


#########################################
# to be run in the root github directory
LABEL='Zoonomia_data'
PROJDIR='figures/exploratory/ldsc_caudate_zoonomia'
DATADIR=here('data/tidy_data/',LABEL)


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
df = df %>% arrange(Clade) %>%
  mutate(Clade = ifelse(grepl('Xen|Afro', as.character(Clade)), 
                             'Xenarthra & Afrotheria', as.character(Clade)),
         Clade = factor(Clade, unique(Clade)))
col_clade = df %>% select(Clade, col_clade)%>% filter(!duplicated(Clade)) %>% deframe()
col_order = df %>% select(Order, col_order) %>% filter(!duplicated(Order)) %>% deframe()

######################################################
# read in the LDSC partitioned heritability estimation
# run the command below to get number of mappable peaks
#  find /projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/data/tidy_data/Zoonomia_data/peaks -type f -name '*.gz' -exec bash -c 'echo -e "$(basename $1)\t$(zcat $1 | wc -l)"' dummy {} \; > /projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/data/tidy_data/Zoonomia_data/tables/zoonomia_peaks.txt
numPeak_fn =here(DATADIR,'tables','zoonomia_peaks.txt')
input = fread(numPeak_fn, col.names = c('file', 'numPeak'), header = F)

# get number of human peaks
humanPeaks_fn = here('data/raw_data/hg38/Corces_2020/peak') %>% 
  list.files(path = ., pattern = '.narrowPeak.gz', full.names = T) %>%
  grep(pattern = 'Corces2020_caudate\\.',value = TRUE)
humanPeaks_fn = humanPeaks_fn[! grepl('Consensus', humanPeaks_fn)]
names(humanPeaks_fn) = ss(basename(humanPeaks_fn), 'narrowPeak.gz')
input2 = lapply(humanPeaks_fn, fread, header = F) %>% 
  rbindlist(fill = T, idcol='file') %>% group_by(file) %>%
  summarise(numPeak = n())
input3 = bind_rows(input, input2)

#########################################
## format groupings and calculate conditional cell type enrichment p-value
pd = input3 %>% 
  mutate(
    Species = ss(file, '\\.', 4),
    Species = ifelse(is.na(Species), 'Homo_sapiens', Species), 
    celltype = ss(file, '\\.', 2) %>%
      factor(c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb',  'Astro', 
               'Microglia', 'OPC', 'Oligo')), 
    cell_group = case_when(
      grepl('MSN|INT', celltype) ~ 'Neuron', 
      TRUE ~ 'Glia'
    )) %>% inner_join(x = df, by = 'Species') %>%
  filter(!is.na(celltype))

pd$celltype %>% table()
to_label = c('Homo_sapiens', 'Mus_musculus', 'Macaca_mulatta')
pd = pd %>% 
  mutate(label = ifelse(Species %in% to_label, Species, NA)) %>%
  mutate(across(ends_with('N50'), ~as.numeric(.x))) %>%
  mutate(across(ends_with('L50'), ~as.numeric(.x)))

pd = pd %>% group_by(celltype) %>%
  mutate(tmp = numPeak[which(Species =='Homo_sapiens')],
         propMappable = numPeak/tmp) %>% select(-tmp) %>%
  ungroup() %>% mutate(label2 = ifelse(propMappable < .4, Species, 
              ifelse(Species %in% to_label, Species, NA)))

pd %>% filter(Species =='Orycteropus_afer') %>% pull(propMappable)

#################################
## make plots for presentation ##
dir.create(here('plots'))
height_ppt = 4; width_ppt = 8;
height_fig = 1.75; width_fig = 2.25; font_fig = 7

## Plot the number of peaks mapped across species ##
plot_fn1 =here(PROJDIR, 'plots','zoonomia_caudate_numPeaks.ppt.pdf')
pdf(width = width_ppt, height = height_ppt, file = plot_fn1, onefile = F)
pp1 = ggplot(data = pd, aes(x = Time.Since.Split.from.Human.TimeTree.median, 
                 y = numPeak/1000)) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se = TRUE) +  
  geom_jitter(aes(fill = Order), pch = 21, position=position_dodge(width=2)) +
  scale_fill_manual(values = col_order) + 
  geom_label_repel(aes(label = label, fill = Order), alpha = .7,direction ='both',
                   size = 2, show.legend = F,na.rm = T, 
                   nudge_y = .7 * mean(pd$numPeak/1000),
                   point.padding = .1, segment.color = 'grey', max.time = 4,
                   label.padding = .1, force_pull = .3, force = 40) +
  facet_wrap(~celltype, nrow = 2) + 
  theme_bw(base_size = 10) + 
  scale_y_continuous( limits = c(0,NA), expand = expansion(mult = c(0, .1))) +
  xlab('MY from Human (TimeTree median est.)') + 
  ylab('Mappable peaks (thousands)') + 
  guides(fill = guide_legend(nrow =3, title.position="top", override.aes = list(size = 3))) + 
  theme(legend.position = "bottom", legend.text=element_text(size=font_fig + 1),
        legend.title=element_text(size=font_fig + 3), legend.key.height=unit(.5,"line"))
print(pp1)
dev.off()


plot_fn2 =here(PROJDIR, 'plots','zoonomia_caudate_propMappable.ppt.pdf')
pdf(width = width_ppt, height = height_ppt, file = plot_fn2, onefile = F)
pp1 = ggplot(data = pd, aes(x = Time.Since.Split.from.Human.TimeTree.median, 
                            y = propMappable)) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se = TRUE) +  
  geom_jitter(aes(fill = Order), pch = 21, position=position_dodge(width=2)) +
  scale_fill_manual(values = col_order) + 
  geom_label_repel(aes(label = label, fill = Order), alpha = .7,direction ='both',
                   size = 2, show.legend = F,na.rm = T, 
                   nudge_y = .7 * mean(pd$propMappable),
                   point.padding = .1, segment.color = 'grey', max.time = 4,
                   label.padding = .1, force_pull = .3, force = 40) +
  facet_wrap(~celltype, nrow = 2) + 
  theme_bw(base_size = 10) + 
  scale_y_continuous( limits = c(0,NA), expand = expansion(mult = c(0, .1))) +
  xlab('MY from Human (TimeTree median est.)') + 
  ylab('Proportion of Mappable Hg Peaks') + 
  guides(fill = guide_legend(nrow =3, title.position="top", override.aes = list(size = 3))) + 
  theme(legend.position = "bottom", legend.text=element_text(size=font_fig + 1),
        legend.title=element_text(size=font_fig + 3), legend.key.height=unit(.5,"line"))
print(pp1)
dev.off()


cell = 'Astro'
for (cell in unique(pd$celltype)){
plot_fn2 =here(PROJDIR, 'plots',paste0('zoonomia_caudate_N50_propMappable_',cell,'.ppt.pdf'))
pdf(width = width_ppt, height = height_ppt, file = plot_fn2, onefile = F)
pp1 = ggplot(data = pd %>% filter(celltype == cell) %>% filter(!is.na(Contig.N50)), 
             aes(x = log10(Contig.N50), y = propMappable)) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1, se = TRUE) +  
  geom_point(aes(fill = Order), alpha = .4, pch = 21) +
  scale_fill_manual(values = col_order) + 
  geom_label_repel(aes(label = label2, fill = Order), alpha = .7,direction ='both',
                   size = 2, show.legend = F,na.rm = T, 
                   point.padding = .1, segment.color = 'grey', max.time = 4,
                   label.padding = .1, force_pull = 1, force = 200) +
  facet_grid(~Clade,  scale = 'free_x', space = 'fixed') + 
  theme_bw(base_size = 11) + 
  scale_y_continuous( limits = c(0,NA), expand = expansion(mult = c(0, .1))) +
  xlab('log10( Genome Contig N50 )') + 
  ylab('Proportion of Mappable Hg Peaks') + 
  guides(fill = guide_legend(nrow =3, title.position="top", 
                override.aes = list(size = 3, alpha = 1))) + 
  theme(legend.position = "bottom", legend.text=element_text(size=font_fig + 1),
        legend.title=element_text(size=font_fig + 3), legend.key.height=unit(.5,"line"))
print(pp1)
dev.off()
}
