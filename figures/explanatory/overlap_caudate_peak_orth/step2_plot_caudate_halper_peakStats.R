# to be run in the root github directory
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
suppressMessages(library(GenomicRanges))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(ggsci))
suppressMessages(library(ggpubr))
suppressMessages(library(ChIPseeker))

LABEL='caudate_conservation_ldsc'
setwd('figures/explanatory/overlap_caudate_peak_orth')
PROJDIR=file.path('../../../data/raw_data/',LABEL)
source('../../../code/raw_code/hal_scripts/gen_enh_ortholog_sets.R')




###############################################################
## load the rData files of caudate peaks in hg38 coordinates ##
save_fn = file.path(PROJDIR, 'rdas', 'human_macaque_mouse_orthologs_peakList.rda')
group_class = c("hgPeaks", "rm2Hg", "hgRmOrth", "mm2Hg", "hgMmOrth")
annot_class = c('Distal.Intergenic','Promoter',"5' UTR", 'Exon', "Intron","3' UTR", 'Downstream')
celltype_class = c('Caudate', 'MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb', 'Astro', 
                   'Microglia', 'OPC', 'Oligo')
if(FALSE){
  load(file.path(PROJDIR, 'rdas', 'human_macaque_orthologs_peakList.rda'))
  load(file.path(PROJDIR, 'rdas', 'human_mouse_orthologs_peakList.rda'))
  
  # split the peakList by genome coordinates for annotation
  all_peaks = c(human_macaque_orthologs, human_mouse_orthologs[-1])
  human_peakList = all_peaks[grepl('hgPeaks|Hg|Orth',names(all_peaks))]
  mouse_peaks = all_peaks[['MmPeaks']]
  macaque_peaks = lapply(all_peaks[['rmPeaks']], convertHalChrName, chrOut = 'UCSC',
                         species = 'Macaca_mulatta') %>% GRangesList()
    
  # annotate peak +/- 5kb for promoter region 
  fromTSS = c(-2000,2000)
  human_peakList = parallel::mclapply(human_peakList, annotatePeaks, 
                                      fromTSS = fromTSS, genome = 'hg38', mc.cores = 4)
  macaque_peaks = annotatePeaks(macaque_peaks, fromTSS = fromTSS, genome = 'rheMac8')
  mouse_peaks = annotatePeaks(mouse_peaks, fromTSS = fromTSS, genome = 'mm10')

  names(human_peakList) = group_class
  save(human_peakList, macaque_peaks, mouse_peaks, file = save_fn)
} else {
  load(save_fn)
}

#######################################
## gather matrix of peak annotations ##
human_peaks = do.call(c,lapply(human_peakList, as.list))
annot_df = lapply(human_peaks, function(gr){
  df = gr %>% as.data.frame() %>% group_by(annot) %>%
    summarise (n = n()) %>%
    mutate(freq = n / sum(n),
           annot = gsub('X','',annot) %>% 
             gsub(pattern = '\\.\\.',replacement = "' ",.) %>%
             factor(., annot_class),
           nPeaks = sum(n))
  return(df)
  }) %>% bind_rows(.id = 'id') %>% 
  mutate(celltype = ss(id,'\\.',2) %>% gsub(pattern = 'Consensus','Caudate',.),
         celltype = factor(celltype, celltype_class), 
         group = factor(ss(id,'\\.', 1), group_class)) %>% 
  dplyr::select(-id)

#################################
## make plots for presentation ##
system(paste('mkdir -p', file.path( 'plots')))
height_ppt = 4; width_ppt = 8
height_fig = 4; width_fig = 2.25; font_fig = 5

plot_fn = file.path('plots','Corces2020_caudate_peaks_barplot_ppt.pdf')
pdf(width = width_ppt, height = height_ppt, file = plot_fn)
ggplot(annot_df, aes(x = group, y = n, fill = annot)) +
  geom_bar(stat = 'identity') + scale_fill_jco() +
  xlab('Peak Group (hg, mm, rm, etc.)') + ylab('Number of Peaks') + 
  scale_x_discrete(position = "top") +
  facet_wrap(~celltype, nrow =2, scales = 'free_y',strip.position ='right') + 
  theme_bw(base_size = 10) + theme(legend.position="bottom") + 
  guides(fill = guide_legend(nrow = 2, override.aes = list(size = 2)))+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 0)) +
  theme(legend.text=element_text(size=10), legend.title=element_text(size=10),
        legend.key.height=unit(.75,"line"), legend.key.size = unit(.3,"line"))
dev.off()

plot_fn = file.path('plots','Corces2020_caudate_percent_barplot_ppt.pdf')
pdf(width = width_ppt, height = height_ppt, file = plot_fn)
ggplot(annot_df, aes(x = group, y = freq, fill = annot)) +
  geom_bar(stat = 'identity') + scale_fill_jco() +
  xlab('Peak Group (hg, mm, rm, etc.)') + ylab('Fraction of Peaks') + 
  scale_x_discrete(position = "top") +
  facet_wrap(~celltype, nrow =2, scales = 'free_y',strip.position ='right') + 
  theme_bw(base_size = 10) + theme(legend.position="bottom") + 
  guides(fill = guide_legend(nrow = 2, override.aes = list(size = 2)))+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 0)) +
  theme(legend.text=element_text(size=10), legend.title=element_text(size=10),
        legend.key.height=unit(.75,"line"), legend.key.size = unit(.3,"line"))
dev.off()

plot_fn = file.path('plots','Corces2020_caudate_peaks_barplot_fig.pdf')
pdf(width = width_fig, height = height_fig, file = plot_fn, onefile = F)
pp1 = ggplot(annot_df, aes(x = group, y = n, fill = annot)) +
  geom_bar(stat = 'identity') + scale_fill_jco() +
  xlab('Peak Annotation') + ylab('Number of Peaks') + 
  scale_x_discrete(position = "top") +
  facet_wrap(~celltype, ncol = 1, scales = 'free_y',strip.position ='right') + 
  theme_bw(base_size = font_fig) + theme(legend.position="bottom") + 
  guides(fill = guide_legend(nrow = 2, override.aes = list(size = 2)))+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 0)) +
  theme(legend.text=element_text(size=font_fig), legend.key.height=unit(.75,"line"), 
        legend.title=element_text(size=font_fig), legend.key.size = unit(.3,"line"))

pp2 = ggplot(annot_df, aes(x = group, y = freq, fill = annot)) +
  geom_bar(stat = 'identity') + scale_fill_jco() +
  xlab('Peak Annotation Fraction') + ylab('Peak Annotation Fraction') + 
  scale_x_discrete(position = "top") +
  facet_wrap(~celltype, ncol = 1, scales = 'free_y',strip.position ='right') + 
  theme_bw(base_size = font_fig) + theme(legend.position="bottom") + 
  guides(fill = guide_legend(nrow = 2, override.aes = list(size = 2)))+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 0)) +
  theme(legend.text=element_text(size=font_fig), legend.key.height=unit(.75,"line"), 
        legend.title=element_text(size=font_fig), legend.key.size = unit(.3,"line"))
ggarrange(pp1, pp2, ncol = 2, common.legend = T, align = 'h', legend="bottom")
dev.off()







