#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F, bitmapType='cairo')
library(tidyverse)
library(tidymodels)
library(corrr)
library(broom)
library(data.table)
library(RColorBrewer)
library(rtracklayer)
library(ggh4x)
library(here)
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# to be run in the root github directory
LABEL='ldsc_atac_and_conservation'
PLOTDIR=file.path('figures/exploratory',LABEL)
sapply(c('plots', 'tables', 'rdas'), function(x)
  dir.create(here(PLOTDIR,x), showWarnings = F, recursive = F))
DATADIR=here('data/raw_data',LABEL)

##############################################################################
## 1) import the mean metric values across metric & Cell TACIT Age
bed_fn1 = list.files(here('data/raw_data/ldsc_caudate_zoonomia/CellTACIT'), 
                     full.names = T, pattern = '.CellTACIT.mean.bed.gz')
bed_fn2 = list.files(file.path(DATADIR, 'peak'), full.names = T) %>% 
  str_subset('summary') %>% str_subset('Consensus', negate =T)

## read in the peaks and 
peak_df = c(bed_fn1, bed_fn2) %>% set_names() %>% 
  lapply(import) %>% lapply(as.data.frame) %>% rbindlist(idcol = 'file')

peak_wide_df = peak_df %>% 
  mutate(file = basename(file), 
         celltype = ss(file, '\\.', 2), 
         metric = ss(file, '\\.', 3), 
         metric = gsub('_summary', '', metric)) %>% 
  dplyr::select(-c(file, name)) %>% 
  pivot_wider(names_from = 'metric', values_from = 'score')

## annotate peaks
annot_df <- annotatePeak(GRanges(peak_wide_df), annoDb='org.Hs.eg.db', 
                         TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, 
                         tssRegion = c(-20000, 20000)) %>% 
  as.GRanges() %>% as.data.frame(row.names = seq(length(.))) %>%
  mutate(annotation = ss(annotation, ' '))

peak_wide_df$annotation = annot_df$annotation
annot_type = c('Intron', 'Distal', 'Promoter', 'Exon', 'UTR')
change_type = c("5'", "3'", 'Downstream')

################################################## 
## 2) compute pairwise correlation across all OCRs 
cols = c( '91mam_GERP', '200mam_PhyloP', '43prim_PhastCons','CellTACIT')
rename= c('91 Mammals GERP', '240 Mammals PhyloP', '43 Primate PhastCons', 
          'Cell-TACIT Age') %>% setNames(cols)
cons_corr_df = peak_wide_df %>% 
  mutate(annotation = ifelse(annotation %in% change_type, 'UTR', annotation)) %>% 
  nest(data = -c(celltype, annotation)) %>% 
  mutate(
    tidied = map(data, ~ correlate(.x %>% dplyr::select(all_of(cols)) %>% 
                                     as.matrix(), method = "spearman", quiet =T))
  ) %>% 
  unnest(cols = tidied) %>% dplyr::select(-data) %>% 
  pivot_longer(cols = -c('celltype', 'annotation', 'term'), names_to= 'term2') %>% 
  mutate(term = rename[term],   term = factor(term, rename), 
         term2 = rename[term2], term2 = factor(term2, rename), 
         annotation = factor(annotation, annot_type),
         celltype = ifelse(celltype== 'MSN_SN', "MSN_D1/D2H", celltype),
         celltype = celltype %>% 
           factor(c('MSN_D1', 'MSN_D2', "MSN_D1/D2H", 'INT_Pvalb', 'Astro', 
                    'Microglia', 'OPC', 'Oligo')))

############################################## 
## 3) make the correlation plot of the scores 
plot_fn = here(PLOTDIR,'plots',paste0('conservation_correlation_matrix.All.pdf'))
pdf(plot_fn, height = 10.5, width = 8)
ggplot(cons_corr_df, aes(x = term, y = term2)) +
  geom_tile(aes(fill = value)) + 
  scale_fill_gradientn(colours=brewer.pal(11,"RdYlBu"), 
                       limits = c(-1, 1),'Spearman Correlation') +
  geom_text(aes(label = signif(value, 2)), size = 2)+
  facet_grid(celltype~annotation) + 
  scale_x_discrete(limits = rev) + theme_classic() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1), 
        axis.title = element_blank(), 
        legend.position = 'bottom')
dev.off()

############################################## 
## 3) make the correlation plot of the scores 
plot_fn = here(PLOTDIR,'plots',paste0('conservation_correlation_barplot.All.pdf'))
pdf(plot_fn, height = 10.5, width = 8)
ggplot(cons_corr_df, aes(x = annotation, y = value)) +
  geom_hline(yintercept = 0, color = 'red') +
  geom_hline(yintercept = 0.5, color ='black', linetype = 'dashed') +
  geom_point(aes(fill = term, shape = term)) + 
  scale_fill_brewer(palette = 'Set1', 'Conservation metric') +
  scale_shape_manual(values = c(21:24), 'Conservation metric') +
  facet_nested_wrap(~term2 + celltype, nrow = 4) + 
  ylab('Spearman Correlation') + theme_classic() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank(), 
        legend.position = 'bottom')
dev.off()


