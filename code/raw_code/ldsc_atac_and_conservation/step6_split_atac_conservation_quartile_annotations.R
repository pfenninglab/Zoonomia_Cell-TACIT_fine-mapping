#######################################
### set up libraries and functions ####
# conda activate r4
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F, bitmapType='cairo')
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
library(rcartocolor)
library(here)
library(rtracklayer)
library(GenomicRanges)
library(ChIPseeker)

suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))

#############################################################################
# 1) import the peaks for each track, CellTACIT-Age, GERP, PhyloP, PhastCons
bed_fn = list.files('data/raw_data/ldsc_caudate_zoonomia/CellTACIT', full.names = T, 
                    pattern = '.CellTACIT.mean.bed.g')
names(bed_fn) = basename(bed_fn) %>% ss('Corces2020.|.mean.bed.gz', 2) %>% 
  str_replace('Cell', 'Cell-') %>% paste0('_Age')

bed_fn2 = list.files('data/raw_data/ldsc_atac_and_conservation/peak', full.names = T, 
                    pattern = '_summary.bed') %>% str_subset('main', negate = T)
names(bed_fn2) = basename(bed_fn2) %>% ss('Corces2020.|_summary.bed', 2)
peaks = lapply(c(bed_fn, bed_fn2), import) 
peaks = peaks[sort(names(peaks))]

############################################
## 2) annotate the peaks by distance from TSS
gr = GRangesList(peaks) %>% unlist()
gr$score = NULL
gr = gr[!duplicated(gr$name)]
names(gr) = NULL

change_type = c("5'" = 'Other', "3'"  = 'Other', 'Downstream' = 'Other', 
                'Exon' = 'Other', 'UTR' = 'Other', 
                'Promoter' = 'Promoter', 
                'Intron' = 'Enhancer', 'Distal' = 'Enhancer')

annot_df <- annotatePeak(gr, annoDb='org.Hs.eg.db', 
                         TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, 
                         tssRegion = c(-5000, 5000)) %>% 
  as.GRanges() %>% as.data.frame(row.names = seq(length(.))) %>%
  mutate(annotation = ss(annotation, ' '))
# keep only non-coding distal/intronic peaks
gr$annotation = annot_df$annotation
gr$annotation = change_type[gr$annotation]
gr_list = split(gr, gr$annotation)

## split peaks by their major annotations of Enhancer, Promoter, Other
peaks_list = lapply(peaks, function(x){
  x %>% data.frame() %>% 
    mutate(annotation = gr$annotation[match(name, gr$name)]) %>% 
    group_by(annotation) %>% 
    mutate(quartile = paste0('quartile', ntile(score, 4))) %>% ungroup() %>% 
    mutate(split = paste(annotation, quartile, sep = '.')) %>% 
    split(f = .$split)
})

DATADIR='data/raw_data/ldsc_atac_and_conservation'
dir.create('data/raw_data/ldsc_atac_and_conservation/peaks_quartile', showWarnings = F)
for (lab in names(peaks_list)){
  ## grab the granges of each cell type, split by quartile
  tmpList = peaks_list[[lab]]
  fn = here(DATADIR, 'peaks_quartile', 
            paste0('Corces2020.',lab,'.',names(tmpList),'.bed.gz'))
  ## export files to bed
  parallel::mcmapply(export, tmpList, fn, mc.cores = 12)
}

peaks_list %>% lapply(function(gr_list) lapply(gr_list, function(x) range(x$score)))





