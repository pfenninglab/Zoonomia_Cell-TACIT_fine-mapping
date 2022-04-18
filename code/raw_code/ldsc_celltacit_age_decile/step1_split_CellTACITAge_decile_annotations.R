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
library(GenomicFeatures)

#################################
# 1) import the Cell Tacit tracks
DATADIR='data/raw_data/ldsc_celltacit_age_decile'
bed_fn = list.files('data/raw_data/ldsc_caudate_zoonomia/CellTACIT', full.names = T, 
                    pattern = '.CellTACIT.mean.bed.g')
names(bed_fn) = basename(bed_fn) %>% ss('Corces2020.|.mean.bed.gz', 2)

#################################
## 2) break up peaks into deciles
cell_tacit_peaks = lapply(bed_fn, import) %>%
  lapply(function(gr){
    # break up Ages into deciles per cell type
    gr$decile = ntile(gr$score, 10)
    gr$decile_rand = ntile(sample(gr$score), 10)
    return(gr)
  })

table_fn = here(DATADIR, 'tables', 'cell_tacit_decile_breaks.xlsx')
cell_tacit_peaks %>% lapply(as.data.frame) %>%
  rbindlist(idcol = 'celltype') %>% 
  group_by(celltype, decile) %>% 
  summarise(age_min = min(score) %>% signif(digits = 3),
            age_max = max(score) %>% signif(digits = 3)) %>% 
  data.frame() %>% writexl::write_xlsx(table_fn)


for (cell in names(cell_tacit_peaks)){
  ## grab the granges of each cell type, split by decile
  tmpList = split(cell_tacit_peaks[[cell]],
                  cell_tacit_peaks[[cell]]$decile)
  fn = here(DATADIR, 'peaks', 
            paste0('Corces2020.',cell,'.decile',names(tmpList),'.bed.gz'))
  ## export files to bed
  mapply(export, tmpList, fn)
  
  ## get the randomized deciles
  tmpList2 = split(cell_tacit_peaks[[cell]],
                  cell_tacit_peaks[[cell]]$decile_rand)
  fn2 = here(DATADIR, 'peaks', 
            paste0('Corces2020.',cell,'.decile_rand',names(tmpList2),'.bed.gz'))
  ## export files to bed
  mapply(export, tmpList2, fn2)
}



#############################################
## 3) break up peaks into quintile
cell_tacit_peaks = lapply(bed_fn, import) %>%
  lapply(function(gr){
  # break up Ages into quintile per cell type
  gr$quintile = ntile(gr$score, 5)
  gr$decile_rand = ntile(sample(gr$score), 5)
  return(gr)
})

table_fn = here(DATADIR, 'tables', 'cell_tacit_quintile_breaks.xlsx')
cell_tacit_peaks %>% lapply(as.data.frame) %>%
  rbindlist(idcol = 'celltype') %>% 
  group_by(celltype, quintile) %>% 
  summarise(age_min = min(score) %>% signif(digits = 3),
            age_max = max(score) %>% signif(digits = 3)) %>% 
  data.frame() %>% writexl::write_xlsx(table_fn)


DATADIR='data/raw_data/ldsc_celltacit_age_decile'
dir.create('data/raw_data/ldsc_celltacit_age_decile/peaks', showWarnings = F)
for (cell in names(cell_tacit_peaks)){
  ## grab the granges of each cell type, split by quintile
  tmpList = split(cell_tacit_peaks[[cell]],
                  cell_tacit_peaks[[cell]]$quintile)
  fn = here(DATADIR, 'peaks', 
            paste0('Corces2020.',cell,'.quintile',names(tmpList),'.bed.gz'))
  ## export files to bed
  mapply(export, tmpList, fn)
  
  
  ## get the randomized quintiles
  tmpList2 = split(cell_tacit_peaks[[cell]],
                   cell_tacit_peaks[[cell]]$decile_rand)
  fn2 = here(DATADIR, 'peaks', 
             paste0('Corces2020.',cell,'.quintile_rand',names(tmpList2),'.bed.gz'))
  ## export files to bed
  mapply(export, tmpList2, fn2)
}






