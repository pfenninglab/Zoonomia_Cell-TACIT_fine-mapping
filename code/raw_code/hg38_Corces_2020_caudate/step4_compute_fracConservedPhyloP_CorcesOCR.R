#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F, bitmapType='cairo')
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
library(here)
library(rtracklayer)

## read in Corces caudate OCR
fn= list.files(here('data/raw_data/hg38/Corces_2020/peak'), 
               pattern = 'Corces2020_caudate\\.', full.names = T) %>%
  grep(pattern = 'Consensus|All', value = T, invert = T)
names(fn) = basename(fn) %>% ss('\\.', 2)

peakList = fn %>% lapply(import) %>% GRangesList()


## read in 241 mammals phyloP >FDR5% threshold
phyloP_FDR5 = '/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/phyloP/phyloP_cutoff_regions/phyloPcons.241mam.fdr.05.hg38.bed.gz'
phyloP_FDR5_gr = phyloP_FDR5 %>% import()
phyloP_FDR5_rle = coverage(phyloP_FDR5_gr)

## get the percent of bases w/ phyloP > 5% FDR cutoff
gr= peakList[[1]]
peakList_cons = lapply(peakList, function(gr){
  tmp = intersect(phyloP_FDR5_gr, gr)
  oo = findOverlaps(subject = gr, query = tmp)
  widthList = split(width(tmp[queryHits(oo)]), subjectHits(oo))

  gr$score = 0
  gr$score[as.numeric(names(widthList))] = 
    sapply(widthList, sum) 
  gr$score = gr$score/width(gr)
  return(gr)
}) %>% GRangesList()

save_fn = here('data/raw_data/hg38/Corces_2020/rdas', 
               'Corces2020_caudate.peakList.phyloPconsFrac.rds')
saveRDS(peakList_cons, file = save_fn)



