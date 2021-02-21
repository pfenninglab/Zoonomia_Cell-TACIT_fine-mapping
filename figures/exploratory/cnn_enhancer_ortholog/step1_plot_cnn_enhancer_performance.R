suppressMessages(library(ArchR))
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
setwd('code/raw_code/hg38_Corces_2020_caudate')
