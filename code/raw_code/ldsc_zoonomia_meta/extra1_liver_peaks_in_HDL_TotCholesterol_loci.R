#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F, bitmapType='cairo')
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
library(here)
library(rtracklayer)
library(GenomicRanges)

DATADIR=here('data/raw_data','ldsc_zoonomia_meta') 
chain <- file.path("/home/bnphan/resources/liftOver_chainz",
                   'hg19ToHg38.over.chain') %>% import.chain()

### 1) read in LD clumped SNPs
snps_fn = list.files(here(DATADIR, 'FUMA'), pattern = 'snps.txt', recursive = T, 
                     full.names = T)
names(snps_fn) = snps_fn %>% ss('/',10) %>% gsub(pattern = 'FUMA_',replacement = '') %>%
  gsub(pattern = '\\.txt',replacement = '') 
snps_df = snps_fn %>% lapply(fread) %>% rbindlist(idcol = 'label') %>% 
  rename('pos' = 'pos_hg19') 
snps_df = snps_df %>% 
  mutate(name = paste0('chr',chr,':',pos_hg19)) %>% 
  select(uniqID, name) %>% deframe() %>% GRanges() %>% 
  rtracklayer::liftOver(chain = chain) %>%
  GenomicRanges::reduce() %>% as.data.frame() %>% 
  rename('group_name' = 'uniqID', 'start' = 'pos_hg38') %>%
  select(c(uniqID:pos_hg38)) %>% select(-seqnames) %>% 
  right_join(snps_df, by = 'uniqID') %>% 
  relocate(pos_hg38, .after = pos_hg19) %>%
  filter(!is.na(pos_hg38))

### 2) read in LD clumped loci
loci_fn = list.files(here(DATADIR, 'FUMA'), pattern = 'GenomicRiskLoci.txt', recursive = T, 
                     full.names = T)
names(loci_fn) = loci_fn %>% ss('/',10) %>% gsub(pattern = 'FUMA_',replacement = '') %>%
  gsub(pattern = '\\.txt',replacement = '') 
loci_df = loci_fn %>% lapply(fread) %>% rbindlist(idcol = 'label') %>% 
  rename('start' = 'start_hg19','end' = 'end_hg19') 
loci_df = loci_df %>% 
  mutate(name = paste0('chr',chr,':',start_hg19, '-',end_hg19)) %>% 
  select(IndSigSNPs, name) %>% deframe() %>% GRanges() %>% 
  rtracklayer::liftOver(chain = chain) %>%
  GenomicRanges::reduce() %>% as.data.frame() %>% 
  rename('group_name' = 'IndSigSNPs', 'start' = 'start_hg38','end' = 'end_hg38') %>%
  select(c(IndSigSNPs:end_hg38)) %>% select(-seqnames) %>% 
  right_join(loci_df, by = 'IndSigSNPs') %>% 
  relocate(start_hg38, end_hg38, .after = end_hg19) %>%
  mutate(loci = paste0('chr',chr,':',start_hg38,'-',end_hg38)) %>%
  filter(!is.na(start_hg38) | is.na(end_hg38))


### 3) read in the liver peaks and overlap w/ SNPs or loci
liverPeak_bed = '/projects/MPRA/Irene/PGLSOutputs/HumanOfRenamedAllNRSummit_liverModified.bed' %>% import()
mcols(liverPeak_bed) %>% data.frame() %>% pull(name) %>% ss('_', 2) %>% table()

ooSNP = findOverlaps( 
  subject = liverPeak_bed, 
  query = snps_df %>% mutate(tmp = paste0('chr',chr,':',pos_hg38)) %>%
    pull(tmp) %>% GRanges()
)
liverPeakSNP_bed = liverPeak_bed[unique(subjectHits(ooSNP))]
snps_df2 = snps_df[unique(queryHits(ooSNP)),]

ooLoci = findOverlaps( 
  subject = liverPeak_bed, 
  query = loci_df %>% 
    mutate(tmp = paste0('chr',chr,':',start_hg38, '-',end_hg38)) %>%
    pull(tmp) %>% GRanges()
)
liverPeakLoci_bed = liverPeak_bed[unique(subjectHits(ooLoci))]
loci_df2 = loci_df[unique(queryHits(ooLoci)),]

## export to bed files
save_rda = here(DATADIR, 'rdas','liverPeaks_overlapGWAS_SNPorLoci.rda')
dir.create(here(DATADIR, 'rdas'), recursive = T, showWarnings = F)
save(liverPeakSNP_bed, liverPeakLoci_bed, snps_df2, loci_df2, file = save_rda)

dir.create(here(DATADIR, 'liverPeaks'), recursive = T, showWarnings = F)
liverPeakSNP_bed %>% export(con = here(DATADIR, 'liverPeaks', 'liverPeaks_overlapGWAS_SNP.bed'))
liverPeakLoci_bed %>% export(con = here(DATADIR, 'liverPeaks', 'liverPeaks_overlapGWAS_Loci.bed'))

snps_df2 %>% write_tsv(here(DATADIR, 'liverPeaks', 'GWAS_SNP_overlap_liverPeaks.txt'))
loci_df2 %>% write_tsv(here(DATADIR, 'liverPeaks', 'GWAS_Loci_overlap_liverPeaks.txt'))



