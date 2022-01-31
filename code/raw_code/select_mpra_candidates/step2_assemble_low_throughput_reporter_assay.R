library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(parallel)
library(here)
library(Biostrings)

options(stringsAsFactors = F)
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
DATADIR='data/raw_data/select_mpra_candidates'
here::i_am('code/raw_code/select_mpra_candidates/step2_assemble_low_throughput_reporter_mpra.R')

################################
## grab the 16-1 freebarcodes
set.seed(20210820)
barcodes = read_tsv('~/src/freebarcodes/barcodes/barcodes16-1.txt', col_names = 'seq') %>%
  mutate(barcode_ID = paste0('free16-1.',seq(n())), 
         tmp = str_count(seq, "CGTCTC") + str_count(seq, 'GAGACG')) %>%
  filter(tmp == 0) %>% dplyr::select(-tmp)
barcodes = barcodes[sample(nrow(barcodes)), ]

barcddes_fa = setNames(DNAStringSet(barcodes$seq), barcodes$barcode_ID)
barcodes_fn = here('~/src/freebarcodes/barcodes/free_barcodes16-1_with_names.fasta')
writeXStringSet(barcddes_fa, barcodes_fn, format = 'fasta',width = 300)


DATADIR2 = '/projects/pfenninggroup/machineLearningForComputationalBiology/Addiction_MPRA/data/tidy_data/assemble_the_array'
xlsx_fn = c(here(DATADIR2, 'tables/Addiction_MPRA_Sequences_Library_20210820.xlsx'), 
         here(DATADIR2, 'tables/PanTissueControl_MPRA_Library_20211115.xlsx'))
used_sequences_df = lapply(xlsx_fn, readxl::read_xlsx) %>% bind_rows()
with(used_sequences_df, table(contributor))
#  AJL   ARP   BNP   IMK 
# 2750  1024 19429  1500
indKeep = which(! barcodes$seq %in% used_sequences_df$barcode_seq)
barcodes = barcodes[indKeep,]

########################################################
## grab the positives and negatives from the first MPRA
parts_fa = here(DATADIR, 'fasta','MPRAii_enhancer_parts.fa') %>% 
  readDNAStringSet() %>% data.frame() %>% rename('.'= 'value') %>% 
  rownames_to_column('name') %>% deframe()

enhancer_fn = here(DATADIR, 'fasta', c("Zoonomia_MSN_low_throughput_reporter_enhancers.fasta"))

enhancer_fa = readDNAStringSet(enhancer_fn)
length(enhancer_fa); mean(width(enhancer_fa))
enhancer_df = enhancer_fa %>% data.frame() %>% rename('.'= 'enhancer_seq') %>% 
  rownames_to_column('enhancer_ID') %>% mutate(count = 1) %>% uncount(count)

## put together the sequences
sequences_df = bind_rows(enhancer_df) %>%
  mutate(barcode_seq = barcodes$seq[seq(n()) + 20], barcode_ID = barcodes$barcode_ID[seq(n()) + 20],
         full_oligo_seq = paste0(parts_fa['prefix'], enhancer_seq, 
                                 parts_fa['linker'], barcode_seq, 
                                 parts_fa['suffix']),
         seq_ID = paste0(enhancer_ID,'#', barcode_ID), 
         contributor = str_sub(enhancer_ID, end = 3),
         RE_match = str_count(full_oligo_seq, "CGTCTC") + str_count(full_oligo_seq, 'GAGACG')) %>%
  relocate(c(seq_ID, full_oligo_seq), .before = everything())

## sanity checks, make sure barcodes not duplicated
sequences_df %>% distinct(enhancer_ID, .keep_all = T) %>% pull(contributor) %>% table()
sum(duplicated(c(sequences_df$barcode_seq, used_sequences_df$barcode_seq)))
sum(duplicated(c(sequences_df$full_oligo_seq, used_sequences_df$full_oligo_seq)))

## sanity checks, make sure barcodes not duplicated
enhancer_fa2 = DNAStringSet(sequences_df$full_oligo_seq)
names(enhancer_fa2) = sequences_df$seq_ID
enhancer_fn = here(DATADIR,'fasta/Zoonomia_MSN_low_throughput_MPRA_Library_20220119.fasta')
writeXStringSet(enhancer_fa2, enhancer_fn, format = 'fasta',width = 300)

sequences_df %>% writexl::write_xlsx(here(DATADIR, 'tables/Zoonomia_MSN_low_throughput_MPRA_Library_20220119.xlsx'))


