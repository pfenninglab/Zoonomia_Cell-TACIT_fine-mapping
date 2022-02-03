library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(parallel)
library(here)
library(data.table)
library(Biostrings)

options(stringsAsFactors = F)
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
DATADIR='data/raw_data/select_mpra_candidates'
here::i_am('code/raw_code/select_mpra_candidates/step5_assemble_Zoonomia_PPRA.R')

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
DATADIR3 = '/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/data/raw_data/select_mpra_candidates'
xlsx_fn = c(here(DATADIR2, 'tables/Addiction_MPRA_Sequences_Library_20210820.xlsx'), 
         here(DATADIR2, 'tables/PanTissueControl_MPRA_Library_20211115.xlsx'), 
         here(DATADIR3, 'tables/Zoonomia_MSN_low_throughput_MPRA_Library_20220119.xlsx'),
         here(DATADIR3, 'tables/AD_MPRA_barcode.xlsx'))
used_sequences_df = lapply(xlsx_fn, readxl::read_xlsx) %>% rbindlist(fill = TRUE)
used_sequences_df$contributor[is.na(used_sequences_df$contributor)] = 'ZC'
with(used_sequences_df, table(contributor))
#  AJL   ARP   BNP   HSA   HSB   IMK   MMA   MMB    ZC 
# 2750  1024 19429     1     1  1500     1     1 24150 

dim(used_sequences_df)
indKeep = which(! barcodes$seq %in% used_sequences_df$barcode_seq)
barcodes = barcodes[indKeep,]


########################################################
## grab the positives and negatives from the first MPRA
parts_fa = here(DATADIR, 'fasta','MPRAii_enhancer_parts.fa') %>% 
  readDNAStringSet() %>% data.frame() %>% rename('.'= 'value') %>% 
  rownames_to_column('name') %>% deframe()

enhancer_fn = here(DATADIR, 'fasta', c("BNP_PPRA_enhancers_20220126.fa",
                                       'IMK_PPRA_SequencesAll.fa', 
                                       'MEW_VL_PPRA2Sequences.fa'))
enhancer_fa = readDNAStringSet(enhancer_fn)
names(enhancer_fa) = str_trim(names(enhancer_fa))
length(enhancer_fa); mean(width(enhancer_fa))
enhancer_df = enhancer_fa %>% data.frame() %>% rename('.'= 'enhancer_seq') %>% 
  rownames_to_column('enhancer_ID') %>% mutate(count = 1) %>% uncount(count)

## put together the sequences
set.seed(22220202)
idx = sample(nrow(barcodes), length(enhancer_fa))
sequences_df = bind_rows(enhancer_df) %>%
  mutate(barcode_seq = barcodes$seq[idx], barcode_ID = barcodes$barcode_ID[idx],
         full_oligo_seq = paste0(parts_fa['prefix'], enhancer_seq, 
                                 parts_fa['linker'], barcode_seq, 
                                 parts_fa['suffix']),
         seq_ID = paste0(enhancer_ID,'#', barcode_ID), 
         contributor = str_sub(enhancer_ID, end = 3),
         RE_match = str_count(full_oligo_seq, "CGTCTC") + str_count(full_oligo_seq, 'GAGACG'),
         RE_match_enh = str_count(enhancer_seq, "CGTCTC") + str_count(enhancer_seq, 'GAGACG'),
         RE_match_bar = str_count(barcode_seq, "CGTCTC") + str_count(barcode_seq, 'GAGACG')) %>%
  relocate(c(seq_ID, full_oligo_seq), .before = everything())

# make sure no RE matches
print(sum(sequences_df$RE_match == 0 ) == nrow(sequences_df))

## sanity checks, make sure barcodes not duplicated
sequences_df %>% distinct(enhancer_ID, .keep_all = T) %>% pull(contributor) %>% table()
sum(duplicated(c(sequences_df$barcode_seq, used_sequences_df$barcode_seq)))
sum(duplicated(c(sequences_df$full_oligo_seq, used_sequences_df$full_oligo_seq)))

## sanity checks, make sure barcodes not duplicated
enhancer_fa2 = DNAStringSet(sequences_df$full_oligo_seq)
names(enhancer_fa2) = sequences_df$seq_ID
enhancer_fn = here(DATADIR,'fasta/Zoonomia_BNP_IMK_MEW_PPRA_Library_20220202.fasta')
writeXStringSet(enhancer_fa2, enhancer_fn, format = 'fasta',width = 600)

sequences_df %>% writexl::write_xlsx(here(DATADIR, 'tables/Zoonomia_BNP_IMK_MEW_PPRA_Library_20220202.xlsx'))

nrep = c(9,9,12)
wells = paste0(mapply(rep, c('A', 'B', 'C'), nrep) %>% unlist(), 
               sapply(nrep, seq)%>% unlist())

sequences_df2 = sequences_df %>% 
  mutate(`Well Position` = wells) %>%
  dplyr::select(`Well Position`, seq_ID, full_oligo_seq) %>%
  dplyr::rename('Name' = 'seq_ID', 'Sequence' = 'full_oligo_seq') %>%
  mutate_if(is.character, ~gsub(':|#', '.', .))

list(`Zoonomia PPRA` = sequences_df2) %>% 
  writexl::write_xlsx(here(DATADIR, 'tables/Zoonomia_BNP_IMK_MEW_PPRA_Library_20220202_eblock96WPtemplate.xlsx'))