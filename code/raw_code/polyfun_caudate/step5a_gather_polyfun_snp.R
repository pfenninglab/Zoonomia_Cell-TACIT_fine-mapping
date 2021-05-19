ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(here)
library(tidyverse)
library(rcartocolor)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(GenomicRanges)
library(rtracklayer)

CODEDIR='code/raw_code/polyfun_caudate'
DATADIR='data/raw_data/polyfun_caudate'
PLOTDIR='figures/exploratory/polyfun_caudate'
i_am(file.path(CODEDIR, 'step5a_gather_polyfun_snp.R'))

###################################
# read in the gwas phenotypes ##
load(here('data/tidy_data/ldsc_gwas','rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% select( -file) %>% 
  mutate(label = ss(as.character(trait), '_')) %>% 
  inner_join(fread(here('data/tidy_data/ldsc_gwas','gwas_list_sumstats_polyfun.tsv')) %>% select(match),
             by = 'match')

###############################
# read in annotated PIP SNPs ##
annot_fn = here(DATADIR, pheno$match) %>%
  lapply(list.files, pattern = '.top_annot.txt.gz', full.names = T) %>%
  unlist()
names(annot_fn) = basename(annot_fn)
input = annot_fn %>% lapply(fread) %>% rbindlist(idcol = 'file', fill = TRUE) %>%
  as_tibble() 

snps_df = input %>% 
  select(c(file:A2, contains('Caud'), contains('Zoonomia'), contains('Zoonomia'),
           contains("ENCODE3"), contains('synonymous'))) %>%
  mutate(match = file %>% basename() %>% ss('\\.', 1),
         value = paste(match, SNP, A1, A2, sep = ':')) %>%
  left_join(pheno, by = 'match') %>% select(-file)

########################################
# liftOver hg19 SNP positions to hg38 ##
chain <- file.path("/home/bnphan/resources/liftOver_chainz",
                   'hg19ToHg38.over.chain') %>% import.chain()
snpRanges_hg19 = snps_df %>% 
  mutate(name = paste0('chr',CHR,':',BP)) %>% 
  select(value, name) %>% deframe() %>% GRanges()

snpRanges_hg38 = snpRanges_hg19 %>% rtracklayer::liftOver(chain = chain) %>%
  GenomicRanges::reduce() %>% as.data.frame() %>% rename('group_name' = 'value') %>%
  select(c(value:start))

snps_df = snps_df %>% left_join(snpRanges_hg38, by = 'value') %>%
  rename('start' = 'POS_hg38', 'BP' = 'POS_hg19') %>% 
  filter( seqnames %>% as.character() %>% ss('chr', 2) == CHR) %>% 
  relocate(POS_hg38, .after = 'CHR') %>% select(-c(seqnames, value))

########################################################################
## annotate finemapped SNPs w/ phyloP, primate PhastCons, HARS, CHARs ##
top_phyloP_lvls = c('Con.top.0-1%', 'Con.top.1-2%', 'Con.top.2-3%', 'Con.top.3-4%',
                    'Acc.top.3-4%', 'Acc.top.2-3%', 'Acc.top.1-2%', 'Acc.top.0-1%',  
                    'Other')
top_phyloP_cols = c(brewer.pal(4,'Blues'),brewer.pal(4,'Reds'),'#bdbdbd')
names(top_phyloP_cols) = c('Con.top.3-4%', 'Con.top.2-3%', 'Con.top.1-2%', 'Con.top.0-1%', 
                           'Acc.top.3-4%', 'Acc.top.2-3%', 'Acc.top.1-2%', 'Acc.top.0-1%',  
                           'Other')

# colors for 43primate phastCons
top_phastCons_lvls = c('PhastCons.top.0-1%', 'PhastCons.top.1-2%', 'PhastCons.top.2-3%', 
                       'PhastCons.top.3-4%', 'PhastCons.top.4-5%', 'Other')
top_phastCons_cols = c(rev(brewer.pal(5,'PuBuGn')),'#bdbdbd')
names(top_phastCons_cols) = top_phastCons_lvls

# colors for 3 ENCODE3 cCRE annotations
ENCODE3_cCRE_lvls = c('PLS', 'pELS', 'dELS', 'Other')
ENCODE3_cCRE_cols = c(brewer.pal(3,'Dark2'),'#bdbdbd')
names(ENCODE3_cCRE_cols) = ENCODE3_cCRE_lvls

# colors for 2 HAR and CHAR annotations
HAR_lvls = c('HAR1500bp', 'CHAR1500bp', 'Other')
HAR_cols = c(brewer.pal(3,'Accent'))
names(HAR_cols) = HAR_lvls

snps_df = snps_df %>% mutate(
  # Group 241mammals phyloP
  top_phyloP = case_when(
    `Zoonomia_phyloPaccl.241mam.top0-1%` == 1 ~ 'Acc.top.0-1%',
    `Zoonomia_phyloPaccl.241mam.top1-2%` == 1 ~ 'Acc.top.1-2%',
    `Zoonomia_phyloPaccl.241mam.top2-3%` == 1 ~ 'Acc.top.2-3%',
    `Zoonomia_phyloPaccl.241mam.top3-4%` == 1 ~ 'Acc.top.3-4%',
    `Zoonomia_phyloPcons.241mam.top0-1%` == 1 ~ 'Con.top.0-1%',
    `Zoonomia_phyloPcons.241mam.top1-2%` == 1 ~ 'Con.top.1-2%',
    `Zoonomia_phyloPcons.241mam.top2-3%` == 1 ~ 'Con.top.2-3%',
    `Zoonomia_phyloPcons.241mam.top3-4%` == 1 ~ 'Con.top.3-4%',
    TRUE ~ 'Other'),
  top_phyloP = factor(top_phyloP, top_phyloP_lvls), 
  # Group primate PhastCons
  top_phastCons = case_when(
    `Zoonomia_phastCons.43prim.top0-1%` == 1 ~ 'PhastCons.top.0-1%',
    `Zoonomia_phastCons.43prim.top1-2%` == 1 ~ 'PhastCons.top.1-2%',
    `Zoonomia_phastCons.43prim.top2-3%` == 1 ~ 'PhastCons.top.2-3%',
    `Zoonomia_phastCons.43prim.top3-4%` == 1 ~ 'PhastCons.top.3-4%',
    `Zoonomia_phastCons.43prim.top4-5%` == 1 ~ 'PhastCons.top.4-5%',
    TRUE ~ 'Other'),
  top_phastCons = factor(top_phastCons, top_phastCons_lvls), 
  # Group ENCODE3 cCREs
  cCRE_group = case_when(
    ENCODE3.dELS == 1 ~ 'dELS',
    ENCODE3.pELS == 1 ~ 'pELS',
    ENCODE3.PLS == 1 ~ 'PLS',
    TRUE ~ 'Other'), 
  cCRE_group = factor(cCRE_group, ENCODE3_cCRE_lvls))

####################################
## save table of fine-mapped SNPs ##
dir.create(here('data/raw_data/polyfun_caudate/rdas'), showWarnings = F)
poly_fn = here('data/raw_data/polyfun_caudate/rdas',
               'polyfun_caudate_finemapped_snps_20210518.rds')
saveRDS(snps_df, file = poly_fn)

#########################
## write SNPs to fasta ##
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
genome <- BSgenome.Hsapiens.UCSC.hg19.masked
dir.create(here('data/raw_data/polyfun_caudate/fasta'), showWarnings = F)

##########################################
# get fasta sequences for effect allele
peakLength = 501
offset <- (peakLength-1)/2 # make for 501 bp sequence
seqEffect <- xscat(getSeq(genome,paste0('chr',snps_df$CHR),snps_df$POS_hg19-offset,snps_df$POS_hg19-1),
                   snps_df$A1,
                   getSeq(genome,paste0('chr',snps_df$CHR),snps_df$POS_hg19+1,snps_df$POS_hg19+offset),
                   sep='')
names(seqEffect) = with(snps_df, paste(SNP, A1, A2, 'Effect', sep = ':'))
effect_fasta = here('data/raw_data/polyfun_caudate/fasta', 
                    paste0('polyfun_caudate_finemapped_snps_',peakLength,'_effect.fasta'))
writeXStringSet(seqEffect, effect_fasta, format = 'fasta',width = peakLength)

#######################################
# get fasta sequences non effect allele
seqNonEffect <- xscat(getSeq(genome,paste0('chr',snps_df$CHR),snps_df$POS_hg19-offset,snps_df$POS_hg19-1),
                      snps_df$A2,
                      getSeq(genome,paste0('chr',snps_df$CHR),snps_df$POS_hg19+1,snps_df$POS_hg19+offset),
                      sep='')
names(seqNonEffect) = with(snps_df, paste(SNP, A1, A2, 'NonEffect', sep = ':'))

nonEff_fasta = here('data/raw_data/polyfun_caudate/fasta', 
                    paste0('polyfun_caudate_finemapped_snps_',peakLength,'_nonEffect.fasta'))
writeXStringSet(seqNonEffect, nonEff_fasta, format = 'fasta',width = peakLength)



#########################################################
## find the finemapped SNPs overlapping CHARs and HARs ##
bed_fn = list.files(path = here('/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/phyloP/phyloP_cutoff_regions'),
                    pattern = '.bed.gz', full.names = T) %>%
  grep(pattern = 'HAR',value = T) %>% grep(pattern = '20210402|20210416',value = T)
names(bed_fn) = basename(bed_fn) %>% ss('\\_20', 1) 
bed_fn = bed_fn[! grepl('flanking', bed_fn)]
bed_gr = bed_fn %>% lapply(import)

overLapsHARexact = bed_gr %>%
  map( ~ countOverlaps(subject = .x, query = GRanges(paste0('chr',snps_df$CHR,':',snps_df$POS_hg38)))) %>%
  bind_cols(snps_df, . ) %>% mutate(group = 'SNP_in_C/HAR') %>%
  relocate(group, .after = A2)

bed_gr2 = bed_gr %>% lapply(function(gr) {
  start(gr) = round((start(gr) + end(gr)) / 2) 
  end(gr) = start(gr)
  return(gr)
})

overLapsHAR750bp = bed_gr2 %>%
  map( ~ countOverlaps(subject = .x, query = GRanges(paste0('chr',snps_df$CHR,':',snps_df$POS_hg38)),
                       maxgap = 750)) %>%
  bind_cols(snps_df, . ) %>% mutate(group = 'SNP750bp_from_C/HAR_center') %>%
  relocate(group, .after = A2)

## combine the two overlapping HARs/CHARs groups
overLapsHAR = overLapsHARexact %>% 
  bind_rows(overLapsHAR750bp) %>% filter_at(vars(contains('HAR')), any_vars(.==1)) %>% 
  mutate(CHAR = paste0('chr',CHR,':',POS_hg38) %>% sapply(function(x){
    tmp = bed_gr2[['zoonomia_CHARs']]
    oo = findOverlaps(subject =  GRanges(x), query = tmp, maxgap = 750)
    return( unique(mcols(tmp)$name[queryHits(oo)])[1])
  }),
  HAR = paste0('chr',CHR,':',POS_hg38) %>% sapply(function(x){
    tmp = bed_gr2[['zoonomia_HARs']]
    oo = findOverlaps(subject =  GRanges(x), query = tmp, maxgap = 750)
    return( unique(mcols(tmp)$name[queryHits(oo)])[1])
  })) %>%
  relocate(c(CHAR, HAR, match:reference), .after = A2) %>%
  mutate_if(is.character,replace_na,'')

overLapsHAR %>% as.data.frame() %>% select(c(CHR, POS_hg38, SNP, CHAR, HAR))

#########################################################
## output the HAR/CHAR overlapped SNPs
dir.create(here('data/raw_data/polyfun_caudate/tables'), showWarnings = F)
har_tsv = here('data/raw_data/polyfun_caudate/tables',
               'polyfun_caudate_finemapped_snps_overlap_HAR_CHAR_20210518.tsv')
write_tsv(overLapsHAR, file = har_tsv)
