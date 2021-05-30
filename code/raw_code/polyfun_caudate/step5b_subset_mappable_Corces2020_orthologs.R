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
source(here('code/raw_code/hal_scripts/gen_enh_ortholog_sets.R'))

i_am(file.path(CODEDIR, 'step5b_subset_mappable_Corces2020_orthologs.R'))

####################################
## 1) read in table of fine-mapped SNPs ##
dir.create(here('data/raw_data/polyfun_caudate/rdas'), showWarnings = F)
poly_fn = here('data/raw_data/polyfun_caudate/rdas',
               'polyfun_caudate_finemapped_snps_20210518.rds')
snps_df = readRDS(file = poly_fn) %>%
  # remove SNPs in coding regions
  filter_at(vars(contains('Coding_UCSC_')), all_vars(. == 0)) %>% 
  relocate(c(trait,match), .before = everything())

snps_gr = snps_df %>% mutate(CHR = paste0('chr',CHR),
                             name = paste(CHR,POS_hg38, sep = ':'),
                             value = paste(CHR,POS_hg38, A1, A2, sep = ':')) %>%
  select(c(value, name)) %>% deframe() %>% GRanges()

# read in the Zoonomia tree and species list
rda_fn = here('data/tidy_data/Zoonomia_data', 
              'rdas','200_Mammals_Genome_Information.rda')
load(file = rda_fn)
col_clade = df %>% select(Clade, col_clade)%>% filter(!duplicated(Clade)) %>% deframe()
col_order = df %>% select(Order, col_order) %>% filter(!duplicated(Order)) %>% deframe()
celltypes = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb', 'Astro', 'Microglia', 'OPC', 'Oligo')


################################################################
# 2) subset OCR regions overlapping finemapped SNP in any OCR peaks
PROJDIR=here('data/raw_data/ldsc_caudate_zoonomia')
MODEL_TYPE = 'hgRmMm_nonCelltypeNonEnhBiasAway10x'
calib_out_fn = here(here('data/raw_data/cnn_enhancer_ortholog'),
                    paste('rdas/Caudate',MODEL_TYPE,'pos_calibration_ecdf.rds', sep = '_'))
model_calibration = readRDS(file = calib_out_fn)

cell = 'MSN_D1'
snpInPeaks_df = data.frame()
for( cell in celltypes){
    outPeaks_rds = here(PROJDIR,'rdas',paste0('Corces2020.',cell, '.allPeaks.avgCNN.predictions.rds'))
    outCalib_tsv = here(PROJDIR,'tables',paste0('Corces2020.',cell, '.allPeaksCalibWithSNPs.avgCNN.predictions.txt.gz'))
    outCalib_rds = here(PROJDIR,'rdas',paste0('Corces2020.',cell, '.allPeaksCalibWithSNPs.avgCNN.predictions.rds'))
    df_peaks = readRDS(outPeaks_rds) %>% filter(grepl('hg38', name)) %>%
      mutate_if(is.numeric, ~ model_calibration[[cell]](.)) %>% 
      mutate_if(is.numeric, ~ ifelse(. == Inf, 1, .))
  
    peaks_gr = df_peaks %>% mutate(name= gsub('^hg38:|:250$', '', name)) %>%
      pull(name) %>% GRanges()

    oo = findOverlaps(query = snps_gr, subject = peaks_gr)

    ## get distance between peak start and SNP
    tmp_gr = snps_gr[unique(queryHits(oo))]
    oo2 = findOverlaps(query = tmp_gr, subject = peaks_gr)
    tmp_df = snps_df[unique(queryHits(oo)),] %>% 
      mutate(celltype = cell, 
             offset = POS_hg38 - start(peaks_gr)[subjectHits(oo2)], 
             peakNames = df_peaks$name[subjectHits(oo2)]) %>% 
      relocate(c(celltype, peakNames), .after = SNP)
    snpInPeaks_df = tmp_df %>% bind_rows(snpInPeaks_df)
    
    # write out the predictions over all mapped peaks
    # if(any(!file.exists(c(outCalib_tsv, outCalib_rds)))){
    write_tsv(df_peaks[unique(subjectHits(oo)),], file = outCalib_tsv)
    saveRDS(df_peaks[unique(subjectHits(oo)),], outCalib_rds)
    # }
}

allPeaksSNPs_fn = here('data/raw_data/polyfun_caudate/rdas',
                       'polyfun_caudate_finemapped_snpsInAllPeaks_20210518.rds')
allPeaksSNPs_tsv = here('data/raw_data/polyfun_caudate/tables',
                       'polyfun_caudate_finemapped_snpsInAllPeaks_20210518.tsv')
if(!all(file.exists(c(allPeaksSNPs_fn, allPeaksSNPs_tsv)))){
snpInPeaks_df %>% count(celltype)
write_tsv(snpInPeaks_df, file = allPeaksSNPs_tsv)
saveRDS(snpInPeaks_df, file = allPeaksSNPs_fn)
} else {
snpInPeaks_df = allPeaksSNPs_fn %>% readRDS()
}


#########################
## write SNPs to fasta ##
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
genome <- BSgenome.Hsapiens.UCSC.hg19.masked
dir.create(here('data/raw_data/polyfun_caudate/fasta'), showWarnings = F)

##########################################
# get fasta sequences for effect allele
peakLength = 501

seqEffectList = snpInPeaks_df %>% mutate(celltype = factor(celltype, celltypes)) %>% 
  split(snpInPeaks_df$celltype) %>% 
  lapply(function(x){
    tmp = with(x, xscat(getSeq(genome,paste0('chr',CHR),POS_hg19-offset,POS_hg19-1),
                A1, getSeq(genome,paste0('chr',CHR),POS_hg19+1,POS_hg19-offset+ peakLength-1),sep=''))
    names(tmp) = with(x, paste(name, 'Effect', sep = ':'))
    return(tmp)
  })

effect_fasta = here('data/raw_data/polyfun_caudate/fasta', 
                    paste0('polyfun_caudate_finemapped_snps_', peakLength, 
                           '_effect_',names(seqEffectList),'.fasta'))
out = map2( .x = seqEffectList, .y = effect_fasta, 
     ~ writeXStringSet(.x, .y, format = 'fasta',width = peakLength))


seqNonEffList = snpInPeaks_df %>% mutate(celltype = factor(celltype, celltypes)) %>% 
  split(snpInPeaks_df$celltype) %>% 
  lapply(function(x){
    tmp = with(x, xscat(getSeq(genome,paste0('chr',CHR),POS_hg19-offset,POS_hg19-1),
                        A2, 
                        getSeq(genome,paste0('chr',CHR),POS_hg19+1,POS_hg19-offset+ peakLength-1),sep=''))
    names(tmp) = with(x, paste(name, 'NonEffect', sep = ':'))
    return(tmp)
  })
nonEff_fasta = here('data/raw_data/polyfun_caudate/fasta', 
                    paste0('polyfun_caudate_finemapped_snps_', peakLength, 
                           '_nonEffect_',names(seqNonEffList),'.fasta'))
out = map2( .x = seqNonEffList, .y = nonEff_fasta, 
      ~ writeXStringSet(.x, .y, format = 'fasta',width = peakLength))




#####################################################################
# 3) subset OCR regions overlapping finemapped SNP in enhancer OCR peaks
snpInPeaks_df = data.frame()
for( cell in celltypes){
  outEnh_rds = here(PROJDIR,'rdas',paste0('Corces2020.',cell, '.enhPeaks.avgCNN.predictions.rds'))
  outCalib_tsv = here(PROJDIR,'tables',paste0('Corces2020.',cell, '.enhPeaksCalibWithSNPs.avgCNN.predictions.txt.gz'))
  outCalib_rds = here(PROJDIR,'rdas',paste0('Corces2020.',cell, '.enhPeaksCalibWithSNPs.avgCNN.predictions.rds'))
    df_peaks = readRDS(outEnh_rds) %>% filter(grepl('hg38', name)) %>%
      mutate_if(is.numeric, ~ model_calibration[[cell]](.)) %>% 
      mutate_if(is.numeric, ~ ifelse(. == Inf, 1, .))
    
    peaks_gr = df_peaks %>% mutate(name= gsub('^hg38:|:250$', '', name)) %>%
      pull(name) %>% GRanges()
    oo = findOverlaps(subject = peaks_gr, query = snps_gr)
    
    snpInPeaks_df = snps_df[unique(queryHits(oo)),] %>% mutate(celltype = cell) %>%
      bind_rows(snpInPeaks_df)
    
  # write out the predictions over all mapped peaks
  if(any(!file.exists(c(outCalib_tsv, outCalib_rds)))){
    write_tsv(df_peaks[unique(subjectHits(oo)),], file = outCalib_tsv)
    saveRDS(df_peaks[unique(subjectHits(oo)),], outCalib_rds)
  }
}

enhPeaksSNPs_fn = here('data/raw_data/polyfun_caudate/rdas',
                       'polyfun_caudate_finemapped_snpsInEnhPeaks_20210518.rds')
if(!file.exists(enhPeaksSNPs_fn)){
  snpInPeaks_df = snpInPeaks_df %>% filter(!is.na(celltype))
  snpInPeaks_df %>% count(celltype)
  saveRDS(snpInPeaks_df, file = enhPeaksSNPs_fn)
}



#####################################################################
# 4) summary of how many SNPs in peaks or peaks overlapping SNPs
### all peaks w/ finemapped SNPs
dfList_allPeaks = 
  here(PROJDIR,'rdas',
       paste0('Corces2020.',celltypes, '.allPeaksCalibWithSNPs.avgCNN.predictions.rds')) %>%
  lapply(readRDS)
names(dfList_allPeaks) = celltypes
sapply(dfList_allPeaks, nrow)

### enhancer peaks w/ finemapped SNPs
dfList_enhPeaks = 
  here(PROJDIR,'rdas',
       paste0('Corces2020.',celltypes,'.enhPeaksCalibWithSNPs.avgCNN.predictions.rds')) %>%
  lapply(readRDS)
names(dfList_enhPeaks) = celltypes
sapply(dfList_enhPeaks, nrow)

## percent of enhancer peaks of all peaks
sapply(dfList_enhPeaks, nrow) / sapply(dfList_allPeaks, nrow) * 100


