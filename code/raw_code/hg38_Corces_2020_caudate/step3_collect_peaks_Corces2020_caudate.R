# to be run in the root github directory
LABEL='cnn_enhancer_ortholog'
setwd('code/raw_code/cnn_enhancer_ortholog')
PROJDIR=file.path('../../../data/raw_data',LABEL)

### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(rtracklayer))
suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))
source('../hal_scripts/narrowPeakFunctions.R')
source('../hal_scripts/gen_enh_ortholog_sets.R')

## parameters for annotation and filtering
fromTSS = c(-20000,20000)
include = c('Distal.Intergenic','Intron')

## chromosomal splits
testSet = c('chr1','chr2')
folds = list(fold1 = c('chr6', 'chr13', 'chr21'),
             fold2 = c('chr7', 'chr14', 'chr18'),
             fold3 = c('chr11','chr17', 'chrX'),
             fold4 = c('chr9', 'chr12'), 
             fold5 = c('chr8', 'chr10'))

########################
## get the hg38 peaks ##
genome = 'hg38'
load(file.path(PROJDIR, 'rdas', paste('caudate_nonEnh_sequences', genome,'.rda', sep = '.')))
bias_fn = file.path(PROJDIR,'fasta',paste(genome, names(negativeSet), 'biasAway10x.fa',sep = '_'))
names(bias_fn) = names(negativeSet)
bias_seq = lapply(bias_fn, readDNAStringSet)
bias_gr = lapply(bias_seq, function(seq){
  gr = paste(ss(names(seq),' ', 1), ss(names(seq),'_', 2), sep = ':') %>% GRanges()
  gr = resize(gr, width = 501)
  names(gr) = mcols(gr)$name = names(seq)
  return(gr)
}) %>% GRangesList()

# annotate the GC-matched negatives, exclude promoters, exonic
bias_gr = annotatePeaks(bias_gr, fromTSS = fromTSS, genome = genome)
biasFilt_gr = filterPeaks(bias_gr, include = include)
biasFilt_gr = mapply(getNonEnhOrthPeaks, inPeaks = biasFilt_gr, 
                          excludePeaks = human_pos_list)
lengths(biasFilt_gr) / 1000
lengths(biasFilt_gr) / lengths(bias_gr)

negativeSplit_list = lapply(folds, function(fold){
  lapply(biasFilt_gr, splitPeakSet, testSet = testSet, validSet = fold)
})


# export postive sequences to summit-centered 501bp fasta file
split = names(negativeSplit_list[[1]][[1]])
system(paste('mkdir -p',  file.path(PROJDIR, 'fasta')))
for(cell in names(negativeSet)){
  for(fold in names(folds)){
    # write the negatives
    neg_fasta_fn = file.path(PROJDIR, 'fasta', 
                             paste(genome, cell, fold, split, 'biasAway10x.fa', sep = '_'))
    negSeq = lapply(negativeSplit_list[[fold]][[cell]], function(gr) bias_seq[[cell]][names(gr)])
    tmp = mapply(writeXStringSet, x = negSeq, filepath = neg_fasta_fn )
  }}






########################
## get the mm10 peaks ##
genome = 'mm10'
load(file.path(PROJDIR, 'rdas', paste('caudate_nonEnh_sequences', genome,'.rda', sep = '.')))
bias_fn = file.path(PROJDIR,'fasta',paste(genome, names(negativeSet), 'biasAway10x.fa',sep = '_'))
names(bias_fn) = names(negativeSet)
bias_seq = lapply(bias_fn, readDNAStringSet)
bias_gr = lapply(bias_seq, function(seq){
  gr = ss(names(seq),':', 1) %>% gsub(pattern = '_',replacement = ':') %>% GRanges()
  gr = resize(gr, width = 501)
  names(gr) = mcols(gr)$name = names(seq)
  return(gr)
}) %>% GRangesList()

# annotate the GC-matched negatives, exclude promoters, exonic
bias_gr = annotatePeaks(bias_gr, fromTSS = fromTSS, genome = genome)
biasFilt_gr = filterPeaks(bias_gr, include = include)
biasFilt_gr = mapply(getNonEnhOrthPeaks, inPeaks = biasFilt_gr, 
                     excludePeaks = mouse_pos_list)

# get the hg38 chromosomes
chainFile =file.path("/home/bnphan/resources/liftOver_chainz", paste0(genome,'ToHg38.over.chain'))
biasFilt_gr = mapply(getLiftedChr, p = biasFilt_gr, chainFile = chainFile) %>% GRangesList()
table(seqnames(biasFilt_gr[[1]]), biasFilt_gr[[1]]$col)
negativeSplit_list = lapply(folds, function(fold){
  lapply(biasFilt_gr, splitPeakSet, testSet = testSet, validSet = fold, useCol = 'col')
})

lengths(biasFilt_gr) / 1000
lengths(biasFilt_gr) / lengths(bias_gr)

# export postive sequences to summit-centered 501bp fasta file
split = names(negativeSplit_list[[1]][[1]])
system(paste('mkdir -p',  file.path(PROJDIR, 'fasta')))
for(cell in names(negativeSet)){
  for(fold in names(folds)){
    # write the negatives
    neg_fasta_fn = file.path(PROJDIR, 'fasta', 
                             paste(genome, cell, fold, split, 'biasAway10x.fa', sep = '_'))
    negSeq = lapply(negativeSplit_list[[fold]][[cell]], function(gr) bias_seq[[cell]][names(gr)])
    tmp = mapply(writeXStringSet, x = negSeq, filepath = neg_fasta_fn )
  }}





############################
## get the rheMac10 peaks ##
genome = 'rheMac10'
load(file.path(PROJDIR, 'rdas', paste('caudate_nonEnh_sequences', genome,'.rda', sep = '.')))
bias_fn = file.path(PROJDIR,'fasta',paste(genome, names(negativeSet), 'biasAway10x.fa',sep = '_'))
names(bias_fn) = names(negativeSet)
bias_seq = lapply(bias_fn, readDNAStringSet)
bias_gr = lapply(bias_seq, function(seq){
  tmp = names(seq)
  tmp = tmp[!is.na(as.numeric(ss(tmp, '_', 2)))]
  gr = tmp %>% gsub(pattern = '_',replacement = ':') %>% GRanges()
  gr = resize(gr, width = 501)
  names(gr) = mcols(gr)$name = tmp
  return(gr)
}) %>% GRangesList()

# annotate the GC-matched negatives, exclude promoters, exonic
bias_gr = annotatePeaks(bias_gr, fromTSS = fromTSS, genome = genome)
biasFilt_gr = filterPeaks(bias_gr, include = include)
biasFilt_gr = mapply(getNonEnhOrthPeaks, inPeaks = biasFilt_gr, 
                     excludePeaks = mouse_pos_list)

# get the hg38 chromosomes
chainFile =file.path("/home/bnphan/resources/liftOver_chainz", paste0(genome,'ToHg38.over.chain'))
biasFilt_gr = mapply(getLiftedChr, p = biasFilt_gr, chainFile = chainFile) %>% GRangesList()
table(seqnames(biasFilt_gr[[1]]), biasFilt_gr[[1]]$col)
negativeSplit_list = lapply(folds, function(fold){
  lapply(biasFilt_gr, splitPeakSet, testSet = testSet, validSet = fold, useCol = 'col')
})

lengths(biasFilt_gr) / 1000
lengths(biasFilt_gr) / lengths(bias_gr)

# export postive sequences to summit-centered 501bp fasta file
split = names(negativeSplit_list[[1]][[1]])
system(paste('mkdir -p',  file.path(PROJDIR, 'fasta')))
for(cell in names(negativeSet)){
  for(fold in names(folds)){
    # write the negatives
    neg_fasta_fn = file.path(PROJDIR, 'fasta', 
                             paste(genome, cell, fold, split, 'biasAway10x.fa', sep = '_'))
    negSeq = lapply(negativeSplit_list[[fold]][[cell]], function(gr) bias_seq[[cell]][names(gr)])
    tmp = mapply(writeXStringSet, x = negSeq, filepath = neg_fasta_fn )
  }}


