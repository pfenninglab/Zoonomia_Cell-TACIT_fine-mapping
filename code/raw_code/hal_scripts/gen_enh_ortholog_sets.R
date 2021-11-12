#############################################
## by BaDoi Phan February 7, 2021
## Sets of functions that do the following: 
##
## Positive sets (OCR detected in tissue of source species):
## 1) use gene annotation of source species to label peaks
## 2) filter out exons, promoter (<20kb from TSS)

## Negative sets halpered to a target species:
## 1) use gene annotation of target species to label peaks
## 2) filter out exons, promoter (<20kb from TSS)
## 3) exclude overlap w/ any OCR peak found for tissue in target species

## current genomes: hg38, mm10, rheMac10, rn6
supported_genomes = c('hg38', 'mm10', 'rheMac10','rheMac8', 'rn6')
require(GenomicRanges) # for working w/ genomic regiosn
require(ChIPseeker) # for annotating regions

# for ease of not in
`%ni%` <- Negate(`%in%`)

annotatePeaks <- function(peaks, genome = 'hg38', fromTSS = c(-5000,5000)){
  txdb = getTxDb(genome)
  annodb = getAnnoDb(genome)
  print(paste('Using TSS boundaries from',fromTSS[1], 'to',fromTSS[2]))
  if(class(peaks) %in% c('CompressedGRangesList', 'GRangesList')){
    ret = lapply(peaks, annotatePeak, TxDb=txdb, annoDb=annodb, tssRegion=fromTSS)
    ret = lapply(ret, function(x){
      ret = as.GRanges(x)
      ret$annot = make.names(ss(ret$annotation,' \\('))
      return(ret)
    })
    ret = GRangesList(ret)
  } else if(class(peaks) == "GRanges"){
    ret <- as.GRanges(annotatePeak(peaks, tssRegion= fromTSS,
                             TxDb=txdb, annoDb=annodb))
    ret$annot = make.names(ss(ret$annotation,' \\('))
  } else{
      stop("Is inputs to peaks a GRanges or GRangesList?")
  }
  return(ret)
}

filterPeaks <-function(peaks, include =c('Distal.Intergenic','Intron')){
  if(class(peaks) %in% c('CompressedGRangesList', 'GRangesList')){
    ret = lapply(peaks, function(x) x[which(x$annot %in% include)])
    ret = GRangesList(ret)
  } else if(class(peaks) == "GRanges"){
    ret = peaks[which(peaks$annot %in% include)]
  } else{
    stop("Is inputs to peaks a GRanges or GRangesList?")
  }
  return(ret)
}


getLiftedChr <- function(p, chainFile){
  require(rtracklayer)
  # import the chain file for liftOver
  chain <- import.chain(chainFile)
  p = nameNarrowPeakRanges(p)
  
  # liftover the peak regions
  p1 = resize(p, width = 1); names(p1) = mcols(p1)$name
  p2 = unlist(GenomicRanges::reduce(rtracklayer::liftOver(p1, chain = chain)))
  p2 = GenomeInfoDb::keepStandardChromosomes(p2,pruning.mode="coarse")

  # add lifted chromosome to col
  mcols(p)$col = NA
  mcols(p)$col = as.character(seqnames(p2))[match(mcols(p)$name, names(p2))]
  p = p %>% sort() %>% as.data.frame() %>% fill(col, .direction = "updown") %>% GRanges()
  
  return(p)
}


getTxDb <- function(genome){
  if(genome =='hg38'){
    suppressMessages(require(TxDb.Hsapiens.UCSC.hg38.knownGene))
    txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
  }else if(genome =='mm10'){
    suppressMessages(require(TxDb.Mmusculus.UCSC.mm10.knownGene))
    txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
  }else if(genome =='rheMac8'){
    print('Using LiftOff of GRCh38 gene annotation -> rheMac8.')
    liftOff_gff= '/home/bnphan/resources/genomes/rheMac8/rheMac8_liftoff_GRCh38.p13_RefSeq.gff3'
    liftOff_txdb= '/home/bnphan/resources/genomes/rheMac8/rheMac8_liftoff_GRCh38.p13_RefSeq.txDb'
    if(!file.exists(liftOff_txdb)){
      txdb = makeTxDbFromGFF(liftOff_gff, dataSource= 'GRCh38 Refseq liftoff to rheMac8',
                             organism = 'Macaca mulatta', dbxrefTag = 'Dbxref')
      AnnotationDbi::saveDb(txdb, file = liftOff_txdb)
      } else{
      txdb = AnnotationDbi::loadDb(liftOff_txdb)
    }
  }else if(genome =='rheMac10'){
    print('Using LiftOff of GRCh38 gene annotation -> rheMac10')
    liftOff_gff= '/home/bnphan/resources/genomes/rheMac10/rheMac10_liftoff_GRCh38.p13_RefSeq.gff3'
    liftOff_txdb= '/home/bnphan/resources/genomes/rheMac10/rheMac10_liftoff_GRCh38.p13_RefSeq.txDb'
    if(!file.exists(liftOff_txdb)){
      txdb = makeTxDbFromGFF(liftOff_gff, dataSource= 'GRCh38 Refseq liftoff to rheMac10',
                             organism = 'Macaca mulatta', dbxrefTag = 'Dbxref')
      AnnotationDbi::saveDb(txdb, file = liftOff_txdb)
    } else{
      txdb = AnnotationDbi::loadDb(liftOff_txdb)
    }
  }else if(genome =='rn6'){
    suppressMessages(require(TxDb.Rnorvegicus.UCSC.rn6.refGene))
    txdb = TxDb.Rnorvegicus.UCSC.rn6.refGene
  }else{
    stop(paste0("Genome is not in current set. Choose from: ",supported_genomes))
  }
  return(txdb)
}

getBSgenome <- function(genome){
  if(genome =='hg38'){
    suppressMessages(require(BSgenome.Hsapiens.UCSC.hg38))
    bsgenome = BSgenome.Hsapiens.UCSC.hg38
  }else if(genome =='mm10'){
    suppressMessages(require(BSgenome.Mmusculus.UCSC.mm10))
    bsgenome = BSgenome.Mmusculus.UCSC.mm10
  }else if(genome =='rheMac8'){
    suppressMessages(require(BSgenome.Mmulatta.UCSC.rheMac8))
    bsgenome = BSgenome.Mmulatta.UCSC.rheMac8
  }else if(genome =='rheMac10'){
    suppressMessages(require(BSgenome.Mmulatta.UCSC.rheMac10))
    bsgenome = BSgenome.Mmulatta.UCSC.rheMac10
  }else if(genome =='rn6'){
    suppressMessages(require(BSgenome.Rnorvegicus.UCSC.rn6))
    bsgenome = BSgenome.Rnorvegicus.UCSC.rn6
  }else{
    stop(paste0("Genome is not in current set. Choose from: ",supported_genomes))
  }
  return(bsgenome)
}


getAnnoDb <- function(genome){
  # check if the genome is supported
  error_msg = paste0("Genome is not in current set. Choose from: ",
                     supported_genomes)
  if(!genome %in% supported_genomes) stop(error_msg)
  
  # get the org.db for each genome
  annodb = case_when(
    genome == 'hg38' ~ 'org.Hs.eg.db',
    genome == 'mm10' ~ 'org.Mm.eg.db',
    grepl('rheMac', genome) ~ 'org.Mmu.eg.db',
    genome == 'rn6' ~ 'org.Rn.eg.db',
  )
  return(annodb)
}

splitPeakSet <- function(gr, testSet = c('chr1','chr2'),
                         validSet = c('chr8','chr9'), useCol = 'seqnames'){
  indTest = gr %>% data.frame() %>% pull(useCol) %in% testSet %>% which()
  indValid= gr %>% data.frame() %>% pull(useCol) %in% validSet %>% which()
  indTrain = gr %>% data.frame() %>% pull(useCol) %ni% c(testSet, validSet) %>% which()
  
  ret = list(test = gr[indTest], train = gr[indTrain], valid = gr[indValid])
  return(ret)
}

getNonEnhOrthPeaks <- function(inPeaks, excludePeaks){
  oo = findOverlaps(query = inPeaks, subject = excludePeaks)
  # return queries NOT IN THE BG SET
  overlaps = unique(queryHits(oo))
  ret = inPeaks[seq_along(inPeaks) %ni% overlaps]
  return(ret)
  }

transferColumn <- function(toPeaks, fromPeaks, colOut = 'col'){
  tmp = fromPeaks %>% data.frame() %>% rename(col = seqnames) %>% 
    dplyr::select(col, name)
  ret = toPeaks %>% data.frame() %>% left_join(y = tmp, by = 'name') %>% 
    arrange(seqnames, start) %>% fill(contains('col')) %>% GRanges()
  return(ret)
}


writeGRangesToFasta <- function(gr, genome, file){
  bsgenome = getBSgenome(genome)
  seq <- getSeq(bsgenome,gr)
  
  ## add names
  names(seq) <- gr$name
  names(seq)[is.na(names(seq))] = 
    paste(genome,seq_len(sum(is.na(names(seq)))), sep = '_')
  
  ## write to fasta file
  if(grepl('.gz$', file)){
    writeXStringSet(seq,file=file, compress= TRUE)
  }else {
    writeXStringSet(seq,file=file)
  }
  return(seq)
}


summitCenter <- function(peaks, width = 501){
  if(class(peaks) %in% c('CompressedGRangesList', 'GRangesList')){
    ret = lapply(peaks, function(x) {
      ret = x
      start(ret) = start(x) + x$peak - (width - 1)/2
      end(ret) = start(ret) + width - 1
      return(ret)
    })
    ret = GRangesList(ret)
  } else if(class(peaks) == "GRanges"){
    ret = peaks
    start(peaks) + peaks$peak - ( width - 1)/2
    end(ret) = start(ret) + width -1
  } else{
    stop("Is inputs to peaks a GRanges or GRangesList?")
  }
  return(ret)
  
  return(ret)
}

convertHalChrName <- function(gr, species = 'Macaca_mulatta', chrOut = 'GenBank',
                              data.dir = here::here('data/tidy_data/Zoonomia_data/')){
  #### gather the zoonomia list of species ####
  zoo_fn = file.path(data.dir,'rdas','200_Mammals_Genome_Information.rds')
  if(!file.exists(zoo_fn)){
    zoo_df = file.path(data.dir,'tables','200_Mammals_Genome_Information.tsv') %>% 
      read_tsv() %>% rename_with(make.names) %>%
        mutate(assembly.report.file = basename(Genome.File.Name) %>% 
                 gsub(pattern = '_genomic.fna',replacement = '_assembly_report.txt.gz'),
               assembly.report.rds = gsub('_assembly_report.txt.gz','_assembly_report.rds', 
                                          assembly.report.file))
    saveRDS(zoo_df, file = zoo_fn)
  }else{
    zoo_df = readRDS(zoo_fn)
  }
  
  #### get the genome assembly report w/ chromosome mapping ####
  ind = match(species, zoo_df$Species)
  assembly_report_fn = file.path(data.dir,'rdas', zoo_df$assembly.report.rds[ind])
  if(!file.exists(zoo_fn)){
    chrmap_fn = file.path(data.dir,'tables', zoo_df$assembly.report.file[ind])
    map = data.table::fread(chrmap_fn, skip = '# Sequence-Name') %>%
      rename_with(gsub, pattern ='# ',replacement =  '') %>%
      rename_with(make.names) %>%
      filter(Assigned.Molecule != 'na') %>%
      filter(GenBank.Accn != 'na')
    # save the map file
    saveRDS(map, file = assembly_report_fn)
  } else {
    map = readRDS(assembly_report_fn)
  }
  
  #### set up dictionary to-from chromosomes ####
  target_chr = dplyr::case_when(chrOut == 'GenBank' ~ 'GenBank.Accn', 
                                chrOut == 'UCSC' ~ 'UCSC.style.name',
                                chrOut == 'RefSeq' ~ 'RefSeq.Accn')

  #### get most likely current naming convention ####
  current_chr = names(map)[which.max(apply(map, 2, function(x)
    sum(seqnames(gr) %in% x, na.rm = T)))]
  dict = map %>% pull(target_chr, current_chr)
  
  #### remap the genomic ranges ####
  gr_renamed = gr %>% as.data.frame() %>%
      filter(seqnames %in% names(dict)) %>%
      mutate(seqnames = dict[as.character(seqnames)]) %>% 
      GRanges()
  return(gr_renamed)
}
