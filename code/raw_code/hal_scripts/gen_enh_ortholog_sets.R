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

annotatePeaks <- function(peaks, genome = 'hg38', fromTSS = c(-5000,5000)){
  txdb = getTxDb(genome)
  annodb = getAnnoDb(genome)
  print(paste('Using TSS boundaries from',fromTSS[1], 'to',fromTSS[2]))
  if(class(peaks) == c('CompressedGRangesList', 'GRangesList')){
    ret = lapply(peaks, annotatePeak, TxDb=txdb, annoDb=annodb, tssRegion=fromTSS)
    ret = lapply(ret, function(x){
      ret = as.GRanges(x)
      ret$annot = make.names(ss(ret$annotation,' \\('))
      return(ret)
    })
    ret = GRangesList(ret)
  } else if(class(peaks) == "GRanges"){
    ret <- as.GRanges(annotatePeak(peaks, tssRegion= fromTSS),
                             TxDb=txdb, annoDb=annodb)
    ret$annot = make.names(ss(ret$annotation,' \\('))
  } else{
      stop("Is inputs to peaks a GRanges or GRangesList?")
  }
  return(ret)
}

getTxDb <- function(genome){
  if(genome =='hg38'){
    suppressMessages(require(TxDb.Hsapiens.UCSC.hg38.knownGene))
    txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
  }else if(genome =='mm10'){
    suppressMessages(require(TxDb.Mmusculus.UCSC.mm10.knownGene))
    txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
  }else if(genome =='rheMac8'){
    suppressMessages(require(TxDb.Mmulatta.UCSC.rheMac8.refGene))
    txdb = TxDb.Mmulatta.UCSC.rheMac10.refGene
  }else if(genome =='rheMac10'){
    suppressMessages(require(TxDb.Mmulatta.UCSC.rheMac10.refGene))
    txdb = TxDb.Mmulatta.UCSC.rheMac10.refGene
  }else if(genome =='rn6'){
    suppressMessages(require(TxDb.Rnorvegicus.UCSC.rn6.refGene))
    txdb = TxDb.Rnorvegicus.UCSC.rn6.refGene
  }else{
    stop(paste0("Genome is not in current set. Choose from: ",supported_genomes))
  }
  return(txdb)
}

getAnnoDb <- function(genome){
  # check if the genome is supported
  error_msg = paste0("Genome is not in current set. Choose from: ",
                     supported_genomes)
  if(!genome %in% supported_genomes) stop(error_msg)
  
  # get the org.db.
  annodb = case_when(
    genome == 'hg38' ~ 'org.Hs.eg.db',
    genome == 'mm10' ~ 'org.Mm.eg.db',
    genome == 'rheMac10' ~ 'org.Mmu.eg.db',
    genome == 'rn6' ~ 'org.Rn.eg.db',
  )
  return(annodb)
}

convertHalChrName <- function(gr, species = 'Macaca_mulatta', chrOut = 'GenBank',
                              data.dir = '../../../data/tidy_data/Zoonomia_data/'){
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
