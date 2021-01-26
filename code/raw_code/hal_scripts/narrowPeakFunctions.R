write_GRangesToNarrowPeak <- function(gr, file, genome = 'hg38', ...) {

  narrowPeakCols = c("seqnames", "start", "end", "name", "score", "strand", 
                     "signalValue", "pValue", "qValue", "peak")

  #give unique names if not already, and summit
  gr = nameNarrowPeakRanges(gr = gr, genome = genome)

  # make sure not using peaks w/ duplicate names
  names(gr) = NULL; x <- as(gr, "data.frame")
  
  # add extra narrowPeak columns if not present
  for (n in c('score','signalValue','pValue','qValue')) {
    if(!n %in% names(x))
      x[,n]  <- -1
  }
  
  # reorder narropeak columns & write narrowPeak file
  x <- x[, narrowPeakCols]
  
  # make sure directory exists
  system(paste('mkdir -p', dirname(file)))
  
  # handle gzipped narrowPeak files
  if(grepl('.gz', file)){
    write.table(x, gzfile(file), sep = "\t", row.names = FALSE, 
                col.names = FALSE, quote = FALSE, ...)
  } else {
    write.table(x, file, sep = "\t", row.names = FALSE, 
                col.names = FALSE, quote = FALSE, ...)
    }
  toReturn <- as(x, "GenomicRanges")
  return(gr)
}

addSummitCenter <- function(gr, force = FALSE){
  # only add if forced or no summits
  if( force | all(is.null(mcols(gr)$peak)))
    # use peak center to make summit
    mcols(gr)$peak = round(abs(start(gr) - end(gr)) / 2)
  return(gr)
}

nameNarrowPeakRanges <- function(gr, genome = 'hg38'){
  # make peak offset summit in the middle of the range
  if(all(is.null(mcols(gr)$peak)))
    mcols(gr)$peak = round(abs(start(gr) - end(gr)) / 2)
  
  # give unique names using genome coordinates
  if(all(is.null(mcols(gr)$name)))
    # make sure right genome is given
    mcols(gr)$name = paste0(genome,':', 
      # use chr:start-end:summit to name peak
      seqnames(gr),':', start(gr),'-',end(gr), ':',mcols(gr)$peak)

  # remove duplicate peaks
  gr = gr[!duplicated(mcols(gr)$name),]
  return(gr)
}

#' make sure p is a GRanges object 
#' with unique peak names and a summit offset from the start position
#' appropriate chain file for liftOver
liftOver_narrowPeak <- function(p, chainFile){
  require(rtracklayer)
  
  #####################################
  # import the chain file for liftOver
  chain <- import.chain(chainFile)
  
  #####################################
  # add summits if missing, make sure peak names are unique
  p = addSummitCenter(p, force = FALSE)
  p = nameNarrowPeakRanges(p)
  names(p) = mcols(p)$name
  
  ############################
  # liftover the peak regions
  p2 = unlist(GenomicRanges::reduce(rtracklayer::liftOver(p, chain = chain)))
  p2 = GenomeInfoDb::keepStandardChromosomes(p2,pruning.mode="coarse")
  
  ####################################################
  # lift over the peak summits, 
  s = p; start(s) = end(s) = start(s) + mcols(s)$peak
  s2 = unlist(GenomicRanges::reduce(rtracklayer::liftOver(s, chain = chain)))
  s2 = GenomeInfoDb::keepStandardChromosomes(s2,pruning.mode="coarse")
  
  ############################
  # keep peaks that survived
  keepIdx = Reduce('intersect', lapply(list(p,p2,s2), names))
  p = p[keepIdx]; p2 = p2[keepIdx];   s2 = s2[keepIdx]
  mcols(p2) = mcols(p)
  mcols(p2)$peak= start(s2) - start(p2)
  
  ##################################################
  # drop the peaks that discordantly mapped summits
  p2 = p2[ start(p2) < start(s2) & end(p2) > end(s2)]
  return(p2)
}

#' given a set of peaks mapped to multiple coordinate systems
#' derived from one peak set with unique mcols(peaks)$name
#' intersect the peaks and return one or all peak sets
keepOrthologs <- function(peakList, idxReturn = NA){
  
  # get peaks intersecting all coordinate systems
  tmpList = lapply(peakList, as, 'data.frame')
  idxKeep = Reduce('intersect',lapply(tmpList, '[[', 'name'))
  
  # subset the 
  peakList = lapply(peakList, function(peak){
    peak[match(idxKeep, mcols(peak)$name)]
  })
  
  # return the desired peakset
  if(idxReturn %in% seq_along(peakList)){
    toReturn = peakList[[idxReturn]]
  } else {
    toReturn = peakList
  }
  return(toReturn)
}

mergeNarrowPeaks <- function(narrowPeakList, width = 501, min.gapwidth = 100){
  # extract the summits
  summits = GRangesList(lapply(narrowPeakList, function(p){
    start(p) = end(p) =  start(p) + mcols(p)$peak
    mcols(p) = NULL
    return(p)
  }))
  
  # merge the summits, keep those not too spread out after merging
  summits = GenomicRanges::reduce(unlist(summits), min.gapwidth=min.gapwidth)
  summits = summits[width(summits) < width]
  
  # merge the peaks
  peaks = GenomicRanges::reduce(unlist(GRangesList(narrowPeakList)))
  
  # take the average of merged summits
  start(summits) = end(summits) = 
    (start(summits) + end(summits))/2
  
  # resize the peaks
  start(summits) = start(summits) - round(width/2)
  end(summits) = end(summits) + round(width/2)
  mcols(summits)$peaks = round(width/2)
  
  # drop non-standard chromosomes
  summits = GenomeInfoDb::keepStandardChromosomes(
    summits,pruning.mode="coarse")
  return(summits)
}


# instead of merging, keep peaks separate but ordered by overlapping
linkNarrowPeaks <- function(narrowPeakList, min.gapwidth = 100){
  # param min.gapwidth to merge peaks
  # find peaks that overlaps other list
  tmpList = lapply(names(narrowPeakList), function(x){
    ooList = suppressWarnings(
      lapply(narrowPeakList[!names(narrowPeakList) %in% x], 
             findOverlaps,subject = narrowPeakList[[x]]))
    idx = Reduce('intersect', lapply(ooList, subjectHits))
    return(narrowPeakList[[x]][idx])
  })
  
  # extract the summits
  summits = GRangesList(lapply(tmpList, function(p){
    start(p) = end(p) =  start(p) + mcols(p)$peak
    mcols(p) = NULL
    return(p)
  }))
  
  # merge the summits, keep those not too spread out after merging
  summits = GenomicRanges::reduce(unlist(summits), min.gapwidth=min.gapwidth)
  summits = sort(summits)
  
  # keep summits overlapping all peak sets
  summits = summits[Reduce('intersect', lapply(tmpList, function(x)
    subjectHits(findOverlaps(x, subject = summits))))]
  
  # take the average of merged summits
  start(summits) = end(summits) = 
    (start(summits) + end(summits))/2
  
  # return the closest peak to summit
  tmpList = lapply(tmpList,function(x) {
    idx = nearest(x = summits, subject = x, select= 'arbitrary')
    return(x[idx])
    })
  
  # remove any peak duplicated in any group
  dupIdx = sapply(tmpList, duplicated)
  keepIdx = which(!apply(dupIdx, 1, all))
  tmpList = lapply(tmpList, function(x) x[keepIdx])
  names(tmpList) = names(narrowPeakList)
  
  return(tmpList)
}


