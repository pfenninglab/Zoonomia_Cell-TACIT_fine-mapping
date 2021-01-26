makeChomAssaySeurat <- function(counts, metadata, ensdb, group.by = NULL, ranges = NULL,
                                min.cells = 10, min.features = 0, verbose = TRUE){
  require('Signac')
  
  # only keep peaks able to be 
  if(!is.null(group.by)){
    print(paste('Dropping peaks with fewer than', 
                min.cells,'cells in each', group.by,'.'))
    idxList = split(rownames(metadata), metadata[, group.by])
    keepIdx = Reduce('intersect', lapply(idxList, function(i) 
      which(rowSums(counts[,i])> min.cells)))
    counts = counts[keepIdx, ]
  }
  
  ### 
  print('Making Seurat object.')
  chrom_assay <- CreateChromatinAssay(
    counts = counts, sep = c(":", "-"),  genome = 'hg38', ranges = ranges,
    fragments = NULL, min.cells = min.cells, min.features = min.features)
  
  obj_seurat <- CreateSeuratObject( counts = chrom_assay, assay = "peaks", 
                                    metadata = metadata)
  obj_seurat@meta.data = cbind(
    obj_seurat@meta.data, metadata[colnames(obj_seurat),])

  print(obj_seurat)
  
  ### extract gene annotations from EnsDb ####
  print('Adding annotation features to Seurat object.')
  annotations <- GetGRangesFromEnsDb(ensdb = ensdb)
  
  # change to UCSC style since the data was mapped to hg38
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "hg38"
  
  # add the gene information to the object
  Annotation(obj_seurat) <- annotations
  
  ## perform TF-IDF on the peak counts ###
  print('Adding TFIDF normalization, top features, and SVD to Seurat object.')
  obj_seurat <- RunTFIDF(obj_seurat, verbose = verbose)
  obj_seurat <- FindTopFeatures(obj_seurat, min.cutoff = 'q0', verbose = verbose)
  obj_seurat <- RunSVD(obj_seurat, verbose = verbose)
  
  ## run UMAP on unintegrated object, ignore first dim, correlates w/ sequencing depth
  print('Adding UMAP normalization to Seurat object.')
  obj_seurat <- RunUMAP(object = obj_seurat, verbose = verbose, 
                        reduction = 'lsi', dims = 2:30)
  obj_seurat <- FindNeighbors(object = obj_seurat, verbose = verbose, 
                              reduction = 'lsi', dims = 2:30)
  
  return(obj_seurat)
}