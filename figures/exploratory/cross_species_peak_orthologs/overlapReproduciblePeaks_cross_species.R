# to be run in the root github directory
setwd('figures/exploratory/cross_species_peak_orthologs')
PROJDIR=file.path('../../../data/raw_data/cross_species_peak_orthologs')

#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F, repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(rtracklayer)); suppressMessages(library(ArchR))

source('../../../code/raw_code/hal_scripts/narrowPeakFunctions.R')

#################################################################
### collect human peak files that mapped to mouse and macaque ###
ortholog_fn = list.files(file.path(PROJDIR, 'peaks'), full.names = T,
                         pattern = '_orthologPeakList.rds')
names(ortholog_fn) = ss(basename(ortholog_fn), '_orthologPeakList',1)
orthologPeakList = lapply(ortholog_fn, readRDS)
celltypes = Reduce('intersect', lapply(orthologPeakList, names))
names(celltypes) = celltypes
orthologPeakList = lapply(celltypes, function(x){
  lapply(orthologPeakList, `[[`, x)
})
orthologPeakList = orthologPeakList[!names(orthologPeakList) =='Consensus']


################################################
# intersect narrowPeak sets with human-mouse and human-macaque orthologs
linkedOrthologPeakList = lapply(orthologPeakList, linkNarrowPeaks, min.gapwidth = 251)
lengths(linkedOrthologPeakList)

# union narrowPeak sets mapped to human coordinates from human, mouse, macaque
mergedOrthologPeakList = lapply(orthologPeakList, mergeNarrowPeaks, 
                                width = 501, min.gapwidth = 251)
lengths(mergedOrthologPeakList)


################################################
# create boolean overlap matrix of unionized peaks w/ orthologous mapped peaks
library(ComplexHeatmap)
library(GenomicRanges)

cellType = c('Astro',"Interneuron",'Microglia','MSN_D1','MSN_D2','Oligo','OPC')
cellTypeCol = RColorBrewer::brewer.pal(length(cellType), 'Paired') 
names(cellTypeCol) = cellType

species = c('Human','Macaque' ,'Mouse')
names(species)= c('Corces_2020_caudate','Stauffer_caudate','BICCN_CP')
speciesCol = RColorBrewer::brewer.pal(3,'Dark2')
names(speciesCol) = species

orthologPeakList = lapply(orthologPeakList, function(m){
  names(m) = species[names(m)]
  return(m)
})

m_list = lapply(names(orthologPeakList), function(n){
  countMatList = suppressWarnings(
    lapply(orthologPeakList[[n]], GenomicRanges::countOverlaps, 
           query = mergedOrthologPeakList[[n]]))
  countMat = do.call('cbind', countMatList)
  countMat[countMat >0] = 1 #binarize
  countMat = countMat[,species]
  ret = ComplexHeatmap::make_comb_mat(
    countMat, mode = 'distinct', remove_empty_comb_set = FALSE, 
    remove_complement_set = FALSE)
  return(ret)
})

names(m_list) = names(orthologPeakList)
m_list = normalize_comb_mat(m_list)
max_set_size = max(sapply(m_list, set_size)) # max peaks per cell type per dataset
max_comb_size = max(sapply(m_list, comb_size)) # max intersection peaks

ht_list = NULL
for(i in seq_along(m_list)) {
  ht_list = ht_list %v%
    UpSet(m_list[[i]], row_title = names(m_list)[i],
          set_order = species, 
          comb_order = NULL,
          top_annotation = HeatmapAnnotation(
            "# Overlap" = anno_barplot( 
              comb_size(m_list[[i]]), ylim = c(0, max_comb_size), annotation_name_rot = 90,
              gp = gpar(fill = cellTypeCol[names(m_list)[i]]),
              height = unit(1.5, "cm")), 
            annotation_name_side = "left"
            ),
          right_annotation = rowAnnotation(
            "# Orthologs" = anno_barplot(
              set_size(m_list[[i]]), ylim = c(0, max_set_size),
              gp = gpar(fill = speciesCol),
              width = unit(2.5, "cm")))
          )
}

upset_fn = file.path('plots', 'cross_species_orthologUpset_allCells.pdf')
pdf(upset_fn, height = 11, width = 4)
draw(ht_list, gap = unit(.7, "cm"))
dev.off()

ht_list2 = NULL
for(i in grep('MSN',names(m_list))) {
  ht_list2 = ht_list2 +
    UpSet(m_list[[i]], column_title = names(m_list)[i],
          set_order = species, 
          comb_order = NULL,
          top_annotation = HeatmapAnnotation(
            "# Overlap" = anno_barplot( 
              comb_size(m_list[[i]]), ylim = c(0, max_comb_size), 
              annotation_name_rot = 90,
              gp = gpar(fill = cellTypeCol[names(m_list)[i]]),
              height = unit(2, "cm")), 
            annotation_name_side = "left"
          ),
          right_annotation = rowAnnotation(
            "# Orthologs" = anno_barplot(
              set_size(m_list[[i]]), ylim = c(0, max_set_size),
              gp = gpar(fill = speciesCol),
              width = unit(2, "cm")))
    )
}

upsetMSN_fn = file.path('plots', 'cross_species_orthologUpset_MSN.pdf')
pdf(upsetMSN_fn, height = 3, width = 8)
draw(ht_list2, gap = unit(1, "cm"))
dev.off()


