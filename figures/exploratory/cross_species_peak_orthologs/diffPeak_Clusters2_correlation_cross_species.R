# to be run in the root github directory
setwd('figures/exploratory/cross_species_peak_orthologs')
PROJDIR=file.path('../../../data/raw_data/cross_species_peak_orthologs')

#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F, repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(rtracklayer)); suppressMessages(library(ArchR))

#########################################################
### get the peaks from experiment and halpered peaks ####
peakListRDS_fn = file.path(PROJDIR, 'rdas', 'multispeciesMergedPeakList.rds')
peakList = readRDS(peakListRDS_fn)

cellTypes = c('Astro','Microglia','Oligo','OPC',
              'MSN_D1','MSN_D2',"Interneuron")


##########################################
### read in Corces2020 ArchR project ####
addArchRGenome("hg38")
PROJDIR2='../../../data/raw_data/hg38/Corces_2020'
ARCHDIR2=file.path(PROJDIR2,paste0('ArchR_Corces2020_caudate_labeled'))
humanProj = loadArchRProject(ARCHDIR2, showLogo = F)

humanProj = humanProj[humanProj$Clusters2 %in% cellTypes, ]
humanMarkersSE = getMarkerFeatures(
  ArchRProj = humanProj, groupBy = "Clusters2", useMatrix = "OrthologPeakMatrix",
  bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon",
  maxCells = 1000, scaleTo = 10^4,  k = 100, binarize = TRUE)

#############################################
### read in mouse BICCN CP ArchR project ####
addArchRGenome("mm10")
PROJDIR3='../../../data/raw_data/mm10/BICCN_mouse_caudoputamen'
ARCHDIR3=file.path(PROJDIR3,paste0('ArchR_BICCN_CP_labeled'))
mouseProj = loadArchRProject(ARCHDIR3, showLogo = F)

mouseProj = mouseProj[mouseProj$Clusters2 %in% cellTypes, ]
mouseMarkersSE = getMarkerFeatures(
  ArchRProj = mouseProj, groupBy = "Clusters2", useMatrix = "OrthologPeakMatrix",
  bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon",
  maxCells = 1000, scaleTo = 10^4,  k = 100, binarize = TRUE)

##############################
### read in ArchR project ####
GENOMEDIR='/home/bnphan/resources/genomes/rheMac10'
load(file.path(GENOMEDIR,'rheMac10_liftoff_GRCh38.p13_ArchR_annotations.rda'))
PROJDIR4='../../../data/raw_data/rheMac10/Stauffer_caudate'
ARCHDIR4=file.path(PROJDIR4,paste0('ArchR_Stauffer_caudate_labeled'))
macaqueProj = loadArchRProject(ARCHDIR4, showLogo = F)

macaqueProj = macaqueProj[macaqueProj$Clusters2 %in% cellTypes, ]
macaqueMarkersSE = getMarkerFeatures(
  ArchRProj = macaqueProj, groupBy = "Clusters2", useMatrix = "OrthologPeakMatrix",
  bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon",
  maxCells = 1000, scaleTo = 10^4,  k = 100, binarize = TRUE)


markerSEList = 
  list(Human = humanMarkersSE, Mouse = mouseMarkersSE, Macaque = macaqueMarkersSE)
system('mkdir -p rdas')
markerRDS_fn = file.path('rdas', 'cross_species_markerSEList.rds')
saveRDS(markerSEList, file = markerRDS_fn)


#####################################
### correlate marker peak Log2FC ####
correlateMarkerSE = function(se1, se2, label1, label2, alpha = 0.10,
                             method = 'spearman', assay = 'Log2FC'){
  # make sure both se have the same rows
  idxName = intersect(rowData(se1)$name, rowData(se2)$name)
  se1 = se1[match(idxName, rowData(se1)$name),]
  se2 = se2[match(idxName, rowData(se2)$name),]
  
  # do for every pairs of cell type shared across SE's
  cells = intersect(colnames(se1), colnames(se2))
  compareCells = expand.grid(cells,cells)
  rownames(compareCells) = apply(compareCells, 1, paste, collapse = '_')
  
  corList = parallel::mclapply(rownames(compareCells), function(rn){
    cell1 = compareCells[rn,1]; cell2 = compareCells[rn,2]
    idx = union(which(assays(se1)$FDR[,cell1]< alpha), 
                    which(assays(se2)$FDR[,cell2]< alpha))
    tmp = cor.test(assays(se1)[[assay]][idx,cell1], 
                   assays(se2)[[assay]][idx,cell2], method = method)
    
    ret = data.frame(cell1 = cell1, cell2 = cell2, 
                     nMarkers1 = length(which(assays(se1)$FDR[,cell1]< alpha)),
                     nMarkers2 = length(which(assays(se2)$FDR[,cell2]< alpha)),
                     nMarkerOverlap = length(idx),
                     label1 = label1, label2 = label2,
                     rho = tmp$estimate, p.value = tmp$p.value,
                     assay = assay, method = method)
    return(ret)
  }, mc.cores =getArchRThreads())
  ret = do.call('rbind', corList)
  return(ret)
}

library(ComplexHeatmap)
library(reshape2)
library(tidyverse)

# perform spearman correlation w/ human vs. mouse and macaque OCR
corSpearman_df = rbind(
  correlateMarkerSE(se1 = humanMarkersSE, se2 = mouseMarkersSE, alpha = 0.05,
                    label1 = 'Human', label2 = 'Mouse', method = 'spearman'),
  correlateMarkerSE(se1 =humanMarkersSE, se2 = macaqueMarkersSE, alpha = 0.05, 
                    label1 = 'Human', label2 = 'Macaque', method = 'spearman'))

# perform pearson correlation 
corPearson_df = rbind(
  correlateMarkerSE(se1 = humanMarkersSE, se2 = mouseMarkersSE, alpha = 0.05,
                    label1 = 'Human', label2 = 'Mouse', method = 'pearson'),
  correlateMarkerSE(se1 =humanMarkersSE, se2 = macaqueMarkersSE, alpha = 0.05,
                    label1 = 'Human', label2 = 'Macaque', method = 'pearson'))

# add new columns
corSpearman_df = corSpearman_df %>%
  mutate(p.bonferroni = p.adjust(p.value, method = 'bonferroni'),
         row.label = paste(cell1, label1, sep = '_'),
         col.label = paste(cell2, label2, sep = '_'))
corPearson_df = corPearson_df %>%
  mutate(p.bonferroni = p.adjust(p.value, method = 'bonferroni'),
         row.label = paste(cell1, label1, sep = '_'),
         col.label = paste(cell2, label2, sep = '_'))

cellType = c('Astro',"Interneuron",'Microglia','MSN_D1','MSN_D2','Oligo','OPC')
cellTypeCol = RColorBrewer::brewer.pal(length(cellType), 'Paired') 
names(cellTypeCol) = cellType
speciesCol = RColorBrewer::brewer.pal(3,'Dark2')
names(speciesCol) = c('Human','Mouse','Macaque')


## transform metrics into matrices
corMatSpearman = acast(corSpearman_df, cell1 + label1 ~ cell2 + label2, 
                       value.var = 'rho')
corMatSpearPBonf = acast(corSpearman_df, cell1 + label1 ~ cell2 + label2, 
                         value.var = 'p.bonferroni')
corMatPearson = acast(corPearson_df, cell1 + label1 ~ cell2 + label2, 
                      value.var = 'rho')
corMatPearPBonf = acast(corPearson_df, cell1 + label1 ~ cell2 + label2, 
                        value.var = 'p.bonferroni')

## set up row & column annotations
idxRow = match(rownames(corMatSpearman), corSpearman_df$row.label)
row_ha = with(corSpearman_df, rowAnnotation(
  nPeaks = anno_barplot(nMarkers1[idxRow], 
                        gp = gpar(fill = cellTypeCol[cell1[idxRow]]),
                        axis_param = list(direction = 'reverse')),
  cellType = cell1[idxRow], col = list(cellType = cellTypeCol),
  show_legend = FALSE))

idxCol = match(colnames(corMatSpearman), corSpearman_df$col.label)
column_ha = with(corSpearman_df, HeatmapAnnotation(
  nPeaks = anno_barplot(nMarkers2[idxCol],
                        gp = gpar(fill = cellTypeCol[cell2[idxCol]])),
  species = label2[idxCol], cellType = cell2[idxCol], 
  col = list(cellType = cellTypeCol, species = speciesCol)))

#######################################################
#### make marker correlation heatmap  plot ###########
system('mkdir -p plots')
heatmap_fn = file.path('plots', 'cross_species_marker_correlation.pdf')

pdf(heatmap_fn, height = 5, width = 10)
hm1 = Heatmap(corMatSpearman,
              # matrix parameters
              heatmap_legend_param = list(
                title = expression(paste("Spearman's ", rho)), 
                heatmap_legend_side = "bottom",
                legend_width = unit(3.5, "in"),
                title_position = "topcenter",
                direction = "horizontal"),
              
              # add text to each grid w/  correlation
              cell_fun = function(j, i, x, y, w, h, col) {
                if(corMatSpearman[i, j] > .2) # by magnitude
                  grid.text(signif(corMatSpearman[i, j],2), x, y)
              },
              
              # row annotations, human
              row_order = sort(rownames(corMatPearson)), 
              row_label = gsub('_Human','',rownames(corMatSpearman)),
              row_title = 'Human Marker Peak',
              left_annotation = row_ha, 
              
              # column annotations, mouse & macaque
              column_order = sort(colnames(corMatPearson)),
              column_title = "Orthologous Marker Peak Log2FC Correlation",
              show_column_names = FALSE,
              top_annotation = column_ha)
draw(hm1, merge_legend = FALSE, heatmap_legend_side = "bottom", 
     annotation_legend_side = "right")

hm2 = Heatmap(corMatPearson,
        # matrix parameters
        heatmap_legend_param = list(
          title = "Pearson's R", heatmap_legend_side = "bottom",
          legend_width = unit(3.5, "in"), title_position = "topcenter",
          direction = "horizontal"),
        
        # add text to each grid w/ signif correlation
        cell_fun = function(j, i, x, y, w, h, col) {
          if(corMatPearson[i, j] > .2) # by magnitude
            grid.text(signif(corMatPearson[i, j],2), x, y)
        },
        
        # row annotations, human
        row_order = sort(rownames(corMatSpearman)), 
        row_label = gsub('_Human','',rownames(corMatSpearman)),
        row_title = 'Human Marker Peak',
        left_annotation = row_ha, 
        
        # column annotations, mouse & macaque
        column_order = sort(colnames(corMatPearson)),
        column_title = "Orthologous Marker Peak Log2FC Correlation",
        show_column_names = FALSE,
        top_annotation = column_ha)
draw(hm2, merge_legend = FALSE, heatmap_legend_side = "bottom", 
     annotation_legend_side = "right")
dev.off()


