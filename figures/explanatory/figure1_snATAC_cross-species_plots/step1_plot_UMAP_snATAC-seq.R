suppressMessages(library(ArchR))
suppressMessages(library(tidyverse))
suppressMessages(library(here))
library(rcartocolor)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

PROJDIR='data/raw_data/'
FIGDIR='figures/explanatory/figure1_snATAC_cross-species_plots'

##################################
### set Arrow File parameters ####
addArchRThreads(threads = 10)
addArchRGenome("hg38")

celltypes = c('MSN_D1', 'MSN_D2', 'MSN_SN', 'INT_Pvalb', 'Astro', 'Microglia', 'OPC', 'Oligo')
cell_cols = carto_pal(8, "Safe")
names(cell_cols) = celltypes
  
##########################################
### load in the human snATAC-seq data ####
PROJDIR_HG=here(PROJDIR,'hg38/Corces_2020','ArchR_Corces2020_caudate_labeled')
projHg = loadArchRProject(PROJDIR_HG)
names(getCellColData(projHg))
table(projHg$Clusters2)

### load in the macaque snATAC-seq data ####
PROJDIR_RM=here(PROJDIR,'rheMac10/Stauffer_caudate','ArchR_Stauffer_caudate_labeled')
projRm = loadArchRProject(PROJDIR_RM)
names(getCellColData(projRm))
table(projRm$Clusters2)

### load in the mouse snATAC-seq data ####
addArchRGenome("mm10")
PROJDIR_MM=here(PROJDIR,'mm10/BICCN_mouse_caudoputamen','ArchR_BICCN_CP_labeled')
projMm = loadArchRProject(PROJDIR_MM)
names(getCellColData(projMm))
table(projMm$Clusters2)
projMm = projMm[projMm$Clusters2%in% c('MSN_SN', 'Astro', 'Microglia', 'OPC', 'Oligo')]
table(projMm$Clusters2)

###########################################
### plot UMAP of human, macaque, mouse ####
pdf('tmp.pdf') 
pHg <- plotEmbedding(ArchRProj = projHg, colorBy = "cellColData", name = "Clusters2", 
                     embedding = "UMAP", pal = cell_cols) + 
  ggtitle('Human snATAC-seq') + theme(legend.position = 'none')

pRm <- plotEmbedding(ArchRProj = projRm, colorBy = "cellColData", name = "Clusters2", 
                     embedding = "UMAP", pal = cell_cols) + 
  ggtitle('Macaque snATAC-seq') + theme(legend.position = 'none')

pMm <- plotEmbedding(ArchRProj = projMm, colorBy = "cellColData", name = "Clusters2", 
                     embedding = "UMAP", pal = cell_cols) + 
  ggtitle('Mouse snATAC-seq') + theme(legend.position = 'none')
dev.off()

hgPlot_fn = here(FIGDIR, 'plots', 'fig1_cross-species_caudate_UMAP.pdf')
pdf(hgPlot_fn, height = .9*4, width = 2.25*4, onefile = T)
ggAlignPlots(pHg, pRm, pMm, type = "h")
dev.off()



### Plot the Rhesus Macaque snATAC track plots ####
pdf('tmp.pdf')
pList2 <- plotBrowserTrack( 
  projRm, groupBy = "Clusters2", pal = carto_pal(8, "Safe"),
  sizes = c(5, 1), plotSummary = c("bulkTrack", "geneTrack"), 
  geneSymbol = markerGenes, tileSize = 100, baseSize = 6, 
  useGroups = celltypes, upstream = 50000, downstream = 50000)
dev.off()

rmPlot_fn = here(FIGDIR, 'plots', 'fig1_rhesusMacaque_caudate_marker_tracks.pdf')
pdf(rmPlot_fn, height = 1*4, width = .75*4, onefile = T)
for (p in pList2) {grid::grid.draw(p); grid::grid.newpage()}
dev.off()


### Plot the BICCN mouse snATAC track plots ####
pdf('tmp.pdf')
pList3 <- plotBrowserTrack( 
  projMm, groupBy = "Clusters2", pal = carto_pal(8, "Safe")[c(3,5:8)],
  sizes = c(5, 1), plotSummary = c("bulkTrack", "geneTrack"), 
  normMethod = c('ReadsInPromoter'),
  geneSymbol = markerGenes2, tileSize = 500, baseSize = 6, 
  useGroups = celltypes, upstream = 50000, downstream = 50000)
dev.off()

mmPlot_fn = here(FIGDIR, 'plots', 'fig1_BICCN_mouse_caudate_marker_tracks.pdf')
pdf(mmPlot_fn, height = .65*4, width = .75*4, onefile = T)
for (p in pList3) {grid::grid.draw(p); grid::grid.newpage()}
dev.off()











