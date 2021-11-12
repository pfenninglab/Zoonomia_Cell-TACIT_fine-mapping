suppressMessages(library(ArchR))
suppressMessages(library(tidyverse))
suppressMessages(library(here))
library(rcartocolor)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)


PROJDIR='data/raw_data/'
FIGDIR='figures/explanatory/figure1_snATAC_cross-species_plots'
celltypes = c('MSN_D1', 'MSN_D2', 'MSN_SN', 'INT_Pvalb', 'Astro', 'Microglia', 'OPC', 'Oligo')
markerGenes <- c("DRD1", "DRD2", "RXFP1",'LHX6', "AQP4", "CX3CR1", "PDGFRA", "MAG")

##################################
### set Arrow File parameters ####
addArchRThreads(threads = 10)
addArchRGenome("hg38")


##########################################
### load in the human snATAC-seq data ####
PROJDIR_HG=here(PROJDIR,'hg38/Corces_2020','ArchR_Corces2020_caudate_labeled')
projHg = loadArchRProject(PROJDIR_HG)
names(getCellColData(projHg))
table(projHg$Clusters2)

### Plot the Human snATAC track plots ####
pdf('tmp.pdf')
pList <- plotBrowserTrack( 
  projHg, groupBy = "Clusters2", pal = carto_pal(8, "Safe"),
  sizes = c(5, 1), plotSummary = c("bulkTrack", "geneTrack"), 
  geneSymbol = markerGenes, tileSize = 100, baseSize = 6, 
  useGroups = celltypes, upstream = 50000, downstream = 50000)
dev.off()


hgPlot_fn = here(FIGDIR, 'plots', 'fig1_Corces2020_caudate_marker_tracks.pdf')
pdf(hgPlot_fn, height = 1*4, width = .75*4, onefile = T)
for (p in pList) {grid::grid.draw(p); grid::grid.newpage()}
dev.off()


############################################
### load in the macaque snATAC-seq data ####
PROJDIR_RM=here(PROJDIR,'rheMac10/Stauffer_caudate','ArchR_Stauffer_caudate_labeled')
projRm = loadArchRProject(PROJDIR_RM)
names(getCellColData(projRm))
table(projRm$Clusters2)

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



############################################
### load in the mouse snATAC-seq data ####
addArchRGenome("mm10")
PROJDIR_MM=here(PROJDIR,'mm10/BICCN_mouse_caudoputamen','ArchR_BICCN_CP_labeled')
projMm = loadArchRProject(PROJDIR_MM)
names(getCellColData(projMm))
table(projMm$Clusters2)
projMm = projMm[projMm$Clusters2%in% c('MSN_SN', 'Astro', 'Microglia', 'OPC', 'Oligo')]
table(projMm$Clusters2)

markerGenes2 <- c("DRD1",'DRD2', "FOXP2", "LHX6", "AQP4", "CX3CR1", "PDGFRA", "MAG")

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











