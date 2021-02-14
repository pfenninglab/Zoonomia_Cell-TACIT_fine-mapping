# to be run in the root github directory
setwd('figures/explanatory/overlap_caudate_peak_orth')

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

### set Arrow File parameters ####
suppressMessages(library(ArchR))
addArchRThreads(threads = 10)
addArchRGenome("hg38")
system('mkdir -p plots')
##############################
### read in ArchR project ####
PROJDIR='../../../data/raw_data/hg38/Corces_2020'
LABEL='Corces2020_caudate'; GENOME = 'hg38'; 
ARCHDIR=file.path(PROJDIR,paste0('ArchR_',LABEL,'_labeled'))
proj = loadArchRProject(ARCHDIR, showLogo = F)

###################################################
## find mouse peaks that can be halpered to hg38 ##
mouse_hal2hg38_fn =file.path('../../../data/raw_data/','mm10') %>%
  list.files(path = ., full.names = T, recursive = T, 
             pattern = '.Mus_musculusToHomo_sapiens.HALPER.narrowPeak.gz')
mouse_hal2hg38_fn = mouse_hal2hg38_fn[
  ss(basename(mouse_hal2hg38_fn), '\\.') %in% c('BICCN_CP', 'Pfenning_Cpu')]
# exclude BICCN_CP.MSN_D1 and BICCN_CP.MSN_D2 & halper peaks for now
mouse_hal2hg38_fn = mouse_hal2hg38_fn[!grepl('BICCN_CP.MSN_D|BICCN_CP.INT', mouse_hal2hg38_fn)]
names(mouse_hal2hg38_fn) = paste(basename(mouse_hal2hg38_fn) %>% ss('\\.', 2),
                                  'Mm',sep = '_')
mouse_hal2hg38_fn = mouse_hal2hg38_fn[order(names(mouse_hal2hg38_fn))]

#####################################################
## find macaque peaks that can be halpered to hg38 ##
macaque_hal2hg38_fn =file.path('../../../data/raw_data/','rheMac10') %>%
  list.files(path = ., full.names = T, recursive = T, 
             pattern = '.Macaca_mulattaToHomo_sapiens.HALPER.narrowPeak.gz')
names(macaque_hal2hg38_fn) = paste(basename(macaque_hal2hg38_fn) %>% ss('\\.', 2),
                                   'Rm',sep = '_')
macaque_hal2hg38_fn = macaque_hal2hg38_fn[order(names(macaque_hal2hg38_fn))]

####################################################################################
## find ArchR marker Peaks and enrichment in Mouse, Macaque, and ChromHMM regions ##
markersPeaks <- getMarkerFeatures( ArchRProj = proj, useMatrix = "PeakMatrix", 
                                   groupBy = "Clusters2", testMethod = "wilcoxon", 
                                   bias = c("TSSEnrichment", "log10(nFrags)"))

# merge the peakSets and annotate ArchR project
regionCustom_fn = c(mouse_hal2hg38_fn, macaque_hal2hg38_fn)
regionCustom_fn = regionCustom_fn[!grepl('Consensus', names(regionCustom_fn))]

proj <- addPeakAnnotations( proj, regions = regionCustom_fn, 
                            name = "Caudate_cons_chrom", force = T)

##################################################################
## make enrichment plots of mouse/macaque cell types w/ human snATAC
system(paste('mkdir -p', file.path( 'plots')))
height_ppt = 4.5; width_ppt = 8
height_fig = 5; width_fig = 4.75; font_fig = 5
height_fig2 = 3; width_fig2 = 2.25

enrichRegions <- peakAnnoEnrichment(seMarker = markersPeaks, ArchRProj = proj, 
                                    peakAnnotation = "Caudate_cons_chrom", 
                                    background= "all",
                                    cutOff = "FDR <= 0.05 & Log2FC >= 1")

plot_fn = 'plots/Corces2020_markerPeak_caudate_conservation_enrichment.ppt.pdf'
pdf(plot_fn, height = height_ppt, width = width_ppt)
heatmapRegions <- plotEnrichHeatmap(enrichRegions, n = 7, transpose = T)
ComplexHeatmap::draw(heatmapRegions, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

plot_fn = 'plots/Corces2020_markerPeak_caudate_conservation_enrichment.fig.pdf'
pdf(plot_fn, height = height_fig2, width = width_fig2)
heatmapRegions <- plotEnrichHeatmap(enrichRegions, n = 7, transpose = F)
ComplexHeatmap::draw(heatmapRegions, heatmap_legend_side = "bot", 
                     annotation_legend_side = "top")
dev.off()



##########################################
## find the caudate ChromHMM annoations ##
epiMap_chromHmm_fn =file.path('../../../data/tidy_data/','epimap') %>%
  list.files(path = ., full.names = T, recursive = T, 
             pattern = 'BRN.CAUD.NUC_ChromHMM') 
names(epiMap_chromHmm_fn) = 
  basename(epiMap_chromHmm_fn) %>% ss('ChromHMM_',2) %>% ss('\\.')

proj <- addPeakAnnotations( proj, regions = epiMap_chromHmm_fn, 
                            name = "Caudate_ChromHMM", force = T)

enrichRegions2 <- peakAnnoEnrichment(seMarker = markersPeaks, ArchRProj = proj, 
                                     peakAnnotation = "Caudate_ChromHMM", 
                                     background= "bgdPeaks",
                                     cutOff = "FDR <= 0.05 & Log2FC >= .5")

plot_fn = 'plots/Corces2020_markerPeak_ChromHMM_enrichment.ppt.pdf'
pdf(plot_fn, height = height_ppt, width = width_ppt)
heatmapRegions2 <- plotEnrichHeatmap(enrichRegions2, pMax = 40,  transpose = T)
ComplexHeatmap::draw(heatmapRegions2, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

plot_fn = 'plots/Corces2020_markerPeak_ChromHMM_enrichment.fig.pdf'
pdf(plot_fn, height = height_fig2, width = width_fig2)
heatmapRegions2 <- plotEnrichHeatmap(enrichRegions2, pMax = 40, transpose = T)
ComplexHeatmap::draw(heatmapRegions2, heatmap_legend_side = "bot", 
                     annotation_legend_side = "top")
dev.off()


