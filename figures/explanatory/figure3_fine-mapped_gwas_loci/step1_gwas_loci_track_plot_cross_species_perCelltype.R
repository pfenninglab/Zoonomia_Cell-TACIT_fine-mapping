#######################################
### set up libraries and functions ####
# conda activate r4
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F, bitmapType='cairo')
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
library(rcartocolor)
library(here)
library(rtracklayer)
library(trackViewer)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)


##############################################
# 1a) read in the Zoonomia tree and group_meta list
rda_fn = here('data/tidy_data/Zoonomia_data', 
              'rdas','200_Mammals_Genome_Information.rda')
load(file = rda_fn)
df_meta = df_meta %>% mutate(Order2 = ifelse(Time.Since.Split.from.Human.TimeTree.median <=94, as.character(Order), as.character(Clade)))

col_clade = df %>% dplyr::select(Clade, col_clade)%>% filter(!duplicated(Clade)) %>% deframe()
col_meta = df_meta %>% mutate(value = Order2, name = col_meta) %>% 
  filter(!duplicated(value)) %>% dplyr::select(value,name) %>% deframe()

col_celltypes = rcartocolor::carto_pal(n = 8, 'Safe')
names(col_celltypes) = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb',  'Astro', 
                         'Microglia', 'OPC', 'Oligo')  

# 1b) read in the GWAS traits
load(here('data/tidy_data/ldsc_gwas','rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% dplyr::select( -file) %>% 
  mutate(label = ss(as.character(trait), '_'))

##################################
#### read the fine-mapped SNPs ###
DATADIR='data/raw_data'
finemap_df = here(DATADIR, 'polyfun_caudate/rdas/polyfun_caudate_finemapped_snps_20210518.rds') %>%
  readRDS() %>% mutate(start = POS_hg38, end = POS_hg38) %>%
  dplyr::select(-c(POS_hg19,POS_hg38, SNPVAR:P, BETA_MEAN:MAF, index:runPolyfun, h2:base))

gwas_group_legend = list(labels = names(group_col)[1:6], col = 'gray80', fill = group_col[1:6])
group_col2 = group_col[1:6]

## read in CellTACIT scores
save_track_fn= here('figures/explanatory/figure3_fine-mapped_gwas_loci/rdas',
                    paste0('fig3_trackViewer.CellTACIT_track.rds'))
if(file.exists(save_track_fn)){
  CellTACIT_track = readRDS(save_track_fn)
} else{
  fn2 = list.files(path = 'figures/explanatory/figure3_fine-mapped_gwas_loci/CellTACIT', full.names = T) %>%
    grep(pattern = 'mean', value = TRUE)
  names(fn2) = ss(basename(fn2), '\\.', 2)
  CellTACIT_track = fn2 %>% sapply(function(ll){
    tr = importScore(ll,ll, format="BED")
    annot_df <- annotatePeak(tr$dat, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, 
                             annoDb='org.Hs.eg.db') %>% 
      as.GRanges() %>% as.data.frame(row.names = seq(length(.))) %>%
      mutate(annotation = ss(annotation, ' ')) %>%
      mutate(peakNames = paste0('hg38:', seqnames,':', start, '-', end, ':250')) %>%
      dplyr::select(peakNames, annotation, SYMBOL, distanceToTSS)
    idx = which(annot_df$annotation %in% c('Distal', 'Intron', 'Promoter' ))
    tr$dat = tr$dat[idx];  tr$dat2 = tr$dat2[idx]
    tr$dat2$score = .1; strand(tr$dat2) <- "-"
    setTrackStyleParam(tr, "color", c(col_celltypes[ss(basename(ll), '\\.', 2)], "black"))
    return(tr) })
  CellTACIT_track = CellTACIT_track[rev(names(col_celltypes))]
  saveRDS(CellTACIT_track, file = save_track_fn)
}

## read in Genetic Map
genetic_map = 'data/tidy_data/DECODE/genetic.map.final.sexavg.gor.gz'
genetic_map_gr = fread(genetic_map) %>% dplyr::select(-cM) %>% 
  rename('cMperMb' = 'score', 'Begin' = 'start', 'End' = 'end', 'Chr' = 'seqnames') %>% GRanges()
genetic_map_track = new("track", dat=genetic_map_gr, type="data", format = 'bedGraph')

## plot for each cell type
celltype = 'MSN_D2'

for (celltype in names(col_celltypes)){
### annotate SNPs w/ CellTACIT scores meta peaks
oo = findOverlaps(query = CellTACIT_track[[celltype]]$dat, subject = finemap_df %>% 
                    mutate(CHR= paste0('chr', CHR)) %>% GRanges())
indList = split(queryHits(oo), subjectHits(oo))
finemap_df$CellTACIT_score = 0
finemap_df$CellTACIT_score[subjectHits(oo)] = mcols(CellTACIT_track[[celltype]]$dat)[queryHits(oo),'score']

## labeling of SNPs across traits
finemap_gr = finemap_df %>% arrange(PIP) %>%
  mutate( color = 'black', fill = group_col2[group], 
          SNPsideID = 'bottom', border = fill,
          height = .1, cex = .6, lwd = .8,
          shape = ifelse(CellTACIT_score >0, 'square', 'circle')) %>% 
  group_by(group, SNP, CHR, start, end, A1, A2) %>%
  mutate(score = length(unique(label))) %>% ungroup() %>%
  distinct(group, SNP, CHR, start, end, A1, A2, score, .keep_all = TRUE) %>%
  GRanges()
seqlevelsStyle(finemap_gr) = 'UCSC'

## group SNPs together across finemapping traits
finemap_df2 = finemap_df %>% arrange(group, subgroup, label, PIP) %>% 
  mutate( color = 'black', fill = group_col2[group], 
          border = group_col2[group], height = .1, cex = .6,lwd = 1,
          label.parameter.rot = 40, SNPsideID = 'top',
          shape = ifelse(CellTACIT_score >0, 'square', 'circle')) %>%
  group_by(group, SNP, CHR, start, end, A1, A2) %>%
  mutate(alpha = 100 * max(PIP), score = round(max(PIP, na.rm = TRUE) * 10),
         PIP = sum(PIP, na.rm = TRUE), 
         numTrait = length(unique(label)), 
         dashline.col = ifelse(score < 1, '#00000000', 'gray80'), 
         label = paste(label, collapse = ', '),
         numSubGroup = length(unique(subgroup)),
         subgroup = paste(unique(subgroup), collapse = ', ')) %>% 
  filter(score >= 1) %>%
  ungroup() %>% distinct(SNP, CHR, start, end, A1, A2, score, .keep_all = TRUE)
summary(finemap_df2$score)
summary(finemap_df2$numSubGroup)

## the fine-mapped SNP track
finemap_gr2 = finemap_df2 %>% 
  dplyr::select(-c (numTrait, numSubGroup, PIP, label)) %>% GRanges()
names(finemap_gr2)= with(mcols(finemap_gr2), paste0(SNP, ':', A1, ':', A2))
names(finemap_gr2)[finemap_gr2$score < 3] = ''
seqlevelsStyle(finemap_gr2) = 'UCSC'
finemap_gr2$feature.height = .3
finemap_track = new("track", dat=finemap_gr2, dat2 = finemap_gr, type="lollipopData")

## save all the files
save_snps_fn= here('figures/explanatory/figure3_fine-mapped_gwas_loci/rdas',
                   paste0('fig3_snpsInd_for_plots.', celltype,'.rds'))
saveRDS(finemap_df, file = save_snps_fn)

save_snps_fn2= here('figures/explanatory/figure3_fine-mapped_gwas_loci/rdas',
                    paste0('fig3_snpsAgg_for_plots.', celltype,'.rds'))
saveRDS(finemap_df2, file = save_snps_fn2)

#########################
## find some good loci ##
loci = finemap_gr %>% resize( 2*5e3, fix="start") %>% 
  resize( 2e4, fix="end") %>% reduce()
oo2 = findOverlaps(query= finemap_gr2, subject = loci)
inds= split(queryHits(oo2), subjectHits(oo2))
loci$score = 0
loci$score[as.numeric(names(inds))] = sapply(inds, function(ind){
  finemap_df2[ind,] %>% mutate(tmp = PIP * CellTACIT_score * numTrait) %>% pull(tmp) %>% sum()
})
loci = loci[loci$score > 40]
start(loci) = start(loci) - 5e4
end(loci) = end(loci) + 5e4
loci = loci[order(loci$score, decreasing = T)]
loci_track = new("track", dat=loci,  type="data", format = 'BED')

save_loci_fn= here('figures/explanatory/figure3_fine-mapped_gwas_loci/rdas',
                   paste0('fig3_loci_for_plots.', celltype,'.rds'))
saveRDS(loci, file = save_loci_fn)

################################
## Zoonomia constraint tracks ##
for(idx in seq_along(loci)) {
  print(paste('Plotting locus' ,idx,'of', length(loci), 'for', celltype,'.'))
  locus = loci[idx]
  constraint_folder = '/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/phyloP'
  
  chr = seqnames(locus) %>% as.character() %>% gsub(pattern = 'chr', replacement = '')
  phyloP_fn = file.path(constraint_folder, 'human-centered-200m-Feb2021',
                        paste0('200m_scoresPhyloP_20210214.', chr,'.bigWig'))
  phastCons_fn = file.path(constraint_folder, 'Primates_PhastCons_scores_hg38', 
                           paste0('43prim_PhastCons.', chr,'.bigWig'))

  zoo_fn = c(phyloP_fn, phastCons_fn)
  names(zoo_fn) = c('PhyloP.200mam', 'PhastCons.43prim')
  
  zoo_track <- zoo_fn %>% sapply(importScore, format="BigWig", ranges = locus) 
  zoo_track[['PhyloP.200mam']]$dat <- zoo_track[['PhyloP.200mam']]$dat[zoo_track[['PhyloP.200mam']]$dat$score > 0]
  
  ##################
  ## gene tracks ##
  ids <- trackViewer::getGeneIDsFromTxDb(locus, TxDb.Hsapiens.UCSC.hg38.knownGene)
  symbols <- mget(ids, org.Hs.egSYMBOL)
  indKeep = which(!grepl('^MIR|^SNOR|-AS[1-9]',symbols))
  genes = geneTrack(ids[indKeep], TxDb.Hsapiens.UCSC.hg38.knownGene, symbols[indKeep], type = 'gene')
  names(genes) = symbols[indKeep][names(genes)]
  
  oo3 = findOverlaps(subject = finemap_gr2, query = locus )
  CellTACIT_scoreMax = finemap_df2$CellTACIT_score[subjectHits(oo3)] %>% max() %>% ceiling()
  
  oo4 = findOverlaps(subject = genetic_map_gr, query = locus )
  recomb_max = mcols(genetic_map_gr)$score[subjectHits(oo4)] %>% max() %>% ceiling()
  
  ##########################
  ## 3) construct tracks  ##
  height_fig = 1 * 2; width_fig = 2.25 * 2
  pdf('tmp2.pdf', height = height_fig, width = width_fig, onefile = F)
  trackList <- trackList( CellTACIT_track[celltype], finemap_track, 
                           loci_track, zoo_track, genetic_map_track, genes, 
                          heightDist = c(.75, 1, 1.5, 1, .5, 1))
  names(trackList)[which(names(trackList)== 'genetic_map_track')]= 'Recombination Map'
  names(trackList)[which(names(trackList)== 'finemap_track')]= 'Polyfun SNPs'
  optSty <- optimizeStyle(trackList)
  trackList <- optSty$tracks
  viewerStyle <- optSty$style
  dev.off()
  
  ## x-axis/scale parameters
  gparList = list(lwd = .6, cex = .6, hjust = 1, fontface = 'bold')
  setTrackViewerStyleParam(viewerStyle, "xaxis", TRUE)
  setTrackViewerStyleParam(viewerStyle, "xgp", gparList)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.13, .2, .01, .06))
  setTrackXscaleParam(trackList[[1]], "draw", FALSE)
  
  ## PhyloP/PhastCons y-axis parameters
  top_phyloP_score = ceiling(max(zoo_track[['PhyloP.200mam']]$dat$score))
  setTrackStyleParam(trackList[['PhyloP.200mam']], "ylim",c(0,top_phyloP_score))
  setTrackStyleParam(trackList[["PhastCons.43prim"]], "ylim",c(0,1))
  setTrackStyleParam(trackList[['Recombination Map']], "ylim",c(0, recomb_max))
  
  ## General Y-axis parameters
  for (i in celltype) { 
    setTrackStyleParam(trackList[[i]], "ylim",c(0, CellTACIT_scoreMax))
    setTrackYaxisParam(trackList[[i]], "main", FALSE)
  }
  for(i in c(names(zoo_fn), 'Recombination Map')) {
    setTrackStyleParam(trackList[[i]], "color", "black")
    setTrackYaxisParam(trackList[[i]], "main", FALSE)
  }
  for (i in names(genes)) {
    setTrackStyleParam(trackList[[i]], "ylabpos", "downstream")
    setTrackStyleParam(trackList[[i]], "color", "black")
  }
  for (i in seq_along(trackList)) {
    setTrackStyleParam(trackList[[i]], "ylabgp", gparList)
    setTrackYaxisParam(trackList[[i]], "gp", gparList)
    setTrackXscaleParam(trackList[[i]], "gp", gparList)
  }
  
  ## loci tracks y-axis parameters, buffer track
  setTrackStyleParam(trackList[['loci_track']], "color",'#00000000')
  setTrackStyleParam(trackList[['loci_track']], "ylabgp", list(col="#00000000"))
  setTrackYaxisParam(trackList[['loci_track']], "gp", list(col="#00000000"))

  regionName = paste0('chr',chr,'_', signif(start(locus)/1e6,4), '-',signif(end(locus)/1e6,4), 'Mb')
  geneNames = paste(sort(names(genes))[1:min(length(genes), 5)], collapse = ',')
  plot_fn = here('figures/explanatory/figure3_fine-mapped_gwas_loci/plots',
                 paste0('fig3_trackPlots.',celltype, '.',regionName, '.', geneNames,'.pdf'))
  pdf(plot_fn, height = height_fig, width = width_fig, onefile = F)
  viewTracks(trackList, gr=locus, viewerStyle=viewerStyle)
  dev.off()
}
}