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

##############################################
## read in PGLS results
DATADIR='data/raw_data/trait_pgls'
alpha = 0.10
pgls_rds = here(DATADIR, 'rdas', 'trait_pgls_Corces2020_finemapped_snps_20210914.rds')
pgls_df = readRDS(pgls_rds) %>% filter(zooTrait != "Brain Residual") %>%
  filter(PGLS_FDR < alpha, grepl('MSN', celltype), grepl('Cross|Morn', label)) %>%
  arrange(PGLS_FDR) %>%
  mutate(name = gsub('^hg38:|:250$', '', peakNames) %>% ss('-', 1),
         name = factor(name, unique(name)))
with(pgls_df, table(celltype, zooTrait))

pgls_gr = pgls_df %>% pull(peakNames) %>% gsub(pattern = '^hg38:|:250$', replacement = '') %>%
  GRanges()
loci = pgls_gr %>% resize( 2e4, fix="start") %>% resize( 2*2e4, fix="end") %>% 
  GenomicRanges::reduce()
mcols(loci)$score = 1
## the PGLS trait-SNP loci
loci_track = new("track", dat=loci,  type="data", format = 'BED')


##################################
#### read the fine-mapped SNPs ###
DATADIR='data/raw_data'
finemap_df = here(DATADIR, 'polyfun_caudate/rdas/polyfun_caudate_finemapped_snps_20210518.rds') %>%
  readRDS() %>% mutate(start = POS_hg38, end = POS_hg38) %>%
  dplyr::select(-c(POS_hg19,POS_hg38, SNPVAR:P, BETA_MEAN:MAF, index:runPolyfun, h2:base))

gwas_group_legend = list(labels = names(group_col)[1:6], col = 'gray80', fill = group_col[1:6])
group_col2 = group_col[1:6]


###########################################################
## Zoo-meta tracks, file1 = predActive, file2 = mappable ##
celltype = 'MSN_D2'
fn= list.files(path = 'figures/explanatory/figure3_fine-mapped_gwas_loci/bed',
               pattern = celltype, full.names = T) %>%
  grep(pattern = 'allPeaks', value = TRUE)

zooMeta_tracks = split(fn, ss(basename(fn), '\\.', 4)) %>% sapply(function(ll){
  tr = importScore(ll[1], ll[1], format="BED")
  strand(tr$dat) <- "+"; strand(tr$dat2) <- "-"
  tr$dat$score = tr$dat$score * 100; tr$dat2$score = 20
  setTrackStyleParam(tr, "color", c(col_celltypes[celltype], "black"))
  return(tr)
})

zooMeta_gr = fn %>% grep(pattern = 'mappable', value = TRUE) %>%
  lapply(import) %>% GRangesList() %>% unlist()

zooMeta_gr[zooMeta_gr$name %in% pgls_df$peakNames] %>% mcols() %>% 
  as.data.frame() %>%pull(score)

reorder = df_meta$group_meta %>% as.character() %>% ss('\\.', 1) %>% rev()
zooMeta_tracks = zooMeta_tracks[reorder]
names(zooMeta_tracks) = gsub('#', ' ', names(zooMeta_tracks))
names(zooMeta_tracks) = gsub('$', 'MYA', names(zooMeta_tracks))

### annotate SNPs w/ Zoo meta peaks
oo = findOverlaps(query = zooMeta_gr, subject = finemap_df %>% 
                    mutate(CHR= paste0('chr', CHR)) %>% GRanges())
indList = split(queryHits(oo), subjectHits(oo))
finemap_df$numHits = 0
finemap_df$numHits[as.numeric(names(indList))] = lengths(indList)

## labeling of SNPs across traits
finemap_gr = finemap_df %>% arrange(PIP) %>%
  mutate( color = 'black', fill = group_col2[group], 
          SNPsideID = 'bottom', border = fill,
          height = .1, cex = .6, lwd = .8,
          shape = ifelse(numHits >0, 'square', 'circle')) %>% 
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
          shape = ifelse(numHits >0, 'square', 'circle')) %>%
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
names(finemap_gr2)[finemap_gr2$score < 1] = ''
seqlevelsStyle(finemap_gr2) = 'UCSC'
finemap_gr2$feature.height = .3
finemap_track = new("track", dat = finemap_gr2, dat2=finemap_gr, type="lollipopData")

################################
## Zoonomia constraint tracks ##
idx = 1
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
  ids <- getGeneIDsFromTxDb(locus, TxDb.Hsapiens.UCSC.hg38.knownGene)
  symbols <- mget(ids, org.Hs.egSYMBOL)
  indKeep = which(!grepl('^MIR|^SNOR|-AS[1-9]',symbols))
  genes = geneTrack(ids[indKeep], TxDb.Hsapiens.UCSC.hg38.knownGene, symbols[indKeep], type = 'gene')
  names(genes) = symbols[indKeep][names(genes)]

  
  ##########################
  ## 3) construct tracks  ##
  pdf('tmp2.pdf', height = .75*2, width = 4.75*2, onefile = F)
  trackList <- trackList(zooMeta_tracks["Primates 0MYA"],  finemap_track, 
                         loci_track, zoo_track, genes, heightDist = c(1, 1.5, 2, 1, 1))
  # names(trackList)[which(names(trackList)== 'finemap_track')]= 'Fine-mapped SNPs'
  optSty <- optimizeStyle(trackList)
  trackList <- optSty$tracks
  viewerStyle <- optSty$style
  dev.off()
  
  ## x-axis/scale parameters
  gparList = list(lwd = .5, cex = .7)
  setTrackViewerStyleParam(viewerStyle, "xaxis", TRUE)
  setTrackViewerStyleParam(viewerStyle, "xgp", gparList)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.2, .15, .01, .05))
  setTrackXscaleParam(trackList[[1]], "draw", FALSE)
  
  ## PhyloP/PhastCons y-axis parameters
  top_phyloP_score = ceiling(max(zoo_track[['PhyloP.200mam']]$dat$score))
  setTrackStyleParam(trackList[['PhyloP.200mam']], "ylim",c(0,top_phyloP_score))
  setTrackStyleParam(trackList[["PhastCons.43prim"]], "ylim",c(0,1))
  
  ## General Y-axis parameters
  for (i in c("Primates 0MYA")) setTrackYaxisParam(trackList[[i]],  "draw", FALSE)
  for(i in names(zoo_fn)) {
    setTrackYaxisParam(trackList[[i]], "main", FALSE)
    setTrackStyleParam(trackList[[i]], "color", "black")
    setTrackYaxisParam(trackList[[2]], "gp", list(lwd = .5, cex = .4)) 
  }
  for (i in names(genes)) {
    setTrackStyleParam(trackList[[i]], "ylabpos", "upstream")
    setTrackStyleParam(trackList[[i]], "color", "black")
    setTrackYaxisParam(trackList[[i]], "gp", list(lwd = .5, cex = .6))
  }
  for (i in seq_along(trackList)) {
    setTrackStyleParam(trackList[[i]], "ylabgp", gparList)
  }
  
  ## SNP tracks y-axis parameters
  setTrackYaxisParam(trackList[['finemap_track']], "gp",gparList)
  setTrackStyleParam(trackList[['finemap_track']], "ylabgp",gparList)
  
  ## loci tracks y-axis parameters, buffer track
  setTrackStyleParam(trackList[['loci_track']], "color",'#00000000')
  setTrackStyleParam(trackList[['loci_track']], "ylabgp", list(col="#00000000"))
  setTrackYaxisParam(trackList[['loci_track']], "gp", list(col="#00000000"))
  
  regionName = paste0('chr',chr,':', signif(start(locus)/1e6,4), '-',signif(end(locus)/1e6,4), 'Mb')
  geneNames = paste(sort(names(genes))[1:min(length(genes), 5)], collapse = ',')
  plot_fn = here('figures/explanatory/figure4_trait_pgls_loci/plots',
                 paste0('fig4_trackPlots.',celltype, '.',regionName, '.', geneNames,'.pdf'))
  pdf(plot_fn, height = .75*2, width = 4.75*2, onefile = F)
  viewTracks(trackList, gr=locus, viewerStyle=viewerStyle)
  dev.off()
  
}