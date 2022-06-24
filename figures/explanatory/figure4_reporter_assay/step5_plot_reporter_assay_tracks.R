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

PLOTDIR='figures/explanatory/figure4_reporter_assay'

## the candidate enhancers
features = GRanges("chr11", IRanges(c(113567061, 113577668), width=c(501, 501),
                                    names=c("Candidate A", 'Candidate B'), 
                                    fill = c('#a6cee3', "#b2df8a"), 
                                    height = c(.1, .1), score = c(1, 1)))

## region around the peaks to plot, chr11:113562061-113587061
locus = GRanges("chr11", IRanges(c(113562061), width=c(25000)))

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

## group SNPs together across finemapping traits
finemap_df2 = finemap_df %>% arrange(group, subgroup, label, PIP) %>% 
  mutate( color = group_col2[group], fill = group_col2[group], 
          border = 'black', height = .1, cex = .6,lwd = 1,
          label.parameter.rot = 40, SNPsideID = 'top',
          shape = 'circle') %>%
  group_by(group, SNP, CHR, start, end, A1, A2) %>%
  mutate(score = round(max(PIP, na.rm = TRUE) * 100),
         PIP = sum(PIP, na.rm = TRUE), 
         numTrait = length(unique(label)), 
         dashline.col = '#00000000', 
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
finemap_track = new("track", dat=finemap_gr2, type="lollipopData")

###################################################
## make a lollipop plot of the SNP and enhancers ##
height_fig = 1 * 1; width_fig = 2.25 * 2
yaxis <- seq(0, 50, 10) ## define the position
names(yaxis) <- yaxis/100 # define the labels
xaxis <- 113565000 + seq(0, 20000, 10000) ## define the position
names(xaxis) <- xaxis # define the labels

pdf(here(PLOTDIR, 'plots', 'track_plot_snps_lollipop.pdf'), height = 2.5, width = width_fig)
lolliplot(finemap_gr2, features, locus, 
          yaxis = yaxis, xaxis = xaxis, 
          ylab = 'PIP', legend = F)
dev.off()


##################################
## make the conservation tracks ##
chr= 'chr11'
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

##########################
## 3) construct tracks  ##
pdf('tmp2.pdf', height = height_fig, width = width_fig, onefile = F)
trackList <- trackList(zoo_track, heightDist = c(1))
# names(trackList) = c('241 Mammal PhyloP', '43 Primate PhastCons')
optSty <- optimizeStyle(trackList)
trackList <- optSty$tracks
viewerStyle <- optSty$style
dev.off()

## x-axis/scale parameters
gparList = list(lwd = .6, cex = .6, hjust = 1, fontface = 'bold')
setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
setTrackViewerStyleParam(viewerStyle, "xgp", gparList)
setTrackViewerStyleParam(viewerStyle, "margin", c(.13, .2, .01, .06))
setTrackXscaleParam(trackList[[1]], "draw", FALSE)

## PhyloP/PhastCons y-axis parameters
top_phyloP_score = ceiling(max(zoo_track[['PhyloP.200mam']]$dat$score))
setTrackStyleParam(trackList[['PhyloP.200mam']], "ylim",c(0,top_phyloP_score))
setTrackStyleParam(trackList[["PhastCons.43prim"]], "ylim",c(0,1))

for (i in seq_along(trackList)) {
  setTrackStyleParam(trackList[[i]], "ylabgp", gparList)
  setTrackYaxisParam(trackList[[i]], "gp", gparList)
  setTrackXscaleParam(trackList[[i]], "gp", gparList)
}

###
plot_fn = here(PLOTDIR,'plots', paste0('track_plot_snps_conservation.pdf'))
pdf(plot_fn, height = height_fig, width = width_fig, onefile = F)
vp = viewTracks(trackList, gr=locus, viewerStyle=viewerStyle)
addGuideLine(c(113567311, 113577918), vp=vp)
dev.off()

