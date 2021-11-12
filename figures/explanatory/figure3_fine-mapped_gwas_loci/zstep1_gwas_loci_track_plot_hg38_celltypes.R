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

############################
#### read the human OCRs ###
hg38_peaks_fn = list.files(path = 'data/raw_data/hg38/Corces_2020/peak/', 
                           pattern = '.narrowPeak.gz', full.names = TRUE) %>%
  grep(pattern = 'Corces2020_caudate\\.', value = TRUE)
names(hg38_peaks_fn) = ss(basename(hg38_peaks_fn), '\\.', 2)

hg38_peaks_tracks = hg38_peaks_fn[names(col_celltypes)] %>% 
  lapply(fread, col.names = c(bedNames, LETTERS[1:4])) %>% 
  rbindlist(idcol = 'celltype') %>% dplyr::select(-all_of(LETTERS[1:4])) %>%
  mutate(featureLayerID = celltype, border = 'black', color = col_celltypes[celltype]) %>%
  split(.$celltype) %>% sapply(function(x) new("track", dat=GRanges(x), type="data", format="BED"))
hg38_peaks_tracks = hg38_peaks_tracks[rev(names(col_celltypes))]

##################################
#### read the fine-mapped SNPs ###
DATADIR='data/raw_data'
finemap_df = here(DATADIR, 'polyfun_caudate/rdas/polyfun_caudate_finemapped_snps_20210518.rds') %>%
  readRDS() %>% mutate(start = POS_hg38, end = POS_hg38) %>%
  dplyr::select(-c(POS_hg19,POS_hg38, SNPVAR:P, BETA_MEAN:MAF, index:runPolyfun, h2:base))

gwas_group_legend = list(labels = names(group_col), col = 'gray80', fill = group_col)
group_col2 = group_col

finemap_df2 = finemap_df %>% filter(PIP > 0) %>% 
  arrange(group, subgroup, label) %>% mutate( color = 'black') %>%
  group_by(SNP, CHR, start, end, A1, A2) %>%
  mutate(score = 10 * sum(PIP, na.rm = TRUE),
         PIP = max(PIP), alpha = 100 * PIP,
         numTrait = length(unique(label)),  
         label = paste(label, collapse = ', '),
         numSubGroup = length(unique(subgroup)),
         subgroup = paste(unique(subgroup), collapse = ', ')) %>% 
  filter(score > 1) %>%
  ungroup() %>% distinct(SNP, CHR, start, end, A1, A2, score, .keep_all = TRUE)
summary(finemap_df2$score)
summary(finemap_df2$numSubGroup)

## the fine-mapped SNP track
finemap_gr2 = finemap_df2 %>% 
  dplyr::select(-c (numTrait, numSubGroup, PIP, label)) %>% GRanges()
names(finemap_gr2)= with(mcols(finemap_gr2), paste0(SNP, ':', A1, ':', A2))
seqlevelsStyle(finemap_gr2) = 'UCSC'
finemap_gr2$feature.height = .05
finemap_track = new("track", dat=finemap_gr2, type="lollipopData")

ind = with(finemap_df2, which(numSubGroup >=4 & score > 10))
oo = findOverlaps(subject = finemap_gr2[ind], query = gr)
loci = finemap_gr2[ind][unique(subjectHits(oo))] %>% 
  resize( 2.5e4, fix="start") %>% resize( 5e4, fix="end") %>% reduce()


################################
## Zoonomia constraint tracks ##
constraint_folder = '/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/phyloP'
phyloP_fn = file.path(constraint_folder, 'human-centered-200m-Feb2021','200m_scoresPhyloP_20210214.11.bigWig')
phastCons_fn = file.path(constraint_folder, 'Primates_PhastCons_scores_hg38', '43prim_PhastCons.11.bigWig')
zoo_fn = c(phyloP_fn, phastCons_fn)
names(zoo_fn) = c('PhyloP', 'PhastCons')

zoo_track <- zoo_fn %>% sapply(importScore, format="BigWig", ranges = loci[1]) 
zoo_track[['PhyloP']]$dat <- zoo_track[['PhyloP']]$dat[zoo_track[['PhyloP']]$dat$score > 0]


##################
## gene tracks ##
ids <- getGeneIDsFromTxDb(loci[1], TxDb.Hsapiens.UCSC.hg38.knownGene)
symbols <- mget(ids, org.Hs.egSYMBOL)
genes = geneTrack(ids, TxDb.Hsapiens.UCSC.hg38.knownGene, symbols, type = 'gene')
names(genes) = symbols[names(genes)]


##########################
## 3) construct tracks  ##
pdf('tmp2.pdf', height = 5, width = 2.25, onefile = F)
trackList <- trackList( zoo_track, hg38_peaks_tracks, finemap_track, genes, 
                        heightDist = c(1,4,1.5,1))
# names(trackList)[which(names(trackList)== 'finemap_track')]= 'Fine-mapped SNPs'
optSty <- optimizeStyle(trackList, theme="bw")
trackList <- optSty$tracks
viewerStyle <- optSty$style

setTrackViewerStyleParam(viewerStyle, "xaxis", TRUE)
setTrackViewerStyleParam(viewerStyle, "margin", c(.05, .05, .01, .05))
top_phyloP_score = ceiling(max(zoo_track[['PhyloP']]$dat$score))
setTrackStyleParam(trackList[['PhyloP']], "ylim",c(0,top_phyloP_score))
setTrackStyleParam(trackList[["PhastCons"]], "ylim",c(0,1))

for (i in seq_along(trackList)) setTrackStyleParam(trackList[[i]], "ylabpos", "top")
for (g in names(genes)) setTrackStyleParam(trackList[[g]], "height", .05/length(genes))
for (n2 in names(col_celltypes)) {
  setTrackStyleParam(trackList[[n2]], "color", col_celltypes[n2])
}

viewTracks(trackList, gr=loci[1], viewerStyle=viewerStyle)
dev.off()




