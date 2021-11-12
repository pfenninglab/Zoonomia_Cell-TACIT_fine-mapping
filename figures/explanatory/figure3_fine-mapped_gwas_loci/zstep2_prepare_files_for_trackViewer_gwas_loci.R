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


##############################################
# 1a) read in the Zoonomia tree and group_meta list
rda_fn = here('data/tidy_data/Zoonomia_data', 
              'rdas','200_Mammals_Genome_Information.rda')
load(file = rda_fn)
df_meta = df_meta %>% mutate(Order2 = ifelse(Time.Since.Split.from.Human.TimeTree.median <=94, as.character(Order), as.character(Clade)))

col_clade = df %>% select(Clade, col_clade)%>% filter(!duplicated(Clade)) %>% deframe()
col_meta = df_meta %>% mutate(value = Order2, name = col_meta) %>% 
  filter(!duplicated(value)) %>% dplyr::select(value,name) %>% deframe()

col_celltypes = rcartocolor::carto_pal(n = 8, 'Safe')
names(col_celltypes) = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb',  'Astro', 
                         'Microglia', 'OPC', 'Oligo')  

############################
# 1b) read in the GWAS traits
load(here('data/tidy_data/ldsc_gwas','rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% select( -file) %>% 
  mutate(label = ss(as.character(trait), '_'))


##############################################
## 2) read in the polyfun, fine-mapped SNPs ##
DATADIR='data/raw_data'
finemap_df = here(DATADIR, 'polyfun_caudate/rdas/polyfun_caudate_finemapped_snps_20210518.rds') %>%
  readRDS() 

finemap_df %>% pull(PIP) %>% summary()
finemap_df2 = finemap_df %>% filter(PIP > 0.1) %>%
  arrange(group, subgroup, label) %>%
  group_by(SNP, CHR, POS_hg38, A1, A2) %>%
  mutate(PIP = max(PIP, na.rm = TRUE),
         numTrait = n(), label = paste(label, collapse = ', '),
         numGroup = length(unique(group)), group = paste(unique(group), collapse = ', '),
         numSubGroup = length(unique(subgroup)),  subgroup = paste(unique(subgroup), collapse = ', ')) %>% 
  ungroup() %>%  arrange(desc(numTrait), desc(PIP)) %>%
  distinct(SNP, CHR, POS_hg38, A1, A2, PIP, numTrait, .keep_all = TRUE)

table(finemap_df2$numTrait)
finemap_df2 %>% group_by(numTrait) %>%
  summarize(mean = mean(PIP), med = median(PIP), max = max(PIP))

finemap_df2 %>% filter(grepl('SUD', group)) %>%
  count(group) %>% arrange(desc(n))

finemap_df2 %>% filter(grepl('Sleep', group)) %>%
  count(group)%>% arrange(desc(n))

finemap_gr2 = finemap_df2 %>% filter(CHR==1) %>%
  dplyr::select(-c(population:base)) %>%
  mutate(start = POS_hg38, end = POS_hg38, score = numTrait , 
         thick = PIP, border = 'black', color = 'black', alpha = round(PIP *100)) %>% GRanges()
names(finemap_gr2)= with(mcols(finemap_gr2), paste0(SNP, ':', A1, ':', A2))
seqlevelsStyle(finemap_gr2) = 'UCSC'


###########################################################
# read in SNPs overlapping enhancer peaks across species ##
finemapEnh_df = here(DATADIR, 'polyfun_caudate/rdas/polyfun_caudate_finemapped_snpsInEnhPeaks_20210518.rds') %>%
  readRDS()  %>% arrange(group, subgroup, label) %>%
  group_by(SNP, CHR, POS_hg38, A1, A2, celltype) %>%
  mutate(PIP = max(PIP, na.rm = TRUE),
         numTrait = n(), label = paste(label, collapse = ', '),
         numGroup = length(unique(group)), group = paste(unique(group), collapse = ', '),
         numSubGroup = length(unique(subgroup)),  subgroup = paste(unique(subgroup), collapse = ', ')) %>%
  ungroup() %>%  arrange(desc(PIP)) %>%
  filter(PIP > 0.10) %>%
  distinct(SNP, CHR, POS_hg38, A1, A2, PIP, numTrait, celltype, .keep_all = TRUE)
finemapEnh_df %>% count(CHR) # CHR3 most SNPs

table(finemapEnh_df$numTrait, finemapEnh_df$celltype)
finemapEnh_df %>% group_by(numTrait, celltype) %>%
  summarize(mean = mean(PIP), med = median(PIP), max = max(PIP)) %>%
  arrange(desc(max), desc(numTrait)) %>% as.data.frame()

finemapEnh_gr = finemapEnh_df %>% filter(celltype == 'MSN_D1', CHR==1) %>%
  dplyr::select(-c(population:base)) %>%
  mutate(start = POS_hg38, end = POS_hg38, score = round(PIP * 10) , 
         thick = PIP, border = 'black', color = 'black', alpha = round(PIP *100)) %>% GRanges()
names(finemapEnh_gr)= with(mcols(finemapEnh_gr), paste0(SNP, ':', A1, ':', A2))
seqlevelsStyle(finemapEnh_gr) = 'UCSC'

finemap_loci = finemapEnh_gr %>% resize( 2e4, fix="start") %>% 
  resize( 4e4, fix="end") %>% reduce()
finemap_loci %>% width() %>% summary()
seqlevelsStyle(finemap_loci) = 'UCSC'


#####################################################
## read the OCRs mappable across the Zoonomia tree ##
peaksMappable_fn = list.files(here('data/raw_data/ldsc_zoonomia_meta/peaks'),
                              pattern = '.predActive.bed.gz', full.names = T) %>%
  grep(pattern = 'Corces2020.MSN_D1.allPeaks', value = T)
names(peaksMappable_fn) = gsub('Corces2020.|.allPeaks|.enhPeaks|.predActive.bed.gz','', basename(peaksMappable_fn))

peaksMappable_fn2= gsub('data/raw_data/ldsc_zoonomia_meta/peaks', 
                        'figures/explanatory/figure3_fine-mapped_gwas_loci/bed', 
                        peaksMappable_fn)
peaksMappable_fn2 = gsub('.bed.gz', '.bed', peaksMappable_fn2)
if(FALSE){
  thecall= paste('zcat', peaksMappable_fn, '>',peaksMappable_fn2)
  sapply(thecall, system)
}
bedNames = c('seqnames', 'start', 'end', 'name', 'score', 'strand')

gr = peaksMappable_fn2 %>% lapply(fread, col.names = bedNames) %>% 
  rbindlist() %>% 
  group_by(name) %>%  mutate(score = sum(score), height = score * .001) %>% ungroup() %>% 
  distinct(name, .keep_all = TRUE) %>% 
  mutate(alpha = round(score/ max(score) *100),
         fill = paste0(col_celltypes[1], stringr::str_pad(alpha, 2, pad = 0))) %>%
  GRanges()


col_pal = scales::colour_ramp(colorspace::sequential_hcl(11, "Plasma"))
gr2 = peaksMappable_fn2 %>% lapply(fread, col.names = bedNames) %>% 
  rbindlist(idcol = 'id') %>% 
  mutate(alpha = round(score*100), 
         featureLayerID = factor(ss(id, '\\.|#',2), names(col_meta)), 
         border = col_meta[featureLayerID], 
         fill = col_pal(score/max(score))) %>%
  GRanges()


pdf('tmp.pdf', width =12,height = 5)
ind = subjectHits(findOverlaps(query = finemapEnh_gr, subject = finemap_loci))
lolliplot(SNP.gr = finemap_gr2, features = gr2, ranges = finemap_loci[2], 
          newpage = FALSE)
dev.off()





peaksMappable_grList = peaksMappable_fn %>% sapply(import) %>% 
  sapply(function(gr){
    gr = gr %>% as.data.frame() %>%
      mutate(fill = col_celltypes[1], 
             alpha = round(score *100)) %>%
      GRanges()
    return(gr)    
  }) %>% GRangesList()

peaksMappable_tracks = peaksMappable_grList %>%  
  sapply(function(gr) new("track", dat=gr, type="data", format="BED")) %>%
  trackList()


##################
## gene tracks ##
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

ids <- getGeneIDsFromTxDb(finemap_loci, TxDb.Hsapiens.UCSC.hg38.knownGene)
symbols <- mget(ids, org.Hs.egSYMBOL)
genes = geneTrack(ids, TxDb.Hsapiens.UCSC.hg38.knownGene, symbols, type = 'gene')
names(genes) = symbols[names(genes)]

## Zoonomia constraint tracks ##
constraint_folder = '/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/phyloP'
phyloP_fn = file.path(constraint_folder, 'human-centered-200m-Feb2021','200m_scoresPhyloP_20210214.2.bigWig')
phastCons_fn = file.path(constraint_folder, 'Primates_PhastCons_scores_hg38', '43prim_PhastCons.2.bigWig')
zoo_fn = c(phyloP_fn, phastCons_fn)
names(zoo_fn) = c('PhyloP.200mam', 'PhastCons.43prim')

zoo_track <- zoo_fn %>% sapply(importScore, format="BigWig", ranges = finemap_loci) 
zoo_track['PhyloP.200mam']$dat <- zoo_track['PhyloP.200mam']$dat[zoo_track['PhyloP.200mam']$dat$score > 0]

##########################
## 3) construct tracks  ##
trackList <- trackList(peaksMappable_tracks, zoo_track, genes)

pdf('tmp2.pdf', width = 12, height = 5, onefile = F)

optSty <- optimizeStyle(trackList)
trackList <- optSty$tracks
viewerStyle <- optSty$style

setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .05, .01, .01))
setTrackXscaleParam(trackList[[1]], "draw", TRUE)
setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.8))

setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .05, .01, .05))
for(i in seq_along(peaksMappable_tracks)){
  setTrackYaxisParam(trackList[[i]], "main", FALSE)
}

viewTracks(trackList, gr=finemap_loci[which.max(table(ind))], viewerStyle=viewerStyle)

dev.off()




