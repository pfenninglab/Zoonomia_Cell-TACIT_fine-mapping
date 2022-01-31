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
library(GenomicRanges)
library(ChIPseeker)
library(GenomicFeatures)

##############################################
# 0) human transcriptome annotation via gencode
sqlite_fn = '/home/bnphan/resources/genomes/GRCh38.p13/gencode.v35.annotation.sqlite'
if(!file.exists(sqlite_fn)){
  gff_fn = '/home/bnphan/resources/genomes/GRCh38.p13/gencode.v35.annotation.gff3'
  txdb = makeTxDbFromGFF(gff_fn)
  saveDb(txdb, file = sqlite_fn)
} else {
  txdb = loadDb(sqlite_fn)
}

##############################################
# 1a) read in the Zoonomia tree and group_meta list
rda_fn = here('data/tidy_data/Zoonomia_data', 
              'rdas','200_Mammals_Genome_Information.rda')
load(file = rda_fn)
labels = c('Dermoptera' = 'Flying Lemurs', 
           'Scandentia' = 'Treeshrews', 
           'Cetartiodactyla' = 'Ungulates', 
           'Eulipotyphla' = 'Hedgehogs')
df_meta = df_meta %>% 
  mutate(Order2 = ifelse(Time.Since.Split.from.Human.TimeTree.median <=94, as.character(Order), as.character(Clade)), label = ss(Species, ',', 1), 
         plotGroup = paste0(plotGroup,' ', ss(as.character(group_meta), '#', 2),'MYA' ))

col_clade = df %>% dplyr::select(Clade, col_clade)%>% filter(!duplicated(Clade)) %>% deframe()
col_meta = df_meta %>% mutate(value = Order2, name = col_meta) %>% 
  filter(!duplicated(value)) %>% dplyr::select(value,name) %>% deframe()
clade_list = with(df, split(Species, Clade))

col_celltypes = rcartocolor::carto_pal(n = 8, 'Safe')
names(col_celltypes) = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb',  'Astro', 
                         'Microglia', 'OPC', 'Oligo')  

#################################################################
# 1b) read in the SNPs w/ CellTACIT scores, phyloP annotations, etc
## read in the CellTACIT Scores per SNP
PLOTDIR='figures/explanatory/figure3_fine-mapped_gwas_loci'
finemap_df = readRDS(here(PLOTDIR,'rdas', 'Data_S6_CellTACIT_Age_caudate_finemapped_snps_20220119.rds'))

## read in the CNN predictions per cell type and allele
poly_fn2 = here('data/raw_data/polyfun_caudate/rdas',
                'polyfun_caudate_finemapped_snps_with_CNN_predictions_20220119.rds')
snps_prediction_df2 = readRDS(file = poly_fn2)
finemap_df = full_join(finemap_df, snps_prediction_df2) %>% 
  filter(!is.na(CHR), !is.na(start)) %>% 
  dplyr::select(-c(population:group_col, POS_hg38:Caud_hgMmOrth.OPC))
finemap_gr = GRanges(finemap_df)
seqlevelsStyle(finemap_gr) = 'UCSC'

## annotate SNPs w/ transcriptome again
annoPeaks = annotatePeak(finemap_gr, tssRegion = c(-3000, 3000), 
                         TxDb = txdb, annoDb="org.Hs.eg.db") %>% as.data.frame()
annoPeaks$annot = ss(annoPeaks$annotation, ' ')
table(ss(annoPeaks$annotation, ' '))
finemap_df$annot = annoPeaks$annot[match(finemap_df$name, annoPeaks$name)] 

######################################################
# 1c) read in the zoo meta predictions for these peaks
cells = c('MSN_D1', 'MSN_D2')
DATADIR= 'data/raw_data/ldsc_zoonomia_meta/peaks'
groups = c('Primates#0', 'Rodentia#89', 'Chiroptera#94')
pred_fn = file.path(DATADIR, paste('Corces2020',rep(cells, each = 3),'allPeaks', 
                                   rep(groups, times = 2),'mappable.bed.gz', sep ='.') )
names(pred_fn) = paste0(rep(cells, each = 3),'.', rep(groups, times = 2))
perGp_pred_df = lapply(pred_fn, fread, col.names = c('chr', 'start', 'end','name',  'score',  'strand')) %>% 
  rbindlist(idcol = 'lab') %>%
  mutate(celltype = ss(lab, '\\.', 1) ,Zoo_meta = ss(lab, '\\.', 2)) %>% 
  dplyr::select(-lab) %>%
  pivot_wider(values_from = score, names_from = Zoo_meta)

#########################################################
# 1d) read in the per species predictions for these peaks
species = c('Homo_sapiens', 'Mus_musculus', 'Rousettus_aegyptiacus')
DATADIR2= 'data/raw_data/ldsc_caudate_zoonomia/predictions'
pred_fn2 = file.path(DATADIR2, paste('Corces2020_caudate',rep(cells, each = 3),
                                     rep(species, times = 2),'avgCNN.predictions.txt.gz', sep ='.') )
names(pred_fn2) = paste0(rep(cells, each = 3),'.', rep(species, times = 2))
perSp_pred_df = lapply(pred_fn2, fread, col.names = c('name', 'score')) %>% 
  rbindlist(idcol = 'lab') %>%
  mutate(celltype = ss(lab, '\\.', 1) ,Species = ss(lab, '\\.', 2)) %>% 
  dplyr::select(-lab) %>%
  pivot_wider(values_from = score, names_from = Species)

## join the two species and group predictions together
## average prediction in species and in clade
pred_df = full_join(perGp_pred_df, perSp_pred_df) %>% 
  rename('name' = 'peakName') %>%
  rowwise() %>%
  mutate(
    Homo_sapiens = mean(c(Homo_sapiens, `Primates#0`), na.rm = T),
    Mus_musculus = mean(c(Mus_musculus, `Rodentia#89`), na.rm = T),
    Rousettus_aegyptiacus = mean(c(Rousettus_aegyptiacus, `Chiroptera#94`), na.rm = T),
  ) %>% dplyr::select(-contains('#'))

##################################################
# 2) read in the loci prioritized w/ CellTACIT age
loci_fn = list.files(here(PLOTDIR, 'rdas'), pattern =  'loci_for_plots', full.names = T)
names(loci_fn) = basename(loci_fn) %>% ss('\\.',2 )
loci_gr = lapply(loci_fn, readRDS) %>% GRangesList() %>% unlist() %>% reduce()
names(loci_gr) = paste0(seqnames(loci_gr),':', start(loci_gr), '-', end(loci_gr))

## add the loci the SNPs map to
finemap_df$loci = NA
gr1 = GRanges(finemap_df %>% dplyr::select('CHR', 'start', 'end'))
seqlevelsStyle(gr1) <- "UCSC" 
oo = findOverlaps(gr1, loci_gr)
finemap_df$loci[queryHits(oo)] = names(loci_gr)[subjectHits(oo)]

## add the peaks in other species predictions that SNPs lie in
finemap_dfList = split(finemap_df, finemap_df$celltype)
finemap_grList = split(finemap_gr, finemap_df$celltype)
pred_dfList = split(pred_df, pred_df$celltype)
pred_grList = split(GRanges(pred_df), pred_df$celltype)
ooList = mapply(findOverlaps, subject = pred_grList, query = finemap_grList[names(pred_grList)])

finemap_df2 = lapply(names(pred_grList), function(cell){
 df = finemap_dfList[[cell]]
 df[queryHits(ooList[[cell]]), c('peakName', species)] = 
   pred_dfList[[cell]][subjectHits(ooList[[cell]]),c('peakName', species)]
 return(df) 
}) %>% rbindlist() %>%
  mutate(group = case_when(
    Rousettus_aegyptiacus > 0.5 & Mus_musculus > 0.5 & Homo_sapiens > 0.5 ~ 'High',
    (Rousettus_aegyptiacus < 0.5 | is.na(Rousettus_aegyptiacus)) & 
      Mus_musculus > 0.5 & Homo_sapiens > 0.5 ~ 'Med',
    (Rousettus_aegyptiacus < 0.5 | is.na(Rousettus_aegyptiacus)) & 
      (Mus_musculus < 0.5 | is.na(Mus_musculus)) & Homo_sapiens > 0.5 ~ 'Low',
    TRUE ~ 'Other'
  ))

table(finemap_df2$group, finemap_df2$celltype)

## SNP selection criteria ##
finemap_df2 = finemap_df2 %>% 
  filter(!is.na(loci), PIP > 0.01, phyloPcons, celltype %in% cells) %>% 
  group_by(loci, celltype, match) %>% 
  mutate(CellTACIT_AgeDiff = max(CellTACIT_Age) - min(CellTACIT_Age),
         scoreDiff = max(Effect + nonEffect) - min(Effect - nonEffect)) %>% 
  ungroup() %>% group_by(SNP) %>%
  mutate(ModelAge_diff = max(CellTACIT_Age) - min(CellTACIT_Age)) %>% 
  ungroup() %>% arrange(-sqrt(scoreDiff*CellTACIT_AgeDiff), start) %>% 
  data.table()

finemap_df2 %>% writexl::write_xlsx('tmp.xlsx')
table(paste0(finemap_df2$loci, finemap_df2$celltype)) %>% sort()


finemap_df %>% filter(SNP == 'rs16918024')
annoPeaks %>% filter(SNP == 'rs16918024')
annoPeaks %>% filter(SNP == 'rs7933981')







