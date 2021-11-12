#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F, bitmapType='cairo')
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
library(rcartocolor)
library(ggrepel)
library(ggforce)
library(ggpubr)
library(ggsci)
library(here)
library(ggh4x)

# to be run in the root github directory
PLOTDIR=here('figures/explanatory/','figure2_ldsc_cross-species_grid')

# read in the GWAS traits
load(here('data/tidy_data/ldsc_gwas','rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% select( -file) %>% 
  mutate(label = ss(as.character(trait), '_'))

################################################
## 1) Gather the OCR orthologs + phyloP alone
keepTraits = c('AD_E', 'SmokInitiation_E')
keepTraits = pheno %>% pull(trait)

## OCR orthologs
keepPeakTypes =  c('Human' = 'Hg peak', 'Human + Macaque' = 'hgRmOrth', 'Human + Mouse' = 'hgMmOrth')
keepPeakTypes2 = names(keepPeakTypes); names(keepPeakTypes2) = keepPeakTypes
enrichments_orthologs = here('figures/exploratory/ldsc_conservation','rdas',
                             'caudate_conservation_prop_herit_hg_rm_mm.rds') %>%
  readRDS() %>% mutate(datatype = 'ATAC-seq') %>%
  filter(trait %in% keepTraits, peaktype %in%keepPeakTypes) %>%
  mutate(peaktype = keepPeakTypes2[as.character(peaktype)],
         celltype = droplevels(celltype), celltype2 = celltype)
# enrichments_orthologs %>% count(celltype, celltype2)

## phylop only
keepPeakTypesPhyloP =  c('PhyloP.accl' = 'PhyloP.accl', 'PhyloP.cons' = 'PhyloP.cons')
keepPeakTypesPhyloP2 = names(keepPeakTypesPhyloP); names(keepPeakTypesPhyloP2) = keepPeakTypesPhyloP
celltypes = enrichments_orthologs %>% pull(celltype) %>% levels()
enrichments_phyloP = here('figures/exploratory/ldsc_conservation','rdas',
                             'caudate_conservation_prop_herit_hg_rm_mm.rds') %>%
  readRDS() %>% mutate(datatype = 'PhyloP') %>%
  filter(trait %in% keepTraits, peaktype %in%keepPeakTypesPhyloP) %>%
  mutate(peaktype = keepPeakTypesPhyloP2[as.character(peaktype)],
         celltype = 'PhyloP', celltype2 = 'PhyloP')
# enrichments_phyloP %>% count(celltype, celltype2)

## mappable to species
keepSpecies = c('Macaque' = 'Macaca_mulatta',
                'Mouse' = 'Mus_musculus',
                'Fruit bat' = 'Rousettus_aegyptiacus')
keepSpecies2 = names(keepSpecies)
names(keepSpecies2) = keepSpecies
enrichments_mappable = here('figures/exploratory/ldsc_caudate_zoonomia',
               'rdas','caudate_conservation_prop_herit_Zoonomia.rds') %>%
  readRDS() %>% mutate(datatype = 'Mappable to') %>%
  filter(trait %in% keepTraits, Species %in% keepSpecies) %>%
  mutate(peaktype = keepSpecies2[as.character(Species)],
         celltype2 = 'Mappable')
# enrichments_mappable %>% count(celltype, celltype2)

## mappable and overlaps + phyloP
enrichments_mapPhylop = here('figures/exploratory/ldsc_caudate_zoonomia',
                             'rdas','caudate_conservation_prop_herit_phyloP_Zoonomia.rds') %>%
  readRDS() %>% mutate(datatype = 'Mappable to') %>%
  filter(trait %in% keepTraits, Species %in% keepSpecies) %>%
  mutate(peaktype = keepSpecies2[as.character(Species)],
         celltype2 = 'PhyloP')
# enrichments_mapPhylop %>% count(celltype, celltype2)

## CellTACIT-ML, Caudate snATAC-seq cell types
keepGroupMeta = c('Human' = 'Primates#0', 'OW Monkey, 20MY' = 'Primates#19.8', 
                  'NW Monkey, 29MY' = 'Primates#28.81', 'Primates, 74MY' = 'Primates#74.1', 
                  'Rodentia, 89MY' = 'Rodentia#89', 'Chiroptera, 94MY' ='Chiroptera#94')
keepGroupMeta2 = names(keepGroupMeta)
names(keepGroupMeta2) = keepGroupMeta
enrichments_meta = here('figures/exploratory/', 'ldsc_zoonomia_meta',
                        'rdas','zoonomia_meta_prop_heritability_Corces2020.rds') %>%
  readRDS() %>% mutate(datatype = 'CellTACIT-ML') %>%
  filter(trait %in% keepTraits, group_meta %in% keepGroupMeta) %>%
  mutate(celltype2 = ifelse(peaktype != 'predActive', 'Mappable', as.character(celltype)),
         peaktype = keepGroupMeta2[as.character(group_meta)])
# enrichments_meta %>% count(celltype, celltype2)


## CellTACIT-ML, bulk M1 and liver ATAC-seq
enrichments_meta_M1ctxLiver = here('figures/exploratory/', 'ldsc_zoonomia_meta',
                        'rdas','zoonomia_meta_prop_heritability_M1ctx_liver.rds') %>%
  readRDS() %>% mutate(datatype = 'CellTACIT-ML') %>%
  filter(trait %in% keepTraits, group_meta %in% keepGroupMeta) %>%
  mutate(celltype2 = ifelse(peaktype != 'predActive', 'Mappable', as.character(celltype)),
         peaktype = keepGroupMeta2[as.character(group_meta)])
# enrichments_meta_M1ctxLiver %>% count(celltype2, celltype, peaktype)

################################################
## 2) Combine dataframes and reset x labels
peaktype2 = c(keepPeakTypesPhyloP2, keepPeakTypes2, keepSpecies2, keepGroupMeta2) %>% unique()
datatype2 = c('PhyloP', 'ATAC-seq', 'Mappable to', 'CellTACIT-ML')
celltypes2 =  c('Mappable', 'PhyloP', celltypes, 'M1ctx', 'Liver')
enrichments = list(enrichments_phyloP, enrichments_orthologs,  enrichments_mappable, 
                   enrichments_mapPhylop, enrichments_meta, enrichments_meta_M1ctxLiver) %>% 
  rbindlist(fill=TRUE) %>% 
  mutate(peaktype = factor(peaktype, peaktype2),
         datatype = factor(datatype, datatype2), 
         celltype2 = factor(celltype2, celltypes2),
         label = trait %>% as.character() %>% ss('_'))

cell_cols =c(carto_pal(8, "Safe"), 'darkgray', 'black','darkseagreen','darkslategray')
names(cell_cols) = c(celltypes,'PhyloP', 'Mappable', 'M1ctx', 'Liver')
enrichments %>% count(datatype, peaktype)
enrichments %>% count(datatype, celltype2)

dir.create(here(PLOTDIR,'rdas'), showWarnings = F)
save_fn = here(PLOTDIR,'rdas','figure2_heritability_enrichments_Corces2020_phyloP_CellTACIT-ML.rds')
saveRDS(enrichments, file = save_fn)