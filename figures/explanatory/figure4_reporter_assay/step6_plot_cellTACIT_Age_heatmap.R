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
library(ggtree)
library(tidytree)
library(Hmisc)
library(treeio)
library(ArchR)
library(aplot)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(reshape)
library(scales)

PLOTDIR='figures/explanatory/figure4_reporter_assay'
DATADIR='data/raw_data'

enhancers = c('hg38:chr11:113567061-113567561:250' = 'A', 
              'hg38:chr11:113577668-113578168:250' = 'B')

## the candidate enhancers
features = GRanges("chr11", IRanges(c(113567061, 113577668), width=c(501, 501),
                                    names=c("Candidate A", 'Candidate B'), 
                                    fill = c('#a6cee3', "#b2df8a"), 
                                    height = c(.1, .1), score = c(1, 1)))

## region around the peaks to plot
locus = GRanges("chr11", IRanges(c(113562061), width=c(25000)))

##################################################
# 0) read in the proportion of each cell types
celltypes = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb','Astro', 'Microglia', 'OPC', 'Oligo') 
## pre-calculate the proportion of each cell type
proj = loadArchRProject(here('data/raw_data/hg38/Corces_2020/ArchR_Corces2020_caudate_labeled'))
proportion_df = getCellColData(proj) %>% as_tibble() %>%
  group_by(Sample, Clusters2) %>% 
  summarise(isNeuronal = case_when(grepl('MSN|INT', Clusters2) ~ 'NeuN+', 
                                   TRUE ~ 'NeuN-'), 
            num = n()) %>% 
  distinct(Sample, Clusters2, isNeuronal, .keep_all = TRUE) %>%
  group_by(Sample, isNeuronal) %>% mutate(prop = num/sum(num)) %>% 
  group_by(isNeuronal, Clusters2) %>% summarise(prop = mean(prop)) %>%
  dplyr::rename('Celltype' = 'Clusters2')


##############################################
# 1a) read in the Zoonomia tree and group_meta list
rda_fn = here('data/tidy_data/Zoonomia_data', 
              'rdas','200_Mammals_Genome_Information.rda')
load(file = rda_fn)
labels = c('Dermoptera' = 'Flying Lemurs', 'Scandentia' = 'Treeshrews', 
           'Cetartiodactyla' = 'Ungulates', 'Eulipotyphla' = 'Hedgehogs')
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
## read in the CellTACIT Scores
save_track_fn= here('figures/explanatory/figure3_fine-mapped_gwas_loci/rdas',
                    paste0('fig3_trackViewer.CellTACIT_track.rds'))
CellTACIT_track = readRDS(file = save_track_fn)


###############################################
## gather individual SNPs annotated w/ CellTACIT Age
save_snps_fn= list.files(here('figures/explanatory/figure3_fine-mapped_gwas_loci/rdas'),
                         pattern = 'fig3_snpsInd_for_plots', full.names = T)
names(save_snps_fn) = basename(save_snps_fn) %>% ss('\\.', 2)

indSNPList = lapply(save_snps_fn, function(x){
  out = readRDS(x)
  out = out %>% dplyr::rename('CellTACIT_Age' = 'CellTACIT_score') %>%
    dplyr::select(-c(population:end)) %>% filter(CellTACIT_Age> 0)
  return(out)
})
sapply(indSNPList, nrow) %>% summary()
sapply(indSNPList, nrow)

sapply(indSNPList, function(x) x %>% filter(PIP > .01) %>% nrow()) %>% summary()
sapply(indSNPList, function(x) x %>% filter(PIP > .01) %>% nrow())


###############################################
## subset the Zoonomia tree to GroupMeta labels
tree2 = read.nhx(textConnection('((Loxodonta_africana:0.1204058581,Tolypeutes_matacus:0.1229192742):0.01145639324,(((Oryctolagus_cuniculus:0.1608919201,Rattus_norvegicus:0.2965409366):0.006031885449,(Tupaia_chinensis:0.1385286196,(Galeopterus_variegatus:0.08750990202,(Eulemur_flavifrons:0.07870638791,(Callithrix_jacchus:0.05610238745,(Macaca_mulatta:0.02760672357,(Nomascus_leucogenys:0.01503824207,(Pongo_abelii:0.01184875793,(Gorilla_gorilla:0.006121144092,(Pan_troglodytes:0.004607178674,Homo_sapiens:0.004410717579):0.00190302907):0.005963570605):0.002256803519):0.00775016763):0.01406957514):0.05984338905):0.01269656686):0.002845408602):0.002347302752):0.01705378614,(Uropsilus_gracilis:0.207400545,(Rousettus_aegyptiacus:0.1237216561,(Sus_scrofa:0.1368474811,(Equus_caballus:0.07947626741,(Canis_lupus_familiaris:0.1226195774,Manis_pentadactyla:0.1145549768):0.006550057057):0.002145230337):0.002284646667):0.005799996479):0.02130276027):0.01145639324);'))

tree_dt = full_join(tree2 %>% as_tibble(), df_meta) %>%
  filter(!is.na(node)) %>% as.treedata()
tree_dt = drop.tip(tree_dt, "Manis_pentadactyla")  
tree_dt =  groupOTU(tree_dt, clade_list, 'Clade2')

# 1b) read in the GWAS traits
load(here('data/tidy_data/ldsc_gwas','rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% dplyr::select( -file) %>% 
  mutate(label = ss(as.character(trait), '_'))

#################################
### get the CellTACITage values ###
CellTACITage_fn = here('data/raw_data/reporter_assay', 'rdas') %>% 
  list.files(pattern = 'CellTACIT.mean.rds', full.names = T)

names(CellTACITage_fn) = basename(CellTACITage_fn) %>% ss('\\.', 2)

pred_df = CellTACITage_fn %>% lapply(readRDS) %>% 
  lapply(function(gr) { gr$name = names(gr); return(gr)}) %>%
  lapply(as.data.frame) %>%
  rbindlist(idcol = 'Celltype') %>% 
  mutate(Enhancer =enhancers[name] %>% factor(),
         isNeuronal = case_when(grepl('MSN|INT', Celltype) ~ 'NeuN+', 
                                TRUE ~ 'NeuN-'), 
         isNeuronal = factor(isNeuronal, c('NeuN+', 'NeuN-'))) %>% 
  inner_join(proportion_df)

pred_df2 = pred_df %>% 
  group_by(Enhancer, isNeuronal) %>%
  mutate(pred_mean = sum(score * prop), 
         pred_sem = sqrt(wtd.var(score, prop*n()))/sqrt(n())) %>%
  dplyr::select(-c(Celltype)) %>% 
  distinct(Enhancer, isNeuronal, .keep_all = T)


###################################
## read in the per clade OCR peaks
fn= list.files(path = here('data/raw_data/reporter_assay','rdas'),
               pattern = '.ZooMeta.rds', full.names = T)
names(fn) =  basename(fn) %>% ss('\\.', 2)

zooMeta_df =  lapply(fn, readRDS) %>% 
  data.table::rbindlist(idcol = 'Celltype') %>%
  inner_join(proportion_df) %>% group_by(group_meta, isNeuronal, name) %>% 
  summarise(score = sum(score * prop)) %>% 
  pivot_wider(names_from = group_meta, values_from = score, values_fill = NA) %>%
  filter(!is.na(name)) %>% 
  pivot_longer(cols = -c(isNeuronal, name), names_to ='group_meta', values_to = 'score') %>%
  mutate(group_meta = factor(group_meta, levels(df_meta$group_meta)),
         name = gsub('^hg38:|:250$','', name) %>% ss('-'),
         name = factor(name, unique(name))) %>%
  inner_join(df_meta %>% dplyr::select(group_meta, plotGroup, label)) %>% 
  filter(!grepl('Pholidota', as.character(group_meta)))
zooMeta_gr = zooMeta_df %>% pull(name)%>% GRanges()

################################
## Zoonomia constraint tracks ##
  chr = unique(seqnames(locus)) %>% as.character()
  oo = findOverlaps(subject = zooMeta_gr, query = locus)
  zooMeta_df2 = zooMeta_df[unique(subjectHits(oo)),]
  zooMeta_gr2 = zooMeta_gr[unique(subjectHits(oo))]
  
  # oo3 = findOverlaps(subject = zooMeta_gr2, query = CellTACIT_track[[celltype]]$dat)
  # zooMeta_df2$CellTACITscore = 0
  # zooMeta_df2$CellTACITscore[subjectHits(oo3)] = CellTACIT_track[[celltype]]$dat$score[queryHits(oo3)]
  # 
  ##########################
  ## 3) construct tracks  ##
  height_fig = 1.2 * 6; width_fig = 1.33 * 6; font_fig = 7
  plot_fn = here(PLOTDIR, 'plots', paste0('track_plot_snps_heatmap.pdf'))
  
  t1 <- ggtree(tree_dt, aes(color = Clade2, fill = Clade2)) + 
    # geom_tiplab(aes(image=uid), geom="phylopic", size = .05,
    #             offset=.7, align=TRUE, linetype = 'blank') +
    geom_tiplab(aes(label=plotGroup),  align=TRUE) + xlim(NA, 1.1) + 
    scale_color_manual(values = col_clade, name = 'Clade') 
  
  p1 <- ggplot(zooMeta_df2 %>% filter(isNeuronal =='NeuN+'), 
               aes(x = name, y=label)) + 
    geom_tile(aes(fill=score)) + 
    scale_fill_viridis_c(option = 'inferno', name ='Pred OCR Act.') + 
    guides(fill = guide_legend(title.position = "right")) +
    scale_x_discrete(position = "bottom")+ theme_tree2() +  
    theme(legend.title = element_text(angle = -90)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  p2 <- ggplot(zooMeta_df2 %>% filter(isNeuronal =='NeuN-'), 
               aes(x = name, y=label)) + 
    geom_tile(aes(fill=score)) + 
    scale_fill_viridis_c(option = 'inferno', name ='Pred OCR Act.') + 
    guides(fill = guide_legend(title.position = "right")) +
    scale_x_discrete(position = "bottom")+ theme_tree2() +  
    theme(legend.title = element_text(angle = -90)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  pg = p1 %>% insert_left(t1, width = .9) %>% 
    insert_right(p2, width = 1)

  pdf(plot_fn, height = height_fig, width = width_fig, onefile = F)
  print(pg) & theme(legend.position='none')
  dev.off()




