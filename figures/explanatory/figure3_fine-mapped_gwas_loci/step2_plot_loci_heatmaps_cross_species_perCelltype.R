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
library(treeio)
library(aplot)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(reshape)
library(scales)

DATADIR='data/raw_data'

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

sapply(indSNPList, function(x) x %>% filter(PIP > .1) %>% nrow()) %>% summary()
sapply(indSNPList, function(x) x %>% filter(PIP > .1) %>% nrow())


###############################################
## subset the Zoonomia tree to GroupMeta labels
tree2 = read.nhx(textConnection('((Loxodonta_africana:0.1204058581,Tolypeutes_matacus:0.1229192742):0.01145639324,(((Oryctolagus_cuniculus:0.1608919201,Rattus_norvegicus:0.2965409366):0.006031885449,(Tupaia_chinensis:0.1385286196,(Galeopterus_variegatus:0.08750990202,(Eulemur_flavifrons:0.07870638791,(Callithrix_jacchus:0.05610238745,(Macaca_mulatta:0.02760672357,(Nomascus_leucogenys:0.01503824207,(Pongo_abelii:0.01184875793,(Gorilla_gorilla:0.006121144092,(Pan_troglodytes:0.004607178674,Homo_sapiens:0.004410717579):0.00190302907):0.005963570605):0.002256803519):0.00775016763):0.01406957514):0.05984338905):0.01269656686):0.002845408602):0.002347302752):0.01705378614,(Uropsilus_gracilis:0.207400545,(Rousettus_aegyptiacus:0.1237216561,(Sus_scrofa:0.1368474811,(Equus_caballus:0.07947626741,(Canis_lupus_familiaris:0.1226195774,Manis_pentadactyla:0.1145549768):0.006550057057):0.002145230337):0.002284646667):0.005799996479):0.02130276027):0.01145639324);'))
tree2@data$node= tree2@data$value

tree_dt = full_join(tree2 %>% as_tibble(), df_meta) %>%
  filter(!is.na(node)) %>% as.treedata()
tree_dt =  groupOTU(tree_dt, clade_list, 'Clade2')


# 1b) read in the GWAS traits
load(here('data/tidy_data/ldsc_gwas','rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% dplyr::select( -file) %>% 
  mutate(label = ss(as.character(trait), '_'))

##################################
## plot for each cell type
celltype = 'MSN_D1'

for (celltype in names(col_celltypes)){
## read in the saved SNP and loci files
save_snps_fn= here('figures/explanatory/figure3_fine-mapped_gwas_loci/rdas',
              paste0('fig3_snpsInd_for_plots.', celltype,'.rds'))
finemap_df = readRDS(file = save_snps_fn)
finemap_gr = finemap_df %>% 
GRanges()
seqlevelsStyle(finemap_gr) = 'UCSC'

save_loci_fn= here('figures/explanatory/figure3_fine-mapped_gwas_loci/rdas',
             paste0('fig3_loci_for_plots.', celltype,'.rds'))
loci = readRDS(file = save_loci_fn)

###################################
## read in the per clade OCR peaks
peaknames = c('seqnames','start', 'end', 'name', 'score' , 'strand')
fn= list.files(path = 'figures/explanatory/figure3_fine-mapped_gwas_loci/bed',
         pattern = celltype, full.names = T) %>%
grep(pattern = 'allPeaks', value = TRUE) %>%
grep(pattern = 'mappable', value = TRUE) 
names(fn) = basename(fn) %>% ss('allPeaks\\.|\\.mappable', 2)

zooMeta_df = lapply(fn, fread, col.names = peaknames, header = FALSE) %>% 
data.table::rbindlist(idcol = 'group_meta') %>% 
arrange(seqnames, start) %>%
pivot_wider(all_of(peaknames[-5]), names_from = group_meta, values_from = score,
         values_fill = NA) %>%
pivot_longer(cols = -all_of(peaknames[-5]), names_to ='group_meta', values_to = 'score') %>%
mutate(group_meta = factor(group_meta, levels(df_meta$group_meta)),
name = gsub('^hg38:|:250$','', name) %>% ss('-'),
name = factor(name, unique(name))) %>%
inner_join(df_meta %>% dplyr::select(group_meta, plotGroup, label))
zooMeta_gr = zooMeta_df %>% GRanges()

################################
## Zoonomia constraint tracks ##
idx = 1
for(idx in seq_along(loci)) {
print(paste('Plotting locus' ,idx,'of', length(loci), 'for', celltype,'.'))
locus = loci[idx]
chr = unique(seqnames(locus)) %>% as.character()
oo = findOverlaps(subject = zooMeta_gr, query = locus)
zooMeta_df2 = zooMeta_df[unique(subjectHits(oo)),]
zooMeta_gr2 = zooMeta_gr[unique(subjectHits(oo))]

zooMeta_df2$numSNP = 0
oo2 = findOverlaps(subject = zooMeta_gr2, query = finemap_gr)
tmp = split(subjectHits(oo2), subjectHits(oo2))
zooMeta_df2$numSNP[as.numeric(names(tmp))] = lengths(tmp)

oo3 = findOverlaps(subject = zooMeta_gr2, query = CellTACIT_track[[celltype]]$dat)
zooMeta_df2$CellTACITscore = 0
zooMeta_df2$CellTACITscore[subjectHits(oo3)] = CellTACIT_track[[celltype]]$dat$score[queryHits(oo3)]

##########################
## 3) construct tracks  ##
height_fig = .9 * 6; width_fig = 2.25 * 6; font_fig = 7
regionName = paste0(chr,'_', signif(start(locus)/1e6,4), '-',signif(end(locus)/1e6,4), 'Mb')
geneNames = paste(sort(names(genes))[1:min(length(genes), 5)], collapse = ',')
plot_fn = here('figures/explanatory/figure3_fine-mapped_gwas_loci/plots',
             paste0('fig3_trackPlots.',celltype, '.',regionName,'.heatMap.pdf'))

t1 <- ggtree(tree_dt, aes(color = Clade2, fill = Clade2)) + 
geom_tiplab(aes(image=uid), geom="phylopic", size = .05,
            offset=.7, align=TRUE, linetype = 'blank') +
geom_tiplab(aes(label=plotGroup),  align=TRUE) + xlim(NA, 1.1) + 
scale_color_manual(values = col_clade, name = 'Clade') 

p1 <- ggplot(zooMeta_df2, aes(x = name, y=label)) + 
geom_tile(aes(fill=score)) + 
scale_fill_viridis_c(option = 'inferno', name ='Pred OCR Act.') + 
guides(fill = guide_legend(title.position = "right")) +
scale_x_discrete(position = "bottom")+ theme_tree2() +  
theme(legend.title = element_text(angle = -90)) +
theme(axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())

p2 <- ggplot(zooMeta_df2 %>% filter(!duplicated(name)), 
             aes(x = name, y='1', fill= CellTACITscore)) + 
geom_tile()  +  
geom_text(aes(label=round(CellTACITscore),
              color = -CellTACITscore), size = 3* font_fig / sqrt(nrow(zooMeta_df2) / 20), 
          check_overlap = TRUE)+ 
scale_fill_viridis_c(option = 'inferno', guide = 'none') + 
scale_color_viridis_c(option = 'mako', guide = 'none') + 
scale_x_discrete(position = "bottom")+ theme_tree2() +  
theme(legend.title = element_text(angle = -90)) +
theme(axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())

p3 <- ggplot(zooMeta_df2 %>% filter(!duplicated(name)), 
             aes(x = name, y='1', fill= numSNP)) + 
geom_tile()  +  geom_text(aes(label=numSNP, color = -numSNP), 
                          size = 3* font_fig / sqrt(nrow(zooMeta_df2) / 20), check_overlap = TRUE)+ 
scale_fill_viridis_c(option = 'plasma', guide = 'none') + 
scale_color_viridis_c(option = 'mako', guide = 'none') +
scale_x_discrete(position = "bottom")+ theme_tree2() +  
theme(legend.position = "none") +
theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

pg = p1 %>% insert_left(t1, width = .3)  %>% 
insert_top(p2, height = .1) %>% 
insert_bottom(p3, height = .1)
pdf(plot_fn, height = height_fig, width = width_fig, onefile = F)
print(pg) & theme(legend.position='bottom')
dev.off()
}}



