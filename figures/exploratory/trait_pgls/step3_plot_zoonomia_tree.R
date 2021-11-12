#######################################
### set up libraries and functions ####
# conda activate r4
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F, bitmapType='cairo')

library(rcartocolor)
library(ggrepel)
library(ggforce)
library(here)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(ggtree)
library(ggtreeExtra)
library(tidyverse)
library(tidytree)
library(aplot)

LABEL='trait_pgls'
CODEDIR='code/raw_code/trait_pgls'
PLOTDIR='figures/exploratory/trait_pgls'

##########################
# read in the GWAS traits
load(here('data/tidy_data/ldsc_gwas','rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% select( -file) %>% 
  mutate(label = ss(as.character(trait), '_'))

##########################
# read in the Zoonomia tree and species list
rda_fn = here('data/tidy_data/Zoonomia_data', 
              'rdas','200_Mammals_Genome_Information.rda')
load(file = rda_fn)
df = df %>% arrange(Clade) %>%
  mutate(Clade = ifelse(grepl('Xen|Afro', as.character(Clade)), 
                        'Xenarthra & Afrotheria', as.character(Clade)),
         Clade = factor(Clade, unique(Clade)))
col_clade = df %>% dplyr::select(c(Clade, col_clade))%>% filter(!duplicated(Clade)) %>% deframe()
col_clade = col_clade[!grepl('Xenarthra', names( col_clade))]
col_order = df %>% dplyr::select(c(Order, col_order)) %>% filter(!duplicated(Order)) %>% deframe()

celltypes = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb',  'Astro', 
              'Microglia', 'OPC', 'Oligo')
col_celltypes = carto_pal(8,'Safe')
names(col_celltypes) = celltypes


######################################################
# read in the LDSC partitioned heritability estimation
DATADIR=here('data/tidy_data/Zoonomia_data')
numPeak_fn =here(DATADIR,'tables','zoonomia_peaks.txt')
input = fread(numPeak_fn, col.names = c('file', 'numPeak'), header = F)

# get number of human peaks
humanPeaks_fn = here('data/raw_data/hg38/Corces_2020/peak') %>% 
  list.files(path = ., pattern = '.narrowPeak.gz', full.names = T) %>%
  grep(pattern = 'Corces2020_caudate\\.',value = TRUE)
humanPeaks_fn = humanPeaks_fn[! grepl('Consensus', humanPeaks_fn)]
names(humanPeaks_fn) = ss(basename(humanPeaks_fn), 'narrowPeak.gz')
input2 = lapply(humanPeaks_fn, fread, header = F) %>% 
  rbindlist(fill = T, idcol='file') %>% group_by(file) %>%
  summarise(numPeak = n())
input3 = bind_rows(input, input2)

#########################################
## format groupings and calculate conditional cell type enrichment p-value
pd = input3 %>%  mutate(
    Species = ss(file, '\\.', 4),
    Species = ifelse(is.na(Species), 'Homo_sapiens', Species), 
    celltype = ss(file, '\\.', 2) %>%
      factor(rev(c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb',  'Astro', 
               'Microglia', 'OPC', 'Oligo'))), 
    cell_group = case_when(
      grepl('MSN|INT', celltype) ~ 'Neuron', 
      TRUE ~ 'Glia'
    )) %>% inner_join(x = df, by = 'Species') %>%
  filter(!is.na(celltype))

pd = pd %>% group_by(celltype) %>%
  mutate(tmp = numPeak[which(Species =='Homo_sapiens')],
         propMappable = numPeak/max(numPeak), 
         label = Species) %>% select(-tmp) 

#########################################
## combine tree data and peaks mappable
clade_list = with(df, split(Species, Clade))

tree2 <- ape::root(tree, outgroup = "Homo_sapiens", edgelabel = TRUE)
tree_dt = full_join(tree2 %>% as_tibble() %>% mutate(Species =label), df) %>%
  filter(!is.na(node)) %>% as.treedata()
tree_dt =  groupOTU(tree_dt, clade_list, 'Clade2')

t1 <- ggtree(tree_dt, aes(color = Clade2, fill = Clade2)) + 
  geom_tiplab(aes(label=plotGroup),  align=TRUE) + xlim(NA, 1.1) + 
  xlim(NA, 1.1) + scale_color_manual(values = col_clade, guide = 'none')
  
p1 <- ggplot(pd, aes(x = Species, y=numPeak)) + 
  geom_bar(stat = 'identity', aes(fill = celltype)) +  coord_flip() +
  scale_x_discrete(position = "bottom")+ theme_tree2() +  
  scale_y_continuous(trans = "reverse")+
  scale_fill_manual(values = col_celltypes, guide = 'none') 


pgls_pdf = here(PLOTDIR, 'plots', 'fig1b_peaks_mappable_zoonomia.pdf')
pdf(pgls_pdf, width = .75*5, height = 2.25*5)
pg = p1 %>% insert_left(t1, width = 1)
print(pg) & theme(legend.position='bottom')
dev.off()


