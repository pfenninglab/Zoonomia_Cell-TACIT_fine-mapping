# to be run in the root github directory
setwd('code/final_code/Zoonomia_data')
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

suppressMessages(library(tidyverse))
suppressMessages(library(ggtree))
suppressMessages(library(ggtreeExtra))
library(wesanderson)
library(rcartocolor)

PROJDIR='../../../data/tidy_data/Zoonomia_data'
OUTDIR=file.path(PROJDIR, 'tables')

######################################
# read in the species in 200m cactus #
tab_fn = file.path(PROJDIR, 'tables','200_Mammals_Genome_Information.tsv')
df = read_tsv(tab_fn, col_types = cols()) %>% rename_with(make.names) %>%
  rename_with(~ gsub(pattern = "([[:punct:]])\\1+",replacement = "\\1",.)) %>%
  rename_with(~ gsub(pattern = "\\.$",replacement = "",.))

df = df %>% arrange(Time.Since.Split.from.Human.TimeTree.median) %>%
  mutate(Order = factor(Order, unique(Order)), 
         Species = factor(Species, unique(Species)), 
         Clade = factor(Clade, unique(Clade)), 
         Family = factor(Family, unique(Family)))

col_order = c(carto_pal(12, "Vivid"),  carto_pal(length(levels(df$Order)) - 12, "Antique"))
names(col_order) = levels(df$Order)

col_clade = c(carto_pal(length(levels(df$Clade)) , "Bold"))
names(col_clade) = levels(df$Clade)

df = df %>% mutate(col_order = col_order[Order], col_clade = col_clade[Clade])

############################
# read in the species tree #
tree_fn = file.path(PROJDIR, 'tree','Zoonomia_ChrX_lessGC40_241species_30Consensus.tree')
tree = read.tree(tree_fn)

######################################
# read in the species in 200m cactus #
rda_fn = file.path(PROJDIR, 'rdas','200_Mammals_Genome_Information.rda')
save(df, tree, file = rda_fn)
