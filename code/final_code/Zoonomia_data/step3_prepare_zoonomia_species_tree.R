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

#####################################################################
# split the species groups up by MYA from human and Order-level LDSC
df = df %>% mutate(
  # group the Xenarthra and Afrotheria Orders together by Clade
  tmp = ifelse(Time.Since.Split.from.Human.TimeTree.median <=94, as.character(Order), as.character(Clade)),
  col_meta = ifelse(Time.Since.Split.from.Human.TimeTree.median <=94, col_order, col_clade),
  group_meta = paste(tmp, Time.Since.Split.from.Human.TimeTree.median, sep = '#'), 
  group_meta = factor(group_meta, unique(group_meta))) %>% select(-tmp)
df %>% count(group_meta) %>% data.frame()

df_meta = df %>% group_by(group_meta) %>% 
  select(-c(Zoonomia.Index,Common.Name, Genome.Download.Command_s:Genome.File.Name,
         starts_with('Time.Since.Split.from.Mouse'))) %>%
  summarise_all(~ if(is.numeric(.)) {mean(., na.rm = TRUE)} 
                else paste(unique(.), collapse = ',')) %>%
  mutate_if(is.character,~ ifelse(.x %in% c('NA', ","), '', .x)) %>%
  mutate_if(is.character,~ gsub('^NA,', '', .x)) %>%
  right_join(df %>% count(group_meta, name = 'numSpecies'), by = 'group_meta')


############################
# read in the species tree #
tree_fn = file.path(PROJDIR, 'tree','Zoonomia_ChrX_lessGC40_241species_30Consensus.tree')
tree = read.tree(tree_fn)

######################################
# read in the species in 200m cactus #
rda_fn = file.path(PROJDIR, 'rdas','200_Mammals_Genome_Information.rda')
save(df,df_meta, tree, file = rda_fn)
