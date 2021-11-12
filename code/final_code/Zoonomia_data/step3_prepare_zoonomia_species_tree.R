# to be run in the root github directory
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

suppressMessages(library(tidyverse))
suppressMessages(library(ggtree))
suppressMessages(library(ggtreeExtra))
library(wesanderson)
library(here)
library(rcartocolor)
library(ggimage)

PROJDIR='data/tidy_data/Zoonomia_data'
OUTDIR=here(PROJDIR, 'tables')

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

df_meta = df %>% arrange(plotGroup != '') %>%
  group_by(group_meta) %>% 
  select(-c(Zoonomia.Index,Common.Name, Genome.Download.Command_s:Genome.File.Name,
         starts_with('Time.Since.Split.from.Mouse'))) %>%
  summarise_all(~ if(is.numeric(.)) {mean(., na.rm = TRUE)} 
                else paste(unique(.), collapse = ',')) %>%
  mutate_if(is.character,~ ifelse(.x %in% c('NA', ","), '', .x)) %>%
  mutate_if(is.character,~ gsub('^NA,', '', .x)) %>%
  mutate_if(is.character,~ gsub(',NA$', '', .x)) %>%
  right_join(df %>% count(group_meta, name = 'numSpecies'), by = 'group_meta')

#######################################################
# grab the phyloPics of as many species as possible
df_meta = df_meta %>%
  mutate(uid = sapply(Species, function(nameConcat){
    if(nameConcat == 'Homo_sapiens') {
      uid = '9fae30cd-fb59-4a81-a39c-e1826a35f612'
      return(uid)
    }
    for (name in unlist(strsplit(nameConcat, ','))){
      x <- gsub("[^a-zA-Z]+", "+", name)
      url <- paste0("http://phylopic.org/api/a/name/search?text=", 
                    x, "&options=scientific+json")
      results <- jsonlite::fromJSON(url)$result   
      if (length(results)==0){
        uid = NA
      } else{
        res <- results[[1]]
        for (id in res$uid) {
          uid <- ggimage:::phylopic_valid_id(id)
          if (!is.na(uid)) 
            return(uid)
            break
        } } } }
  ))

names(df_meta$uid) = NULL

############################
# read in the species tree #
tree_fn = file.path(PROJDIR, 'tree','Zoonomia_ChrX_lessGC40_241species_30Consensus.tree')
tree = read.tree(tree_fn)

############################
# read in the phenotypes to investigated
trait_fn = file.path(PROJDIR, 'tables','phenotypes_to_investigate.xlsx')
trait_df = readxl::read_excel(trait_fn, sheet = '2. phenotypes_list', na = c('', "NA", 'na')) %>%
  select(-contains(c("DELETE", 'MaTrics','Sex','Number', 'juven', 'Diet', 'beha',
                     'restraint', 'temperature', 'lab_condition', 'Hypselodonty', 'Mean_',
                     'Basal', 'Bmr','animals_sampled', 'Solitary', 'Social', 
                     'Group', 'Patern','ageabsolute', 
                     'Adaptation', 'Telemetry', '24HR', 'Light', 'EEG','...'))) %>% 
  rename_with(make.names) %>% rename_with(~gsub('...C', '', .x) %>% str_to_sentence()) %>% 
  mutate(
    Species.binomial = gsub(' ', '_', Species.binomial) %>% 
      stringr::str_to_sentence() %>% stringr::str_extract('^([^_]*_[^_]*)'),
    Species_synonyms = gsub(' ', '_', Species_synonyms) %>% stringr::str_to_sentence(),
    Species.binomial = ifelse(Species_synonyms %in% df$Species, Species_synonyms,
                              ifelse(Species.binomial %in% df$Species, Species.binomial, NA))) %>%
  filter(!is.na(Species.binomial), !duplicated(Species.binomial)) %>%
  rename('Species' = 'Species.binomial')

######################################
# read in the species in 200m cactus #
rda_fn = file.path(PROJDIR, 'rdas','200_Mammals_Genome_Information.rda')
save(df,df_meta, tree, trait_df, file = rda_fn)

zoo_meta_table = file.path(PROJDIR, 'tables','Table_S7_Zoo_meta_Genome_Information.xlsx')
writexl::write_xlsx(df_meta, path = zoo_meta_table)

zoo_meta_table2 = file.path(PROJDIR, 'tables','Table_S7_Zoonomia_Genome_Information.xlsx')
writexl::write_xlsx(df_meta, path = zoo_meta_table)

zoo_meta_table3 = file.path(PROJDIR, 'tables','Table_S7_brain_sleep_phenotypes.xlsx')
writexl::write_xlsx(df_meta, path = zoo_meta_table3)

trait_df %>% filter(grepl('Boreo', Taxonomic.lineage)) %>%
  summarise(sum(!is.na(Brain.resid)))

trait_df %>% filter(grepl('Euarc|Glires', Taxonomic.lineage)) %>%
  summarise(sum(!is.na(Brain.resid)))



