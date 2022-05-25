#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
library(tidyverse)
library(here)

# to be run in the root github directory
DATADIR='data/tidy_data/track_files'
bigbed_fn = list.files(here(DATADIR,'bigBed'))

df = tibble("Directory Name" = 'bigBed', "File Name" = bigbed_fn) %>%
  mutate("Species" = ss(`File Name`, '\\.', 3), 
         'Celltype' = ss(`File Name`, '\\.', 2), 
    "Shortened File Name" = paste(Celltype, Species, sep = '.'), 
    "Trackdb File Name" = paste0('trackdb/trackdb_ortholog_', Celltype,'.txt'), 
    "Predictions in Original Open Chromatin Regions?" = 
      case_when(Species == 'Homo_sapiens' ~ 'Yes', TRUE ~ 'No'),
    "Track Description" = paste0("Predicted activity for human caudate ", 
                                 Celltype," open chromatin region orthologs in ",Species))

zoo_df = here("data/tidy_data/Zoonomia_data/tables/200_Mammals_Genome_Information.tsv") %>%
  read_tsv() %>%
  dplyr::rename("Chromosome Naming Convention" = "Chromosome Naming Convention in Cactus Alignment") %>%
  dplyr::select(Species, `Chromosome Naming Convention`)

df = df %>% inner_join(zoo_df) %>%
  relocate(`Chromosome Naming Convention`, .after = `Trackdb File Name`)

dfList = split(df, df$Celltype) %>% lapply(function(df) df %>% dplyr::select(-Celltype))
writexl::write_xlsx(dfList, path = here(DATADIR,'trackdb', 'Prediction_Browser_Tracks_data_description.xlsx'))
