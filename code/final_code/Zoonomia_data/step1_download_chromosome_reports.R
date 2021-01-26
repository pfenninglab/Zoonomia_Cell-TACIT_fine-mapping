# to be run in the root github directory
setwd('code/final_code/Zoonomia_data')

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

suppressMessages(library(tidyverse))

PROJDIR='../../../data/tidy_data/Zoonomia_data'
OUTDIR=file.path(PROJDIR, 'tables')

######################################
# read in the species in 200m cactus #
tab_fn = file.path(PROJDIR, 'tables','200_Mammals_Genome_Information.tsv')
df = read_tsv(tab_fn)
names(df) = make.names(names(df))

df= df %>% filter(Genome.Download.Command_s. != 'Could not find genome') %>%
  filter(grepl('_genomic.fna.gz', Genome.Download.Command_s.))
df = df %>% mutate(Assembly.Report = 
                     paste0(basename(ss(gsub('_genomic.fna.gz', '_assembly_report.txt', 
                              Genome.Download.Command_s.), ' ',2)),'.gz'))

#################################
# download the assembly reports #
thecall = paste(gsub('_genomic.fna.gz', '_assembly_report.txt', 
                     df$`Genome.Download.Command_s.`), '-P', OUTDIR)
# parallel::mclapply(thecall, system, mc.cores = parallel::detectCores() - 3)

##########################
# gzip all the txt files #
chrtable_fn = list.files(path = file.path(PROJDIR, 'tables'), 
                         pattern = '_assembly_report.txt', full.names = T)
thecall = paste('gzip', chrtable_fn)
# parallel::mclapply(thecall, system, mc.cores = parallel::detectCores() - 3)

gziptable_fn = file.path(PROJDIR, 'tables', df$Assembly.Report)
all(file.exists(gziptable_fn))
names(gziptable_fn) = ss(basename(gziptable_fn),'_', 3)
gziptable_fn = gziptable_fn[order(names(gziptable_fn))]
library(data.table)

chrTableList = parallel::mclapply(gziptable_fn, function(file){
  tab = fread(file, skip = '# Sequence-Name')
  names(tab) = gsub('# ', '', names(tab)) %>% make.names()
  return(tab)
}, mc.cores = parallel::detectCores() - 3)


