suppressMessages(library(ArchR))
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
setwd('code/raw_code/hg38_Corces_2020_caudate')

##################################
### set Arrow File parameters ####
addArchRGenome("hg38")
PROJDIR='../../../data/raw_data/hg38/Corces_2020'
ARCHDIR=file.path(PROJDIR,'ArchR_Corces2020_caudate')

### set Arrow File parameters ####
proj = loadArchRProject(ARCHDIR, showLogo = F)
barcodes = getCellNames(proj)
head(barcodes)

#############################################
### read in excel table w/ cell clusters ####
library(tidyverse)
library(readxl)
table_fn = file.path(PROJDIR, 'tables', 
                     'SupplementaryDataSet2_scATAC-QC-And-Cluster-Residence_v1.xlsx')
tab = readxl::read_excel(table_fn, sheet = 'scATAC-seq QC Metadata', skip = 21)
tab = tab %>% mutate(
  cellNames = paste0(Donor_ID,'.',Region,'#', `10x_SingleCell_Barcode`,'-1'))

# get main cell types in caudate only, filter by clusters > 40 cells
tab2 = tab %>% filter(Region=='CAUD' & cellNames %in% barcodes) %>% 
  filter(! FinalClusters %in% c('Doublets', 'Excluded', 'Nonneuronal')) %>% 
  group_by(FinalClusters) %>% filter(n() > 40) %>% ungroup()
tab2 %>% group_by(FinalClusters) %>% summarise(count = n())

#####################################
# subset to cells from final analyses
ARCHDIR2=file.path(PROJDIR,'ArchR_Corces2020_caudate_labeled2')
indKeep = match(tab2$cellNames, getCellNames(proj))
indKeep = indKeep[!is.na(indKeep)]
proj2 = subsetArchRProject( ArchRProj = proj, cells = getCellNames(proj)[indKeep],
                           outputDirectory = ARCHDIR2, force = TRUE)

#####################################
# add cell metadata to ArchR project
tab2 = tab2[match(getCellNames(proj2), tab2$cellNames),]
tab2 = as.data.frame(tab2)
iter = names(tab2)[! names(tab2) %in% c(names(getCellColData(proj2)),'cellNames')]
iter = iter[-c(1:9)]
for (name in iter){
proj2 = addCellColData(
  proj2, data = tab2[,name], name = name, 
  cells = proj2$cellNames, force = TRUE)
}
proj2$FinalClusters = make.names(proj2$FinalClusters)
proj2 = saveArchRProject(proj2)
