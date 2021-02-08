ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(ArchR))

setwd('code/raw_code/cross_species_peak_orthologs')

#####################################
### load Corces2018 ArchR project ###
PROJDIR='../../../data/raw_data/hg38/Corces_2020'
ARCHDIR=file.path(PROJDIR,'ArchR_Corces2020_caudate_labeled')
addArchRGenome(genome = 'hg38')
proj = loadArchRProject(ARCHDIR)

# add motif enrichment matrix, add motif deviations matrix
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "JASPAR2018", 
                            species = 'Homo sapiens', name = "Motif_JASPAR2018")
proj <- addDeviationsMatrix(proj,  peakAnnotation = "Motif_JASPAR2018", force = TRUE)
proj = saveArchRProject(ArchRProj = proj)


#########################################
### load Mouse BICCN CP ArchR project ###
PROJDIR2='../../../data/raw_data/mm10/BICCN_mouse_caudoputamen'
ARCHDIR2=file.path(PROJDIR2,'ArchR_BICCN_CP_labeled')
addArchRGenome(genome = 'mm10')

# add motif enrichment matrix, add motif deviations matrix
proj2 = loadArchRProject(ARCHDIR2)
proj2 <- addMotifAnnotations(ArchRProj = proj2, motifSet = "JASPAR2018", 
                            species = 'Mus musculus', name = "Motif_JASPAR2018")
proj2 <- addDeviationsMatrix(proj2,  peakAnnotation = "Motif_JASPAR2018", force = TRUE)
proj2 = saveArchRProject(ArchRProj = proj2)


#########################################
### load Mouse BICCN CP ArchR project ###
PROJDIR3='../../../data/raw_data/rheMac10/Stauffer_caudate'
ARCHDIR3=file.path(PROJDIR3,'ArchR_Stauffer_caudate_labeled')
GENOMEDIR='/home/bnphan/resources/genomes/rheMac10'
load(file.path(GENOMEDIR,'rheMac10_liftoff_GRCh38.p13_ArchR_annotations.rda'))

# add motif enrichment matrix, add motif deviations matrix
proj3 = loadArchRProject(ARCHDIR3)
proj3 <- addMotifAnnotations(ArchRProj = proj3, motifSet = "JASPAR2018", 
                             species = 'Homo sapiens', name = "Motif_JASPAR2018")
proj3 <- addDeviationsMatrix(proj3,  peakAnnotation = "Motif_JASPAR2018", force = TRUE)
proj3 = saveArchRProject(ArchRProj = proj3)




