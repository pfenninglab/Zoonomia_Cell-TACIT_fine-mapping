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


##############################################
# 1a) read in the Zoonomia tree and group_meta list
rda_fn = here('data/tidy_data/Zoonomia_data', 
              'rdas','200_Mammals_Genome_Information.rda')
load(file = rda_fn)
df_meta = df_meta %>% mutate(Order2 = ifelse(Time.Since.Split.from.Human.TimeTree.median <=94, as.character(Order), as.character(Clade)))

col_clade = df %>% dplyr::select(Clade, col_clade)%>% filter(!duplicated(Clade)) %>% deframe()
col_meta = df_meta %>% mutate(value = Order2, name = col_meta) %>% 
  filter(!duplicated(value)) %>% dplyr::select(value,name) %>% deframe()

col_celltypes = rcartocolor::carto_pal(n = 8, 'Safe')
names(col_celltypes) = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb',  'Astro', 
                         'Microglia', 'OPC', 'Oligo')  

############################
# 1b) read in the GWAS traits
load(here('data/tidy_data/ldsc_gwas','rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% dplyr::select( -file) %>% 
  mutate(label = ss(as.character(trait), '_'))


##############################################
## 2) read in the polyfun, fine-mapped SNPs ##
DATADIR='data/raw_data'
finemap_df = here(DATADIR, 'polyfun_caudate/rdas/polyfun_caudate_finemapped_snps_20210518.rds') %>%
  readRDS() %>% mutate(start = POS_hg38, end = POS_hg38) %>%
  dplyr::select(-c(POS_hg19,POS_hg38, SNPVAR:P, BETA_MEAN:MAF, index:runPolyfun, h2:base))

gwas_group_legend = list(labels = names(group_col), col = 'gray80', fill = group_col)
group_col2 = group_col

finemap_df2 = finemap_df %>% filter(PIP > 0.01) %>% 
  arrange(group, subgroup, label) %>%
  mutate(stack.factor = as.character(group), color = sapply(group, function(x) list(c(group_col2[x], 'white')))) %>%
  group_by(group, SNP, CHR, start, end, A1, A2) %>%
  mutate(value1 = 100 * max(PIP, na.rm = TRUE), value2 = 100-value1,
         numTrait = length(unique(label)),  
         label = paste(label, collapse = ', '),
        numSubGroup = length(unique(subgroup)),  
        subgroup = paste(unique(subgroup), collapse = ', ')) %>% 
  ungroup() %>% distinct(group, SNP, CHR, start, end, A1, A2, PIP, .keep_all = TRUE)

## the SNPs
finemap_gr2 = finemap_df2 %>% 
  dplyr::select(-c (numTrait, numSubGroup, PIP)) %>% GRanges()
names(finemap_gr2)= with(mcols(finemap_gr2), paste0(SNP, ':', A1, ':', A2))
seqlevelsStyle(finemap_gr2) = 'UCSC'


#############################################
############# read the human OCRs ###########
hg38_peaks_fn = list.files(path = 'data/raw_data/hg38/Corces_2020/peak/', 
                      pattern = '.narrowPeak.gz', full.names = TRUE) %>%
  grep(pattern = 'Corces2020_caudate\\.', value = TRUE)
names(hg38_peaks_fn) = ss(basename(hg38_peaks_fn), '\\.', 2)
hg38_peaks_gr = hg38_peaks_fn[names(col_celltypes)] %>% 
  lapply(fread, col.names = c(bedNames, LETTERS[1:4])) %>% 
  rbindlist(idcol = 'celltype') %>% dplyr::select(-all_of(LETTERS[1:4])) %>%
  mutate(featureLayerID = celltype, border = 'black', fill = col_celltypes[celltype]) %>%
  GRanges()

oo = findOverlaps(subject = finemap_gr2[ finemap_df2$numTrait >=3], query = gr)
loci = finemap_gr2[ finemap_df2$numTrait >= 3][unique(subjectHits(oo))] %>% 
  resize( 2e4, fix="start") %>% resize( 4e4, fix="end") %>% reduce()

pdf('tmp.pdf', width =12,height = 5)
ind = subjectHits(findOverlaps(query = finemap_gr2, subject = loci))
lolliplot(finemap_gr2, hg38_peaks_gr, ranges = loci[20], type="pie.stack", 
          legend=gwas_group_legend, dashline.col="gray")
dev.off()