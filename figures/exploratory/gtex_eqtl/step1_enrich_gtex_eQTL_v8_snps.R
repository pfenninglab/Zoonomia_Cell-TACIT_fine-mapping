library(tidyverse)
library(ggplot2)
library(data.table)
library(here)
library(broom)
library(rtracklayer)
library(GenomicRanges)
library(speedglm)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

CODEDIR='code/final_code/gtex_eqtl'
DATADIR='data/tidy_data/gtex_eqtl'
PLOTDIR='figures/exploratory/gtex_eqtl'

dir.create(here(PLOTDIR, 'plots'), showWarnings = F, recursive = T)
dir.create(here(PLOTDIR, 'tables'), showWarnings = F, recursive = T)
dir.create(here(PLOTDIR, 'rdas'), showWarnings = F, recursive = T)

quartiles = c('0-25%', '25-50%', '50-75%', '75-100%')
names(quartiles) = paste0('quartile', 1:4); 

celltypes = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb',  'Astro', 
              'Microglia', 'OPC', 'Oligo')

##########################################
## 1) read in the Cell TACIT Age quartiles
annotated_eQTL_fn = here(DATADIR, 'rdas', 'GTEx_v8_Caudate_eQTL.celltacitAge.rds')

if(file.exists(annotated_eQTL_fn)){
  snp_dfList = readRDS(annotated_eQTL_fn) %>% split(.$celltype)
} else{
  peak_fn = list.files('data/raw_data/ldsc_celltacit_age_decile/peaks', 
                       pattern = '.bed.gz', full.names = T) %>% str_subset('quartile')
  peak_fn = peak_fn[!grepl('rand', peak_fn)]
  names(peak_fn) = basename(peak_fn) %>% gsub('Corces2020.|CellTACIT.|.bed.gz', '', .)
  peak_list = peak_fn %>% lapply(import)
  
  ## read in the GTEx and variant_id SNPs
  snp_gr = readRDS(here(DATADIR, 'rdas', 'Brain_Caudate_basal_ganglia_hg38.v8.gr.rds'))
  snp_df = as.data.frame(snp_gr)  %>% rownames_to_column('variant_id') %>% 
    as_tibble() %>% dplyr::rename('chr' = 'seqnames') %>% 
    mutate(pos = start) %>% dplyr::select(-c(start, end, strand, width))
  
  ## annotate each SNP with the overlapping peak
  snp_df2 = lapply(peak_list, function(gr){
    countOverlaps(gr, query = snp_gr) %>% 
      enframe() %>% filter(value > 0)
  }) %>% rbindlist(idcol = 'peak') %>% 
    mutate(celltype = ss(peak, '\\.', 1), quartile = ss(peak, '\\.', 2)) %>%
    pivot_wider(id_cols = c('name'), names_from = celltype, 
                values_from = quartile, values_fill = 'non-overlap') %>% 
    dplyr::rename('variant_id' = 'name') %>% 
    left_join(x = snp_df, y = .) %>% 
    pivot_longer(cols = all_of(celltypes), names_to = "celltype", values_to = "quantile") %>% 
    replace(is.na(.), 'non-overlap') %>% 
    mutate(quantile = factor(quantile),
           quantile = relevel(quantile, ref = 'non-overlap'),
           celltype = factor(celltype, celltypes))
  
  table(snp_df2$celltype)
  table(snp_df2$quantile)
  
  ## write out the annotations
  snp_df2 %>% saveRDS(annotated_eQTL_fn)  
  snp_dfList = snp_df2 %>% split(.$celltype)
  rm(snp_df2, snp_df, snp_gr)
}


###########################################
## 3) read in the GTEx fine-mapped variant_id SNPs
alpha = 0.05
save_enrichments_fn = here(DATADIR,'rdas',paste0('GTEx_v8_Caudate_eQTL.Caudate.glm.rds'))
if(!file.exists(save_enrichments_fn)){
  variant_id_df = here(DATADIR,'rdas', 'Brain_Caudate_basal_ganglia.v8.df.rds') %>% 
    readRDS() %>% mutate(signif_FDR = FDR < alpha) %>% 
    arrange(FDR) %>% # ascending 
    filter(!duplicated(variant_id))# keep SNP variant_id w/ the lowest FDR value
  
  table(variant_id_df$FDR < alpha)
  
  ## test enrichment
  snp_enrichment_df = lapply(snp_dfList, function(df){
      ## assemble SNP annotation w/ sn-variant_id
      ret = variant_id_df %>% inner_join(df) # add the SNP-CellTACIT age annotation
      ## calculate the enrichment w/ GLM
      test = speedglm(signif_FDR ~ quantile, data = ret,family=binomial())
      tidied = tidy(test)
      return(tidied)
    }) %>% rbindlist(idcol = 'CellTACIT_celltype') %>% 
    filter(term != '(Intercept)') %>%
    mutate(OR = exp(estimate), 
           OR_max = exp(estimate + std.error), 
           OR_min = exp(estimate - std.error), 
           term = gsub('quantile', '', term), 
           FDR = p.adjust(p.value, 'fdr')) %>% 
    arrange(p.value)
  
  summary(snp_enrichment_df$p.value)
  ## save enrichment
  snp_enrichment_df %>% saveRDS(save_enrichments_fn)
}







