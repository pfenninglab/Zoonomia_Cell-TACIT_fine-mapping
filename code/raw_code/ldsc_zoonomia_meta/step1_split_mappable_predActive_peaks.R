#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F, bitmapType='cairo')
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(here))
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))

## hi this is BaDoi
DATADIR = here('data/raw_data/ldsc_zoonomia_meta'); 
dir.create(DATADIR,showWarnings = F); dir.create(here(DATADIR, 'peaks'),showWarnings = F)
PREDDIR=here('data/raw_data/ldsc_caudate_zoonomia/rdas')
PROJDIR=here('data/raw_data/ldsc_caudate_zoonomia')

##############################################
# 1) read in the Zoonomia tree and species list
rda_fn = here('data/tidy_data/Zoonomia_data', 
              'rdas','200_Mammals_Genome_Information.rda')
load(file = rda_fn)
col_meta = df_meta %>% select(col_meta)%>% filter(!duplicated(col_meta)) %>% deframe()

#####################################################################
# 2) read in the caudate cell types w/ prediction matrics of allPeaks
celltypes = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb', 'Astro', 'Microglia', 'OPC', 'Oligo')

cell = 'MSN_D1'
for (cell in celltypes){
  print(paste('Working on:', cell))
  for (type in c('allPeaks', 'enhPeaks')){
    print(paste('Getting',type,'for:', cell))
    
    ## read in the peak x species mtx -> peak x group long
    fn = here(PREDDIR,paste('Corces2020',cell,type, 'avgCNN.predictions.rds', sep = '.'))
    df_allMeta = readRDS(fn) %>% 
      filter(grepl('hg38', name)) %>% 
      pivot_longer(cols = !name, names_to = 'Species', values_to = 'score') %>%
      right_join(df %>% select(c(Species, group_meta)), by = 'Species') %>% 
      filter(!is.na(score)) %>% group_by(group_meta, name) %>%
      summarise(score = mean(score, na.rm = T)) %>%
      ungroup()
    
    ## get the peaks mappable to any of the groups & make bed file
    print(paste('Extracting mappable and predActive',type,'for:', cell))
    gr_mappable = df_allMeta %>%
      mutate(seqnames = name %>% ss(':', 2), 
             start = name %>% ss(':', 3) %>% ss('-', 1), 
             end = name %>% ss(':', 3) %>% ss('-', 2)) %>%
      filter(!grepl('-', seqnames)) %>%
      data.frame() %>% split(.$group_meta) %>% 
      lapply(GRanges)
    out_mappable_fn =  here(DATADIR, 'peaks', paste('Corces2020',cell,type, names(gr_mappable), 
                                                    'mappable.bed.gz', sep = '.'))
    tmp = map2(.x = gr_mappable, .y = out_mappable_fn, ~export(.x, .y))
    
    ## get the peaks pred active on average across group and make bed file
    out_predActive_fn = here(DATADIR, 'peaks', paste('Corces2020',cell,type, names(gr_predActive), 
                                            'predActive.bed.gz', sep = '.'))
    gr_predActive = gr_mappable %>% lapply(function(gr) gr[ gr$score > 0.5])
    tmp = map2(.x = gr_predActive, .y = out_predActive_fn, ~export(.x, .y))
    }
}



########################################################################
# CellTACIT score: mcompute calibrated CNN score then compute cross-species 
MODEL_TYPE = 'hgRmMm_nonCelltypeNonEnhBiasAway10x'
DATADIR=here('data/raw_data/cnn_enhancer_ortholog')
calib_out_fn = here(DATADIR,paste('rdas/Caudate',MODEL_TYPE,'pos_calibration_ecdf.rds', sep = '_'))
model_calibration = readRDS(file = calib_out_fn)

cell = 'MSN_D1'
for( cell in celltypes){
  outCellTACIT_rds = here(PROJDIR,'rdas',paste0('Corces2020.',cell, '.CellTACIT.mean.rds'))
  outCellTACIT_bed = here(PROJDIR,'CellTACIT',paste0('Corces2020.',cell, '.CellTACIT.mean.bed.gz'))
  if(file.exists(outCellTACIT_bed)){
 
  ## read in the peak x species mtx -> peak x group long
  fn = here(PREDDIR,paste('Corces2020',cell, 'allPeaks.avgCNN.predictions.rds', sep = '.'))
  df_allMeta = readRDS(fn) %>% 
    filter(grepl('hg38', name)) %>% 
    # for peaks w/ predicted open, get the calibrated positive score
    mutate_if(is.numeric, ~ model_calibration[[cell]](.)) %>% 
    pivot_longer(cols = !name, names_to = 'Species', values_to = 'score') %>%
    right_join(df %>% select(c(Species, group_meta)), by = 'Species') %>% 
    group_by(group_meta, name) %>% summarise(score = mean(score, na.rm = T)) %>%
    mutate(score = ifelse(is.na(score),0, score)) %>% 
    ungroup()

  df_allMeta2 = df_allMeta %>%
    mutate(MYA = group_meta %>% as.character() %>% ss('#', 2), 
           MYA = as.numeric(MYA) + 1) %>%
    group_by(name) %>%
    ## weighted geometric mean, MYA weighted by mean clade-wise enhancer activity score
    # summarise(num = sum( score * log(MYA), na.rm = T), 
    #           den = sum( score, na.rm = T ), 
    #           score = exp(num / den)) %>%
    # dplyr::select(-c(num, den)) %>%
    ## weighted ari
    summarise(score = sum(score * MYA) / n()) %>%
    mutate( score = ifelse(is.na(score), 1, score))
  
  gr_allMeta2 = df_allMeta2 %>%
    mutate(seqnames = ss(name, ":", 2), 
           start = ss(name, ":", 3) %>% ss('-', 1), 
           end = ss(name, ":", 3) %>% ss('-', 2)) %>%
    column_to_rownames(var ='name') %>% GRanges()
  
  gr_allMeta2 %>% saveRDS(outCellTACIT_rds)
  export(gr_allMeta2, outCellTACIT_bed)
}
}


##############################################################################
# 3) read in the mouse cortical interneuron cell types w/ prediction matrics of allPeaks
celltypes = c('PV', 'SST', "VIP")
PREDDIR2='/projects/pfenninggroup/mouseCxStr/NeuronSubtypeATAC/Zoonomia_CNN/predictions/'
HALDIR='/projects/pfenninggroup/mouseCxStr/NeuronSubtypeATAC/Zoonomia_CNN'
cell = 'VIP'; type = 'enhPeaks'

if(FALSE){
  for (cell in celltypes){
    for (type in c('enhPeaks')){
      print(paste('Getting',type,'for:', cell))
      ## read in the peak x species mtx -> peak x group long
      pred_fn = here(PREDDIR2,paste0('mouse', cell,'peaks_AvgActivityPredictions_boreoeutheria.txt'))
      orth_fn = here(HALDIR, paste0('mouse_',cell),paste0('mouseReproducible', cell,'_hg38_HALPERout.bed'))
      df_orth = orth_fn %>% 
        fread(header = FALSE, col.names = c('chr','start','end', 'peak','V1', LETTERS[1:4])) %>%
        select(c(chr:end, V1)) %>% filter(!grepl('Un|alt|\\.', chr)) %>%
        mutate(name = paste0('hg38:',chr,':', start,'-', end)) %>%
        select(c(name, V1))
      
      df_allMeta = fread(pred_fn) %>% 
        pivot_longer(cols = !V1, names_to = 'Species', values_to = 'score') %>%
        inner_join(df_orth, by = 'V1') %>% 
        filter(grepl('hg38', name)) %>% 
        right_join(df %>% select(c(Species, group_meta)), by = 'Species') %>% 
        filter(!is.na(score)) %>% group_by(group_meta, name) %>%
        summarise(score = mean(score, na.rm = T)) %>%
        ungroup()
      
      ## get the peaks mappable to any of the groups & make bed file
      print(paste('Extracting mappable and predActive',type,'for:', cell))
      gr_mappable = df_allMeta %>%
        mutate(group_meta = droplevels(group_meta)) %>%
        mutate(seqnames = name %>% ss(':', 2), 
               start = name %>% ss(':', 3) %>% ss('-', 1), 
               end = name %>% ss(':', 3) %>% ss('-', 2)) %>%
        filter(!grepl('-', seqnames)) %>% 
        data.frame() %>% split(.$group_meta) %>% 
        lapply(GRanges)
      out_mappable_fn =  here(DATADIR, 'peaks', paste('BICCN_mouse_CATlas',cell,type, names(gr_mappable), 
                                                      'mappable.bed.gz', sep = '.'))
      tmp = map2(.x = gr_mappable, .y = out_mappable_fn, ~export(.x, .y))
      
      ## get the peaks pred active on average across group and make bed file
      gr_predActive = gr_mappable %>% lapply(function(gr) gr[ gr$score > 0.5])
      out_predActive_fn = here(DATADIR, 'peaks', paste('BICCN_mouse_CATlas',cell,type, names(gr_predActive), 
                                                       'predActive.bed.gz', sep = '.'))
      tmp = map2(.x = gr_predActive, .y = out_predActive_fn, ~export(.x, .y))
    }
  }
}


##############################################################################
# 3a) read in the mouse cortical interneuron cell types w/ prediction matrics of allPeaks
celltypes = c('PV', 'SST', "VIP")
PREDDIR2='/projects/pfenninggroup/mouseCxStr/NeuronSubtypeATAC/Zoonomia_CNN/predictions/'
HALDIR='/projects/pfenninggroup/mouseCxStr/NeuronSubtypeATAC/Zoonomia_CNN'
cell = 'VIP'; type = 'enhPeaks'

if(TRUE){
  for (cell in celltypes){
    for (type in c('enhPeaks')){
      print(paste('Getting',type,'for:', cell))
      ## read in the peak x species mtx -> peak x group long
      pred_fn = here(PREDDIR2,paste0('human', cell,'peaks_AvgActivityPredictions_boreoeutheria.txt'))

      df_allMeta = fread(pred_fn) %>% 
        rename('V1' = 'name') %>%
        pivot_longer(cols = !name, names_to = 'Species', values_to = 'score') %>%
        filter(grepl('hg38', name)) %>% 
        mutate(name = gsub('hg38_','hg38.', name)) %>%
        right_join(df %>% select(c(Species, group_meta)), by = 'Species') %>% 
        filter(!is.na(score)) %>% group_by(group_meta, name) %>%
        summarise(score = mean(score, na.rm = T)) %>%
        ungroup()
      
      ## get the peaks mappable to any of the groups & make bed file
      print(paste('Extracting mappable and predActive',type,'for:', cell))
      gr_mappable = df_allMeta %>%
        mutate(group_meta = droplevels(group_meta)) %>%
        mutate(seqnames = name %>% ss('\\.', 2), 
               start = name %>% ss('\\.', 3), end = name %>% ss('\\.', 4)) %>%
        filter(!grepl('-', seqnames)) %>% 
        data.frame() %>% split(.$group_meta) %>% 
        lapply(GRanges)
      out_mappable_fn =  here(DATADIR, 'peaks', paste('BICCN_human_M1_SNAREseq_v1',cell,type, names(gr_mappable), 
                                                      'mappable.bed.gz', sep = '.'))
      tmp = map2(.x = gr_mappable, .y = out_mappable_fn, ~export(.x, .y))
      
      ## get the peaks pred active on average across group and make bed file
      gr_predActive = gr_mappable %>% lapply(function(gr) gr[ gr$score > 0.5])
      lengths(gr_predActive)/ lengths(gr_mappable) * 100
      out_predActive_fn = here(DATADIR, 'peaks', paste('BICCN_human_M1_SNAREseq_v1',cell,type, names(gr_predActive), 
                                                       'predActive.bed.gz', sep = '.'))
      tmp = map2(.x = gr_predActive, .y = out_predActive_fn, ~export(.x, .y))
    }
  }
}


#########################################
# 4) read in the bulk M1 and liver peaks
celltypes = c('M1ctx', "Liver")
pred_fns = c('M1ctx' = '/home/dschaffe/enhancer-clustering/data/cortex/MoMaRaBaEnhancersNRSummit.txt',
           'Liver' = '/home/dschaffe/enhancer-clustering/data/liver/MoRaMaCoPiEnhancersNRSummit.txt')
orth_fns = c('M1ctx' = '/projects/MPRA/Irene/PGLSOutputs/HumanOfRenamedAllNRSummit_cortexModified.bed',
           'Liver' = '/projects/MPRA/Irene/PGLSOutputs/HumanOfRenamedAllNRSummit_liverModified.bed')
colnames = c('peak', '/home/dschaffe/enhancer-clustering/data/BoreoeutheriaSpeciesNamesNew.txt' %>% fread(sep = '\t',header = F) %>% pull(V1))
cell = 'M1ctx'; type = 'enhPeaks'

if(FALSE){
  for (cell in celltypes){
    for (type in c('enhPeaks')){
      print(paste('Getting',type,'for:', cell))
      ## read in the peak x species mtx -> peak x group long
      pred_fn = pred_fns[cell]; orth_fn = orth_fns[cell]
      df_orth = orth_fn %>% 
        fread(header = FALSE, col.names = c('chr','start','end','peak')) %>%
        select(c(chr:end, peak)) %>%
        mutate(name = paste0('hg38:',chr,':', start,'-', end)) %>%
        select(c(name, peak))
      
      df_allMeta = fread(pred_fn, header = FALSE, col.names = colnames) %>% 
        pivot_longer(cols = !peak, names_to = 'Species', values_to = 'score') %>%
        inner_join(df_orth, by = 'peak') %>% 
        filter(grepl('hg38', name)) %>% 
        mutate(Species = gsub(' ', '_', Species), 
               score = ifelse( score == -1, NA, score)) %>%
        right_join(df %>% select(c(Species, group_meta)), by = 'Species') %>% 
        filter(!is.na(score)) %>% group_by(group_meta, name) %>%
        summarise(score = mean(score, na.rm = T)) %>%
        ungroup()
      
      ## get the peaks mappable to any of the groups & make bed file
      print(paste('Extracting mappable and predActive',type,'for:', cell))
      gr_mappable = df_allMeta %>%
        mutate(group_meta = droplevels(group_meta)) %>%
        mutate(seqnames = name %>% ss(':', 2), 
               start = name %>% ss(':', 3) %>% ss('-', 1), 
               end = name %>% ss(':', 3) %>% ss('-', 2)) %>%
        filter(!grepl('-', seqnames)) %>% 
        data.frame() %>% split(.$group_meta) %>% 
        lapply(GRanges)
      out_mappable_fn =  here(DATADIR, 'peaks', paste('BulkATACseq',cell,type, names(gr_mappable), 
                                                      'mappable.bed.gz', sep = '.'))
      tmp = map2(.x = gr_mappable, .y = out_mappable_fn, ~export(.x, .y))
      
      
      ## get the peaks pred active on average across group and make bed file
      out_predActive_fn = here(DATADIR, 'peaks', paste('BulkATACseq',cell,type, names(gr_predActive), 
                                                       'predActive.bed.gz', sep = '.'))
      gr_predActive = gr_mappable %>% lapply(function(gr) gr[ gr$score > 0.5])
      tmp = map2(.x = gr_predActive, .y = out_predActive_fn, ~export(.x, .y))
    }
  }
}
