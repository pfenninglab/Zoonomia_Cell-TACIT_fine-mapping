#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F, bitmapType='cairo')
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
library(rcartocolor)
library(ggrepel)
library(ggforce)
library(ggpubr)
library(ggsci)
library(here)

#########################################
# to be run in the root github directory
LABEL='ldsc_atac_and_conservation'
PLOTDIR=here('figures/exploratory/', 'ldsc_atac_and_conservation')
DATADIR=here('data/raw_data/',LABEL)

##########################
# read in the GWAS traits
load(here('data/tidy_data/ldsc_gwas','rdas','gwas_list_sumstats.rda'))
pheno = pheno %>% select( -file) %>% 
  mutate(label = ss(as.character(trait), '_'))

##############################################
# 1) read in the Zoonomia tree and group_meta list
rda_fn = here('data/tidy_data/Zoonomia_data', 
              'rdas','200_Mammals_Genome_Information.rda')
load(file = rda_fn)
df_meta = df_meta %>% mutate(Order2 = ifelse(Time.Since.Split.from.Human.TimeTree.median <=94, as.character(Order), as.character(Clade)))
col_meta = df_meta %>% mutate(value = Order2, name = col_meta) %>% 
  filter(!duplicated(value)) %>% dplyr::select(value,name) %>% deframe()

#################################
## save enrichments table to RDS
dir.create(here(PLOTDIR,'rdas'), showWarnings = F)
enrichments1 = readRDS(here(PLOTDIR,'rdas','zoo_meta_atac_and_conservation_prop_heritability_Corces2020.rds'))
enrichments2 = readRDS(here('figures/exploratory', 'ldsc_zoonomia_meta',
                            'rdas','zoonomia_meta_prop_heritability_Corces2020.rds'))
enrichments = bind_rows(enrichments1, enrichments2) %>% distinct()
enrichments$peaktype %>% table()
enrichments$celltype %>% table()

#################################
## make plots for presentation ##
dir.create(here(PLOTDIR,'plots'))
dir.create(here(PLOTDIR,'plots','prop_herit'))
height_ppt = 4; width_ppt = 8;
height_fig = 1.75; width_fig = 2.25; font_fig = 7
plot_traits = sort(unique(enrichments$trait))

###############
#### make plots
for(cell in unique(enrichments$celltype)){
  plot_fn = here(PLOTDIR,'plots',paste('zoonomia_caudate', 
                                       cell, 'prop_herit_enrichments.ppt.pdf', sep = '_'))
  pdf(width = width_ppt, height = height_ppt, file = plot_fn, onefile = T)
  for(lab in plot_traits){
    ## subset data to cell type and GWAS
    tmp = enrichments %>% filter(celltype %in% cell, trait %in% lab)
    if(nrow(tmp) > 0){
      ## the SNP enrichments ##
      pp1 = ggplot(data = tmp, aes(x = Time.Since.Split.from.Human.TimeTree.median, y = Enrichment)) +
        # # add phyloP enrichments lines
        # geom_hline(aes(yintercept = mean(tmp$phyloPcon.Enr), 
        #                linetype = "phyloP.Con"),color = "#882255") + 
        # geom_hline(aes(yintercept = mean(tmp$phyloPacc.Enr), 
        #                linetype = "phyloP.Acc"),color = "#888888")  +
        # add line fit
        geom_smooth(aes(color = peaktype), method = "lm", formula = y ~ x + I(x^2), size = 1, se = TRUE) +  
        geom_errorbar(aes(fill = Order2, shape = peaktype, ymin = Enrichment_min, ymax = Enrichment_max),
                      width = 0.1, position=position_dodge(width=2)) + 
        geom_jitter(aes(fill = Order2, shape = peaktype), position=position_dodge(width=2)) +
        scale_fill_manual(values = col_meta) + 
        scale_shape_manual(values = c(21:25)) + 
        scale_linetype_manual(name = "Zoonomia", values = c(1, 1), 
                              guide = guide_legend(title.position="top", nrow =3,override.aes = 
                                                     list(size = 1,color = c("#882255", "#888888")))) + 
        theme_bw(base_size = 10) + ggtitle(paste(cell, lab, sep = ', ')) +
        xlab(paste0('MY from Human (TimeTree median est.)')) + ylab(paste0('Heritability Enrichment')) + 
        guides(fill = guide_legend(nrow =3, title.position="top", override.aes = list(pch = 21 ,size = 3)),
               shape = guide_legend(nrow =3, title.position="top", override.aes = list(size = 3))) + 
        theme(legend.position = "bottom", legend.text=element_text(size=font_fig),
              legend.title=element_text(size=font_fig), legend.key.height=unit(.5,"line"))
      
      ## the proportion of SNP heritability ##
      pp2 = ggplot(data = tmp, aes(x = Time.Since.Split.from.Human.TimeTree.median, 
                                   y = Proportion_of_h2g)) +
        # # add phyloP enrichments lines
        # geom_hline(aes(yintercept = mean(tmp$phyloPcon.Prop), 
        #                linetype = "phyloP.Con"),color = "#882255") + 
        # geom_hline(aes(yintercept= mean(tmp$phyloPacc.Prop), 
        #                linetype = "phyloP.Acc"),color = "#888888")  +
        # add line fit
        geom_smooth(aes(color = peaktype), method = "lm", formula = y ~ x + I(x^2), size = 1, se = TRUE) +  
        geom_errorbar(aes(fill = Order2, shape = peaktype, ymin = Proportion_of_h2g_min, ymax = Proportion_of_h2g_max),
                      width = 0.1, position=position_dodge(width=2)) + 
        geom_jitter(aes(fill = Order2, shape = peaktype), position=position_dodge(width=2)) +
        scale_fill_manual(values = col_meta) + 
        scale_shape_manual(values = c(21:25)) + 
        scale_linetype_manual(name = "Zoonomia", values = c(1, 1), 
                              guide = guide_legend(title.position="top", nrow =3,override.aes = 
                                                     list(size = 1,color = c("#882255", "#888888")))) + 
        theme_bw(base_size = 10) + ggtitle(paste(cell, lab, sep = ', ')) +
        xlab(paste0('MY from Human (TimeTree median est.)')) + ylab(paste0('Proportion SNP Heritability')) + 
        guides(fill = guide_legend(nrow =3, title.position="top", override.aes = list(pch = 21 ,size = 3)),
               shape = guide_legend(nrow =3, title.position="top", override.aes = list(size = 3))) + 
        theme(legend.position = "bottom", legend.text=element_text(size=font_fig),
              legend.title=element_text(size=font_fig), legend.key.height=unit(.5,"line"))
      
      ## Plot the proportion of SNP heritability ##
      plot_fn1 = here(PLOTDIR,'plots', 'prop_herit',
                      paste('zoonomia_caudate',lab, cell, 'prop_herit_enrichments.ppt.pdf', sep = '_'))
      pdf(width = width_ppt, height = height_ppt, file = plot_fn1, onefile = F)
      pp = ggarrange(pp1, pp2, nrow = 1, common.legend = TRUE, 
                     align = 'h', legend="bottom")
      print(pp)
      dev.off()
      print(pp)
    }}
  dev.off()
}



