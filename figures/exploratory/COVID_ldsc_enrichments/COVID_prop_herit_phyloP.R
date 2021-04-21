#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(repr.plot.width=11, repr.plot.height=8.5)
options(stringsAsFactors = F, bitmapType='cairo')
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
library(rcartocolor)
library(ggrepel)
library(ggforce)
library(ggpubr)
library(ggsci)
library(here)

# to be run in the root github directory
LABEL='COVID_ldsc_enrichments'
SETDIR='figures/exploratory/COVID_ldsc_enrichments'
PROJDIR=file.path('data/raw_data/',LABEL)

##########################################
# read in the SNP heritability estimates 
pheno = file.path('data/tidy_data/ldsc_gwas/gwas_list_sumstats.tsv') %>% 
  read_tsv(col_type = cols()) %>%
  filter(!grepl('XX', fig_group)) %>%
  mutate(trait = factor(trait, unique(trait)), 
         group = factor(group, unique(group)))
pheno = pheno %>% mutate(group_col = group_col[group]) %>% filter(group == 'Immune')

##########################################
# read in the SNP heritability estimates 
h2_fn = file.path('data/tidy_data/ldsc_gwas/heritability') %>% 
  list.files(path = ., pattern = 'log',full.names = T)
names(h2_fn) = ss(basename(h2_fn),'\\.')
h2 = h2_fn %>% lapply(fread, skip = 'two-step', col.names ='log', nrows = 5, sep = '\n') %>% 
  bind_rows(.id = 'match') %>%
  separate(log, sep = ':', into = c('var', 'val')) %>% 
  mutate(var = make.names(var), val = str_trim(val)) %>%
  pivot_wider(names_from = var, values_from = val) %>%
  separate(Total.Observed.scale.h2, sep = ' ', into = c('h2', 'h2_se')) %>%
  separate(Intercept, sep = ' ', into = c('Inter', 'Inter_se')) %>% 
  mutate_if(is.character, str_replace_all, pattern = "\\(", replacement = "") %>% 
  mutate_if(is.character, str_replace_all, pattern = "\\)", replacement = "") %>%
  select(match:Inter_se) %>% mutate(across(!match, as.numeric)) %>% 
  arrange(desc(h2)) %>% mutate(
    h2_Z = h2 / h2_se, h2_P = pnorm(-h2_Z), h2_logP = -log10(h2_P))


#####################################################
# read in how many SNPs used to estimate heritability 
M = h2_fn %>% lapply(fread, skip = 'After merging with reference panel LD', 
                     col.names ='log', nrows = 1, sep = '\n') %>% 
  bind_rows(.id = 'match') %>%
  separate(log, sep = ',', into = c('var', 'num_SNP')) %>% 
  mutate(var = make.names(var), num_SNP = ss(str_trim(num_SNP),' ')) %>%
  select(match,num_SNP) %>% mutate(across(!match, as.numeric)) %>%
  arrange(desc(num_SNP))

##############################################
# merge gwas phenotypes and SNP heritabilities
pheno = left_join(h2, M, by = 'match') %>% left_join(pheno, ., by = 'match') %>%  
  mutate(h2_perSNP = h2 / num_SNP) %>% filter(h2 > 0, h2 < 1) 

out = pheno$h2
names(out) = pheno$trait
print(out)

out2 = pheno$h2_se
names(out2) = pheno$trait
print(out2)

#########################################
# read in the LDSC partitioned heritability estimation
enrich_fn =here('data/raw_data/caudate_conservation_ldsc/prop_herit') %>% 
  list.files(path = ., pattern = '-HGI_2021.EUR.results.gz', full.names = T)
names(enrich_fn) = ss(basename(enrich_fn), '.results.gz')
enrich_fn = enrich_fn[grepl('200m_scoresPhyloP', names(enrich_fn))]
input = lapply(enrich_fn, read_tsv, n_max = 1, col_type = cols()) %>% 
  rbindlist(fill = T, idcol='file') %>% select(file:Coefficient_SE)
input %>% data.frame() %>% head()

#########################################
## format groupings and calculate conditional cell type enrichment p-value
enrichments = input %>% 
  mutate(
    match = ss(file, '\\.', 4), 
    Categories = ss(file,'\\.',3), 
    peaktype = case_when(
      grepl('accl', Categories) ~ 'Accelerated', 
      grepl('cons', Categories) ~ 'Conserved',
      TRUE ~ 'NULL'),
    peaktype = factor(peaktype, c('Accelerated','Conserved')),
    ) %>% inner_join(x = pheno, by = 'match') %>%
  mutate(trait = factor(trait, rev(levels(trait))))

enrichments$peaktype %>% table()

enrichments = enrichments %>% mutate(
  ## Truncate heritability predictions to be nonnative
  Observed_scale_h2 = pmax(Observed_scale_h2, 0),
  Observed_scale_h2_min = pmax(Observed_scale_h2 - Observed_scale_h2_SE, 0), 
  Observed_scale_h2_max = pmax(Observed_scale_h2 + Observed_scale_h2_SE, 0), 

  # percent of total h2g heritability
  Proportion_of_h2g = Observed_scale_h2 / h2,
  Proportion_of_h2g_min = Observed_scale_h2_min / h2,
  Proportion_of_h2g_max = Observed_scale_h2_max / h2,
  
  # SNP enrichment
  Enrichment = Proportion_of_h2g / Proportion_of_SNPs,
  Enrichment_min = Proportion_of_h2g_min / Proportion_of_SNPs,
  Enrichment_max = Proportion_of_h2g_max / Proportion_of_SNPs,
)

enrichments %>% data.frame() %>% head()
summary(enrichments$Enrichment)

## normalize the coefficients by per SNP heritability
# compute Padj w/ bonferroni family mutliple hypothesis correction
alpha = 0.1;
enrichments = enrichments %>% group_by(peaktype) %>%
  mutate(
  # normalize conditional coefficients
  Coef_norm = Coefficients / h2_perSNP,
  Coef_norm_se = Coefficient_SE / h2_perSNP, 
  # calculate 2-sided p-value for non-zero effect
  Coef_z = Observed_scale_h2 / Observed_scale_h2_SE,
  Coefficient_P_value = pnorm(-Coef_z ),
  Padj = p.adjust(Coefficient_P_value, 'fdr'), 
  logPadj = -log10(Padj),
  p.signif = ifelse(Padj < alpha, paste('FDR <',alpha), 'NS'),
  p.signif = factor(p.signif, levels = c('NS', paste('FDR <',alpha)))
  ) %>% ungroup() %>% filter(complete.cases(Coefficient_P_value))

# look at top signif PhyloP
enrichments %>% filter(Padj < alpha) %>% 
  select(match, peaktype, Padj) %>%
  group_by(match) %>% top_n(1, -Padj) %>% ungroup() %>% 
  data.frame() %>% arrange(Padj)

#################################
## make plots for presentation ##
system(paste('mkdir -p', file.path( 'plots')))
height_ppt = 4; width_ppt = 8;
height_fig = 1.75; width_fig = 2.25; font_fig = 7
plot_traits = sort(unique(enrichments$trait))

# make plots
dir.create(file.path(SETDIR, 'plots'))
plot_fn = file.path(SETDIR, 'plots', paste('COVID_prop_herit_enrichments.ppt.pdf'))
pdf(width = width_ppt, height = height_ppt, file = plot_fn, onefile = F)
## first row - proportion
pp1 = ggplot(data = enrichments, 
             aes(y = Proportion_of_h2g, x = peaktype, fill = peaktype, color = p.signif)) +
  geom_bar( stat = 'identity') + 
  geom_errorbar(aes(ymin=Proportion_of_h2g_min, ymax=Proportion_of_h2g_max), width=.4) +
  scale_fill_manual(name = "PhyloP:", values = c("#882255", "#888888")) +
  scale_color_manual(values = c(alpha('black', .4), alpha('black',1))) + 
  facet_grid( ~ trait , scales = 'free_x', space = 'free_x' ) +  
  xlab('') +  ylab('Prop. Heritability') + 
  theme_bw(base_size = 12) +
  guides(fill = FALSE, colour = guide_legend(override.aes = list(fill = 'white')),
         pattern = guide_legend(
           override.aes = list(fill = 'white', color = 'black', pattern_size = .2,
                               pattern_spacing = 0.01, pattern_density = 0.2))) + 
  theme(legend.position = "bottom", 
        legend.text=element_text(size=7), legend.title=element_text(size=7),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank())

## 2nd row - enrichments 
pp2 = ggplot(data = enrichments, 
             aes(y = Enrichment, x = peaktype, fill = peaktype, color = p.signif)) +
  geom_bar( stat = 'identity') + 
  geom_errorbar(aes(ymin=Enrichment_min, ymax=Enrichment_max), width=.4) +
  scale_fill_manual(name = "PhyloP:", values = c("#882255", "#888888")) +
  scale_color_manual(values = c(alpha('black', .4), alpha('black',1)), guide = 'none') + 
  facet_grid( ~ trait , scales = 'free', space = 'free_x' ) +  
  xlab('') +  ylab('Enrichment') + 
  theme_bw(base_size = 12) + 
  guides(fill = FALSE, colour = guide_legend(override.aes = list(fill = 'white')),
         pattern = guide_legend(
           override.aes = list(fill = 'white', color = 'black', pattern_size = .2,
                               pattern_spacing = 0.01, pattern_density = 0.2))) + 
  theme(legend.position = "bottom",
        legend.text=element_text(size=7), legend.title=element_text(size=7),
        axis.title.x = element_blank())

## make plots
pp = ggarrange(pp1, pp2, nrow = 2, common.legend = TRUE, 
               align = 'h', legend="right")
print(pp)
dev.off()
