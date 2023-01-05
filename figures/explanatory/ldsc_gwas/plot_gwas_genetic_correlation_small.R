#######################################
### set up libraries and functions ####
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
library(ggrepel)
library(RColorBrewer)
library(rcartocolor)
library(wesanderson)
library(ComplexHeatmap)
library(ggcorrplot)
library(reshape2)
library(here)

########################################
# to be run in the root github directory
LABEL='ldsc_gwas'; DATADIR=here('data/tidy_data',LABEL)
PROJDIR=here('figures/explanatory/ldsc_gwas')
load(file = here(DATADIR,'rdas','gwas_list_sumstats.rda'))

#################################
# 1) sample size vs. effect size
height_fig = 2.5; width_fig = 4; font_fig = 7
dir.create('figures/explanatory/ldsc_gwas/plots', showWarnings = F, recursive = T)

label = c('EduAttain_E','Morningness_E', 'Depression_E', 'Schizophrenia_E',
          'SmokInitiation_E', 'OpioidDep_E','PTSD_E')
shape_signif = c('Larger_GWAS' = 23, 'Smaller_GWAS' = 21)

tmp = pheno %>% mutate(label2 = ifelse(trait %in% label, as.character(trait), NA))
plot_fn= here('figures/explanatory/ldsc_gwas/plots',
                   'ldsc_h2_magnitude_significance.pdf')
pdf(plot_fn, height = height_fig, width = width_fig)
ggplot(data = tmp, aes(x = h2_Z, y = h2, fill = group, shape = signif_group)) + 
  geom_vline(xintercept = 15, linetype= 'dashed', color = 'grey') +
  geom_point(size =1.5) + theme_bw(base_size = 8) + 
  scale_fill_manual(values = group_col, guide = 'none') + 
  scale_shape_manual(values = shape_signif, name = '') + 
  geom_label_repel(aes(label = label2), box.padding = .1, label.size = .01,
                   max.overlaps = 40, size = 2, show.legend = F,na.rm = T,
                   point.padding = .1, segment.color = 'grey50', max.time = 2,
                   min.segment.length = .1, alpha = .5, nudge_y = .25,
                   label.padding = .1, force_pull = 1, force = 20) +
  xlab('h2 confidence (Z-score) ') + ylab('SNP h2 heritability') + 
  guides( shape = guide_legend(title.position="right", nrow =1, 
                               override.aes = list(size =3))) + 
  theme(legend.position = "bottom", legend.text=element_text(size=font_fig-1),
        legend.title=element_text(size=font_fig), 
        legend.key.size = unit(.25, 'line'),
        legend.key.height=unit(.25,"line"), legend.key.width=unit(.35,"line"))
dev.off()

####################################
# 2) read in gwas correlation files
alpha = .05; n_tests = choose(nrow(pheno),2) # 1891
rg_fn = here(DATADIR,'heritability') %>% 
  list.files(path = ., pattern = '_gwasCorrelation.log',full.names = T)
names(rg_fn) = ss(basename(rg_fn),'\\.')
rg_fn = rg_fn[names(rg_fn) %in% pheno$match]
pheno$match[which(!pheno$match %in% names(rg_fn))]

# read in the gwas correlations, subset to tested GWAS
rg_df = rg_fn %>% lapply(fread, skip = 'Summary of Genetic Correlation Results', nrows = 67) %>% 
  bind_rows(.id = 'match') %>% left_join(pheno, by = 'match') %>%
  replace_na(list(rg = 0, p = 1))%>%
  mutate(p1 = basename(p1) %>% ss('\\.'),  p2 = basename(p2) %>% ss('\\.'),
         trait2 = pheno$trait[match(p2, pheno$match)]) %>%
  filter(p2 %in% p1) %>%
  mutate( rg = pmin(1, abs(rg)) * sign(rg),
          Padj = pmin(1, p * n_tests),
          logPadj = p.adjust(Padj, 'bonferroni'), 
          p.signif = case_when(
            Padj < alpha ~ paste('Pbonf <', alpha), 
            TRUE ~ "NS")) %>% arrange( abs(rg) * logPadj)

pheno = pheno %>% filter(group %in% c('SUD', "Psych")) %>% 
  filter(!trait %in% c('Cross-Psych_E')) %>% 
  filter(!grepl('_A|AUDIT|Anor|OCD', trait))
rg_df = rg_df %>% filter(trait %in% pheno$trait) %>% filter(trait2 %in% pheno$trait)

rg_df %>% pull(p.signif) %>% table()
rg_df %>% pull(trait) %>% unique() %>% length()
rg_df %>% pull(trait2) %>% unique() %>% length()

#######################################
# 3) plot the gwas correlation matrix
corMat = acast(rg_df, trait~trait2, fun.aggregate = mean, value.var = 'rg')
fMat = acast(rg_df, trait~trait2, fun.aggregate = mean, value.var = 'Padj')
diag(corMat) = 1; diag(fMat) = 0
corMat[is.na(corMat)] = 0; fMat[is.na(fMat)] = 1

# cluster the gwas correlations
corMat2 = corMat; corMat2[fMat > .5] = 0
clust = hclust(dist(1-corMat2), method = 'ward.D2')

rg_plot_fn = here(PROJDIR, 'plots/ldsc_rg_gwas_correlation_small.pdf')
pdf(rg_plot_fn, height = 4, width = 5)
pp = Heatmap(corMat, na_col = "gray", cluster_rows = clust,
        cluster_columns = clust, name = 'Correlation',
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 7),
        row_split = 3, column_split = 3,
        col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        row_title = NULL
        )

pp
dev.off()






