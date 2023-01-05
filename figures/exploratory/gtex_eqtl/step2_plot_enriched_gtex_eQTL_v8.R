library(tidyverse)
library(ggplot2)
library(data.table)
library(here)
library(rcartocolor)
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

PLOTDIR='figures/exploratory/gtex_eqtl'
DATADIR='data/tidy_data/gtex_eqtl'

dir.create(here(PLOTDIR, 'plots'), showWarnings = F, recursive = T)
dir.create(here(PLOTDIR, 'tables'), showWarnings = F, recursive = T)
dir.create(here(PLOTDIR, 'rdas'), showWarnings = F, recursive = T)

quartiles = c('0-25%', '25-50%', '50-75%', '75-100%')
names(quartiles) = paste0('quartile', 1:4); 

celltypes = c('MSN_D1', 'MSN_D2', "MSN_SN", 'INT_Pvalb',  'Astro', 
              'Microglia', 'OPC', 'Oligo')

enrichment_fn = list.files(here(DATADIR,'rdas'), full.names = T,pattern = '.glm.rds') %>% 
  str_subset('v7')
names(enrichment_fn) = basename(enrichment_fn) %>% ss('\\.', 2)

## read in the enrichment tests and correct by p.bonferroni
alpha = 0.05
snp_enrichment_df = lapply(enrichment_fn, readRDS) %>% rbindlist(idcol = 'brainQTL_celltype') %>% 
  dplyr::rename('quantile' = 'term') %>% 
  mutate(p.bonferroni = p.adjust(p.value, 'bonferroni'), 
         isSignif = p.bonferroni < alpha,
         brainQTL_celltype = 'GTEx Caudate eQTL',
         CellTACIT_celltype = factor(CellTACIT_celltype, celltypes), 
         quantile = quartiles[quantile], 
         quantile = factor(quantile, quartiles)) %>% 
  arrange(desc(estimate * -log10(p.value)))

snp_enrichment_df %>% writexl::write_xlsx(here(PLOTDIR, 'tables', 'table_SXX_GTExQTL_celltype_SNP_enrichmentByGroup_quartiles_main.xlsx'))


################################
## 4) plot the log odds ratio
plot_fn = here(PLOTDIR, 'plots', 'GTExQTL_celltype_SNP_enrichmentByGroup_quartiles_main.pdf')
pdf(plot_fn, height = .9*2, width = 0.55*2)
ggplot(snp_enrichment_df, 
       aes(x = quantile, y = OR, fill = CellTACIT_celltype, alpha = p.bonferroni < 0.05)) + 
  geom_hline(yintercept = 1, color = 'red', linetype = 'dashed') + 
  geom_errorbar(aes(ymin=OR_min, ymax=OR_max), 
                width=1, position=position_dodge(width=.75)) + 
  geom_point(pch =21, position=position_dodge(width=0.75), size = 2) +
  scale_fill_carto_d(palette = 'Safe') + 
  scale_alpha_manual(values = c(0.3, 1)) + 
  facet_wrap(~brainQTL_celltype, nrow = 1, scales = 'fixed') + 
  theme_bw(base_size = 6) + ylim(c(0.5, 3.2)) +
  xlab('Cell-TACIT Age Quantile') + 
  ylab('Odds Ratio (eQTL vs. non-eQTL)') +
  guides(fill = guide_legend(override.aes = list(size = 1),nrow = 2, title=''),
         alpha = guide_legend(nrow = 2, title='P.Bonf<0.05')) + 
  theme(plot.margin = margin(2, 2, 2, 2), legend.position = 'none', 
        legend.spacing.x = unit(.07, 'cm'), legend.spacing.y = unit(.5, 'cm'), 
        legend.key.size = unit(.2, "cm"), 
        legend.box.margin=margin(-5,-10,-5,-10),
        axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))
dev.off()

