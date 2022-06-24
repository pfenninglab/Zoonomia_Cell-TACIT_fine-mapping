PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
DATADIR=$PROJDIR/data/tidy_data/gtex_eqtl
mkdir -p $DATADIR; cd $DATADIR

wget https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_v8_finemapping_CAVIAR.tar
wget https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar

tar xvf GTEx_v8_finemapping_CAVIAR.tar
rm GTEx_v8_finemapping_CAVIAR.tar

## extract only the caudate specific files
tar xvf GTEx_Analysis_v8_eQTL.tar
tar -zxvf GTEx_Analysis_v8_eQTL.tar \
GTEx_Analysis_v8_eQTL/Brain_Caudate_basal_ganglia.v8.egenes.txt.gz

tar -zxvf GTEx_Analysis_v8_eQTL.tar \
GTEx_Analysis_v8_eQTL/Brain_Caudate_basal_ganglia.v8.signif_variant_gene_pairs.txt.gz
rm GTEx_Analysis_v8_eQTL.tar

## download the all SNP-gene associations from release v7
wget https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL_all_associations.tar.gz

tar -zxvf GTEx_Analysis_v7_eQTL_all_associations.tar.gz \
GTEx_Analysis_v7_eQTL_all_associations/Brain_Caudate_basal_ganglia.allpairs.txt.gz
rm GTEx_Analysis_v7_eQTL_all_associations.tar.gz


## download the all SNP-gene associations from release v8
gsutil -u assign555 -m cp \
"gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_eQTL_all_associations/Brain_Caudate_basal_ganglia.allpairs.txt.gz" \
"gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_eQTL_all_associations/Brain_Caudate_basal_ganglia.v8.allpairs.parquet" \
${DATADIR}/GTEx_Analysis_v8_eQTL