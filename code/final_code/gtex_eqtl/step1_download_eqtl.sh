PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
DATADIR=$PROJDIR/data/tidy_data/gtex_eqtl
mkdir -p $DATADIR; cd $DATADIR
wget https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_v8_finemapping_CAVIAR.tar
wget https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar

tar xvf GTEx_v8_finemapping_CAVIAR.tar
rm GTEx_v8_finemapping_CAVIAR.tar
cd GTEx_v8_finemapping_CAVIAR

