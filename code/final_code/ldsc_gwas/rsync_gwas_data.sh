DATADIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
TODIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/data/tidy_data/ldsc_gwas

rsync -Pav $DATADIR/gwas/munged $TODIR


