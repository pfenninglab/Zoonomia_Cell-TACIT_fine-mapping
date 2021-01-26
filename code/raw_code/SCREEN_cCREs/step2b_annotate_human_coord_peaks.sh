SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/SCREEN_cCREs
DATADIR=${SETDIR}/data/raw_data/SCREEN_cCREs
mkdir -p $DATADIR/peaks $CODEDIR/logs
cd $CODEDIR/logs

#####################
# annotate for LDSC
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
ANNOTDIR=${DATADIR}/annotation
mkdir -p $ANNOTDIR
for BED in ${DATADIR}/peaks/*g38*; do
	NAME=$(basename $BED | sed 's/.narrowPeak.gz//g')
	sbatch --mem 3G --job-name=${NAME}.AFR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh -i ${BED} -n ${NAME} -g hg38 -p AFR -o $ANNOTDIR
	sbatch --mem 2G --job-name=${NAME}.EUR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh -i ${BED} -n ${NAME} -g hg38 -p EUR -o $ANNOTDIR
done
