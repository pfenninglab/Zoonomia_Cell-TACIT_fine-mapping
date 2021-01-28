SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/meuleman_dhs_index
DATADIR=${SETDIR}/data/raw_data/meuleman_dhs_index
mv $CODEDIR/annotate* $CODEDIR/logs; mv $DATADIR/annotate* $DATADIR/logs; cd $CODEDIR

#####################
# annotate for LDSC
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
ANNOTDIR=${DATADIR}/annotation
mkdir -p $ANNOTDIR

for BED in ${DATADIR}/peaks/*.narrowPeak.gz; do
	NAME=$(basename $BED | sed 's/.narrowPeak.gz//g')
	if [ ! -f ${ANNOTDIR}/${NAME}.AFR.1.l2.ldscore.gz ]; then echo
		sbatch --mem 5G --job-name=${NAME}.AFR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
				-i ${BED} -n ${NAME} -g hg38 -p AFR -o $ANNOTDIR
		sbatch --mem 5G --job-name=${NAME}.EUR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
				-i ${BED} -n ${NAME} -g hg38 -p EUR -o $ANNOTDIR
	fi
done
