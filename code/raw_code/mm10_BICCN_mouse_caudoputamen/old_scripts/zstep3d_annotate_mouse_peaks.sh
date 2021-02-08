#####################
## directories 
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/mm10_BICCN_mouse_caudoputamen
DATADIR=${SETDIR}/data/raw_data/mm10/Mouse_cSNAIL_D1D2
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
ANNOTDIR=${SETDIR}/data/raw_data/caudate_conservation_ldsc/annotation
mkdir -p $ANNOTDIR
cd $CODEDIR

#####################
# annotate for LDSC
for BED in ${DATADIR}/halper/*Mus_musculusToHomo_sapiens.HALPER.narrowPeak.gz; do
	NAME=$(basename $BED | sed 's/.narrowPeak.gz//g'| sed 's/.Mus_musculusToHomo_sapiens.HALPER//g')
	sbatch --mem 15G -p pool1 --job-name=${NAME}.AFR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
			-i ${BED} -n ${NAME} -g hg38 -p AFR -o $ANNOTDIR
	sbatch --mem 15G -p pool1 --job-name=${NAME}.EUR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
			-i ${BED} -n ${NAME} -g hg38 -p EUR -o $ANNOTDIR
done

#####################
# annotate for LDSC
CODEDIR=${SETDIR}/code/raw_code/mm10_BICCN_mouse_caudoputamen
DATADIR=${SETDIR}/data/raw_data/mm10/BICCN_mouse_caudoputamen
for BED in ${DATADIR}/halper/*Mus_musculusToHomo_sapiens.HALPER.narrowPeak.gz; do
	NAME=$(basename $BED | sed 's/.narrowPeak.gz//g'| sed 's/.Mus_musculusToHomo_sapiens.HALPER//g')
	sbatch --mem 5G -p pool1 --job-name=${NAME}.AFR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
			-i ${BED} -n ${NAME} -g hg38 -p AFR -o $ANNOTDIR
	sbatch --mem 5G -p pool1 --job-name=${NAME}.EUR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
			-i ${BED} -n ${NAME} -g hg38 -p EUR -o $ANNOTDIR
done

