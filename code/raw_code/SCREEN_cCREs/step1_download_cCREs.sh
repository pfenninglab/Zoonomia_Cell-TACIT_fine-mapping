SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/SCREEN_cCREs
DATADIR=${SETDIR}/data/raw_data/SCREEN_cCREs
mkdir -p $DATADIR/peaks $CODEDIR/logs
cd $CODEDIR; mv $CODEDIR/annotate* $CODEDIR/logs; mv $DATADIR/annotate* $DATADIR/logs

#################################################################
# download the Meuleman et al. 2020 hg38 human index DHS regions
# wget http://gcp.wenglab.org/GRCh38-ccREs.bed -P $DATADIR/peaks
# wget http://gcp.wenglab.org/mm10-ccREs.bed -P $DATADIR/peaks
gzip $DATADIR/peaks/*.bed

# extract list of tissues, change tissue name to lower
HG_PEAKS=${DATADIR}/peaks/GRCh38-ccREs.bed.gz
MM_PEAKS=${DATADIR}/peaks/mm10-ccREs.bed.gz
zcat $HG_PEAKS $MM_PEAKS | cut -f6 | sed 's/,/\n/g' | sort | uniq | \
	awk -F '\t' 'OFS = "\t" {$2 = $1; gsub("-", "_", $2);  print $1, $2}' \
	> ${DATADIR}/cCRE_annotations.txt

################################################
# use the middle of cCRE to estimate the summit
# split cCRE by annotation, pELS, CTCF-only, etc.
NUM=$(wc -l ${DATADIR}/cCRE_annotations.txt | cut -d ' ' -f1)
for i in `seq $NUM` ; do
	MATCH=$(awk -v IND="$i" 'NR==IND{print $1}' ${DATADIR}/cCRE_annotations.txt)
	ANNOT=$(awk -F'\t' -v IND="$i" 'NR==IND{print $2}' ${DATADIR}/cCRE_annotations.txt)

	## print out the human peans
	echo "Extracting peaks for $ANNOT for human."
	zcat $HG_PEAKS | grep $MATCH | awk -F '\t' -v OFS='\t' \
	'{print $1, $2, $3, $4":"$6, 0, ".", 0, -1, -1, int(($3 - $2)/2)}' |
	gzip > ${DATADIR}/peaks/hg38_cCRES.${ANNOT}.narrowPeak.gz

	## print out the mouse peans
	echo -e "Extracting peaks for $ANNOT for mouse.\n"
	zcat $MM_PEAKS | grep $MATCH | awk -F '\t' -v OFS='\t' \
	'{print $1, $2, $3, $4":"$6, 0, ".", 0, -1, -1, int(($3 - $2)/2)}' |
	gzip > ${DATADIR}/peaks/mm10_cCRES.${ANNOT}.narrowPeak.gz
done

## print out the human peans
echo "Extracting peaks for all peaks for human."
zcat $HG_PEAKS | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $4":"$6, 0, ".", 0, -1, -1, int(($3 - $2)/2)}' |
gzip > ${DATADIR}/peaks/hg38_cCRES.All.narrowPeak.gz

## print out the mouse peans
echo -e "Extracting peaks for all peaks for mouse.\n"
zcat $MM_PEAKS | awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $4":"$6, 0, ".", 0, -1, -1, int(($3 - $2)/2)}' |
gzip > ${DATADIR}/peaks/mm10_cCRES.All.narrowPeak.gz

########################################
# lift halper peaks to rheMac8 and mm10 
HALDIR=${DATADIR}/halper; mkdir -p $HALDIR
for PEAK in ${DATADIR}/peaks/mm10_cCRES*.narrowPeak.gz; do
	sbatch --mem 20G -p pfen1 -w compute-1-40 ${SETDIR}/code/raw_code/hal_scripts/halper_map_peak_orthologs.sh \
	-b $PEAK -o $HALDIR -s Mus_musculus -t Homo_sapiens
done

for PEAK in ${DATADIR}/peaks/hg38_cCRES*.narrowPeak.gz; do
	sbatch --mem 20G -p pfen1 -w compute-1-40 ${SETDIR}/code/raw_code/hal_scripts/halper_map_peak_orthologs.sh \
	-b $PEAK -o $HALDIR -s Homo_sapiens -t Mus_musculus
done


#####################
# annotate for LDSC
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
ANNOTDIR=${DATADIR}/annotation
mkdir -p $ANNOTDIR
for BED in ${DATADIR}/peaks/hg38_cCRES*.narrowPeak.gz; do
	NAME=$(basename $BED | sed 's/.narrowPeak.gz//g')
	sbatch --job-name=${NAME}.AFR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh -i ${BED} -n ${NAME} -g hg38 -p AFR -o $ANNOTDIR
	sbatch --job-name=${NAME}.EUR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh -i ${BED} -n ${NAME} -g hg38 -p EUR -o $ANNOTDIR
done


