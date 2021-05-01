SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
SNPLIST=${GWASDIR}/listHM3.noMHC.txt
CODEDIR=${SETDIR}/code/raw_code/ldsc_caudate_zoonomia; 
DATADIR=${SETDIR}/data/raw_data/ldsc_caudate_zoonomia
ZOONOMIADIR=${SETDIR}/data/tidy_data/Zoonomia_data
PEAKDIR=${SETDIR}/data/raw_data/hg38/Corces_2020/halper_zoo
ANNOTDIR=$DATADIR/annotations
cd $CODEDIR; mkdir -p $ANNOTDIR

#################################################################
# get the species and cell types that need to be scored by CNNs
CELLS="MSN_D1 MSN_D2 MSN_SN INT_Pvalb Astro Oligo OPC Microglia"
BGANNOTDIR=${SETDIR}/data/raw_data/caudate_conservation_ldsc/annotation
BGNAME=Corces2020_caudate.All_Roadmap_DHS_cCREs_BICCN_Stauffer_phyloConsAccl.BG

##############################################
# Predict SPECIES ortholog CELLTYPE activity
for CELL in $CELLS; do
	CTS_AFR_FN=${ZOONOMIADIR}/caudate_zoonomia_mapped_predActive.${CELL}.AFR.hg38.celltypes.ldcts; > $CTS_AFR_FN
	CTS_EUR_FN=${ZOONOMIADIR}/caudate_zoonomia_mapped_predActive.${CELL}.EUR.hg38.celltypes.ldcts; > $CTS_EUR_FN
	for SLURM_ARRAY_TASK_ID in {1..240}; do
	SPECIES=$(awk -F'\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $2}' ${ZOONOMIADIR}/tables/200_Mammals_Genome_Information.tsv )
	# 5) get the hg38 coorindates of predicted active peaks in this species
	BED=${DATADIR}/peaks/Corces2020_caudate.${CELL}.Homo_sapiensTo${SPECIES}.predActive.bed.gz
	NAME=Corces2020_caudate.${CELL}.PredActiveIn${SPECIES}
	if [[ -f $BED ]]; then
		## Add the annotations from base human cell types
		echo -e "${NAME}\t${ANNOTDIR}/${NAME}.AFR.,${BGANNOTDIR}/${BGNAME}.AFR." >> $CTS_AFR_FN
		echo -e "${NAME}\t${ANNOTDIR}/${NAME}.EUR.,${BGANNOTDIR}/${BGNAME}.EUR." >> $CTS_EUR_FN
	else echo "The file ${BED} does not exist."; 
	fi
done
done


