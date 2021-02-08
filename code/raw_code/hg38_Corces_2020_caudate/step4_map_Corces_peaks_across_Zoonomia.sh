# directories 
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/hg38_Corces_2020_caudate
DATADIR=${SETDIR}/data/raw_data/hg38/Corces_2020
ZOONOMIADIR=${SETDIR}/data/tidy_data/Zoonomia_data
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
cd $CODEDIR

##################
# for foreground! $2 ~ /Homo_sapiens/ 
SOURCE='Homo_sapiens'
OUTDIR=${DATADIR}/haper_zoo
TARGETS=$(awk -F'\t' 'FNR >1 {print $2}' ${ZOONOMIADIR}/tables/200_Mammals_Genome_Information.tsv | \
	sed '/Homo_sapiens/d' | sed '/Carlito_syrichta/d' | sed '/Manis_tricuspis/d'| paste -s -d ',')

for BEDFILE in ${DATADIR}/peak/Corces2020_caudate.*.narrowPeak.gz ${DATADIR}/peak/Corces2020_caudate.All_Roadmap_DHS_cCREs_BICCN_Stauffer.BG.bed.gz; do
	NAME=$(basename $BEDFILE .narrowPeak.gz)
	OUTFILE=${OUTDIR}/${NAME}.${SOURCE}To${TARGET}.HALPER.narrowPeak.gz
	if [[ ! -f "$OUTFILE" ]]; then 
		# sbatch --mem 5G -p pfen1 -w compute-1-40 ${CODEDIR}/../hal_scripts/halper_map_peak_orthologs.sh \
		# -s ${SOURCE} -t ${TARGETS} -o ${OUTDIR} -b ${BEDFILE}
	fi
done


################################
# recover any unmmapped regions
SOURCE='Homo_sapiens'
OUTDIR=${DATADIR}/haper_zoo

for BEDFILE in ${DATADIR}/peak/Corces2020_caudate.MSN_D2.narrowPeak.gz; do
REMAINING=''
NAME=$(basename $BEDFILE .narrowPeak.gz)
TARGETS=$(awk -F'\t' 'FNR >1 {print $2}' ${ZOONOMIADIR}/tables/200_Mammals_Genome_Information.tsv | \
sed '/Homo_sapiens/d' | sed '/Carlito_syrichta/d' | sed '/Manis_tricuspis/d' | shuf)
for TARGET in ${TARGETS}; do
OUTFILE=${OUTDIR}/${NAME}.${SOURCE}To${TARGET}.HALPER.narrowPeak.gz
if [[ ! -f "$OUTFILE" ]]; then 
echo $TARGET
REMAINING+=",${TARGET}"		
fi
done
REMAINING=$(echo $REMAINING | sed 's/^,//g') #strip leading comma
# sbatch --mem 5G -p pfen1 -w compute-1-40 ${CODEDIR}/../hal_scripts/halper_map_peak_orthologs.sh \
# 	-s ${SOURCE} -t ${REMAINING} -o ${OUTDIR} -b ${BEDFILE}
done


