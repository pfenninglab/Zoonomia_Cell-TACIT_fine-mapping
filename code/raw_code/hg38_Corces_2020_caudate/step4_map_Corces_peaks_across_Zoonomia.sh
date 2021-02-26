# directories 
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/hg38_Corces_2020_caudate
DATADIR=${SETDIR}/data/raw_data/hg38/Corces_2020
ZOONOMIADIR=${SETDIR}/data/tidy_data/Zoonomia_data
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
cd $CODEDIR

SOURCE='Homo_sapiens'
OUTDIR=${DATADIR}/haper_zoo
TARGETS=$(awk -F'\t' 'FNR >1 {print $2}' ${ZOONOMIADIR}/tables/200_Mammals_Genome_Information.tsv | \
	sed '/Homo_sapiens/d' | sed '/Carlito_syrichta/d' | sed '/Manis_tricuspis/d'| paste -s -d ',')


# ############################
# # initial halLiftover call
# for BEDFILE in ${DATADIR}/peak/Corces2020_caudate.*.narrowPeak.gz ${DATADIR}/peak/Corces2020_caudate.All_Roadmap_DHS_cCREs_BICCN_Stauffer.BG.bed.gz; do
# 	NAME=$(basename $BEDFILE .narrowPeak.gz)
# 	OUTFILE=${OUTDIR}/${NAME}.${SOURCE}To${TARGET}.HALPER.narrowPeak.gz
# 	if [[ ! -f "$OUTFILE" ]]; then 
# 		# sbatch --mem 5G -p pfen1 -w compute-1-40 ${CODEDIR}/../hal_scripts/halper_map_peak_orthologs.sh \
# 		# -s ${SOURCE} -t ${TARGETS} -o ${OUTDIR} -b ${BEDFILE}
# 	fi
# done
 

############################################
# Remove early terminated halLiftover jobs #
function checkFile(){
NUMLIFT=$(zcat $FILE | cut -f4 | cut -f1-2 -d':'| sort | uniq | wc -l | cut -d ' ' -f1)
gzip -t $FILE
if [[ $(echo $?) == '1' ]]; then 
# if gzipped file is corrupt
echo "gzipped file is corrupt: $(basename $FILE)."
rm $FILE
elif [[ $NUMLIFT < 23 ]]; then 
echo "too few source chromosomes in $FILE."
rm $FILE
fi
}

for FILE in $(ls -t ${OUTDIR}/*.gz); do
	echo $(basename $FILE)
	checkFile
done


##########################################
# recover any unmmapped Corces 2020 peaks
SOURCE='Homo_sapiens'
OUTDIR=${DATADIR}/haper_zoo
for BEDFILE in ${DATADIR}/peak/Corces2020_caudate.*.narrowPeak.gz; do
REMAINING=''
NAME=$(basename $BEDFILE .narrowPeak.gz)
TARGETS=$(awk -F'\t' 'FNR >1 {print $2}' ${ZOONOMIADIR}/tables/200_Mammals_Genome_Information.tsv | \
sed '/Homo_sapiens/d' | sed '/Carlito_syrichta/d' | sed '/Manis_tricuspis/d' | shuf)
for TARGET in ${TARGETS}; do 
# check if the narrowPeak file exists
FILE=${OUTDIR}/${NAME}.${SOURCE}To${TARGET}.HALPER.narrowPeak.gz
if [[ ! -f "$FILE" ]]; then 
	REMAINING+=",${TARGET}"; 
fi
done
if [[ $REMAINING != '' ]]; then
REMAINING=$(echo $REMAINING | sed 's/^,//g') #strip leading comma
sbatch --mem 4G -p pfen1 -w compute-1-39 ${CODEDIR}/../hal_scripts/halper_map_peak_orthologs.sh \
	-s ${SOURCE} -t ${REMAINING} -o ${OUTDIR} -b ${BEDFILE}
fi
done



########################################
# recover any unmmapped background peaks
SOURCE='Homo_sapiens'
OUTDIR=${DATADIR}/haper_zoo
BEDFILE=${DATADIR}/peak/Corces2020_caudate.All_Roadmap_DHS_cCREs_BICCN_Stauffer.BG.bed.gz
NAME=Corces2020_caudate.All_Roadmap_DHS_cCREs_BICCN_Stauffer.BG
TARGETS=$(awk -F'\t' 'FNR >1 {print $2}' ${ZOONOMIADIR}/tables/200_Mammals_Genome_Information.tsv | \
sed '/Homo_sapiens/d' | sed '/Carlito_syrichta/d' | sed '/Manis_tricuspis/d' | shuf)
REMAINING=''
for TARGET in ${TARGETS}; do OUTFILE=${OUTDIR}/${NAME}.${SOURCE}To${TARGET}.HALPER.narrowPeak.gz
if [[ ! -f "$OUTFILE" ]]; then REMAINING+=",${TARGET}"; fi
done
REMAINING=$(echo $REMAINING | sed 's/^,//g') #strip leading comma
sbatch --mem 12G -p pool3-bigmem -t 4-4 -w compute-1-35 ${CODEDIR}/../hal_scripts/halper_map_peak_orthologs.sh -s ${SOURCE} -t ${REMAINING} -o ${OUTDIR} -b ${BEDFILE}

