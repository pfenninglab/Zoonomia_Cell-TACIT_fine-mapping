# directories
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/ldsc_caudate_zoonomia
DATADIR=${SETDIR}/data/raw_data/hg38/Corces_2020
ZOONOMIADIR=${SETDIR}/data/tidy_data/Zoonomia_data
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
ANNOTDIR=$ZOONOMIADIR/annotation

cd $CODEDIR

SOURCE='Homo_sapiens'; OUTDIR=${DATADIR}/halper_zoo
TARGETS=$(awk -F'\t' 'FNR >1 {print $2}' ${ZOONOMIADIR}/tables/200_Mammals_Genome_Information.tsv | \
sed '/Homo_sapiens/d' | sed '/Carlito_syrichta/d' | sed '/Manis_tricuspis/d'| paste -s -d ',')


############################
# initial halLiftover call
for BEDFILE in ${DATADIR}/peak/Corces2020_caudate.*.narrowPeak.gz ${DATADIR}/peak/Corces2020_caudate.All_Roadmap_DHS_cCREs_BICCN_Stauffer.BG.bed.gz; do
	NAME=$(basename $BEDFILE .narrowPeak.gz)
	OUTFILE=${OUTDIR}/${NAME}.${SOURCE}To${TARGET}.HALPER.narrowPeak.gz
	if [[ ! -f "$OUTFILE" ]]; then 
		sbatch --mem 5G -p pfen1 -w compute-1-40 ${CODEDIR}/../hal_scripts/halper_map_peak_orthologs.sh \
		-s ${SOURCE} -t ${TARGETS} -o ${OUTDIR} -b ${BEDFILE}
	fi
done


###############################################
# 2) Remove early terminated halLiftover jobs #
function checkFile(){
HASFILE=TRUE
NUMLIFT=$(zcat $FILE | cut -f4 | cut -f1-2 -d':'| sort | uniq | wc -l | cut -d ' ' -f1)
gzip -t $FILE
if [[ $(echo $?) != '0' ]]; then 
# if gzipped file is corrupt
echo "gzipped file is corrupt: $(basename $FILE)."
HASFILE=FALSE; rm $FILE
elif [[ $NUMLIFT < 23 ]]; then 
echo "too few source chromosomes in $FILE."
HASFILE=FALSE; rm $FILE; fi
}

checkFile2(){
HASFILE=TRUE
for CHR in {1..22}; do 
F1=$(echo $FILE2 |sed "s/#/$CHR/g")
if [[ ! -f $F1 ]]; then 
# if no annotation found for chr
HASFILE=FALSE
elif [[ $(cat $F1) == "0" ]]; then
# if annotation crashed and 0 SNPs annotated
HASFILE=FALSE
elif [[ $F1 -ot $BED ]]; then
# if newer foreground file
HASFILE=FALSE
fi; done
}

#### remove files that didn't have complete annotations
TOLOOK=$(ls ${OUTDIR}/*.HALPER.narrowPeak.gz | sed '/Consensus/d;/All/d')
for FILE in $TOLOOK; do
NAME=$(basename $FILE | sed 's/.HALPER.narrowPeak.gz//g' |sed 's/Homo_sapiensTo/mappableTo./g')
# check if mapping is right
checkFile
checkFile2
if [[ $HASFILE == "FALSE" ]]; then
echo "Removing $(basename $FILE .HALPER.narrowPeak.gz) files."
rm $(echo $FILE | sed 's/.HALPER.narrowPeak.gz/.halLiftover.sFile.bed.gz/g') $(echo $FILE | sed 's/.HALPER.narrowPeak.gz/.halLiftover.tFile.bed.gz/g')
fi

# # check if annotations are right
# if [[ $HASFILE == "TRUE" ]]; then
# FILE2=${ANNOTDIR}/${NAME}.AFR.#.l2.M; checkFile2
# echo "Checking on $(basename $FILE .HALPER.narrowPeak.gz) annotations."
# if [[ $HASFILE == "FALSE" ]]; then
# echo "Removing $(basename $FILE .HALPER.narrowPeak.gz) files."
# rm $FILE $(echo $FILE | sed 's/.HALPER.narrowPeak.gz/.halLiftover.sFile.bed.gz/g') $(echo $FILE | sed 's/.HALPER.narrowPeak.gz/.halLiftover.tFile.bed.gz/g')
# fi; fi; done

######
# TOLOOK=$(ls -r ${OUTDIR}/*.HALPER.narrowPeak.gz | sed '/Consensus/d;/All/d')
# for FILE in $TOLOOK; do
# echo $(basename $FILE)
# checkFile
# done

##########################################
# 3) recover any unmmapped Corces 2020 peaks
SOURCE='Homo_sapiens'
OUTDIR=${DATADIR}/halper_zoo
PEAKS=$(ls ${DATADIR}/peak/Corces2020_caudate.*.narrowPeak.gz | sed '/Consensus/d;/All/d')
for BEDFILE in $PEAKS ; do
REMAINING=''
NAME=$(basename $BEDFILE .narrowPeak.gz)
TARGETS=$(awk -F'\t' 'FNR >1 {print $2}' ${ZOONOMIADIR}/tables/200_Mammals_Genome_Information.tsv | \
sed '/Homo_sapiens/d' | sed '/Carlito_syrichta/d' | sed '/Manis_tricuspis/d' | shuf)
for TARGET in ${TARGETS}; do 
# check if the narrowPeak file exists
FILE=${OUTDIR}/${NAME}.${SOURCE}To${TARGET}.HALPER.narrowPeak.gz
if [[ ! -f "$FILE" ]]; then 
echo "No ${TARGET} found for ${NAME}."
REMAINING+=",${TARGET}";
fi
done
if [[ $REMAINING != '' ]]; then
REMAINING=$(echo $REMAINING | sed 's/^,//g') #strip leading comma
sbatch --mem 10G -p pfen1 -w compute-1-39,compute-1-40 ${CODEDIR}/../hal_scripts/halper_map_peak_orthologs.sh \
	-s ${SOURCE} -t ${REMAINING} -o ${OUTDIR} -b ${BEDFILE}
fi
done

find ${SETDIR}/data/tidy_data/Zoonomia_data/peaks -type f -name '*.gz' -exec bash -c 'echo -e "$(basename $1)\t$(zcat $1 | wc -l)"' dummy {} \; > ${SETDIR}/data/tidy_data/Zoonomia_data/tables/zoonomia_peaks.txt


########################################
# recover any unmmapped background peaks
SOURCE='Homo_sapiens'
OUTDIR=${DATADIR}/halper_zoo
BEDFILE=${DATADIR}/peak/Corces2020_caudate.All_Roadmap_DHS_cCREs_BICCN_Stauffer_phyloConsAccl_phastCons.BG.bed.gz
NAME=Corces2020_caudate.All_Roadmap_DHS_cCREs_BICCN_Stauffer_phyloConsAccl_phastCons.BG
TARGETS=$(awk -F'\t' 'FNR >1 {print $2}' ${ZOONOMIADIR}/tables/200_Mammals_Genome_Information.tsv | \
sed '/Homo_sapiens/d' | sed '/Carlito_syrichta/d' | sed '/Manis_tricuspis/d' | shuf)
REMAINING=''
for TARGET in ${TARGETS}; do OUTFILE=${OUTDIR}/${NAME}.${SOURCE}To${TARGET}.HALPER.narrowPeak.gz
if [[ ! -f "$OUTFILE" ]]; then REMAINING+=",${TARGET}"; fi
done
REMAINING=$(echo $REMAINING | sed 's/^,//g') #strip leading comma
sbatch --mem 12G -p pool3-bigmem -t 4-4 -w compute-1-35 ${CODEDIR}/../hal_scripts/halper_map_peak_orthologs.sh -s ${SOURCE} -t ${REMAINING} -o ${OUTDIR} -b ${BEDFILE}

