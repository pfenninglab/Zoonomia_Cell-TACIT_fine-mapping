# directories 
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
ZOONOMIADIR=${SETDIR}/data/tidy_data/Zoonomia_data
CODEDIR=${SETDIR}/code/raw_code/snp_annotations
DATADIR=${SETDIR}/data/raw_data/snp_annotations
OUTDIR=${DATADIR}/halper_zoo
cd $CODEDIR

##########################################
# Remove early terminated halLiftover jobs
for FILE in ${OUTDIR}/*.halLiftover.tFile.bed.gz; do
gzip -t $FILE
# if gzipped file is corrupt
if [[ $(echo $?) == '1' ]]; then 
echo "Removing corrupt file: $(basename $FILE)."
rm $FILE
fi
done

####################################
# recover the unmapped grch38 snps
SOURCE='Homo_sapiens'
TARGETS=$(awk -F'\t' 'FNR >1 {print $2}' ${ZOONOMIADIR}/tables/200_Mammals_Genome_Information.tsv | \
	sed '/Homo_sapiens/d' | sed '/Carlito_syrichta/d' | sed '/Manis_tricuspis/d' | shuf)
BEDFILE=${DATADIR}/bed/UKB_imputed_snps.GRCh38.unmappedlifted.bed.gz
for TARGET in $TARGETS; do
NAME=$(basename $BEDFILE .bed.gz)
OUTFILE=${OUTDIR}/${NAME}.${SOURCE}To${TARGET}.halLiftover.tFile.bed.gz
if [[ ! -f "$OUTFILE" ]]; then 
echo "not found $(basename $OUTFILE)"
sbatch --mem 3G -p pfen1 -w compute-1-39 ${CODEDIR}/../hal_scripts/halper_map_peak_orthologs.sh \
-s ${SOURCE} -t ${TARGET} -o ${OUTDIR} -b ${BEDFILE} --snp
fi
done


###############################
# recover the chromosome jobs
for TARGET in $TARGETS; do
for CHR in {1..22}; do
BEDFILE=${DATADIR}/bed/UKB_imputed_snps.GRCh38.${CHR}.bed.gz
NAME=UKB_imputed_snps.GRCh38.${CHR}
OUTFILE=${OUTDIR}/${NAME}.${SOURCE}To${TARGET}.halLiftover.tFile.bed.gz
if [[ ! -f "$OUTFILE" ]]; then 
echo "not found $(basename $OUTFILE)"
sbatch --mem 3G -p pfen1 -w compute-1-39 ${CODEDIR}/../hal_scripts/halper_map_peak_orthologs.sh \
-s ${SOURCE} -t ${TARGET} -o ${OUTDIR} -b ${BEDFILE} --snp
fi
done
done
