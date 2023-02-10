# directories 
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
CODEDIR=${SETDIR}/code/raw_code/ldsc_atac_and_conservation
OUTDIR=${SETDIR}/data/raw_data/ldsc_atac_and_conservation
ANNOTDIR=$OUTDIR/annotation
mkdir -p $OUTDIR/peak_meta $ANNOTDIR
cd $CODEDIR; mv annot* logs

checkFile(){
HASFILE=TRUE
for CHR in {1..22}; do 
F1=$(echo $FILE2 |sed "s/@/$CHR/g")
if [[ ! -f $F1 ]]; then 
echo "no annotation found for chr${CHR}."
HASFILE=FALSE
# elif [[ $(cat $F1) == "0" ]]; then
# echo "annotation crashed and 0 SNPs annotated for chr${CHR}."
# HASFILE=FALSE
elif [[ $F1 -ot $BED ]]; then
echo "newer foreground file than annotations for chr${CHR}."
HASFILE=FALSE
fi
done
}

# use same background for all of the peaks mappable across Zoonomia
BGANNOTDIR=${SETDIR}/data/raw_data/caudate_conservation_ldsc/annotation
BGNAME=Corces2020_caudate.All_Roadmap_DHS_cCREs_BICCN_Stauffer_phyloConsAccl_phastCons.BG
PEAKTYPES="Enhancer Promoter Other"
PEAKTYPES="Other"

for PEAKTYPE in $PEAKTYPES; do
CTS_AFR_FN=${OUTDIR}/conservation_peak_meta.${PEAKTYPE}.AFR.hg38.celltypes.ldcts; > $CTS_AFR_FN
CTS_EUR_FN=${OUTDIR}/conservation_peak_meta.${PEAKTYPE}.EUR.hg38.celltypes.ldcts; > $CTS_EUR_FN

############################################
## annotate cell types mapped across species
for BED in ${OUTDIR}/peak/Corces2020.*.${PEAKTYPE}.predActive.bed.gz ${OUTDIR}/peak_meta/Corces2020.*.${PEAKTYPE}.*.bed.gz; do
NAME=$(basename $BED .bed.gz)
echo "Working on ${NAME}."

## check that AFR annotations are complete
FILE2=${ANNOTDIR}/${NAME}.AFR.@.l2.M; checkFile
if [[ $HASFILE == "FALSE" ]]; then 
sbatch --mem 10G -p pfen1 --array=1-22 \
--output=/dev/null --error=/dev/null --job-name=AFR.${NAME} \
${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BED} -n ${NAME} -g hg38 -p AFR -o $ANNOTDIR
sleep 5s
fi

## check that EUR annotations are complete
FILE2=${ANNOTDIR}/${NAME}.EUR.@.l2.M; checkFile
if [[ $HASFILE == "FALSE" ]]; then 
sbatch --mem 10G -p pfen1 --array=1-22 \
--output=/dev/null --error=/dev/null --job-name=EUR.${NAME} \
${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BED} -n ${NAME} -g hg38 -p EUR -o $ANNOTDIR
sleep 5s
fi

echo -e "${NAME}\t${ANNOTDIR}/${NAME}.AFR.,${BGANNOTDIR}/${BGNAME}.AFR." >> $CTS_AFR_FN
echo -e "${NAME}\t${ANNOTDIR}/${NAME}.EUR.,${BGANNOTDIR}/${BGNAME}.EUR." >> $CTS_EUR_FN
done
done
