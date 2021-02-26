# directories 
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
CODEDIR=${SETDIR}/code/raw_code/hg38_Corces_2020_caudate
DATADIR=${SETDIR}/data/raw_data/hg38/Corces_2020

# place to put peaks mappable to Zoonomia species
ZOONOMIADIR=${SETDIR}/data/tidy_data/Zoonomia_data
ANNOTDIR=$ZOONOMIADIR/annotation
ANNOTDIR2=${SETDIR}/data/raw_data/caudate_conservation_ldsc/annotation/
OUTDIR=${DATADIR}/haper_zoo
mkdir -p $ZOONOMIADIR/peaks 

checkFile(){
HASFILE=TRUE
for CHR in {1..22}; do 
F1=$(echo $FILE |sed "s/#/$CHR/g")
if [[ ! -f $F1 ]]; then 
	# if no annotation found for chr
	HASFILE=FALSE
elif [[ $(cat $F1) == "0" ]]; then
	# if annotation crashed and 0 SNPs annotated
	HASFILE=FALSE
fi
done
}

# use same background for all of the peaks mappable across Zoonomia
BGANNOTDIR=${SETDIR}/data/raw_data/caudate_conservation_ldsc/annotation
BGNAME=Corces2020_caudate.All_Roadmap_DHS_cCREs_BICCN_Stauffer_phyloConsAccl.BG
cd $CODEDIR; mv annot* logs

CELLTYPES=$(ls ${OUTDIR}/*.HALPER.narrowPeak.gz| sed 's/.*\///g' | cut -d '.' -f2 | sort | uniq)
CELLTYPES="MSN_D1 MSN_D2 MSN_SN INT_Pvalb Astro Oligo OPC Microglia"
# CELLTYPES="MSN_D2 MSN_D1 MSN_SN"

## annotate cell types mapped across species
for CELL in $CELLTYPES; do
CTS_AFR_FN=${ZOONOMIADIR}/caudate_zoonomia_mapped.${CELL}.AFR.hg38.celltypes.ldcts; > $CTS_AFR_FN
CTS_EUR_FN=${ZOONOMIADIR}/caudate_zoonomia_mapped.${CELL}.EUR.hg38.celltypes.ldcts; > $CTS_EUR_FN
## Add the annotations from base human cell types
echo -e "Corces2020_caudate.${CELL}.mappableTo.Homo_sapiens\t${ANNOTDIR2}/Corces2020_caudate.${CELL}.AFR.,${BGANNOTDIR}/${BGNAME}.AFR." >> $CTS_AFR_FN
echo -e "Corces2020_caudate.${CELL}.mappableTo.Homo_sapiens\t${ANNOTDIR2}/Corces2020_caudate.${CELL}.EUR.,${BGANNOTDIR}/${BGNAME}.EUR." >> $CTS_EUR_FN
for FILE in ${OUTDIR}/*$CELL*.HALPER.narrowPeak.gz; do
NAME=$(basename $FILE | sed 's/.HALPER.narrowPeak.gz//g' |sed 's/Homo_sapiensTo/mappableTo./g')
echo "Annotating $NAME"
BED=${ZOONOMIADIR}/peaks/${NAME}.bed.gz
if [[ ! -f "${BED}" ]]; then 
zcat $FILE | sed 's/hg38://g'| awk '{print $4}'| tr ':' '\t'| \
	tr '-' '\t'| cut -f1-3 | sort -k1,1 -k2,2n | gzip > $BED
fi

# for AFR annotations
FILE=${ANNOTDIR}/${NAME}.AFR.#.l2.M; checkFile
if [[ $HASFILE == "FALSE" ]]; then 
sbatch --mem 4G -p short1,gpu,interactive,pool1,pool3-bigmem --time=2:00:00 --array=1-22%1 \
--output=/dev/null --error=/dev/null \
--job-name=${NAME}.AFR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BED} -n ${NAME} -g hg38 -p AFR -o $ANNOTDIR
sleep 1m
fi

# for EUR annotations
FILE=${ANNOTDIR}/${NAME}.EUR.#.l2.M; checkFile
if [[ $HASFILE == "FALSE" ]]; then 
sbatch --mem 4G -p short1,gpu,interactive,pool1,pool3-bigmem --time=2:00:00 --array=1-22%1 \
--output=/dev/null --error=/dev/null \
--job-name=${NAME}.AFR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BED} -n ${NAME} -g hg38 -p EUR -o $ANNOTDIR
sleep 1m
fi

echo -e "${NAME}\t${ANNOTDIR}/${NAME}.AFR.,${BGANNOTDIR}/${BGNAME}.AFR." >> $CTS_AFR_FN
echo -e "${NAME}\t${ANNOTDIR}/${NAME}.EUR.,${BGANNOTDIR}/${BGNAME}.EUR." >> $CTS_EUR_FN
done
done

