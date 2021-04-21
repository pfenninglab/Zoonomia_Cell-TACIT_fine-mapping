# directories 
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
CODEDIR=${SETDIR}/code/raw_code/ldsc_caudate_zoonomia
DATADIR=${SETDIR}/data/raw_data/hg38/Corces_2020

# place to put peaks mappable to Zoonomia species
ZOONOMIADIR=${SETDIR}/data/tidy_data/Zoonomia_data
ANNOTDIR=$ZOONOMIADIR/annotation
ANNOTDIR2=${SETDIR}/data/raw_data/caudate_conservation_ldsc/annotation/
OUTDIR=${DATADIR}/halper_zoo
mkdir -p $ZOONOMIADIR/peaks 

checkFile(){
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
fi
done
}

# use same background for all of the peaks mappable across Zoonomia
BGANNOTDIR=${SETDIR}/data/raw_data/caudate_conservation_ldsc/annotation
BGNAME=Corces2020_caudate.All_Roadmap_DHS_cCREs_BICCN_Stauffer_phyloConsAccl_phastCons.BG
cd $CODEDIR; mv annot* logs

CELLTYPES=$(ls ${OUTDIR}/*.HALPER.narrowPeak.gz| sed 's/.*\///g;/Consensus/d;/Roadmap/d' | cut -d '.' -f2 | sort | uniq)
# CELLTYPES='Astro'

############################################
## annotate cell types mapped across species
for CELL in $CELLTYPES; do
CTS_AFR_FN=${ZOONOMIADIR}/caudate_zoonomia_mapped.${CELL}.AFR.hg38.celltypes.ldcts; > $CTS_AFR_FN
CTS_EUR_FN=${ZOONOMIADIR}/caudate_zoonomia_mapped.${CELL}.EUR.hg38.celltypes.ldcts; > $CTS_EUR_FN
## Add the annotations from base human cell types
echo -e "Corces2020_caudate.${CELL}.mappableTo.Homo_sapiens\t${ANNOTDIR2}/Corces2020_caudate.${CELL}.AFR.,${BGANNOTDIR}/${BGNAME}.AFR." >> $CTS_AFR_FN
echo -e "Corces2020_caudate.${CELL}.mappableTo.Homo_sapiens\t${ANNOTDIR2}/Corces2020_caudate.${CELL}.EUR.,${BGANNOTDIR}/${BGNAME}.EUR." >> $CTS_EUR_FN
## Get the annotations from the zoonomia species
for FILE in ${OUTDIR}/*$CELL*.HALPER.narrowPeak.gz; do
NAME=$(basename $FILE | sed 's/.HALPER.narrowPeak.gz//g' |sed 's/Homo_sapiensTo/mappableTo./g')
echo "Annotating $NAME"
BED=${ZOONOMIADIR}/peaks/${NAME}.bed.gz
if [[ ! -f "${BED}" || "${BED}" -ot $FILE ]]; then
echo "Getting human peaks for ${NAME}."
zcat $FILE | sed 's/hg38://g'| awk '{print $4}'| tr ':' '\t'| \
tr '-' '\t'| cut -f1-3 | sort -k1,1 -k2,2n | gzip > $BED
fi

## for AFR annotations
FILE2=${ANNOTDIR}/${NAME}.AFR.#.l2.M; checkFile
if [[ $HASFILE == "FALSE" ]]; then 
sbatch --mem 4G -p pfen1 --array=1-22 \
--output=/dev/null --error=/dev/null --job-name=AFR.${NAME} \
${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BED} -n ${NAME} -g hg38 -p AFR -o $ANNOTDIR
# sleep 1m
fi

# for EUR annotations
FILE2=${ANNOTDIR}/${NAME}.EUR.#.l2.M; checkFile
if [[ $HASFILE == "FALSE" ]]; then 
sbatch --mem 4G -p pfen1 --array=1-22 \
--output=/dev/null --error=/dev/null --job-name=EUR.${NAME} \
${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BED} -n ${NAME} -g hg38 -p EUR -o $ANNOTDIR
# sleep 1m
fi

echo -e "${NAME}\t${ANNOTDIR}/${NAME}.AFR.,${BGANNOTDIR}/${BGNAME}.AFR." >> $CTS_AFR_FN
echo -e "${NAME}\t${ANNOTDIR}/${NAME}.EUR.,${BGANNOTDIR}/${BGNAME}.EUR." >> $CTS_EUR_FN
done
done

