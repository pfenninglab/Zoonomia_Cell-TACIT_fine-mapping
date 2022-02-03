# directories 
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
CODEDIR=${SETDIR}/code/raw_code/ldsc_celltacit_age_decile
DATADIR=${SETDIR}/data/raw_data/ldsc_celltacit_age_decile

# place to put peaks mappable to Zoonomia species
ANNOTDIR=$DATADIR/annotation
mkdir -p $ANNOTDIR $CODEDIR/logs

checkFile(){
HASFILE=TRUE
for CHR in {1..22}; do 
F1=$(echo $FILE |sed "s/@/$CHR/g")
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

CELLTYPES="MSN_D1 MSN_D2 MSN_SN INT_Pvalb Astro Microglia OPC Oligo"
# CELL='Astro'

############################################
## annotate cell types mapped across species
for CELL in $CELLTYPES; do
CTS_AFR_FN=${DATADIR}/CellTACITage_quintile.${CELL}.AFR.hg38.celltypes.ldcts; > $CTS_AFR_FN
CTS_EUR_FN=${DATADIR}/CellTACITage_quintile.${CELL}.EUR.hg38.celltypes.ldcts; > $CTS_EUR_FN
## Get the annotations from the zoonomia species
for BED in ${DATADIR}/peaks/*$CELL*quintile*.bed.gz; do
NAME=$(basename $BED .bed.gz)
echo "Annotating $NAME"

## for AFR annotations
FILE=${ANNOTDIR}/${NAME}.AFR.@.l2.M; checkFile
if [[ $HASFILE == "FALSE" ]]; then 
sbatch --mem 4G -p short1,interactive,pool1,pool3-bigmem,pfen1 --array=1-22 \
--time 2:00:00 --output=/dev/null --error=/dev/null --job-name=AFR.${NAME} \
${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BED} -n ${NAME} -g hg38 -p AFR -o $ANNOTDIR
# sleep 1m
fi

# for EUR annotations
FILE=${ANNOTDIR}/${NAME}.EUR.@.l2.M; checkFile
if [[ $HASFILE == "FALSE" ]]; then 
sbatch --mem 4G -p short1,interactive,pool1,pool3-bigmem,pfen1 --array=1-22 \
--time 2:00:00 --output=/dev/null --error=/dev/null --job-name=EUR.${NAME} \
${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BED} -n ${NAME} -g hg38 -p EUR -o $ANNOTDIR
# sleep 1m
fi

echo -e "${NAME}\t${ANNOTDIR}/${NAME}.AFR.,${BGANNOTDIR}/${BGNAME}.AFR." >> $CTS_AFR_FN
echo -e "${NAME}\t${ANNOTDIR}/${NAME}.EUR.,${BGANNOTDIR}/${BGNAME}.EUR." >> $CTS_EUR_FN
done
done

