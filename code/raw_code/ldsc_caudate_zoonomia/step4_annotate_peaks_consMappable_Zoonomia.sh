# directories 
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
CODEDIR=${SETDIR}/code/raw_code/ldsc_caudate_zoonomia
DATADIR=${SETDIR}/data/raw_data/ldsc_caudate_zoonomia

# place to put peaks mappable to Zoonomia species
ZOONOMIADIR=${SETDIR}/data/tidy_data/Zoonomia_data
ANNOTDIR=$DATADIR/annotation
HALPERDIR=${SETDIR}/data/raw_data/hg38/Corces_2020/halper_zoo
mkdir -p $ZOONOMIADIR/peaks ${DATADIR}/peaks_consMappable ${DATADIR}/peaks $DATADIR/annotation

CONS_BED=${DATADIR}/peaks/phyloPam_cons.fdr.05.hg38.bed
if [[ ! -f $CONS_BED ]]; then zcat ${GWASDIR}/phyloP/phyloP_cutoff_regions/phyloPam_cons.fdr.05.hg38.bed.gz > $CONS_BED; fi

# use same background for all of the peaks mappable across Zoonomia
BGANNOTDIR=${SETDIR}/data/raw_data/caudate_conservation_ldsc/annotation
BGNAME=Corces2020_caudate.All_Roadmap_DHS_cCREs_BICCN_Stauffer_phyloConsAccl.BG
cd $CODEDIR;
CELLTYPES="MSN_D1 MSN_D2 MSN_SN INT_Pvalb Astro Microglia OPC Oligo"


##############################################################
## 1) get the cell types mappable peaks overlapping w/ phyloP
for CELL in $CELLTYPES; do
CTS_AFR_FN=${ZOONOMIADIR}/caudate_zoonomia_mapped_phyloP.${CELL}.AFR.hg38.celltypes.ldcts; > $CTS_AFR_FN
CTS_EUR_FN=${ZOONOMIADIR}/caudate_zoonomia_mapped_phyloP.${CELL}.EUR.hg38.celltypes.ldcts; > $CTS_EUR_FN
## Add the annotations from base human cell types
echo "Annotating human $CELL overlapping phyloPam_cons."
BED=${SETDIR}/data/raw_data/hg38/Corces_2020/peak/Corces2020_caudate.${CELL}.narrowPeak.gz
NAME2=Corces2020_caudate.${CELL}.consMappable.Homo_sapiens
BED2=${DATADIR}/peaks_consMappable/${NAME2}.bed.gz
if [[ ! -f "${BED2}" || "${BED2}" -ot $FILE || "${BED2}" -ot $BED ]]; then
echo "Getting human peaks in ${CELL} overlapping phyloPam_cons."
bedtools intersect -u -a <(zcat $BED | sort -k1,1 -k2,2n) -b $CONS_BED | gzip > $BED2; fi
echo -e "${NAME2}\t${ANNOTDIR}/${NAME2}.AFR.,${BGANNOTDIR}/${BGNAME}.AFR." >> $CTS_AFR_FN
echo -e "${NAME2}\t${ANNOTDIR}/${NAME2}.EUR.,${BGANNOTDIR}/${BGNAME}.EUR." >> $CTS_EUR_FN

## get the other species cell types
for FILE in $HALPERDIR/*$CELL*.HALPER.narrowPeak.gz; do
NAME=$(basename $FILE | sed 's/.HALPER.narrowPeak.gz//g' |sed 's/Homo_sapiensTo/mappableTo./g')
NAME2=$(basename $FILE | sed 's/.HALPER.narrowPeak.gz//g' |sed 's/Homo_sapiensTo/consMappable./g')
BED=${ZOONOMIADIR}/peaks/${NAME}.bed.gz
BED2=${DATADIR}/peaks_consMappable/${NAME2}.bed.gz
if [[ ! -f "${BED2}" || "${BED2}" -ot $FILE || "${BED2}" -ot $BED ]]; then
echo "Getting peaks in ${NAME} overlapping phyloPam_cons."
bedtools intersect -u -a <(zcat $BED | sort -k1,1 -k2,2n) -b $CONS_BED | gzip > $BED2; fi
echo -e "${NAME2}\t${ANNOTDIR}/${NAME2}.AFR.,${BGANNOTDIR}/${BGNAME}.AFR." >> $CTS_AFR_FN
echo -e "${NAME2}\t${ANNOTDIR}/${NAME2}.EUR.,${BGANNOTDIR}/${BGNAME}.EUR." >> $CTS_EUR_FN
done
done


####################################################
## 2) annotate the peaks mappable & overlap PhyloP
mv annot* logs
checkFile(){ 
HASFILE=TRUE
for CHR in {1..22}; do F1=$(echo $FILE2 |sed "s/#/$CHR/g")
if [[ ! -f $F1 ]]; then HASFILE=FALSE # if no annotation found for chr
elif [[ $(cat $F1) == "0" ]]; then HASFILE=FALSE # if annotation crashed and 0 SNPs annotated
elif [[ $F1 -ot $BED2 ]]; then HASFILE=FALSE; fi; done
} # if newer foreground file

for CELL in $CELLTYPES; do
CTS_AFR_FN=${ZOONOMIADIR}/caudate_zoonomia_mapped_phyloP.${CELL}.AFR.hg38.celltypes.ldcts
for NAME2 in $(cut $CTS_AFR_FN -f1); do
BED2=${DATADIR}/peaks_consMappable/${NAME2}.bed.gz

### for AFR annotations
FILE2=${ANNOTDIR}/${NAME2}.AFR.#.l2.M; checkFile
if [[ $HASFILE == "FALSE" ]]; then 
sbatch --mem 4G -p short1,interactive,pool1,pfen1 --array=1-22 --time 2:00:00 \
--output=/dev/null --error=/dev/null --job-name=AFR.${NAME2} \
${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BED2} -n ${NAME2} -g hg38 -p AFR -o $ANNOTDIR
sleep 15s; fi

### for EUR annotations
FILE2=${ANNOTDIR}/${NAME2}.EUR.#.l2.M; checkFile
if [[ $HASFILE == "FALSE" ]]; then 
sbatch --mem 4G -p short1,interactive,pool1,pfen1 --array=1-22 --time 2:00:00 \
--output=/dev/null --error=/dev/null --job-name=EUR.${NAME2} \
${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BED2} -n ${NAME2} -g hg38 -p EUR -o $ANNOTDIR
sleep 15s; fi
done; done


