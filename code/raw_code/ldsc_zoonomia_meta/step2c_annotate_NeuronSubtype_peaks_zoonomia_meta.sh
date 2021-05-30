# directories 
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
CODEDIR=${SETDIR}/code/raw_code/ldsc_zoonomia_meta
DATADIR=${SETDIR}/data/raw_data/ldsc_zoonomia_meta
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
PHYLODIR=${GWASDIR}/phyloP/phyloP_cutoff_regions
ANNOTDIR=$DATADIR/annotation
mkdir -p $ANNOTDIR $CODEDIR/logs
cd $CODEDIR; mv annot* logs

#######################################
## merge human epigenome backgrounds ##
BG1=${SETDIR}/data/raw_data/SCREEN_cCREs/peaks/GRCh38-ccREs.bed.gz
BG2=/home/bnphan/src/atac-seq-pipeline/genome/hg38/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.hg19_to_hg38.bed.gz
PHYLOCONS=${PHYLODIR}/phyloPcons.241mam.fdr.05.hg38.bed.gz
PHYLOACCL=${PHYLODIR}/phyloPaccl.241mam.fdr.05.hg38.bed.gz
PHASTCONS=${GWASDIR}/phyloP/phastCons_cutoff_regions/phastCons.43prim.fdr.05.hg38.bed.gz
BGNAME=BICCN_mouse_CATlas.All_Roadmap_cCREs_phyloConsAccl_phastCons.BG
BGFILE=${DATADIR}/peaks/${BGNAME}.bed.gz

# cat ${BG1} ${BG2} ${PHYLOCONS} ${PHYLOACCL} ${PHASTCONS} \
# 	${DATADIR}/peaks/BICCN_mouse_CATlas*Primates#0.mappable.bed.gz | \
# 	zcat | cut -f 1-3 | sort --parallel=10 -k1,1 -k2,2n | gzip > ${DATADIR}/${BGNAME}.tmp.bed.gz 
# bedtools merge -i ${DATADIR}/${BGNAME}.tmp.bed.gz | gzip > ${BGFILE}
# rm ${DATADIR}/${BGNAME}.tmp.bed.gz


checkFile(){
HASFILE=TRUE
for CHR in {1..22}; do 
F1=$(echo $FILE |sed "s/@/$CHR/g")
if [[ ! -f $F1 ]]; then 
# if no annotation found for chr
HASFILE=FALSE; echo "File not found ${F1}."; break ;
elif [[ $(cat $F1) == "0" ]]; then
# if annotation crashed and 0 SNPs annotated
HASFILE=FALSE; echo "Has 0 annotation column in ${F1}."; break ;
elif [[ $F1 -ot $BED ]]; then
# if newer foreground file
HASFILE=FALSE; echo "Newer foreground"; break ;
fi
done
}

#######################################
# for background for binary annotations
FILE=${ANNOTDIR}/${BGNAME}.AFR.@.l2.M; checkFile
if [[ $HASFILE == "FALSE" ]]; then 
	sbatch --mem 50G -p gpu,pfen1,pfen3 --time 12:00:00 --job-name=${BGNAME}.AFR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
	-i ${BGFILE} -n ${BGNAME} -g hg38 -p AFR -o $ANNOTDIR
fi

FILE=${ANNOTDIR}/${BGNAME}.EUR.@.l2.M; checkFile
if [[ $HASFILE == "FALSE" ]]; then 
	sbatch --mem 50G -p gpu,pfen1,pfen3 --time 12:00:00 --job-name=${BGNAME}.EUR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
	-i ${BGFILE} -n ${BGNAME} -g hg38 -p EUR -o $ANNOTDIR
fi

CELLTYPES="PV SST VIP"

############################################
## annotate cell types mapped across species
for CELL in $CELLTYPES; do
CTS_AFR_FN=${DATADIR}/zoonomia_meta.${CELL}.AFR.hg38.celltypes.ldcts; > $CTS_AFR_FN
CTS_EUR_FN=${DATADIR}/zoonomia_meta.${CELL}.EUR.hg38.celltypes.ldcts; > $CTS_EUR_FN
## Get the annotations from the zoonomia species
for BED in $(ls ${DATADIR}/peaks/*$CELL*.bed.gz | sed '/BG/d'); do
NAME=$(basename $BED .bed.gz)
echo "Annotating $NAME"

## for AFR annotations
FILE=${ANNOTDIR}/${NAME}.AFR.@.l2.M; checkFile
if [[ $HASFILE == "FALSE" ]]; then #--output=/dev/null --error=/dev/null 
sbatch --mem 10G -p short1,interactive,pool1,pool3-bigmem,gpu,pfen1 --array=1-22 \
--time 2:00:00 --output=/dev/null --error=/dev/null --job-name=AFR.${NAME} \
${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BED} -n ${NAME} -g hg38 -p AFR -o $ANNOTDIR
# sleep 1m
fi

# for EUR annotations
FILE=${ANNOTDIR}/${NAME}.EUR.@.l2.M; checkFile
if [[ $HASFILE == "FALSE" ]]; then # 
sbatch --mem 10G -p short1,interactive,pool1,pool3-bigmem,gpu,pfen1 --array=1-22 \
--time 2:00:00 --output=/dev/null --error=/dev/null --job-name=EUR.${NAME} \
${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BED} -n ${NAME} -g hg38 -p EUR -o $ANNOTDIR
# sleep 1m
fi

echo -e "${NAME}\t${ANNOTDIR}/${NAME}.AFR.,${ANNOTDIR}/${BGNAME}.AFR." >> $CTS_AFR_FN
echo -e "${NAME}\t${ANNOTDIR}/${NAME}.EUR.,${ANNOTDIR}/${BGNAME}.EUR." >> $CTS_EUR_FN
done
done

