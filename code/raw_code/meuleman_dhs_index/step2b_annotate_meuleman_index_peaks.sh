SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/meuleman_dhs_index
DATADIR=${SETDIR}/data/raw_data/meuleman_dhs_index
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
mv $CODEDIR/annotate* $CODEDIR/logs; mv $DATADIR/annotate* $DATADIR/logs; cd $CODEDIR

#######################################
## merge human epigenome backgrounds ##
BG1=/home/bnphan/src/atac-seq-pipeline/genome/hg38/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.hg19_to_hg38.bed.gz
BG2=${SETDIR}/data/raw_data/meuleman_dhs_index/peaks/DHS_Index_and_Vocabulary_mappedToMm10.All.full.narrowPeak.gz
BG3=${SETDIR}/data/raw_data/SCREEN_cCREs/peaks/GRCh38-ccREs.bed.gz
BGNAME=meuleman_dhs_index.DHS_Roadmap_cCREs.BG
BGFILE=${DATADIR}/${BGNAME}.bed.gz

# cat ${BG1} ${BG2} ${BG3} | zcat | cut -f 1-3 | \
# 	sort --parallel=10 -k1,1 -k2,2n | gzip > ${DATADIR}/${BGNAME}.tmp.bed.gz 
# bedtools merge -i ${DATADIR}/${BGNAME}.tmp.bed.gz | gzip > ${BGFILE}
# rm ${DATADIR}/${BGNAME}.tmp.bed.gz

#############################################
## annotate for LDSC w/ binary annotations ##
ANNOTDIR=${DATADIR}/annotation; mkdir -p $ANNOTDIR

## for background
sbatch --mem 5G -p pfen3 --job-name=${BGNAME}.AFR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BGFILE} -n ${BGNAME} -g hg38 -p AFR -o $ANNOTDIR
sbatch --mem 5G -p pfen3 --job-name=${BGNAME}.EUR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BGFILE} -n ${BGNAME} -g hg38 -p EUR -o $ANNOTDIR

# for foreground
CTS_AFR_FN=${DATADIR}/Meuleman_DHS_binary_AFR_hg38_celltypes.ldcts; > $CTS_AFR_FN
CTS_EUR_FN=${DATADIR}/Meuleman_DHS_binary_EUR_hg38_celltypes.ldcts; > $CTS_EUR_FN
for BED in $(ls ${DATADIR}/peaks/*.narrowPeak.gz | sed '/All/d'); do
NAME=$(basename $BED | sed 's/.narrowPeak.gz//g')
if [[ ! -f "${ANNOTDIR}/${NAME}.AFR.1.l2.ldscore.gz" ]]; then 
echo "Annotations for ${NAME} not found."
sbatch --mem 5G -p pfen3 --job-name=${NAME}.AFR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BED} -n ${NAME} -g hg38 -p AFR -o $ANNOTDIR
fi
if [[ ! -f "${ANNOTDIR}/${NAME}.EUR.1.l2.ldscore.gz" ]]; then 
echo "Annotations for ${NAME} not found."
sbatch --mem 5G -p pfen3 --job-name=${NAME}.EUR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BED} -n ${NAME} -g hg38 -p EUR -o $ANNOTDIR
fi
echo -e "${NAME}\t${ANNOTDIR}/${NAME}.AFR.,${ANNOTDIR}/${BGNAME}.AFR." >> $CTS_AFR_FN
echo -e "${NAME}\t${ANNOTDIR}/${NAME}.EUR.,${ANNOTDIR}/${BGNAME}.EUR." >> $CTS_EUR_FN
done



#####################################
## annotate for phyloP traits LDSC ##
ANNOTDIR2=${DATADIR}/annot_phyloP; mkdir -p $ANNOTDIR2
BWFILE=${GWASDIR}/phyloP/200m_scoresPhyloP_20201221.@.bigWig

## for background
sbatch --mem 5G -p pfen3 --job-name=${BGNAME}.AFR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BGFILE} -w ${BWFILE} -n ${BGNAME} -g hg38 -p AFR -o $ANNOTDIR2
sbatch --mem 5G -p pfen3 --job-name=${BGNAME}.EUR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BGFILE} -w ${BWFILE} -n ${BGNAME} -g hg38 -p EUR -o $ANNOTDIR2

# for foreground
CTS_AFR_FN2=${DATADIR}/Meuleman_DHS_phyloP_AFR_hg38_celltypes.ldcts; > $CTS_AFR_FN2
CTS_EUR_FN2=${DATADIR}/Meuleman_DHS_phyloP_EUR_hg38_celltypes.ldcts; > $CTS_EUR_FN2
for BED in $(ls ${DATADIR}/peaks/*.narrowPeak.gz | sed '/All/d'); do
NAME=$(basename $BED | sed 's/.narrowPeak.gz//g')
if [[ ! -f "${ANNOTDIR2}/${NAME}.AFR.1.l2.ldscore.gz" ]]; then
echo "Annotations for ${NAME} not found."
sbatch --mem 5G -p pfen3 --job-name=${NAME}.AFR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BED} -n ${NAME} -g hg38 -p AFR -o $ANNOTDIR2 -w ${BWFILE}
fi
if [[ ! -f "${ANNOTDIR2}/${NAME}.EUR.1.l2.ldscore.gz" ]]; then 
echo "Annotations for ${NAME} not found."
sbatch --mem 5G -p pfen3 --job-name=${NAME}.EUR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BED} -n ${NAME} -g hg38 -p EUR -o $ANNOTDIR2 -w ${BWFILE}
fi
echo -e "${NAME}\t${ANNOTDIR2}/${NAME}.AFR.,${ANNOTDIR2}/${BGNAME}.AFR." >> $CTS_AFR_FN2
echo -e "${NAME}\t${ANNOTDIR2}/${NAME}.EUR.,${ANNOTDIR2}/${BGNAME}.EUR." >> $CTS_EUR_FN2
done





