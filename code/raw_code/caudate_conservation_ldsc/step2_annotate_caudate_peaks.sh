#################
## directories ##
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/caudate_conservation_ldsc
DATADIR=${SETDIR}/data/raw_data/hg38/Corces_2020
MACDIR=${SETDIR}/data/raw_data/rheMac10/Stauffer_caudate/halper
MUSDIR1=${SETDIR}/data/raw_data/mm10/BICCN_mouse_caudoputamen/halper
MUSDIR2=${SETDIR}/data/raw_data/mm10/Mouse_cSNAIL_D1D2/halper
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
PHYLODIR=${GWASDIR}/phyloP/phyloP_cutoff_regions
ANNOTDIR=${SETDIR}/data/raw_data/caudate_conservation_ldsc/annotation
mkdir -p $ANNOTDIR

cd $CODEDIR; mv annot* logs

#######################################
## merge human epigenome backgrounds ##
BG1=${DATADIR}/peak/Corces2020_caudate.Consensus.narrowPeak.gz
BG2=${MUSDIR1}/BICCN_CP.Consensus.Mus_musculusToHomo_sapiens.HALPER.narrowPeak.gz
BG3=${MACDIR}/Stauffer_caudate.Consensus.GenBankRheMac8.Macaca_mulattaToHomo_sapiens.HALPER.narrowPeak.gz
BG4=${SETDIR}/data/raw_data/SCREEN_cCREs/peaks/GRCh38-ccREs.bed.gz
BG5=${SETDIR}/data/raw_data/meuleman_dhs_index/peaks/DHS_Index_and_Vocabulary_mappedToMm10.All.full.narrowPeak.gz
BG6=/home/bnphan/src/atac-seq-pipeline/genome/hg38/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.hg19_to_hg38.bed.gz
PHYLOCONS=${PHYLODIR}/phyloPam_cons.fdr.05.hg38.bed.gz
PHYLOACCL=${PHYLODIR}/phyloPam_accl.fdr.05.hg38.bed.gz
BGNAME=Corces2020_caudate.All_Roadmap_DHS_cCREs_BICCN_Stauffer_phyloConsAccl_phastCons.BG
BGFILE=${DATADIR}/peak/${BGNAME}.bed.gz

# cat ${BG1} ${BG2} ${BG3} ${BG4} ${BG5} ${BG6} ${PHYLOCONS} ${PHYLOACCL} ${PHYLODIR}/PhastCons* \
# 	${MUSDIR1}/*Mus_musculusToHomo_sapiens.HALPER.narrowPeak.gz \
# 	${MUSDIR2}/*Mus_musculusToHomo_sapiens.HALPER.narrowPeak.gz | \
# 	zcat | cut -f 1-3 | sort --parallel=10 -k1,1 -k2,2n | gzip > ${DATADIR}/${BGNAME}.tmp.bed.gz 
# bedtools merge -i ${DATADIR}/${BGNAME}.tmp.bed.gz | gzip > ${BGFILE}
# rm ${DATADIR}/${BGNAME}.tmp.bed.gz

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

#######################################
# for background for binary annotations
FILE=${ANNOTDIR}/${BGNAME}.AFR.#.l2.M; checkFile
if [[ $HASFILE == "FALSE" ]]; then 
	sbatch --mem 24G -p gpu --time 12:00:00 --job-name=${BGNAME}.AFR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
	-i ${BGFILE} -n ${BGNAME} -g hg38 -p AFR -o $ANNOTDIR
fi

FILE=${ANNOTDIR}/${BGNAME}.EUR.#.l2.M; checkFile
if [[ $HASFILE == "FALSE" ]]; then 
	sbatch --mem 24G -p gpu --time 12:00:00 --job-name=${BGNAME}.EUR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
	-i ${BGFILE} -n ${BGNAME} -g hg38 -p EUR -o $ANNOTDIR
fi

#############################################
## annotate forgrounds for LDSC w/ binary annotations ##
# for foreground
CTS_AFR_FN=${SETDIR}/data/raw_data/caudate_conservation_ldsc/caudate_conservation_binary_AFR_hg38_celltypes.ldcts; > $CTS_AFR_FN
CTS_EUR_FN=${SETDIR}/data/raw_data/caudate_conservation_ldsc/caudate_conservation_binary_EUR_hg38_celltypes.ldcts; > $CTS_EUR_FN
FILES=$(ls ${DATADIR}/peak/*.narrowPeak.gz ${MACDIR}/*ToHomo_sapiens.HALPER.narrowPeak.gz \
${MUSDIR1}/*ToHomo_sapiens.HALPER.narrowPeak.gz ${MUSDIR2}/*Mus_musculusToHomo_sapiens.HALPER.narrowPeak.gz \
${PHYLODIR}/phylo*.fdr.05.hg38.bed.gz ${PHYLODIR}/PhastCons*.bed.gz )
for BED in $FILES; do
# for AFR annotations
NAME=$(basename $BED | sed 's/.narrowPeak.gz//g;s/_viterbi_conservedRegions//g;s/.fdr.05.hg38//g')
FILE=${ANNOTDIR}/${NAME}.AFR.#.l2.M; checkFile
if [[ $HASFILE == "FALSE" ]]; then 
echo "Annotations for ${NAME} not found."
sbatch --mem 24G -p gpu --time 12:00:00 --job-name=${NAME}.AFR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BED} -n ${NAME} -g hg38 -p AFR -o $ANNOTDIR
fi

# for EUR annotations
FILE=${ANNOTDIR}/${NAME}.EUR.#.l2.M; checkFile
if [[ $HASFILE == "FALSE" ]]; then 
echo "Annotations for ${NAME} not found."
sbatch --mem 24G -p gpu --time 12:00:00 --job-name=${NAME}.EUR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BED} -n ${NAME} -g hg38 -p EUR -o $ANNOTDIR
fi
echo -e "${NAME}\t${ANNOTDIR}/${NAME}.AFR.,${ANNOTDIR}/${BGNAME}.AFR." >> $CTS_AFR_FN
echo -e "${NAME}\t${ANNOTDIR}/${NAME}.EUR.,${ANNOTDIR}/${BGNAME}.EUR." >> $CTS_EUR_FN
done














#####################################
## annotate for phyloP traits LDSC ##
ANNOTDIR2=${SETDIR}/data/raw_data/caudate_conservation_ldsc/annot_phyloP
BWFILE=${GWASDIR}/phyloP/200m_scoresPhyloP_20210214.@.bigWig
mkdir -p $ANNOTDIR2

## for foreground
CTS_AFR_FN2=${SETDIR}/data/raw_data/caudate_conservation_ldsc/caudate_conservation_phyloP_AFR_hg38_celltypes.ldcts; > $CTS_AFR_FN2
CTS_EUR_FN2=${SETDIR}/data/raw_data/caudate_conservation_ldsc/caudate_conservation_phyloP_EUR_hg38_celltypes.ldcts; > $CTS_EUR_FN2
for BED in ${DATADIR}/peak/*.narrowPeak.gz ${MUSDIR}/*ToHomo_sapiens.HALPER.narrowPeak.gz ${MACDIR}/*ToHomo_sapiens.HALPER.narrowPeak.gz ${MUSDIR2}/*Mus_musculusToHomo_sapiens.HALPER.narrowPeak.gz; do
NAME=$(basename $BED | sed 's/.narrowPeak.gz//g')
if [[ ! -f "${ANNOTDIR2}/${NAME}.AFR.1.l2.M" ]]; then
echo "Annotations for ${NAME} not found."
sbatch --mem 2G -p gpu --time 12:00:00 --job-name=${NAME}.AFR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BED} -n ${NAME} -g hg38 -p AFR -o $ANNOTDIR2 -w ${BWFILE}
fi
if [[ ! -f "${ANNOTDIR2}/${NAME}.EUR.1.l2.M" ]]; then 
echo "Annotations for ${NAME} not found."
sbatch --mem 2G -p gpu --time 12:00:00 --job-name=${NAME}.EUR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BED} -n ${NAME} -g hg38 -p EUR -o $ANNOTDIR2 -w ${BWFILE}
fi
echo -e "${NAME}\t${ANNOTDIR2}/${NAME}.AFR.,${ANNOTDIR2}/${BGNAME}.AFR." >> $CTS_AFR_FN2
echo -e "${NAME}\t${ANNOTDIR2}/${NAME}.EUR.,${ANNOTDIR2}/${BGNAME}.EUR." >> $CTS_EUR_FN2
done

## for background phyloP annotations
if [[ ! -f "${ANNOTDIR2}/${BGNAME}.AFR.1.l2.M" ]]; then 
sbatch --mem 24G -p gpu --time 12:00:00 --job-name=${BGNAME}.AFR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BGFILE} -w ${BWFILE} -n ${BGNAME} -g hg38 -p AFR -o $ANNOTDIR2
fi
if [[ ! -f "${ANNOTDIR2}/${BGNAME}.EUR.1.l2.M" ]]; then 
sbatch --mem 24G -p gpu --time 12:00:00 --job-name=${BGNAME}.EUR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BGFILE} -w ${BWFILE} -n ${BGNAME} -g hg38 -p EUR -o $ANNOTDIR2
fi


