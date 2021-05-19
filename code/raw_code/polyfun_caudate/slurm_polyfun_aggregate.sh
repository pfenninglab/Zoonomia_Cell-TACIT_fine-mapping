#!/bin/bash
#SBATCH --partition=pfen1
#SBATCH --job-name=aggregate
#SBATCH --mem=30G
#SBATCH --error=logs/runAggregate_%A.txt
#SBATCH --output=logs/runAggregate_%A.txt

SETWD='/projects/pfenninggroup/machineLearningForComputationalBiology/zoonomia_finemapping'
DATADIR=${SETWD}/data/raw_data/zoonomia_annotations
ANNOTDIR=${DATADIR}/annotation
POLYFUNDIR='/home/bnphan/src/polyfun'
PIP_CUTOFF=0
source ~/.bashrc; conda activate polyfun

###########################################################
# aggregate finemapping results over different SNP blocks
echo "Aggregating ${PREFIX}."
FINESNP=${OUTDIR}/../${PREFIX}.aggregate.txt.gz
if [[ ! -f $FINESNP ]]; then # --allow-missing-jobs 
python ${POLYFUNDIR}/aggregate_finemapper_results.py \
--out-prefix ${OUTDIR}/polyfun_all --sumstats ${SUMSTATS} \
--pvalue-cutoff ${CUTOFF} --out ${FINESNP}
fi

################################################
# filter out SNPs that are in the credible set
PIPSNPS=${OUTDIR}/../${PREFIX}.causal_set.txt.gz
echo "Extracting top PIP SNPs on ${CHR} with PIP greater than ${PIP_CUTOFF}."
if [[ ! -f $PIPSNPS && -f $FINESNP ]]; then
# credible set is in last column, 3rd to last column is PIP
zcat $FINESNP | awk -v VAR=${PIP_CUTOFF} '{if( $(NF - 3 ) >= VAR && $NF != 0 ) print }' | gzip > $PIPSNPS
fi

##############################################
# filter and annotate top PIP snps per chunks
ANONTSNP=${OUTDIR}/../${PREFIX}.top_annot.txt.gz
if [[ -f $FINESNP ]]; then
if [[ ! -f $ANONTSNP || $(zcat $ANONTSNP | wc -l) != $(zcat $PIPSNPS | wc -l) ]]; then
for CHR in {1..22}; do
echo "Annotating SNPs on ${CHR} with PIP greater than ${PIP_CUTOFF}."
FILE=${OUTDIR}/polyfun_all.chr${CHR}.gz
OUT=$(echo $FILE | sed 's/polyfun_all/top_annot/g')
rm $OUT
if [[ ! -f $OUT ]]; then 
# find which column is CHR
COL=$(zcat $PIPSNPS | awk -v col=CHR 'NR==1{for(i=1;i<=NF;i++) {if($i==col) {print i ;exit}} }')
zcat $PIPSNPS | awk  -v COL=$COL -v VAR=${CHR} '{if( NR == 1 || $COL == VAR ) print }' | gzip > $FILE
python ${POLYFUNDIR}/extract_annotations.py \
--annot ${ANNOTDIR}/Caudate_Zoonomia_annot_baselineLF.${CHR}.annot.parquet \
--pips $FILE --pip-cutoff $PIP_CUTOFF --out $OUT
rm $FILE; fi; done
cat ${OUTDIR}/top_annot.*.gz | zcat | awk '{if( NR == 1 || $2 ~ /^[0-9]+$/ ) print }' | gzip > $ANONTSNP
# rm ${OUTDIR}/top_annot.*.gz
fi; fi
