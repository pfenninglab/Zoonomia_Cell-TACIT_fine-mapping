#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=interactive,gpu
#SBATCH --time=8:00:00
#SBATCH --job-name=zooM1Liver
#SBATCH --dependency=afterok:1690518
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --error=logs/est_herit_predActive_M1ctx_liver_%A_%a.txt
#SBATCH --output=logs/est_herit_predActive_M1ctx_liver_%A_%a.txt
#SBATCH --array=1-64%10

log2results() {
awk -F '\t' '/Total Observed scale h2*/{flag=1;next}/Lambda GC/{flag=0}flag' ${OUT}.log | \
# append line with leading white space to previous
sed ':a;N;$!ba;s/\n //g'| \
# parse Partitioned heritability rownames, starts w/ colon
awk -F":" -v OFS='\t' '{gsub("[[:space:]]", "_", $1);gsub(":", "", $0); print}' | \
# transpose row to columns
awk -v OFS='\t' '
{ 
for (i=1; i<=NF; i++){
a[NR,i] = $i
}
}
NF>p { p = NF }
END {    
for(j=1; j<=p; j++){
str=a[1,j]
for(i=2; i<=NR; i++){
str=str"\t"a[i,j];
}
print str
}
}' | awk -v OFS='\t' '{gsub("_[0-9]+$", "", $1); gsub("L2$", "", $1); print}' | \
gzip > ${OUT}.results.gz
}

# get the GWAS for this array job
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
GWASDIR=${SETDIR}/data/tidy_data/ldsc_gwas; SNPLIST=${GWASDIR}/listHM3.noMHC.txt
CODEDIR=${SETDIR}/code/raw_code/ldsc_zoonomia_meta; 
DATADIR=${SETDIR}/data/raw_data/ldsc_zoonomia_meta
OUTDIR=${DATADIR}/prop_herit_M1ctx_liver; mkdir -p $OUTDIR
cd $CODEDIR; 

source ~/.bashrc; conda activate ldsc 

# get the GWAS and reference population
GWAS=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $1}' ${GWASDIR}/gwas_list_sumstats.tsv)
POP=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $2}' ${GWASDIR}/gwas_list_sumstats.tsv)
GWAS_Label=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND+ 1) {print $3}' ${GWASDIR}/gwas_list_sumstats.tsv)
CELLS="M1ctx Liver"

###################################################################
# run LD score regression over the Corces caudate binary annotations
for CELL in $CELLS; do
echo "Working on ${CELL}."
CTS_FN=${DATADIR}/zoonomia_meta.${CELL}.${POP}.hg38.celltypes.ldcts
for IND in $(wc -l $CTS_FN | cut -d ' ' -f1 | xargs seq ); do
## use DHS/roadmap/cCRE/Corces2020 merged peaks as bg instead of All
LAB=$(awk -v IND=$IND 'FNR == IND {print $1}' $CTS_FN)
CELLTYPES=$(awk -v IND=$IND 'FNR == IND {print $2}' $CTS_FN)
CELL2=$(echo $CELLTYPES | cut -d ',' -f1)
OUT=${OUTDIR}/zoonomia_meta.M1ctx_liver.${GWAS_Label}.${POP}.${CELL}.${LAB}
if [[ ! -f "${OUT}.results.gz" || "${OUT}.results.gz" -ot ${CELL2}.1.l2.M ]]; then
ldsc.py --h2 $GWAS --w-ld-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/weights/1000G.${POP}.weights.hm3_noMHC. \
--ref-ld-chr ${CELLTYPES},${GWASDIR}/1000G_ALL_Phase3_hg38_files/baselineLD_v2.2/baselineLD_v2.2.${POP}. \
--print-coefficients --out ${OUT} 
log2results; # rm ${OUT}.log
if [[ $(zcat ${OUT}.results.gz) == '' ]]; then rm ${OUT}.results.gz; fi
fi
done

## extract only the top cell type enrichment
OUTFILE=${OUTDIR}/zoonomia_meta.M1ctx_liver.${GWAS_Label}.${POP}.${CELL}.agg.gz
# if [[ ! -f $OUTFILE || "$OUTFILE" -ot ${CELL2}.1.l2.M ]]; then
echo "Aggregating for ${CELL}."
gunzip ${OUTDIR}/zoonomia_meta.M1ctx_liver.${GWAS_Label}.${POP}.${CELL}*.results.gz
awk ' FNR == 2 || NR==1 {print}' \
	${OUTDIR}/zoonomia_meta.M1ctx_liver.${GWAS_Label}.${POP}.${CELL}*.results |
	gzip >  ${OUTFILE}
gzip ${OUTDIR}/zoonomia_meta.M1ctx_liver.${GWAS_Label}.${POP}.${CELL}*.results
# fi
done


