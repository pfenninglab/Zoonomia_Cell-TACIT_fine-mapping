#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1,pool3-bigmem,pfen1,pfen_bigmem
#SBATCH --time=24:00:00
#SBATCH --job-name=est_herit
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=15G
#SBATCH --error=logs/estim_heritability_%A_%a.txt
#SBATCH --output=logs/estim_heritability_%A_%a.txt
#SBATCH --array=1-64%32

log2results() {
awk -F '\t' '/Total Observed scale h2*/{flag=1;next}/Lambda GC/{flag=0}flag' ${OUT}.log | \
# append line with leading white space to previous
sed ':a;N;$!ba;s/\n //g'| \
# parse Partitioned heritability rownames, starts w/ colon
awk -F":" -v OFS='\t' '{gsub("[[:space:]]", "_", $1);gsub(":", "", $0); print}' | \
# transpose row to columns
awk -v OFS='\t' '
{ 
for (i=1; i<=NF; i++)  {
a[NR,i] = $i
}
}
NF>p { p = NF }
END {    
for(j=1; j<=p; j++) {
str=a[1,j]
for(i=2; i<=NR; i++){
str=str"\t"a[i,j];
}
print str
}
}' | awk -v OFS='\t' '{gsub("_[0-9]+$", "", $1); gsub("L2$", "", $1); print}' | \
gzip > ${OUT}.results.gz
}

## set up the different directors of where files are
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
GWASDIR=${SETDIR}/data/tidy_data/ldsc_gwas; SNPLIST=${GWASDIR}/listHM3.noMHC.txt
CODEDIR=${SETDIR}/code/raw_code/caudate_conservation_ldsc; DATADIR=${SETDIR}/data/raw_data/caudate_conservation_ldsc
cd $CODEDIR; mkdir -p ${DATADIR}/prop_herit
source ~/.bashrc; conda activate ldsc

## get the GWAS and reference population
GWAS=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $1}' ${GWASDIR}/gwas_list_sumstats.tsv)
POP=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $2}' ${GWASDIR}/gwas_list_sumstats.tsv)
GWAS_Label=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND+ 1) {print $3}' ${GWASDIR}/gwas_list_sumstats.tsv)

#####################################################
## Estimate SNP heritability on caudate cell types ##
CTS_FN1=${DATADIR}/caudate_conservation_binary_${POP}_hg38_celltypes.ldcts
for IND in $(wc -l $CTS_FN1 | cut -d ' ' -f1 | xargs seq ); do
## use DHS/roadmap/cCRE/Corces2020 merged peaks as bg instead of All
LAB=$(awk -v IND=$IND 'FNR == IND {print $1}' $CTS_FN1)
CELLTYPES=$(awk -v IND=$IND 'FNR == IND {print $2}' $CTS_FN1)
OUT=${DATADIR}/prop_herit/caudate_conservation_binary.${LAB}.${GWAS_Label}.${POP}
if [[ ! -f "${OUT}.results.gz" ]]; then # --print-snps ${SNPLIST} 
ldsc.py --h2 $GWAS --print-coefficients \
--w-ld-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/weights/1000G.${POP}.weights.hm3_noMHC. \
--ref-ld-chr ${CELLTYPES},${GWASDIR}/1000G_ALL_Phase3_hg38_files/baselineLD_v2.2/baselineLD_v2.2.${POP}. \
--out ${OUT} 
log2results
fi
done

######################################################################
## merge together results file of the single heritability estimates ##
gunzip ${DATADIR}/prop_herit/caudate_conservation_binary.*.${GWAS_Label}.${POP}.results.gz
awk ' FNR == 2 || NR==1 {print}' \
	${DATADIR}/prop_herit/caudate_conservation_binary.*.${GWAS_Label}.${POP}.results |
	gzip > ${DATADIR}/prop_herit/caudate_conservation_binary.${GWAS_Label}.${POP}.agg.gz
gzip ${DATADIR}/prop_herit/caudate_conservation_binary.*.${GWAS_Label}.${POP}.results

