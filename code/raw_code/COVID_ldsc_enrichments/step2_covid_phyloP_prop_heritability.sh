#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen1
#SBATCH --job-name=covid
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --error=logs/estim_heritability_%A_%a.txt
#SBATCH --output=logs/estim_heritability_%A_%a.txt
#SBATCH --array=65-67

log2results() {
awk -F '\t' '/Total Observed scale h2*/{flag=1;next}/Lambda GC/{flag=0}flag' ${OUT}.log | \
# append line with leading white space to previous
sed ':a;N;$!ba;s/\n //g'| \
# parse Partitioned heritability rownames, starts w/ colon
awk -F":" -v OFS='\t' '{gsub("[[:space:]]", "_", $1);gsub(":", "", $0); print}' | \
# transpose row to columns
awk -v OFS='\t' '
{ for (i=1; i<=NF; i++) { a[NR,i] = $i } }
NF>p { p = NF }
END { for(j=1; j<=p; j++) {
str=a[1,j];
for(i=2; i<=NR; i++){ str=str"\t"a[i,j]; }
print str } }' | \
awk -v OFS='\t' '{gsub("_[0-9]+$", "", $1); gsub("L2$", "", $1); print}' | \
gzip > ${OUT}.results.gz
}

## set up the different directors of where files are
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments; 
SNPLIST=${GWASDIR}/1000G_EUR_Phase3_GRCh38_files/hapmap3_snps/listHM3.noMHC.txt
CODEDIR=${SETDIR}/code/raw_code/COVID_ldsc_enrichments
DATADIR=${SETDIR}/data/raw_data/COVID_ldsc_enrichments
cd $CODEDIR
source activate ldsc; mkdir -p ${DATADIR}/prop_herit

# for SLURM_ARRAY_TASK_ID in {65..67}; do
## get the GWAS and reference population
GWAS=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $1}' /projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/data/tidy_data/ldsc_gwas/gwas_list_sumstats.tsv)
POP=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $2}' /projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/data/tidy_data/ldsc_gwas/gwas_list_sumstats.tsv)
GWAS_Label=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND+ 1) {print $3}' /projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/data/tidy_data/ldsc_gwas/gwas_list_sumstats.tsv)

#####################################################
## Estimate SNP heritability on caudate cell types ##
for LAB in accl.qval0.05 accl.qval0.01 cons.qval0.05 cons.qval0.01; do
## use DHS/roadmap/cCRE/Corces2020 merged peaks as bg instead of All
CELLTYPES=${DATADIR}/annotations/phyloP.am.${LAB}.
OUT=${DATADIR}/prop_herit/phyloP.am.${LAB}.${GWAS_Label}.${POP}
if [[ ! -f "${OUT}.results.gz" ]]; then # 
ldsc.py --h2 $GWAS --print-coefficients --print-snps ${SNPLIST} \
--w-ld-chr ${GWASDIR}/1000G_EUR_Phase3_GRCh38_files/weights/weights.hm3_noMHC. \
--ref-ld-chr ${CELLTYPES},${GWASDIR}/1000G_EUR_Phase3_GRCh38_files/baselineLD_v2.2/baselineLD. \
--out ${OUT} 
log2results
fi
done