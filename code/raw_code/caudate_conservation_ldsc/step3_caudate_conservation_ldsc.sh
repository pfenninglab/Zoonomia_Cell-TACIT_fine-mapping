#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1
#SBATCH --time 1-0
#SBATCH --job-name=ldsc_caud
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --error=logs/ldsc_caud_%A_%a.txt
#SBATCH --output=logs/ldsc_caud_%A_%a.txt
#SBATCH --array=65-68

# get the GWAS for this array job
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
GWASDIR=${SETDIR}/data/tidy_data/ldsc_gwas; SNPLIST=${GWASDIR}/listHM3.noMHC.txt
CODEDIR=${SETDIR}/code/raw_code/caudate_conservation_ldsc; DATADIR=${SETDIR}/data/raw_data/caudate_conservation_ldsc
cd $CODEDIR; source activate ldsc

# get the GWAS and reference population
GWAS=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $1}' ${GWASDIR}/gwas_list_sumstats.tsv)
POP=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $2}' ${GWASDIR}/gwas_list_sumstats.tsv)
GWAS_Label=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND+ 1) {print $3}' ${GWASDIR}/gwas_list_sumstats.tsv)

OUTDIR=${DATADIR}/enrichments; mkdir -p $OUTDIR

###################################################################
# run LD score regression over the Corces caudate binary annotations
CTS_FN1=${DATADIR}/caudate_conservation_binary_${POP}_hg38_celltypes.ldcts
if [[ ! -f "$OUTDIR/caudate_conservation_binary.${GWAS_Label}.${POP}.cell_type_results.txt" ]]; then
	ldsc.py --ref-ld-chr-cts $CTS_FN1 \
	--ref-ld-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/baseline_v1.1/baseline_v1.1.${POP}. \
	--w-ld-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/weights/1000G.${POP}.weights.hm3_noMHC. \
	--h2-cts $GWAS --out $OUTDIR/caudate_conservation_binary.${GWAS_Label}.${POP}
fi

# #################################################################################
# # run LD score regression over the Corces caudate peaks intersect Zoonomia phyloP
# CTS_FN2=${DATADIR}/caudate_conservation_phyloP_${POP}_hg38_celltypes.ldcts
# if [[ ! -f "$OUTDIR/caudate_conservation_phyloP.${GWAS_Label}.${POP}.cell_type_results.txt" ]]; then
# 	ldsc.py --ref-ld-chr-cts $CTS_FN2 \
# 	--ref-ld-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/baseline_v1.1/baseline_v1.1.${POP}. \
# 	--w-ld-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/weights/1000G.${POP}.weights.hm3_noMHC. \
# 	--h2-cts $GWAS --out $OUTDIR/caudate_conservation_phyloP.${GWAS_Label}.${POP}
# fi

