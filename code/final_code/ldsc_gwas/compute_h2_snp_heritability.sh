#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1
#SBATCH --time 0-8
#SBATCH --job-name=est_h2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --error=/dev/null
#SBATCH --output=/dev/null
#SBATCH --array=1-67

# get the GWAS for this array job
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
GWASDIR=${SETDIR}/data/tidy_data/ldsc_gwas; SNPLIST=${GWASDIR}/listHM3.noMHC.txt
cd $SETDIR; source activate ldsc

# for SLURM_ARRAY_TASK_ID in {1..67}; do
# get the GWAS and reference population
GWAS=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $1}' ${GWASDIR}/gwas_list_sumstats.tsv)
POP=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $2}' ${GWASDIR}/gwas_list_sumstats.tsv)
GWAS_Label=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND+ 1) {print $3}' ${GWASDIR}/gwas_list_sumstats.tsv)

OUTDIR=$GWASDIR/heritability; mkdir -p $OUTDIR

#####################################
# run compute h2 SNP heritability
ldsc.py --h2 $GWAS \
--ref-ld-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/weights/1000G.${POP}.weights.hm3_noMHC. \
--w-ld-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/weights/1000G.${POP}.weights.hm3_noMHC.  \
--out ${OUTDIR}/${GWAS_Label}.${POP}
# done
