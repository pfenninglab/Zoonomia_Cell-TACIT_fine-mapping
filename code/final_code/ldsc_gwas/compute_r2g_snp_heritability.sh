#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=interactive
#SBATCH --time 0-8
#SBATCH --job-name=est_rg
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --error=/dev/null
#SBATCH --output=/dev/null
#SBATCH --array=1-68

# get the GWAS for this array job
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
GWASDIR=${SETDIR}/data/tidy_data/ldsc_gwas; SNPLIST=${GWASDIR}/listHM3.noMHC.txt
cd $SETDIR/code/final_code/ldsc_gwas; source activate ldsc

# for SLURM_ARRAY_TASK_ID in {1..67}; do
# get the GWAS and reference population
GWAS=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $1}' ${GWASDIR}/gwas_list_sumstats.tsv)
POP=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $2}' ${GWASDIR}/gwas_list_sumstats.tsv)
GWAS_Label=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND+ 1) {print $3}' ${GWASDIR}/gwas_list_sumstats.tsv)
GWAS2=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR!=(IND+ 1) && NR > 1 {print $1}' ${GWASDIR}/gwas_list_sumstats.tsv | tr '\n' ',' | sed 's/,$//g')
OUTDIR=$GWASDIR/heritability; mkdir -p $OUTDIR

#####################################
# run compute h2 SNP heritability
if [[ ! -f ${OUTDIR}/${GWAS_Label}.${POP}_gwasCorrelation.log ]]; then
ldsc.py --rg ${GWAS},${GWAS2} \
--ref-ld-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/weights/1000G.${POP}.weights.hm3_noMHC. \
--w-ld-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/weights/1000G.${POP}.weights.hm3_noMHC.  \
--out ${OUTDIR}/${GWAS_Label}.${POP}_gwasCorrelation
fi
# done
