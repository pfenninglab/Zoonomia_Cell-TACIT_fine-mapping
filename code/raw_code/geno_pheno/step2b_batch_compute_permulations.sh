#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1,pool3-bigmem
##SBATCH --partition=pfen1
#SBATCH --time 3-00:00:00
#SBATCH --job-name=(°o°)
#SBATCH --mem=23G
#SBATCH --array=1-113
#SBATCH --error=logs/permulations_%A_%a.txt
#SBATCH --output=logs/permulations_%A_%a.txt

source activate r4

PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
cd $PROJDIR/code/raw_code/geno_pheno/

(( IDX_CELL = ($SLURM_ARRAY_TASK_ID -1) / 36 + 1))
(( NUM2 = ($SLURM_ARRAY_TASK_ID -1) % 36 + 1))
(( IDX_TRAIT = ($NUM2 -1) / 6 + 1))
(( IDX_PERM = ($NUM2 -1) % 6 + 1))

Rscript --vanilla step2_calculate_permulation_pvalue_striatum_peaks.R \
-c $IDX_CELL -t $IDX_TRAIT -p $IDX_PERM
