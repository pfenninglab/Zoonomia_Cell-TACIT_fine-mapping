#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=short1,pool1,interactive
##SBATCH --partition=pfen3
##SBATCH --time 0-02:00:00
#SBATCH --job-name=(-_-)
#SBATCH --mem=20G
#SBATCH --array=1-30
#SBATCH --error=logs/perm_p_%A_%a.txt
#SBATCH --output=logs/perm_p_%A_%a.txt

source activate r4

PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
cd $PROJDIR/code/raw_code/geno_pheno/

(( IDX_CELL = ($SLURM_ARRAY_TASK_ID -1) / 6 + 1))
(( IDX_TRAIT = ($SLURM_ARRAY_TASK_ID -1) % 6 + 1))

Rscript --vanilla step3_add_pgls_empiricalPvalue_striatum_peaks.R \
-c $IDX_CELL -t $IDX_TRAIT
