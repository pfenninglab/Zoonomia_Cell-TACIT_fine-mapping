#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen1
#SBATCH --time 1-0
#SBATCH --job-name=biasAway
#SBATCH --mem=10G
#SBATCH --error=logs/biasAway_%A_%a_out.txt
#SBATCH --output=logs/biasAway_%A_%a_out.txt
#SBATCH --array=1-8

SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/cnn_enhancer_ortholog
DATADIR=${SETDIR}/data/raw_data/cnn_enhancer_ortholog
cd $CODEDIR

#############################################
# get the cell type to be used for training #
CELLS=( NULL MSN_D1 MSN_D2 MSN_SN INT_Pvalb Astro Microglia Oligo OPC )
CELLTYPE=${CELLS[$SLURM_ARRAY_TASK_ID]}
MULTIPLIER=10

##################################################### 
# generate negatives from hg38, mm10, and rheMac10 ##
for GENOME in hg38 mm10 rheMac10; do
BGDIR=$HOME/resources/biasaway/${GENOME}/501bp
FGFASTA=$DATADIR/fasta/${GENOME}_${CELLTYPE}_positive.fa
BGFASTA=$DATADIR/fasta/${GENOME}_${CELLTYPE}_biasAway${MULTIPLIER}x.fa
biasaway c --foreground $FGFASTA --nfold $MULTIPLIER --deviation 2.6 \
--step 50 --seed 1 --winlen 100 --bgdirectory $BGDIR > $BGFASTA
done
