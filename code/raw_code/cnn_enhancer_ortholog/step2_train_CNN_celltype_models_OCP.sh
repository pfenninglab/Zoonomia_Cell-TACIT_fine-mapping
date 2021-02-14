#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=gpu
#SBATCH --time=0-12
#SBATCH --job-name=cnn_ocp
#SBATCH --gres=gpu:1
#SBATCH --mem=62G
#SBATCH --array=5-8
#SBATCH --error=logs/cnn_ocp_%A_%a.txt
#SBATCH --output=logs/cnn_ocp_%A_%a.txt

SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/cnn_enhancer_ortholog
DATADIR=${SETDIR}/data/raw_data/cnn_enhancer_ortholog
cd $CODEDIR

#############################################
# get the cell type to be used for training #
CELLS=( NULL MSN_D1 MSN_D2 MSN_SN INT_Pvalb Astro Microglia Oligo OPC )
CELLTYPE=${CELLS[$SLURM_ARRAY_TASK_ID]}
LABEL=${CELLTYPE}_hgRmMm_enhVsNonEnhOrth

source activate tf2

#############################################
### merge training set from each genome #####
TRAINPOSFILE=$DATADIR/fasta/${CELLTYPE}_trainPos.fa
cat $DATADIR/fasta/hg38_${CELLTYPE}_train_positive.fa $DATADIR/fasta/mm10_${CELLTYPE}_train_positive.fa \
	$DATADIR/fasta/rheMac10_${CELLTYPE}_train_positive.fa > $TRAINPOSFILE

TRAINNEGFILE=$DATADIR/fasta/${CELLTYPE}_trainNeg.fa
cat $DATADIR/fasta/hg38_${CELLTYPE}_train_negative.fa $DATADIR/fasta/mm10_${CELLTYPE}_train_negative.fa \
	$DATADIR/fasta/rheMac10_${CELLTYPE}_train_negative.fa > $TRAINNEGFILE

### merge validation set from each genome #####
VALIDPOSFILE=$DATADIR/fasta/${CELLTYPE}_validPos.fa
cat $DATADIR/fasta/hg38_${CELLTYPE}_valid_positive.fa $DATADIR/fasta/mm10_${CELLTYPE}_valid_positive.fa \
	$DATADIR/fasta/rheMac10_${CELLTYPE}_valid_positive.fa > $VALIDPOSFILE

VALIDNEGFILE=$DATADIR/fasta/${CELLTYPE}_validNeg.fa
cat $DATADIR/fasta/hg38_${CELLTYPE}_valid_negative.fa $DATADIR/fasta/mm10_${CELLTYPE}_valid_negative.fa \
	$DATADIR/fasta/rheMac10_${CELLTYPE}_valid_negative.fa > $VALIDNEGFILE

#### cyclical learning rate parameters
BATCH_SIZE=1000
EPOCHS=23
NUM_CYCLES=2.35
BASE_LR=1e-2
MAX_LR=1e-1
BASE_M=.85
MAX_M=.99
L2_REG=1e-10
# DROPOUT=.2

#############################################
## train CNN models with default parameters
for DROPOUT in .05 .1 .15 .2 .25 .3; do
	python train_singleTask_CNN_classifier_OCP.py --mode 'train' --out_dir $DATADIR \
		--conv_width 11 --l2_reg $L2_REG --dropout $DROPOUT --batch_size $BATCH_SIZE \
		--epochs $EPOCHS --numCycles $NUM_CYCLES --base_lr $BASE_LR --max_lr $MAX_LR \
		--base_m $BASE_M --max_m $MAX_M --verbose 2 --cyclical_momentum \
		$LABEL $TRAINPOSFILE $TRAINNEGFILE $VALIDPOSFILE $VALIDNEGFILE 

	## score validation sequences with default CNN models with default parameters
	python train_singleTask_CNN_classifier_OCP.py --mode 'evaluate' --out_dir $DATADIR \
		--conv_width 11 --l2_reg $L2_REG --dropout $DROPOUT --batch_size $BATCH_SIZE \
		--epochs $EPOCHS --numCycles $NUM_CYCLES --base_lr $BASE_LR --max_lr $MAX_LR \
		--base_m $BASE_M --max_m $MAX_M --verbose 2 --cyclical_momentum \
		$LABEL $TRAINPOSFILE $TRAINNEGFILE $VALIDPOSFILE $VALIDNEGFILE
done
