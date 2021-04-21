#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen3
#SBATCH --job-name=vBoth
#SBATCH --gres=gpu:1
#SBATCH --mem=63G
#SBATCH --array=1-40%2
#SBATCH --error=logs/cnn_ocp_c_nonEnhBiasAway10x_%A_%a.txt
#SBATCH --output=logs/cnn_ocp_c_nonEnhBiasAway10x_%A_%a.txt

SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/cnn_enhancer_ortholog
DATADIR=${SETDIR}/data/raw_data/cnn_enhancer_ortholog
cd $CODEDIR; source activate tf2

NEGSET=nonEnhBiasAway10x

#############################################
# get the cell type to be used for training #
CELLS=( NULL MSN_D1 MSN_D2 MSN_SN INT_Pvalb Astro Microglia Oligo OPC )
(( FOLD = ($SLURM_ARRAY_TASK_ID -1) % 5 + 1))
(( CELL_IND = ($SLURM_ARRAY_TASK_ID -1) / 5 + 1))
CELLTYPE=${CELLS[$CELL_IND]}
PREFIX=${CELLTYPE}_fold${FOLD}_hgRmMm_${NEGSET}

#####################################################################
### merge positive training and validation set from each genome #####
TRAINPOSFILE=$DATADIR/fasta/${PREFIX}_${NEGSET}_trainPos.fa
if [[ ! -f  $TRAINPOSFILE ]]; then
cat $DATADIR/fasta/hg38_${CELLTYPE}_fold${FOLD}_train_positive.fa \
	$DATADIR/fasta/mm10_${CELLTYPE}_fold${FOLD}_train_positive.fa \
	$DATADIR/fasta/rheMac10_${CELLTYPE}_fold${FOLD}_train_positive.fa > $TRAINPOSFILE
fi

VALIDPOSFILE=$DATADIR/fasta/${PREFIX}_${NEGSET}_validPos.fa
if [[ ! -f  $VALIDPOSFILE ]]; then
cat $DATADIR/fasta/hg38_${CELLTYPE}_fold${FOLD}_valid_positive.fa \
	$DATADIR/fasta/mm10_${CELLTYPE}_fold${FOLD}_valid_positive.fa \
	$DATADIR/fasta/rheMac10_${CELLTYPE}_fold${FOLD}_valid_positive.fa > $VALIDPOSFILE
fi


#####################################################################
### merge negative training and validation set from each genome #####
TRAINNEGFILE=$DATADIR/fasta/${PREFIX}_${NEGSET}_trainNeg.fa
if [[ ! -f  $TRAINNEGFILE ]]; then
cat $DATADIR/fasta/hg38_${CELLTYPE}_fold${FOLD}_train_nonEnhNeg.fa \
	$DATADIR/fasta/mm10_${CELLTYPE}_fold${FOLD}_train_nonEnhNeg.fa \
	$DATADIR/fasta/rheMac10_${CELLTYPE}_fold${FOLD}_train_nonEnhNeg.fa \
	$DATADIR/fasta/hg38_${CELLTYPE}_fold${FOLD}_train_biasAway10x.fa \
	$DATADIR/fasta/mm10_${CELLTYPE}_fold${FOLD}_train_biasAway10x.fa \
	$DATADIR/fasta/rheMac10_${CELLTYPE}_fold${FOLD}_train_biasAway10x.fa > $TRAINNEGFILE
fi

VALIDNEGFILE=$DATADIR/fasta/${PREFIX}_${NEGSET}_validNeg.fa
if [[ ! -f  $VALIDNEGFILE ]]; then
cat $DATADIR/fasta/hg38_${CELLTYPE}_fold${FOLD}_valid_nonEnhNeg.fa \
	$DATADIR/fasta/mm10_${CELLTYPE}_fold${FOLD}_valid_nonEnhNeg.fa \
	$DATADIR/fasta/rheMac10_${CELLTYPE}_fold${FOLD}_valid_nonEnhNeg.fa \
	$DATADIR/fasta/hg38_${CELLTYPE}_fold${FOLD}_valid_biasAway10x.fa \
	$DATADIR/fasta/mm10_${CELLTYPE}_fold${FOLD}_valid_biasAway10x.fa \
	$DATADIR/fasta/rheMac10_${CELLTYPE}_fold${FOLD}_valid_biasAway10x.fa > $VALIDNEGFILE
fi

#### cyclical learning rate parameters
BATCH_SIZE=1000
EPOCHS=23
NUM_CYCLES=2.35
BASE_LR=1e-2
MAX_LR=1e-1
BASE_M=.85
MAX_M=.99
L2_REG=1e-10
DROPOUT=.2

#############################################
## train CNN models with default parameters
for DROPOUT in .25; do
python train_singleTask_CNN_classifier_OCP.py --mode 'train' --out_dir $DATADIR \
	--conv_width 11 --l2_reg $L2_REG --dropout $DROPOUT --batch_size $BATCH_SIZE \
	--epochs $EPOCHS --numCycles $NUM_CYCLES --base_lr $BASE_LR --max_lr $MAX_LR \
	--base_m $BASE_M --max_m $MAX_M --verbose 2 --cyclical_momentum \
	--train_fasta_pos $TRAINPOSFILE --train_fasta_neg $TRAINNEGFILE \
	--valid_fasta_pos $VALIDPOSFILE --valid_fasta_neg $VALIDNEGFILE \
	--prefix $PREFIX

## score validation sequences with default CNN models with default parameters
python train_singleTask_CNN_classifier_OCP.py --mode 'evaluate' --out_dir $DATADIR \
	--conv_width 11 --l2_reg $L2_REG --dropout $DROPOUT --batch_size $BATCH_SIZE \
	--epochs $EPOCHS --numCycles $NUM_CYCLES --base_lr $BASE_LR --max_lr $MAX_LR \
	--base_m $BASE_M --max_m $MAX_M --verbose 2 --cyclical_momentum \
	--train_fasta_pos $TRAINPOSFILE --train_fasta_neg $TRAINNEGFILE \
	--valid_fasta_pos $VALIDPOSFILE --valid_fasta_neg $VALIDNEGFILE \
	--prefix $PREFIX
done


