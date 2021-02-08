#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen3
#SBATCH --time=0-12
#SBATCH --job-name=cnn_train
#SBATCH --gres=gpu:1
#SBATCH --mem=95G
#SBATCH --array=52-55%4
#SBATCH --error=logs/cnn_ocp_v2_training_%A_%a_output.txt
#SBATCH --output=logs/cnn_ocp_v2_training_%A_%a_output.txt

PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/addiction_gwas_enrichment/celltype_specific_ml_models
cd $PROJDIR

# cell type and fold index
CELLTYPE=`awk -v IND="${SLURM_ARRAY_TASK_ID}" 'BEGIN { FS = " " } ;NR==IND + 1 {print $1}' tables/celltype_model_reproducible_vs_gcNegatives10x_5-fold-CV.txt`
FOLD=`awk -v IND="${SLURM_ARRAY_TASK_ID}" 'BEGIN { FS = " " } ;NR==IND + 1 {print $2}' tables/celltype_model_reproducible_vs_gcNegatives10x_5-fold-CV.txt`

# input training files 
TRAINPOSFILE=$PROJDIR/FASTA_CV/${CELLTYPE}_${FOLD}_trainPos.fa
TRAINNEGFILE=$PROJDIR/FASTA_CV/${CELLTYPE}_${FOLD}_trainNeg10x.fa

# input validation files for 10-fold cross evalutation
VALIDPOSFILE=$PROJDIR/FASTA_CV/${CELLTYPE}_${FOLD}_validPos.fa
VALIDNEGFILE=$PROJDIR/FASTA_CV/${CELLTYPE}_${FOLD}_validNeg10x.fa

LABEL=${CELLTYPE}_${FOLD}

# cyclical learning rate parameters
BATCH_SIZE=1000
EPOCHS=23
NUM_CYCLES=2.35
BASE_LR=1e-2
MAX_LR=1e-1
BASE_M=.85
MAX_M=.99
L2_REG=1e-10

source activate tf2
### train CNN models with default parameters
python train_singleTask_CNN_classifier_OCP_v2.py \
	--conv_width 11 --l2_reg $L2_REG \
	--batch_size $BATCH_SIZE --cyclical_momentum \
	--epochs $EPOCHS --numCycles $NUM_CYCLES \
	--base_lr $BASE_LR --max_lr $MAX_LR \
	--base_m $BASE_M --max_m $MAX_M \
	--mode 'train' --verbose 0 \
	$LABEL $TRAINPOSFILE $TRAINNEGFILE $VALIDPOSFILE $VALIDNEGFILE 


### score validation sequences with default CNN models with default parameters
python train_singleTask_CNN_classifier_OCP_v2.py \
	--conv_width 11 --l2_reg $L2_REG \
	--batch_size $BATCH_SIZE --cyclical_momentum \
	--epochs $EPOCHS --numCycles $NUM_CYCLES \
	--base_lr $BASE_LR --max_lr $MAX_LR \
	--base_m $BASE_M --max_m $MAX_M \
	--mode 'evaluate' --verbose 0 \
	$LABEL $TRAINPOSFILE $TRAINNEGFILE $VALIDPOSFILE $VALIDNEGFILE 

## input validation files for 10-fold cross evalutation
# TESTPOSFILE=$PROJDIR/FASTA_CV/${CELLTYPE}_testPos.fa
# TESTNEGFILE=$PROJDIR/FASTA_CV/${CELLTYPE}_testNeg10x.fa
# LABEL=${CELLTYPE}_${FOLD}_training

# # score test sequences with default CNN models with default parameters
# python train_singleTask_CNN_classifier_OCP_v2.py \
# 	--conv_width 11 --l2_reg $L2_REG \
# 	--batch_size $BATCH_SIZE --cyclical_momentum \
# 	--epochs $EPOCHS --numCycles $NUM_CYCLES \
# 	--base_lr $BASE_LR --max_lr $MAX_LR \
# 	--base_m $BASE_M --max_m $MAX_M \
# 	--mode 'evaluate' --verbose 1 \
# 	$LABEL $TRAINPOSFILE $TRAINNEGFILE $TESTPOSFILE $TESTNEGFILE 



