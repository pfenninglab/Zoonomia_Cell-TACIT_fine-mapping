#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=gpu,pfen3
#SBATCH --job-name=evalCell
#SBATCH --time 12:00:00
#SBATCH --gres=gpu:1
#SBATCH --mem=93G
#SBATCH --array=1-40
#SBATCH --error=logs/calibration_d_nonCelltypeNonEnhBiasAway10x_%A_%a.txt
#SBATCH --output=logs/calibration_d_nonCelltypeNonEnhBiasAway10x_%A_%a.txt

SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/cnn_enhancer_ortholog
DATADIR=${SETDIR}/data/raw_data/cnn_enhancer_ortholog
cd $CODEDIR; 
source activate tf2

NEGSET=nonCelltypeNonEnhBiasAway10x

#############################################
# get the cell type to be used for training #
CELLS=( NULL MSN_D1 MSN_D2 MSN_SN INT_Pvalb Astro Microglia Oligo OPC )
(( FOLD = ($SLURM_ARRAY_TASK_ID -1) % 5 + 1))
(( CELL_IND = ($SLURM_ARRAY_TASK_ID -1) / 5 + 1))
CELLTYPE=${CELLS[$CELL_IND]}
PREFIX=${CELLTYPE}_fold${FOLD}_hgRmMm_${NEGSET}
VALIDPOSFILE=$DATADIR/fasta/${PREFIX}_${NEGSET}_validPos.fa.gz
VALIDNEGFILE=$DATADIR/fasta/${PREFIX}_${NEGSET}_validNeg.fa.gz

#############################################
## score species specific evaluation metrics 
for MODEL in $(ls $DATADIR/models/${CELLTYPE}_fold${FOLD}_hgRmMm*/*.h5); do
# evaluate celltype specificity
python train_singleTask_CNN_classifier_OCP.py --mode 'predict' --out_dir $DATADIR \
	--model_name $MODEL --prefix $PREFIX --predict_out validPositiveForCalibration \
	--predict_fasta $VALIDPOSFILE

python train_singleTask_CNN_classifier_OCP.py --mode 'predict' --out_dir $DATADIR \
	--model_name $MODEL --prefix $PREFIX --predict_out validNegativeForCalibration \
	--predict_fasta $VALIDNEGFILE
done