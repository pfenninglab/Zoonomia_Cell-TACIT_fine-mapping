#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=gpu,pfen3
#SBATCH --job-name=testEval
#SBATCH --time 12:00:00
#SBATCH --gres=gpu:1
#SBATCH --mem=93G
#SBATCH --array=1-40
#SBATCH --error=logs/testEval_nonCelltypeNonEnhBiasAway10x_%A_%a.txt
#SBATCH --output=logs/testEval_nonCelltypeNonEnhBiasAway10x_%A_%a.txt

SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/cnn_enhancer_ortholog
DATADIR=${SETDIR}/data/raw_data/cnn_enhancer_ortholog
cd $CODEDIR; source activate tf2

NEGSET=nonCelltypeNonEnhBiasAway10x

#############################################
# get the cell type to be used for training #
CELLS=( NULL MSN_D1 MSN_D2 MSN_SN INT_Pvalb Astro Microglia Oligo OPC )
(( FOLD = ($SLURM_ARRAY_TASK_ID -1) % 5 + 1))
(( CELL_IND = ($SLURM_ARRAY_TASK_ID -1) / 5 + 1))
CELLTYPE=${CELLS[$CELL_IND]}
PREFIX=${CELLTYPE}_fold${FOLD}_hgRmMm_${NEGSET}


#####################################################################
### merge positive training and validation set from each genome #####
TESTPOSFILE=$DATADIR/fasta/${PREFIX}_${NEGSET}_testPos.fa.gz
if [[ ! -f  $TESTPOSFILE ]]; then
cat $DATADIR/fasta/hg38_${CELLTYPE}_fold${FOLD}_test_positive.fa.gz \
	$DATADIR/fasta/mm10_${CELLTYPE}_fold${FOLD}_test_positive.fa.gz \
	$DATADIR/fasta/rheMac10_${CELLTYPE}_fold${FOLD}_test_positive.fa.gz > $TESTPOSFILE
fi

TESTNEGFILE=$DATADIR/fasta/${PREFIX}_${NEGSET}_testNeg.fa.gz
if [[ ! -f  $TESTNEGFILE ]]; then
cat $DATADIR/fasta/hg38_${CELLTYPE}_fold${FOLD}_test_nonEnhNeg.fa.gz \
	$DATADIR/fasta/mm10_${CELLTYPE}_fold${FOLD}_test_nonEnhNeg.fa.gz \
	$DATADIR/fasta/rheMac10_${CELLTYPE}_fold${FOLD}_test_nonEnhNeg.fa.gz \
	$DATADIR/fasta/hg38_${CELLTYPE}_fold${FOLD}_test_biasAway10x.fa.gz \
	$DATADIR/fasta/mm10_${CELLTYPE}_fold${FOLD}_test_biasAway10x.fa.gz \
	$DATADIR/fasta/rheMac10_${CELLTYPE}_fold${FOLD}_test_biasAway10x.fa.gz \
	$DATADIR/fasta/hg38_${CELLTYPE}_fold${FOLD}_test_nonCelltypeNeg.fa.gz \
	$DATADIR/fasta/mm10_${CELLTYPE}_fold${FOLD}_test_nonCelltypeNeg.fa.gz \
	$DATADIR/fasta/rheMac10_${CELLTYPE}_fold${FOLD}_test_nonCelltypeNeg.fa.gz > $TESTNEGFILE
fi

################################################################
## score test sequences with default CNN models with default parameters
for MODEL in $(ls $DATADIR/models/${CELLTYPE}_fold${FOLD}_hgRmMm*_${NEGSET}/*.h5); do
python train_singleTask_CNN_classifier_OCP.py --mode 'evaluate' --out_dir $DATADIR \
--model_name $MODEL --prefix $PREFIX --predict_out "${NEGSET}_test" \
--valid_fasta_pos $TESTPOSFILE --valid_fasta_neg $TESTNEGFILE --verbose 0
done


###################################################
## score cell type and species specific evalutation
for GENOME in 'hg38' 'mm10' 'rheMac10'; do
for MODEL in $(ls $DATADIR/models/${CELLTYPE}_fold${FOLD}_hg*_${NEGSET}/*.h5); do
# # evaluate species specificity
TESTPOSFILE=$DATADIR/fasta/${GENOME}Only_${CELLTYPE}_fold${FOLD}_test_positive.fa.gz
TESTNEGFILE=$DATADIR/fasta/${GENOME}_${CELLTYPE}_fold${FOLD}_test_nonEnhNeg.fa.gz
python train_singleTask_CNN_classifier_OCP.py --mode 'evaluate' --out_dir $DATADIR \
--model_name $MODEL --prefix $PREFIX --predict_out ${GENOME}Only_nonEnh_test \
--valid_fasta_pos $TESTPOSFILE --valid_fasta_neg $TESTNEGFILE --verbose 0

# evaluate celltype specificity
TESTPOSFILE=$DATADIR/fasta/${GENOME}CelltypeOnly_${CELLTYPE}_fold${FOLD}_test_positive.fa.gz
TESTNEGFILE=$DATADIR/fasta/${GENOME}NonCelltype_${CELLTYPE}_fold${FOLD}_test_negative.fa.gz
python train_singleTask_CNN_classifier_OCP.py --mode 'evaluate' --out_dir $DATADIR \
--model_name $MODEL --prefix $PREFIX --predict_out ${GENOME}CelltypeOnly_NonCelltype_test \
--valid_fasta_pos $TESTPOSFILE --valid_fasta_neg $TESTNEGFILE --verbose 0
done
done


#####################################################################
### merge positive training and validation set for only humans #####
PREFIX=${CELLTYPE}_fold${FOLD}_hg_${NEGSET}
TESTNEGFILE=$DATADIR/fasta/${PREFIX}_${NEGSET}_testNeg.fa.gz
if [[ ! -f  $TESTNEGFILE ]]; then
cat $DATADIR/fasta/hg38_${CELLTYPE}_fold${FOLD}_test_nonEnhNeg.fa.gz \
	$DATADIR/fasta/hg38_${CELLTYPE}_fold${FOLD}_test_biasAway10x.fa.gz \
	$DATADIR/fasta/hg38_${CELLTYPE}_fold${FOLD}_test_nonCelltypeNeg.fa.gz > $TESTNEGFILE
fi

################################################################
## score test sequences with default CNN models with default parameters
for MODEL in $(ls $DATADIR/models/${CELLTYPE}_fold${FOLD}_hg_${NEGSET}/*.h5); do
python train_singleTask_CNN_classifier_OCP.py --mode 'evaluate' --out_dir $DATADIR \
--model_name $MODEL --prefix $PREFIX --predict_out "${NEGSET}_test" \
--valid_fasta_pos $DATADIR/fasta/hg38_${CELLTYPE}_fold${FOLD}_test_positive.fa.gz \
--valid_fasta_neg $TESTNEGFILE --verbose 0 --force
done

