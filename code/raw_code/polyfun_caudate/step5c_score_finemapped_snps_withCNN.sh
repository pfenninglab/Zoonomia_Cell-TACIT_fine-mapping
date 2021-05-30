#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen3
#SBATCH --job-name=scoreSNPs
#SBATCH --time 12:00:00
#SBATCH --gres=gpu:1
#SBATCH --mem=60G
#SBATCH --array=1-40%6
#SBATCH --error=logs/scoreSNPs_%A_%a.txt
#SBATCH --output=logs/scoreSNPs_%A_%a.txt

SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/cnn_enhancer_ortholog
DATADIR=${SETDIR}/data/raw_data/cnn_enhancer_ortholog
CODEDIR2=${SETDIR}/code/raw_code/polyfun_caudate
DATADIR2=${SETDIR}/data/raw_data/polyfun_caudate
cd $CODEDIR2; 
source activate tf2

NEGSET=nonCelltypeNonEnhBiasAway10x

#############################################
# get the cell type to be used for training #
CELLS=( NULL MSN_D1 MSN_D2 MSN_SN INT_Pvalb Astro Microglia Oligo OPC )
(( FOLD = ($SLURM_ARRAY_TASK_ID -1) % 5 + 1))
(( CELL_IND = ($SLURM_ARRAY_TASK_ID -1) / 5 + 1))
CELLTYPE=${CELLS[$CELL_IND]}

EFFECT_FASTA=$DATADIR2/fasta/polyfun_caudate_finemapped_snps_501_effect.fasta
NONEFF_FASTA=$DATADIR2/fasta/polyfun_caudate_finemapped_snps_501_nonEffect.fasta

EFFECT_FASTA2=$DATADIR2/fasta/polyfun_caudate_finemapped_snps_501_effect_${CELLTYPE}.fasta
NONEFF_FASTA2=$DATADIR2/fasta/polyfun_caudate_finemapped_snps_501_nonEffect_${CELLTYPE}.fasta

###############################
## score all SNPs with models
for MODEL in $(ls $DATADIR/models/${CELLTYPE}_fold${FOLD}_hgRmMm*/*${NEGSET}*.h5); do
# effect allele
python ${CODEDIR}/train_singleTask_CNN_classifier_OCP.py --mode 'predict' --out_dir $DATADIR2 \
--model_name $MODEL --predict_out 'effectSNPs' --predict_fasta $EFFECT_FASTA

# noneffect allele
python ${CODEDIR}/train_singleTask_CNN_classifier_OCP.py --mode 'predict' --out_dir $DATADIR2 \
--model_name $MODEL --predict_out 'nonEffSNPs' --predict_fasta $NONEFF_FASTA

# effect allele in celltypes
python ${CODEDIR}/train_singleTask_CNN_classifier_OCP.py --mode 'predict' --out_dir $DATADIR2 \
--model_name $MODEL --predict_out "effectSNP_${CELLTYPE}" --predict_fasta $EFFECT_FASTA2 --force

# noneffect allele in celltypes
python ${CODEDIR}/train_singleTask_CNN_classifier_OCP.py --mode 'predict' --out_dir $DATADIR2 \
--model_name $MODEL --predict_out "nonEffSNP_${CELLTYPE}" --predict_fasta $NONEFF_FASTA2 --force
done



