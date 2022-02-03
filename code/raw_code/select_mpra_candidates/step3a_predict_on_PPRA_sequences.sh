#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --job-name=evalSeq
#SBATCH --time 12:00:00
#SBATCH --job-name=scoreMappable
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --error=logs/score_ppra_sequences_%A_%a.txt
#SBATCH --output=logs/score_ppra_sequences_%A_%a.txt

### Set up the directories
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${SETDIR}/code/raw_code/select_mpra_candidates; 
DATADIR=${SETDIR}/data/raw_data/select_mpra_candidates
cd $CODEDIR; mkdir -p $DATADIR/peaks $DATADIR/fasta
SIZE=501; CUTOFF=0.5; source ~/.bashrc
conda activate tf2

CELLS="MSN_D1 MSN_D2 MSN_SN INT_Pvalb Astro Oligo OPC Microglia"

##############################################
# Predict SPECIES ortholog CELLTYPE activity
for FASTA in $(ls $DATADIR/fasta/*.fa | sed '/MPRAii_/d'); do
OUT_PREF=$(basename $FASTA .fa)
## for every cell type, score
for CELL in $CELLS; do
AVG_PREDICTIONS=$DATADIR/predictions/${OUT_PREF}.${CELL}.predictions.txt.gz
if [[ ! -f $AVG_PREDICTIONS ]]; then
## for every fold, score
for FOLD in {1..5}; do
PREFIX=${CELL}_fold${FOLD}_hgRmMm_nonCelltypeNonEnhBiasAway10x
for MODEL in $(ls ${SETDIR}/data/raw_data/cnn_enhancer_ortholog/models/${PREFIX}/*.h5); do
python ${SETDIR}/code/raw_code/cnn_enhancer_ortholog/train_singleTask_CNN_classifier_OCP.py \
--mode 'predict' --out_dir $DATADIR --predict_fasta $FASTA \
--model_name $MODEL --predict_out ${OUT_PREF}.${CELL} --verbose 2 --force
done; done
## average score across 5 folds, column 3 is prediction score
echo "Averaging predictions across folds."
paste <(zcat $DATADIR/predictions/${CELL}_fold1_hgRmMm_nonCelltypeNonEnhBiasAway10x/*${OUT_PREF}.${CELL}* | cut -f1,3|sed '1d' ) \
<(zcat $DATADIR/predictions/${CELL}_fold2_hgRmMm_nonCelltypeNonEnhBiasAway10x/*${OUT_PREF}.${CELL}* | cut -f3|sed '1d' ) \
<(zcat $DATADIR/predictions/${CELL}_fold3_hgRmMm_nonCelltypeNonEnhBiasAway10x/*${OUT_PREF}.${CELL}* | cut -f3|sed '1d' ) \
<(zcat $DATADIR/predictions/${CELL}_fold4_hgRmMm_nonCelltypeNonEnhBiasAway10x/*${OUT_PREF}.${CELL}* | cut -f3|sed '1d' ) \
<(zcat $DATADIR/predictions/${CELL}_fold5_hgRmMm_nonCelltypeNonEnhBiasAway10x/*${OUT_PREF}.${CELL}* | cut -f3|sed '1d' ) | \
awk -v OFS='\t' '{sum = 0; for (i = 2; i <= NF; i++) sum += $i; sum /= (NF - 1); print $1, sum}' | \
gzip > $AVG_PREDICTIONS
else echo "Found $(basename $AVG_PREDICTIONS). Moving on."; fi
done; done

