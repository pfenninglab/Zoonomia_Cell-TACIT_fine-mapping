#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=gpu,pfen3
#SBATCH --gres=gpu:1
#SBATCH --job-name=evalCell
#SBATCH --time 12:00:00
#SBATCH --job-name=scoreMappable
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --error=logs/score_mappable_peaks_%A_%a.txt
#SBATCH --output=logs/score_mappable_peaks_%A_%a.txt
#SBATCH --array=100-240%10

### Set up the directories
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
GWASDIR=${SETDIR}/data/tidy_data/ldsc_gwas; SNPLIST=${GWASDIR}/listHM3.noMHC.txt
CODEDIR=${SETDIR}/code/raw_code/ldsc_caudate_zoonomia; 
DATADIR=${SETDIR}/data/raw_data/ldsc_caudate_zoonomia
ZOONOMIADIR=${SETDIR}/data/tidy_data/Zoonomia_data
FASTADIR=/data/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/200MammalsFastas
PEAKDIR=${SETDIR}/data/raw_data/hg38/Corces_2020/halper_zoo
cd $CODEDIR; mkdir -p $DATADIR/peaks $DATADIR/fasta
SIZE=501; CUTOFF=0.5; source ~/.bashrc

#################################################################
# get the species and cell types that need to be scored by CNNs
CELLS="MSN_D1 MSN_D2 MSN_SN INT_Pvalb Astro Oligo OPC Microglia"
SPECIES=$(awk -F'\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $2}' ${ZOONOMIADIR}/tables/200_Mammals_Genome_Information.tsv )
FASTA_FN=$(awk -F'\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $23}' ${ZOONOMIADIR}/tables/200_Mammals_Genome_Information.tsv | sed $'s/[^[:print:]\t]//g')
if [[ ! -f $FASTA_FN ]]; then exit 1; fi  # one genome doesn't have fasta
IS_CHR_SAME=$(awk -F'\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $20 == $21 }' ${ZOONOMIADIR}/tables/200_Mammals_Genome_Information.tsv)
DICT=$(ls ${FASTADIR}/*${SPECIES}.txt | head -1)

##############################################
# Predict SPECIES ortholog CELLTYPE activity
for CELL in $CELLS; do
echo "Scoring peaks mapped to ${SPECIES} for ${CELL}."
BED=${PEAKDIR}/Corces2020_caudate.${CELL}.Homo_sapiensTo${SPECIES}.HALPER.narrowPeak.gz
if [[ ! -f $BED ]]; then echo "Cannot find ${BED}."; return; fi

# 1) get fasta sequences from mapped peaks
TMPBED=${DATADIR}/peaks/Corces2020_caudate.${CELL}.Homo_sapiensTo${SPECIES}.HALPER.narrowPeak.renamedBed
if [[ $IS_CHR_SAME == 0 ]]; then
echo "Converting chr names with $(basename $DICT)."
zcat $BED | awk -v OFS='\t' 'FNR==NR {dict[$1]=$2; next} {$1=($1 in dict) ? dict[$1] : $1}1' $DICT - > $TMPBED
else zcat $BED > $TMPBED; fi

# 2) get fasta sequences from mapped peaks
echo "Getting fasta for $(basename $BED) file from $(basename $FASTA_FN)."
TMPFASTA=${DATADIR}/fasta/Corces2020_caudate.${CELL}.Homo_sapiensTo${SPECIES}.HALPER.narrowPeak.tmpfasta.gz
bedtools getfasta -nameOnly -fi $FASTA_FN -bed $TMPBED | gzip > $TMPFASTA

# 3) score with CNNs of this cell type across 5 folds
conda activate tf2
for FOLD in {1..5}; do
PREFIX=${CELL}_fold${FOLD}_hgRmMm_nonCelltypeNonEnhBiasAway10x
for MODEL in $(ls ${SETDIR}/data/raw_data/cnn_enhancer_ortholog/models/${PREFIX}/*.h5); do
python ${SETDIR}/code/raw_code/cnn_enhancer_ortholog/train_singleTask_CNN_classifier_OCP.py \
--mode 'predict' --out_dir $DATADIR --predict_fasta $TMPFASTA \
--model_name $MODEL --predict_out ${SPECIES}.${CELL} --verbose 2
done; done
conda deactivate

# 4) average score across 5 folds, column 3 is prediction score
AVG_PREDICTIONS=$DATADIR/predictions/Corces2020_caudate.${CELL}.Homo_sapiensTo${SPECIES}.avgCNN.predictions.txt.gz
paste <(zcat $DATADIR/predictions/${CELL}_fold1_hgRmMm_nonCelltypeNonEnhBiasAway10x/*${SPECIES}.${CELL}* | cut -f1,3|sed '1d;$d' ) \
<(zcat $DATADIR/predictions/${CELL}_fold2_hgRmMm_nonCelltypeNonEnhBiasAway10x/*${SPECIES}.${CELL}* | cut -f3|sed '1d;$d' ) \
<(zcat $DATADIR/predictions/${CELL}_fold3_hgRmMm_nonCelltypeNonEnhBiasAway10x/*${SPECIES}.${CELL}* | cut -f3|sed '1d;$d' ) \
<(zcat $DATADIR/predictions/${CELL}_fold4_hgRmMm_nonCelltypeNonEnhBiasAway10x/*${SPECIES}.${CELL}* | cut -f3|sed '1d;$d' ) \
<(zcat $DATADIR/predictions/${CELL}_fold5_hgRmMm_nonCelltypeNonEnhBiasAway10x/*${SPECIES}.${CELL}* | cut -f3|sed '1d;$d' ) | \
awk -v OFS='\t' '{sum = 0; for (i = 2; i <= NF; i++) sum += $i; sum /= (NF - 1); print $1, sum}' | \
gzip > $AVG_PREDICTIONS

# 5) get the hg38 coorindates of predicted active peaks in this species
PRED_ACTIVE_BED=${DATADIR}/peaks/Corces2020_caudate.${CELL}.Homo_sapiensTo${SPECIES}.predActive.bed.gz
zcat $AVG_PREDICTIONS | awk -v CUTOFF=$CUTOFF '$2 > CUTOFF {print $1}' | sed 's/hg38://g'| \
awk '{print $1}'| tr ':' '\t'| tr '-' '\t'| cut -f1-3 | sort -k1,1 -k2,2n | gzip > $PRED_ACTIVE_BED

rm $TMPFASTA $TMPBED
done




