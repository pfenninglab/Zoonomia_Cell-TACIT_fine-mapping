#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=gpu,pfen3
#SBATCH --gres=gpu:1
#SBATCH --time 12:00:00
#SBATCH --job-name=scoreRA
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=45G
#SBATCH --error=logs/score_reporter_enhancer_zoonomia_%A_%a.txt
#SBATCH --output=logs/score_reporter_enhancer_zoonomia_%A_%a.txt
#SBATCH --array=1-241%10

### Set up the directories
SETDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
GWASDIR=${SETDIR}/data/tidy_data/ldsc_gwas; SNPLIST=${GWASDIR}/listHM3.noMHC.txt
ZOONOMIADIR=${SETDIR}/data/tidy_data/Zoonomia_data
FASTADIR=/data/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/200MammalsFastas
PEAKDIR=${SETDIR}/data/raw_data/hg38/Corces_2020/halper_zoo
CODEDIR=${SETDIR}/code/raw_code/reporter_assay; 
DATADIR=${SETDIR}/data/raw_data/reporter_assay
cd $CODEDIR; mkdir -p $DATADIR/peaks $DATADIR/fasta
SIZE=501; CUTOFF=0.5; source ~/.bashrc


## the HSA and HSB peaks from MSN_D2
HSA='hg38:chr11:113567061-113567561:250'
HSB='hg38:chr11:113577668-113578168:250'
CELL1='MSN_D2'

CELLS="MSN_D1 MSN_D2 MSN_SN INT_Pvalb Astro Oligo OPC Microglia"
SPECIES=$(awk -F'\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $2}' ${ZOONOMIADIR}/tables/200_Mammals_Genome_Information.tsv )
FASTA_FN=/projects/pfenninggroup/machineLearningForComputationalBiology/halLiftover_chains/data/raw_data/fasta/${SPECIES}.fa
if [[ ! -f $FASTA_FN ]]; then  echo "No genome fasta found. Quiting";exit 1; fi  # one genome doesn't have fasta


##############################################
# Predict SPECIES ortholog CELLTYPE activity
if [[ $SPECIES == 'Homo_sapiens' ]]; then 
BED=${SETDIR}/data/raw_data/hg38/Corces_2020/peak/Corces2020_caudate.${CELL1}.narrowPeak.gz 
else BED=${PEAKDIR}/Corces2020_caudate.${CELL1}.Homo_sapiensTo${SPECIES}.HALPER.narrowPeak.gz; fi
if [[ ! -f $BED ]]; then echo "Cannot find ${BED}."; return; fi

# 1) get temporary bed with the right chr names, filter for just the HSA and HSB enhancers
TMPBED=${DATADIR}/peaks/Zoonomia_reporter_assay.${SPECIES}.HALPER.narrowPeak.renamedBed
zcat $BED | awk -v CHAR1=${HSA} -v CHAR2=${HSB} '$4==CHAR1 ||$4==CHAR2 ' > $TMPBED

# 2) get fasta sequences from mapped peaks
echo "Getting fasta for $(basename $BED) file from $(basename $FASTA_FN)."
TMPFASTA=${DATADIR}/fasta/Zoonomia_reporter_assay.${SPECIES}.HALPER.narrowPeak.tmp.fasta.gz
bedtools getfasta -nameOnly -fi $FASTA_FN -bed $TMPBED | gzip > $TMPFASTA

if [[ $(zcat $TMPFASTA| wc -l) -gt 0 ]]; then 
for CELL in $CELLS; do
AVG_PREDICTIONS=${DATADIR}/predictions/Zoonomia_reporter_assay.${CELL}.${SPECIES}.avgCNN.predictions.txt.gz
# if [[ ! -f $AVG_PREDICTIONS ]]; then
echo "$(basename $AVG_PREDICTIONS) not found. Scoring peaks now."
# 3) score with CNNs of this cell type across 5 folds
conda activate tf2
for FOLD in {1..5}; do
PREFIX=${CELL}_fold${FOLD}_hgRmMm_nonCelltypeNonEnhBiasAway10x
for MODEL in $(ls ${SETDIR}/data/raw_data/cnn_enhancer_ortholog/models/${PREFIX}/*.h5); do
python ${SETDIR}/code/raw_code/cnn_enhancer_ortholog/train_singleTask_CNN_classifier_OCP.py \
--mode 'predict' --out_dir $DATADIR --predict_fasta $TMPFASTA \
--model_name $MODEL --predict_out Zoonomia_reporter_assay.${SPECIES}.${CELL} --verbose 2 --force
done
done
conda deactivate

# 4) average score across 5 folds, column 3 is prediction score
echo "Averaging predictions across folds."
paste <(zcat $DATADIR/predictions/${CELL}_fold1_hgRmMm_nonCelltypeNonEnhBiasAway10x/*Zoonomia_reporter_assay.${SPECIES}.${CELL}* | cut -f1,3|sed '1d' ) \
<(zcat $DATADIR/predictions/${CELL}_fold2_hgRmMm_nonCelltypeNonEnhBiasAway10x/*Zoonomia_reporter_assay.${SPECIES}.${CELL}* | cut -f3|sed '1d' ) \
<(zcat $DATADIR/predictions/${CELL}_fold3_hgRmMm_nonCelltypeNonEnhBiasAway10x/*Zoonomia_reporter_assay.${SPECIES}.${CELL}* | cut -f3|sed '1d' ) \
<(zcat $DATADIR/predictions/${CELL}_fold4_hgRmMm_nonCelltypeNonEnhBiasAway10x/*Zoonomia_reporter_assay.${SPECIES}.${CELL}* | cut -f3|sed '1d' ) \
<(zcat $DATADIR/predictions/${CELL}_fold5_hgRmMm_nonCelltypeNonEnhBiasAway10x/*Zoonomia_reporter_assay.${SPECIES}.${CELL}* | cut -f3|sed '1d' ) | \
awk -v OFS='\t' '{sum = 0; for (i = 2; i <= NF; i++) sum += $i; sum /= (NF - 1); print $1, sum}' | \
gzip > $AVG_PREDICTIONS
# else echo "Found $(basename $AVG_PREDICTIONS). Moving on."
# fi 

PRED_ACTIVE_BED=${DATADIR}/peaks/Zoonomia_reporter_assay.${CELL}.${SPECIES}.predActive.bed.gz
if [[ ! -f $PRED_ACTIVE_BED ]]; then
# 5) get the hg38 coorindates of predicted active peaks in this species
echo "Extract the species predicted active peaks."
zcat $AVG_PREDICTIONS | awk -v CUTOFF=$CUTOFF '$2 > CUTOFF {print $1}' | sed 's/hg38://g'| \
awk '{print $1}'| tr ':' '\t'| tr '-' '\t'| cut -f1-3 | sort -k1,1 -k2,2n | gzip > $PRED_ACTIVE_BED
else echo "Found $(basename $PRED_ACTIVE_BED). Moving on."
fi
done
else
echo "No mappable peaks found for ${SPECIES}."
fi

rm $TMPFASTA $TMPBED

