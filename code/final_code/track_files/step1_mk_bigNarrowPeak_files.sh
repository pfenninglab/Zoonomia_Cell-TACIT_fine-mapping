#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1
#SBATCH --job-name=mk_trackfile
#SBATCH --time 12:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --error=logs/make_bigBed_tracks_%A_%a.txt
#SBATCH --output=logs/make_bigBed_tracks_%A_%a.txt
#SBATCH --array=1-241%50

### Set up the directories
PROJDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate
CODEDIR=${PROJDIR}/code/final_code/track_files
DATADIR=${PROJDIR}/data/tidy_data/track_files
ZOONOMIADIR=${PROJDIR}/data/tidy_data/Zoonomia_data
FASTADIR=/data/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/200MammalsFastas
PEAKDIR=${PROJDIR}/data/raw_data/hg38/Corces_2020/halper_zoo
PREDDIR=${PROJDIR}/data/raw_data/ldsc_caudate_zoonomia/predictions
cd $CODEDIR; mkdir -p $DATADIR/peaks $DATADIR/bigBed

source ~/.bashrc


#################################################################
# get the species and cell types that need to be scored by CNNs
CELLS="MSN_D1 MSN_D2 MSN_SN INT_Pvalb Astro Oligo OPC Microglia"
SPECIES=$(awk -F'\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $2}' ${ZOONOMIADIR}/tables/200_Mammals_Genome_Information.tsv )
FASTA_FN=/projects/pfenninggroup/machineLearningForComputationalBiology/halLiftover_chains/data/raw_data/2bit/fasta/${SPECIES}.fa
CHROMSIZE=/projects/pfenninggroup/machineLearningForComputationalBiology/halLiftover_chains/data/raw_data/2bit/fasta/sizes.${SPECIES}
if [[ ! -f $FASTA_FN ]]; then  echo "No genome fasta found. Quiting";exit 1; fi  # one genome doesn't have fasta
if [[ ! -f $CHROMSIZE ]]; then 
echo "Chromosome sizes for ${SPECIES} not found. Computing now."
faidx $FASTA_FN -i chromsizes > $CHROMSIZE; 
sort -k 1,1 $CHROMSIZE > ${SPECIES}.tmp.txt
cat ${SPECIES}.tmp.txt > $CHROMSIZE; rm ${SPECIES}.tmp.txt
fi # make the chromosome sizes file

##############################################
# Predict SPECIES ortholog CELLTYPE activity
for CELL in $CELLS; do
echo "Grabbing the HALPER peaks mapped to ${SPECIES} for ${CELL}."
BIGBED=${DATADIR}/bigBed/Corces2020_caudate.${CELL}.${SPECIES}.scaledPredictions.narrowPeak.bb
if [[ ! -f $BIGBED ]]; then
OUT_PEAK=${DATADIR}/peaks/Corces2020_caudate.${CELL}.${SPECIES}.scaledPredictions.narrowPeak
if [[ $SPECIES == 'Homo_sapiens' ]]; then 
BED=${PROJDIR}/data/raw_data/hg38/Corces_2020/peak/Corces2020_caudate.${CELL}.narrowPeak.gz 
else BED=${PEAKDIR}/Corces2020_caudate.${CELL}.Homo_sapiensTo${SPECIES}.HALPER.narrowPeak.gz; fi

echo "Grabbing the Cell TACIT preditions of ${SPECIES} for ${CELL}."
PREDICT_SCORES=${PREDDIR}/Corces2020_caudate.${CELL}.${SPECIES}.avgCNN.predictions.txt.gz

echo "Adding Cell TACIT scores to narrowPeak file of ${SPECIES} for ${CELL}."
## note SCORE 5th column is 0-1 score scaled by 1000 for bigBed format
awk -F"\t" 'BEGIN { OFS="\t" }; FNR==NR{var[$1]=$2;next;} {print $1, $2, $3, $4, int(var[$4] * 1000), "." , $7, $8, $9, $10}' \
<(zcat $PREDICT_SCORES) <(zcat $BED | LC_COLLATE=C sort -k 1,1 -k2,2n ) > $OUT_PEAK

## filter to only chromosomes in the chrom.size file
OUT_PEAK2=${DATADIR}/peaks/Corces2020_caudate.${CELL}.${SPECIES}.filtered.narrowPeak
awk 'NR==FNR {a[$1]++; next} $1 in a' $CHROMSIZE $OUT_PEAK > $OUT_PEAK2

echo "$(basename $BIGBED) not found. Making track bigNarrowPeak from peaks now."
bedToBigBed -as=bigNarrowPeak.as -type=bed6+4 $OUT_PEAK2 $CHROMSIZE $BIGBED
else echo "Found $(basename $BIGBED). Moving on."
fi

## clean up and compress the predictions
gzip -f ${OUT_PEAK2}; rm $OUT_PEAK
printf "\n"
done

